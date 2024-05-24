import convert_functions
import numpy as np


class HydrogenProductionSystem:

    def __init__(self, time_left, delivery_period, H2MassRemaining, battery_left,
                ltp_PF_wind, ltp_PF_solar, ltp_price, ltp_CO2Intensity):

        # The delivery period is in hours
        # The times series are supposed to be dimensioned correctly

        import pypsa
        import convert_functions
        self.time_left = time_left
        self.delivery_period = delivery_period
        self.H2MassRemaining = H2MassRemaining
        self.battery_left = battery_left

        self.ltp_PF_wind = ltp_PF_wind
        self.ltp_PF_solar = ltp_PF_solar
        self.ltp_price = ltp_price
        self.ltp_CO2Intensity = ltp_CO2Intensity

        # Creation of the network
        n = pypsa.Network()


        # Buses creation
        n.add("Bus", "ACBus")
        n.add("Bus", "DCBus")
        n.add("Bus", "H2Bus")
        n.add("Bus", "BatteryBus")
        n.add("Bus", "NetworkBus")

        # Generators creation
        self.installed_power = 1000
        n.add("Generator", "Wind", bus="ACBus", p_nom=self.installed_power)
        n.add("Generator", "Solar", bus="DCBus", p_nom=self.installed_power)
        n.add("Generator", "GlobalNetwork", bus="NetworkBus", p_nom_extendable=True, p_max_pu=1000, p_min_pu=-1000)

        # Hydrolyzer creation
        self.hydrolyzer_capacity = 1000
        n.add("Store", "H2gen", bus="H2Bus")
        n.add("Load", "H2gen", bus="H2Bus")


        # Battery creation
        battery_capacity = 1000
        n.add("Store", "Battery", bus="BatteryBus", e_nom_extendable=True, e_nom_max=battery_capacity)

        # Link values initialization
        self.eff_charge = 0.9
        self.eff_discharge = 0.9
        self.eff_electrolysis = 0.6
        self.eff_converter = 0.95

        # Links creation
        n.add("Link", "H2Link", bus0="ACBus", bus1="H2Bus", efficiency=self.eff_electrolysis, p_nom_extendable=True)
        n.add("Link", "ChargeLink", bus0="DCBus", bus1="BatteryBus", efficiency=self.eff_charge, p_nom_extendable=True)
        n.add("Link", "DischargeLink", bus0="BatteryBus", bus1="DCBus", efficiency=self.eff_discharge, p_nom_extendable=True)
        n.add("Link", "ACtoDCLink", bus0="ACBus", bus1="DCBus", efficiency=self.eff_converter, p_nom_extendable=True)
        n.add("Link", "DCtoACLink", bus0="DCBus", bus1="ACBus", efficiency=self.eff_converter, p_nom_extendable=True)
        n.add("Link", "BuyLink", bus0="NetworkBus", bus1="ACBus", p_nom_extendable=True)
        n.add("Link", "SellLink", bus0="ACBus", bus1="NetworkBus", p_nom_extendable=True)

        # Optimization parameters initialization

        self.alpha = 0.5
        self.LHV = 33.3    # kWh/kg


        self.n = n

        # initializing global data lists
        self.ltp_pf_series = []
        self.ltp_H2_series = []
        self.ltp_H2_target = []
        self.dp_pf_series = []
        self.dp_Hydro_plan = []
        self.dr_pf_series = []
        self.dr_H2_prod = []

        self.ltp_target_mass = 0
        self.dp_hydro_hourly_schedule = 0
        self.dp_battery_hourly_schedule = 0
        self.total_production = 0
        self.CO2_emissions = 0
        self.electricity_cost = 0
        self.electricity_revenue = 0
        self.electricity_balance = 0
        self.opportunity_cost = 0



    def lt_planner(self):

        n = self.n
        time_left = self.time_left
        delivery_period = self.delivery_period
        n.set_snapshots(range(delivery_period+24))
        n.stores.e_nom_extendable["H2gen"] = True
        n.loads.p_set["H2gen"] = 0
        n.stores.e_initial["Battery"] = self.battery_left

        n.generators_t.p_max_pu["Solar"] = self.ltp_PF_solar
        n.generators_t.p_max_pu["Wind"] = self.ltp_PF_wind
        n.generators_t.p_min_pu["Solar"] = self.ltp_PF_solar
        n.generators_t.p_min_pu["Wind"] = self.ltp_PF_wind

        m = n.optimize.create_model()

        # total_hydrolyzer_energy = m.variables["Store-e"].loc[time_left+23, "H2gen"]
        # total_H2_required_energy = convert_functions.H2_to_P(self.H2MassRemaining, self.LHV)
        # cstr = total_hydrolyzer_energy == total_H2_required_energy
        # m.add_constraints(cstr, name="target match")

        total_hydrolyzer_energy = -m.variables["Store-p"].loc[delivery_period+24-(10+time_left), "H2gen"]
        for i in range(delivery_period+24-(10+time_left)+1, delivery_period+24-10):
            total_hydrolyzer_energy = total_hydrolyzer_energy - m.variables["Store-p"].loc[i, "H2gen"]
        total_H2_required_energy = convert_functions.H2_to_P(self.H2MassRemaining, self.LHV)
        cstr = total_hydrolyzer_energy == total_H2_required_energy
        m.add_constraints(cstr, name="target match")

        # total_hydrolyzer_energy = m.variables["Link-p"].loc[0, "H2Link"]*self.eff_electrolysis
        # for i in range(1, time_left+23):
        #     total_hydrolyzer_energy = total_hydrolyzer_energy + m.variables["Link-p"].loc[i, "H2Link"]*self.eff_electrolysis
        # total_H2_required_energy = convert_functions.H2_to_P(self.H2MassRemaining, self.LHV)
        # cstr = total_hydrolyzer_energy == total_H2_required_energy
        # m.add_constraints(cstr, name="target match")

        for i in range(delivery_period+24):
            cstr = m.variables["Link-p"].loc[i, "H2Link"] <= self.hydrolyzer_capacity
            m.add_constraints(cstr, name="limited capacity at time {}".format(i))

        initial_energy = m.variables["Store-e"].loc[0, "Battery"]
        final_energy = m.variables["Store-e"].loc[delivery_period+23, "Battery"]
        cstr = initial_energy == final_energy
        m.add_constraints(cstr, name="battery sustainability")

        expr = m.variables["Link-p"].loc[0, "BuyLink"] * (self.alpha * self.ltp_CO2Intensity[0] + (1 - self.alpha) * self.ltp_price[0])
        for i in range(1, delivery_period+24):
            expr = expr + m.variables["Link-p"].loc[i, "BuyLink"] * (self.alpha * self.ltp_CO2Intensity[i] + (1 - self.alpha) * self.ltp_price[i])

        m.add_objective(expr, overwrite=True, sense="min")

        n.optimize.solve_model(solver_name="gurobi")

        self.ltp_target_mass = convert_functions.P_to_H2((n.stores_t.p["H2gen"][delivery_period+24-(10+time_left):delivery_period+24-(10+time_left)+24].sum()), self.LHV)
        self.ltp_H2_target = self.ltp_H2_target + [self.ltp_target_mass]

        self.ltp_H2_series = self.ltp_H2_series + [convert_functions.P_to_H2(n.stores_t.p["H2gen"], self.LHV)]

        self.ltp_pf_series = self.ltp_pf_series + [[n.generators_t.p, n.loads_t.p, n.stores_t.p, n.stores_t.e]]


    def daily_planner(self):
        # Time series must be 34h long

        n = self.n
        n.set_snapshots(range(34))
        n.stores.e_nom_extendable["H2gen"] = True
        n.loads.p_set["H2gen"] = 0
        n.stores.e_initial["Battery"] = self.battery_left
        delivery_period = self.delivery_period
        time_left = self.time_left

        dp_PF_wind = self.ltp_PF_wind[delivery_period+24-(10+time_left):delivery_period+24-(10+time_left)+34]
        dp_PF_solar = self.ltp_PF_solar[delivery_period+24-(10+time_left):delivery_period+24-(10+time_left)+34]
        dp_price = self.ltp_price[delivery_period+24-(10+time_left):delivery_period+24-(10+time_left)+34]
        dp_CO2Intensity = self.ltp_CO2Intensity[delivery_period+24-(10+time_left):delivery_period+24-(10+time_left)+34]

        n.generators_t.p_max_pu["Solar"] = dp_PF_solar
        n.generators_t.p_max_pu["Wind"] = dp_PF_wind
        n.generators_t.p_min_pu["Solar"] = dp_PF_solar
        n.generators_t.p_min_pu["Wind"] = dp_PF_wind


        m = n.optimize.create_model()

        total_hydrolyzer_energy = m.variables["Store-p"].loc[0, "H2gen"]
        for i in range(1, 24):
            total_hydrolyzer_energy = total_hydrolyzer_energy + m.variables["Store-p"].loc[i, "H2gen"]
        total_H2_required_energy = convert_functions.H2_to_P(self.ltp_target_mass, self.LHV)
        cstr = total_hydrolyzer_energy == total_H2_required_energy
        m.add_constraints(cstr, name="target match")

        initial_energy = m.variables["Store-e"].loc[0, "Battery"]
        final_energy = m.variables["Store-e"].loc[33, "Battery"]
        cstr = initial_energy == final_energy
        m.add_constraints(cstr, name="battery sustainability")


        # total_hydrolyzer_energy = m.variables["Link-p"].loc[0, "H2Link"]*self.eff_electrolysis
        # for i in range(1, 23):
        #     total_hydrolyzer_energy = total_hydrolyzer_energy + m.variables["Link-p"].loc[i, "H2Link"]*self.eff_electrolysis
        # total_H2_required_energy = convert_functions.H2_to_P(self.ltp_target_mass, self.LHV)
        # cstr = total_hydrolyzer_energy == total_H2_required_energy
        # m.add_constraints(cstr, name="target match")

        for i in range(34):
            cstr = m.variables["Link-p"].loc[i, "H2Link"] <= self.hydrolyzer_capacity
            m.add_constraints(cstr, name="limited capacity at time {}".format(i))

        expr = m.variables["Link-p"].loc[0, "BuyLink"] * (self.alpha * dp_CO2Intensity[0] + (1 - self.alpha) * dp_price[0])
        for i in range(1, 34):
            expr = expr + m.variables["Link-p"].loc[i, "BuyLink"] * (self.alpha * dp_CO2Intensity[i] + (1 - self.alpha) * dp_price[i])

        # total_energy = m.variables["Link-p"].loc[0, "BuyLink"]
        # for i in range(1, 33):
        #     total_energy = total_energy + m.variables["Link-p"].loc[i, "BuyLink"]
        # expr = total_energy * (self.alpha * dp_CO2Intensity + (1 - self.alpha) * dp_price)

        m.add_objective(expr, overwrite=True, sense="min")

        n.optimize.solve_model(solver_name="gurobi")

        self.dp_hydro_hourly_schedule = -n.stores_t.p["H2gen"][0:24]
        # self.dp_battery_hourly_schedule = n.stores_t.p["Battery"][0:24]

        self.dp_pf_series = self.dp_pf_series + [[n.generators_t.p, n.loads_t.p, n.stores_t.p, n.stores_t.e]]

        self.dp_Hydro_plan = self.dp_Hydro_plan +[self.dp_hydro_hourly_schedule]



    # Time series must be 24h long
    def realisation(self):

        n = self.n
        delivery_period = self.delivery_period
        time_left = self.time_left
        n.set_snapshots(range(24))
        n.stores.e_nom_extendable["H2gen"] = False
        n.stores.e_initial["Battery"] = self.battery_left

        dr_PF_wind = self.ltp_PF_wind[delivery_period+24-(10+time_left):delivery_period+24-(10+time_left)+24]
        dr_PF_solar = self.ltp_PF_solar[delivery_period+24-(10+time_left):delivery_period+24-(10+time_left)+24]
        dr_price = self.ltp_price[delivery_period+24-(10+time_left):delivery_period+24-(10+time_left)+24]
        dr_CO2Intensity = self.ltp_CO2Intensity[delivery_period+24-(10+time_left):delivery_period+24-(10+time_left)+24]

        n.loads_t.p_set["H2gen"] = self.dp_hydro_hourly_schedule
        # n.stores_t.p_set["Battery"] = self.dp_battery_hourly_schedule

        n.generators_t.p_max_pu["Solar"] = dr_PF_solar
        n.generators_t.p_max_pu["Wind"] = dr_PF_wind
        n.generators_t.p_min_pu["Solar"] = dr_PF_solar
        n.generators_t.p_min_pu["Wind"] = dr_PF_wind

        m = n.optimize.create_model()

        initial_energy = m.variables["Store-e"].loc[0, "Battery"]
        final_energy = m.variables["Store-e"].loc[23, "Battery"]
        cstr = initial_energy == final_energy
        m.add_constraints(cstr, name="battery sustainability")

        for i in range(0, 24):
            cstr = m.variables["Link-p"].loc[i, "H2Link"] <= self.hydrolyzer_capacity
            m.add_constraints(cstr, name="limited capacity at time {}".format(i))

        expr = m.variables["Link-p"].loc[0, "BuyLink"] * (self.alpha * dr_CO2Intensity[0] + (1 - self.alpha) * dr_price[0])
        for i in range(1, 24):
            expr = expr + m.variables["Link-p"].loc[i, "BuyLink"] * (self.alpha * dr_CO2Intensity[i] + (1 - self.alpha) * dr_price[i])

        # total_energy = m.variables["Link-p"].loc[0, "BuyLink"]
        # for i in range(1, 23):
        #     total_energy = total_energy + m.variables["Link-p"].loc[i, "BuyLink"]
        # expr = total_energy * (self.alpha * dr_CO2Intensity + (1 - self.alpha) * dr_price)
        m.add_objective(expr, overwrite=True, sense="min")

        self.n.optimize.solve_model(solver_name="gurobi")

        self.battery_left = n.stores_t.e["Battery"][23]

        self.dr_pf_series = self.dr_pf_series + [[n.generators_t.p, n.loads_t.p, n.stores_t.p, n.stores_t.e]]

        self.dr_H2_prod = self.dr_H2_prod + [convert_functions.P_to_H2(n.loads_t.p["H2gen"], self.LHV)]

        # print(n.loads_t.p.sum())
        self.total_production = convert_functions.P_to_H2(n.loads_t.p["H2gen"].sum(), self.LHV)

        self.CO2_emissions = np.multiply(n.links_t.p0["BuyLink"], dr_CO2Intensity).sum()
        self.electricity_cost = np.multiply(n.links_t.p0["BuyLink"], dr_price).sum()
        self.electricity_revenue = np.multiply(n.links_t.p0["SellLink"], dr_price).sum()
        self.opportunity_cost = self.opportunity_cost + np.multiply(dr_PF_wind, dr_price).sum()*self.installed_power+np.multiply(dr_PF_solar, dr_price).sum()*self.installed_power
        self.electricity_balance = self.electricity_cost - self.electricity_revenue + self.opportunity_cost


if __name__ == '__main__':
    hydrogen_plant = HydrogenProductionSystem(40, 1000, 1, 1, 1, 1)
    hydrogen_plant.lt_planner()
    hydrogen_plant.daily_planner()
    hydrogen_plant.realisation()
    print(hydrogen_plant.ltp_pf_series)
    print(hydrogen_plant.dp_pf_series)
    print(hydrogen_plant.dr_pf_series)
