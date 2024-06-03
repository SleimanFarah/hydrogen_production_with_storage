import convert_functions
import numpy as np


class HydrogenProductionSystem:

    def __init__(self, time_left, delivery_period, H2MassRemaining, battery_left,
                operative_PF_wind, operative_PF_solar, operative_price, operative_CO2Intensity):

        # The delivery period is in hours
        # The times series are supposed to be dimensioned correctly

        import pypsa
        import convert_functions
        self.time_left = time_left
        self.delivery_period = delivery_period
        self.H2MassRemaining = H2MassRemaining
        self.battery_left = battery_left

        self.operative_PF_wind = operative_PF_wind
        self.operative_PF_solar = operative_PF_solar
        self.operative_price = operative_price
        self.operative_CO2Intensity = operative_CO2Intensity

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
        self.hydrolyzer_ramp_limit = 0.5
        n.add("Store", "H2gen", bus="H2Bus")
        n.add("Load", "H2gen", bus="H2Bus")


        # Battery creation
        battery_capacity = 0
        n.add("Store", "Battery", bus="BatteryBus", e_nom_extendable=True, e_nom_max=battery_capacity)

        # Link values initialization
        self.eff_charge = 0.9
        self.eff_discharge = 0.9
        self.eff_electrolysis = 0.6
        self.eff_converter = 0.95

        # Links creation
        n.add("Link", "H2Link", bus0="ACBus", bus1="H2Bus", efficiency=self.eff_electrolysis, p_nom_extendable=True, p_nom_max=self.hydrolyzer_capacity)
        n.add("Link", "ChargeLink", bus0="DCBus", bus1="BatteryBus", efficiency=self.eff_charge, p_nom_extendable=True, p_nom_max=1000)
        n.add("Link", "DischargeLink", bus0="BatteryBus", bus1="DCBus", efficiency=self.eff_discharge, p_nom_extendable=True, p_nom_max=1000)
        n.add("Link", "ACtoDCLink", bus0="ACBus", bus1="DCBus", efficiency=self.eff_converter, p_nom_extendable=True)
        n.add("Link", "DCtoACLink", bus0="DCBus", bus1="ACBus", efficiency=self.eff_converter, p_nom_extendable=True)
        n.add("Link", "BuyLink", bus0="NetworkBus", bus1="ACBus", p_nom_extendable=True, p_nom_max=1000)
        n.add("Link", "SellLink", bus0="ACBus", bus1="NetworkBus", p_nom_extendable=True, p_nom_max=1000)

        # Optimization parameters initialization

        self.alpha = 0
        self.CO2_price = 0.006  # price of emitting 1 kg of C02
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
        if self.time_left == 24:
            self.ltp_target_mass = self.H2MassRemaining
        else:
            n = self.n
            time_left = self.time_left
            delivery_period = self.delivery_period
            n.set_snapshots(range(time_left))
            n.stores.e_nom_extendable["H2gen"] = True
            n.loads.p_set["H2gen"] = 0
            n.stores.e_initial["Battery"] = self.battery_left

            ltp_PF_wind = self.operative_PF_wind[0:time_left]
            ltp_PF_solar = self.operative_PF_solar[0:time_left]
            ltp_price = self.operative_price[0:time_left]
            ltp_CO2Intensity = self.operative_CO2Intensity[0:time_left]

            n.generators_t.p_max_pu["Solar"] = ltp_PF_solar
            n.generators_t.p_max_pu["Wind"] = ltp_PF_wind
            n.generators_t.p_min_pu["Solar"] = ltp_PF_solar
            n.generators_t.p_min_pu["Wind"] = ltp_PF_wind

            m = n.optimize.create_model()

                # total_hydrolyzer_energy = m.variables["Store-e"].loc[time_left+23, "H2gen"]
                # total_H2_required_energy = convert_functions.H2_to_P(self.H2MassRemaining, self.LHV)
                # cstr = total_hydrolyzer_energy == total_H2_required_energy
                # m.add_constraints(cstr, name="target match")

            total_hydrolyzer_energy = -m.variables["Store-p"].loc[0, "H2gen"]
            for i in range(1, time_left):
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

            for i in range(time_left-1):
                cstr = -self.hydrolyzer_capacity*self.hydrolyzer_ramp_limit*self.eff_electrolysis <= m.variables["Store-p"].loc[i, "H2gen"] - m.variables["Store-p"].loc[i+1, "H2gen"] <= self.hydrolyzer_capacity*self.hydrolyzer_ramp_limit*self.eff_electrolysis
                m.add_constraints(cstr, name="ramp limit at time {}".format(i))

            initial_energy = m.variables["Store-e"].loc[0, "Battery"]
            final_energy = m.variables["Store-e"].loc[time_left-1, "Battery"]
            cstr = initial_energy == final_energy
            m.add_constraints(cstr, name="battery sustainability")

            expr = m.variables["Link-p"].loc[0, "BuyLink"] * self.alpha * ltp_CO2Intensity[0] + (1 - self.alpha) * ltp_price[0] * (m.variables["Link-p"].loc[0, "BuyLink"]-m.variables["Link-p"].loc[0, "SellLink"])
            for i in range(1, time_left):
                expr = expr + m.variables["Link-p"].loc[i, "BuyLink"] * self.alpha * ltp_CO2Intensity[i] + (1 - self.alpha) * ltp_price[i] * (m.variables["Link-p"].loc[i, "BuyLink"]-m.variables["Link-p"].loc[i, "SellLink"])

            m.add_objective(expr, overwrite=True, sense="min")

            n.optimize.solve_model(solver_name="gurobi")

            self.ltp_target_mass = -convert_functions.P_to_H2((n.stores_t.p["H2gen"][time_left-34:time_left-10].sum()), self.LHV)
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
        time_left = self.time_left

        dp_PF_wind = self.operative_PF_wind[time_left-10:time_left+24]
        dp_PF_solar = self.operative_PF_solar[time_left-10:time_left+24]
        dp_price = self.operative_price[time_left-10:time_left+24]
        dp_CO2Intensity = self.operative_CO2Intensity[time_left-10:time_left+24]

        n.generators_t.p_max_pu["Solar"] = dp_PF_solar
        n.generators_t.p_max_pu["Wind"] = dp_PF_wind
        n.generators_t.p_min_pu["Solar"] = dp_PF_solar
        n.generators_t.p_min_pu["Wind"] = dp_PF_wind


        m = n.optimize.create_model()

        total_hydrolyzer_energy = m.variables["Store-p"].loc[0, "H2gen"]
        for i in range(1, 24):
            total_hydrolyzer_energy = total_hydrolyzer_energy + m.variables["Store-p"].loc[i, "H2gen"]
        total_H2_required_energy = -convert_functions.H2_to_P(self.ltp_target_mass, self.LHV)
        # total_H2_required_energy = convert_functions.H2_to_P(-10, self.LHV)
        cstr = total_hydrolyzer_energy == total_H2_required_energy
        m.add_constraints(cstr, name="target match")

        initial_energy = m.variables["Store-e"].loc[0, "Battery"]
        final_energy = m.variables["Store-e"].loc[33, "Battery"]
        cstr = initial_energy == final_energy
        m.add_constraints(cstr, name="battery sustainability")

        for i in range(33):
            cstr = -self.hydrolyzer_capacity*self.hydrolyzer_ramp_limit*self.eff_electrolysis <= m.variables["Store-p"].loc[i, "H2gen"] - m.variables["Store-p"].loc[i+1, "H2gen"] <= self.hydrolyzer_capacity*self.hydrolyzer_ramp_limit*self.eff_electrolysis
            m.add_constraints(cstr, name="ramp limit at time {}".format(i))


        # total_hydrolyzer_energy = m.variables["Link-p"].loc[0, "H2Link"]*self.eff_electrolysis
        # for i in range(1, 23):
        #     total_hydrolyzer_energy = total_hydrolyzer_energy + m.variables["Link-p"].loc[i, "H2Link"]*self.eff_electrolysis
        # total_H2_required_energy = convert_functions.H2_to_P(self.ltp_target_mass, self.LHV)
        # cstr = total_hydrolyzer_energy == total_H2_required_energy
        # m.add_constraints(cstr, name="target match")

        # for i in range(34):
        #     cstr = m.variables["Link-p"].loc[i, "H2Link"] <= self.hydrolyzer_capacity
        #     m.add_constraints(cstr, name="limited capacity at time {}".format(i))

        expr = m.variables["Link-p"].loc[0, "BuyLink"] * self.alpha * dp_CO2Intensity[0] + (1 - self.alpha) * dp_price[0] * (m.variables["Link-p"].loc[0, "BuyLink"]-m.variables["Link-p"].loc[0, "SellLink"])
        for i in range(1, 34):
            expr = expr + m.variables["Link-p"].loc[i, "BuyLink"] * self.alpha * dp_CO2Intensity[i] + (1 - self.alpha) * dp_price[i] * (m.variables["Link-p"].loc[i, "BuyLink"]-m.variables["Link-p"].loc[i, "SellLink"])

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

        dr_PF_wind = self.operative_PF_wind[time_left-10:time_left+14]
        dr_PF_solar = self.operative_PF_solar[time_left-10:time_left+14]
        dr_price = self.operative_price[time_left-10:time_left+14]
        dr_CO2Intensity = self.operative_CO2Intensity[time_left-10:time_left+14]

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

        # for i in range(0, 24):
        #     cstr = m.variables["Link-p"].loc[i, "H2Link"] <= self.hydrolyzer_capacity
        #     m.add_constraints(cstr, name="limited capacity at time {}".format(i))

        expr = m.variables["Link-p"].loc[0, "BuyLink"] * self.alpha * dr_CO2Intensity[0] + (1 - self.alpha) * dr_price[0] * (m.variables["Link-p"].loc[0, "BuyLink"]-m.variables["Link-p"].loc[0, "SellLink"])
        for i in range(1, 24):
            expr = expr + m.variables["Link-p"].loc[i, "BuyLink"] * self.alpha * dr_CO2Intensity[i] + (1 - self.alpha) * dr_price[i] * (m.variables["Link-p"].loc[i, "BuyLink"]-m.variables["Link-p"].loc[i, "SellLink"])

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

    def benchmark(self, total_time, time_period, delivery_mass):
        n = self.n
        n.set_snapshots(range(total_time))
        n.stores.e_nom_extendable["H2gen"] = False
        n.stores.e_initial["Battery"] = self.battery_left

        n.generators_t.p_max_pu["Solar"] = self.operative_PF_solar
        n.generators_t.p_max_pu["Wind"] = self.operative_PF_wind
        n.generators_t.p_min_pu["Solar"] = self.operative_PF_solar
        n.generators_t.p_min_pu["Wind"] = self.operative_PF_wind

        m = n.optimize.create_model()

        initial_energy = m.variables["Store-e"].loc[0, "Battery"]
        final_energy = m.variables["Store-e"].loc[23, "Battery"]
        cstr = initial_energy == final_energy
        m.add_constraints(cstr, name="battery sustainability")

        for i in range(total_time):
            cstr = -self.hydrolyzer_capacity*self.hydrolyzer_ramp_limit*self.eff_electrolysis <= m.variables["Store-p"].loc[i, "H2gen"] - m.variables["Store-p"].loc[i+1, "H2gen"] <= self.hydrolyzer_capacity*self.hydrolyzer_ramp_limit*self.eff_electrolysis
            m.add_constraints(cstr, name="ramp limit at time {}".format(i))

        for j in range(total_time // time_period):
            total_hydrolyzer_energy = m.variables["Store-p"].loc[0, "H2gen"]
            for i in range(1, time_period):
                total_hydrolyzer_energy = total_hydrolyzer_energy + m.variables["Store-p"].loc[i, "H2gen"]
            total_H2_required_energy = -convert_functions.H2_to_P(delivery_mass, self.LHV)
            cstr = total_hydrolyzer_energy == total_H2_required_energy
            m.add_constraints(cstr, name="target match {}".format(j))

        expr = m.variables["Link-p"].loc[0, "BuyLink"] * self.alpha * ltp_CO2Intensity[0] + (1 - self.alpha) * ltp_price[0] * (m.variables["Link-p"].loc[0, "BuyLink"]-m.variables["Link-p"].loc[0, "SellLink"])
        for i in range(1, total_time):
            expr = expr + m.variables["Link-p"].loc[i, "BuyLink"] * self.alpha * ltp_CO2Intensity[i] + (1 - self.alpha) * ltp_price[i] * (m.variables["Link-p"].loc[i, "BuyLink"]-m.variables["Link-p"].loc[i, "SellLink"])


if __name__ == '__main__':
    import random
    PF_wind = [random.uniform(0, 1) for _ in range(48)]
    PF_solar = [random.uniform(0, 1) for _ in range(48)]
    price = [random.uniform(1, 10) for _ in range(48)]
    CO2int = [random.uniform(1, 10) for _ in range(48)]
    hydrogen_plant = HydrogenProductionSystem(24, 24, 200, 0, PF_wind, PF_solar, price, CO2int)
    # hydrogen_plant.lt_planner()
    hydrogen_plant.daily_planner()
    hydrogen_plant.realisation()
    # print(hydrogen_plant.ltp_pf_series)
    # print(hydrogen_plant.ltp_H2_target)
    # print(dp_Hydro_plan)
    print(hydrogen_plant.dp_pf_series)
    print(hydrogen_plant.dr_pf_series)
    print(hydrogen_plant.total_production)
    # print(dr_pf_series[0][1].sum())
    print(hydrogen_plant.electricity_balance/hydrogen_plant.total_production)
    print(hydrogen_plant.CO2_emissions/hydrogen_plant.total_production)
