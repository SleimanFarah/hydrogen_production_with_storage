import convert_functions
import numpy as np


class BenchmarkOptimizer:

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

