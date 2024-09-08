import os
import pandas as pd

import convert_functions
import numpy as np
import pypsatopo
from function_positive_only import *
from function_positive_only_float import *
import pypsa
import xarray as xr
import matplotlib.pyplot as plt

outputs_folder = "outputs"

def create_array_value_every_nth_element(length, n, x):
    # Initialize the array with zeros
    array = np.zeros(length)

    # Set every nth element to x
    for ii in range(n - 1, length, n):
        array[ii] = (ii+1)/n*x

    return array


def annualised_cost(discount_rate, life_expectancy, cost, maintenance):
    if discount_rate != 0.0:
        crf = discount_rate / (1.0 - (1.0 + discount_rate) ** (-life_expectancy))
    else:
        crf = 1

    annualised_cost = crf * cost + maintenance
    return annualised_cost

# For parallel simulation on PRIME
def run_system_simulation(year, alpha, time_period):
    hps = HydrogenProductionSystem(year=year, delivery_period=time_period, h2_kg_annual_target=108000, alpha_co2=alpha,
                                   wind_capacity=1.0, solar_capacity=1.0, electrolyser_capacity=1.0,
                                   battery_capacity=1.0, battery_initial_soc=0.5)
    hps.create_network()
    hps.simulate_benchmark()
    # hps.simulate_day_to_day()


class HydrogenProductionSystem:

    def __init__(self, year, delivery_period, alpha_co2, h2_kg_annual_target, wind_capacity, solar_capacity, electrolyser_capacity, battery_capacity, battery_initial_soc):

        _years = range(2017, 2023)

        self.year = year
        self.delivery_period = delivery_period
        self.alpha_co2 = alpha_co2
        self.wind_p_nom = wind_capacity

        self.electrolyser_p_nom = electrolyser_capacity
        self.battery_e_nom = battery_capacity
        self.battery_initial_soc = battery_initial_soc
        self.grid_initial_soc = 0.5
        self.electrolyser_initial_soc = 0.0

        self.grid_e_nom = 10e6

        self.wind_cf_t = self.read_data_from_csv(file_name_pattern="data/wind/wind_capacity_factor", years=_years)
        self.solar_cf_t = self.read_data_from_csv(file_name_pattern="data/solar/solar_capacity_factor", years=_years)
        self.el_price_t = self.read_data_from_csv(file_name_pattern="data/electricity_price/electricity_price",
                                                  years=_years)
        self.co2_int_t = self.read_data_from_csv(file_name_pattern="data/co2_intensity/co2_intensity", years=_years)

        # Link values initialisation
        self.eff_battery_charge = 0.985
        self.eff_battery_discharge = 0.975
        self.eff_electrolysis = 0.6
        self.eff_inverter = 0.9
        self.eff_rectifier = 0.9

        # Optimization parameters initialization
        self.co2_price = 0.075  # price of emitting 1 kg of C02

        # Solar
        self.solar_p_nom = solar_capacity
        self.solar_fom_euro_per_mw = 9500
        self.solar_vom_euro_per_mwh = 0.0

        self.discount_rate = 0.07

        # Electrolyser cost
        self.electrolyser_eff = 0.6
        self.electrolyser_life_years = 10.0
        self.electrolyser_investment_cost = 700000.0
        self.electrolyser_fom_euro_per_mw = 14000.0
        self.electrolyser_vom_euro_per_mwh = 0.0
        self.electrolyser_annualised_cost = annualised_cost(self.discount_rate, self.electrolyser_life_years,
                                                                 self.electrolyser_investment_cost,
                                                                 self.electrolyser_fom_euro_per_mw)
        # Hydrogen
        # ********
        self.h2_lhv = 120  # MJ/kg
        self.h2_kg_annual_target = h2_kg_annual_target

        self.get_electrolyser_full_load_hours_annually()
        self.get_n_devliery_and_n_day_per_delivery()
        self.get_electrolyser_full_load_hours_per_delivery()


        # Battery
        self.battery_life_years = 10.0
        self.battery_investment_cost = 700000.0
        self.battery_fom_euro_per_mw = 14000.0
        self.battery_vom_euro_per_mwh = 1.8
        self.battery_annualised_cost = annualised_cost(self.discount_rate, self.electrolyser_life_years,
                                                            self.electrolyser_investment_cost,
                                                            self.electrolyser_fom_euro_per_mw)

        # operation cost of components

        self.wind_vom_euro_per_mwh = 1.35

        self.benchmark_pf_series = pd.DataFrame()
        self.benchmark_H2_prod = pd.DataFrame()
        self.benchmark_h2h_CO2 = pd.DataFrame()

        self.remaining_deliveries = self.n_delivery
        self.remaining_days_in_delivery = self.n_days_per_delivery
        self.remaining_electrolyser_full_load_hours_per_delivery = self.electrolyser_full_load_hours_per_delivery
        self.produced_electrolyser_full_load_hours_per_day = 0.0
        self.electrolyser_e_nom = self.electrolyser_full_load_hours_per_delivery*self.n_delivery

        self.all_snapshots = self.wind_cf_t.index
        self.now = self.wind_cf_t.index[self.all_snapshots.year == self.year][0] - pd.Timedelta(hours=14)
        self.now_idx = np.argwhere(self.all_snapshots == self.now)[0][0]

        self.network_time_series_output_file_name = f"network_time_series_history_df_{self.year}_{self.delivery_period}_h2_{self.h2_kg_annual_target}_alpha{self.alpha_co2}_wc{self.wind_p_nom}_sc{self.solar_p_nom}_ec{self.electrolyser_p_nom}_bc{self.battery_e_nom}_bisoc{self.battery_initial_soc}"

        self.network_time_series_history_df = pd.DataFrame()

    def update_snapshots(self):

        if self.remaining_days_in_delivery > 1:
            past_hours = (self.remaining_days_in_delivery - 2)*24
            future_hours = 2*24
        else:
            past_hours = -14
            future_hours = 34 - past_hours

        self.snapshots_ltp = self.all_snapshots[self.now_idx-past_hours:self.now_idx+future_hours]
        self.snapshots_dp = self.snapshots_ltp[-34:]
        self.snapshots_hpp = self.snapshots_dp[0:24]


    def get_electrolyser_full_load_hours_annually(self):

        self.electrolyser_full_load_hours_annually = self.h2_lhv * self.h2_kg_annual_target/(3600.0*self.electrolyser_eff)

        return self.electrolyser_full_load_hours_annually

    def get_electrolyser_full_load_hours_per_delivery(self):

        self.electrolyser_full_load_hours_per_delivery = self.electrolyser_full_load_hours_annually*self.n_days_per_delivery/365.0

        return self.electrolyser_full_load_hours_per_delivery

    def get_n_devliery_and_n_day_per_delivery(self):

        if self.delivery_period == "day":
            self.n_days_per_delivery = 1  # number of days in a delivery period
            self.n_delivery = 365  # number of delivery periods in the simulation
        elif self.delivery_period == "week":
            self.n_days_per_delivery = 7
            self.n_delivery = 52
        elif self.delivery_period == "month":
            self.n_days_per_delivery = 30
            self.n_delivery = 12
        elif self.delivery_period == "year":
            self.n_days_per_delivery = 365
            self.n_delivery = 1
        elif self.delivery_period == "test":
            self.n_days_per_delivery = 1
            self.n_delivery = 1
        self.n_total_hours = 24*self.n_days_per_delivery*self.n_delivery

    def create_network(self):
        # Creation of the network
        self.network = pypsa.Network()
        self.network.set_snapshots(self.wind_cf_t.index)

        self.network.add("Carrier", "ac")
        self.network.add("Carrier", "dc")

        # Buses creation
        self.network.add("Bus", "ac_bus", carrier="ac")
        self.network.add("Bus", "dc_bus", carrier="dc")
        self.network.add("Bus", "grid_bus", carrier="ac")
        self.network.add("Bus", "h2_bus", carrier="ac")

        # Add wind
        # ********
        self.network.add(
            "Generator",
            "wind_turbine",
            bus="ac_bus",
            carrier="ac",
            p_nom=self.wind_p_nom,
            # p_nom_extendable=self.wind_p_nom_extendable,
            p_max_pu=self.wind_cf_t.values.ravel(),
            # p_nom_max=self.wind_p_nom_max,
            # p_nom_min=self.wind_p_nom_min,
            # capital_cost=self.wind_capital_cost_mw,
            marginal_cost=self.wind_vom_euro_per_mwh*(1.0 - self.alpha_co2)
        )
        # Add solar
        # *********
        self.network.add(
            "Generator",
            "solar_pv",
            bus="dc_bus",
            carrier="dc",
            p_nom=self.solar_p_nom,
            # p_nom_extendable=self.solar_p_nom_extendable,
            p_max_pu=self.solar_cf_t.values.ravel(),
            # p_nom_max=self.solar_p_nom_max,
            # p_nom_min=self.solar_p_nom_min,
            # capital_cost=self.solar_capital_cost_mw,
            marginal_cost=self.solar_vom_euro_per_mwh*(1.0 - self.alpha_co2)
        )
        # Add grid
        # ********
        self.network.add(
            "Store",
            "grid",
            carrier="ac",
            bus="grid_bus",
            e_cyclic=False,
            e_nom=self.grid_e_nom,
            e_initial=self.grid_initial_soc*self.grid_e_nom,
            capital_cost=0.0,
            marginal_cost=0.0
        )

        # Add electrolyser
        # ****************
        self.network.add(
            "Store",
            "electrolyser",
            carrier="ac",
            bus="h2_bus",
            e_nom=self.electrolyser_e_nom,
            e_initial = self.electrolyser_initial_soc*self.electrolyser_e_nom,
            e_cyclic=False,
            # e_min_pu=np.ones(len(self.all_snapshots))*self.remaining_electrolyser_full_load_hours_per_delivery,
            capital_cost=self.electrolyser_annualised_cost,
            marginal_cost=self.electrolyser_vom_euro_per_mwh
        )

        # Battery creation
        if self.battery_e_nom > 0:
            self.network.add(
                "Store",
                "battery",
                carrier="dc",
                bus="dc_bus",
                e_nom=self.battery_e_nom,
                e_initial=self.battery_initial_soc*self.battery_e_nom,
                e_cyclic=True,  # This flag is necessary here for the ltp.
                capital_cost=self.battery_annualised_cost,
                marginal_cost=self.battery_vom_euro_per_mwh
            )

        self.network.add(
            "Link",
            "export_to_grid_link",
            bus0="ac_bus",
            bus1="grid_bus",
            carrier="ac",
            p_nom=999,
            capital_cost=0.0,
            marginal_cost=-self.el_price_t.values.ravel()*(1.0 - self.alpha_co2)
        )

        self.network.add(
            "Link",
            "import_from_grid_link",
            bus0="grid_bus",
            bus1="ac_bus",
            carrier="ac",
            p_nom=999,
            marginal_cost=self.el_price_t.values.ravel() * (1.0 - self.alpha_co2) + self.co2_int_t.values.ravel() * self.co2_price * self.alpha_co2
        )
        self.network.add(
            "Link",
            "inverter",
            bus0="dc_bus",
            bus1="ac_bus",
            carrier="dc",
            p_nom=999,
            efficiency=self.eff_inverter,
            marginal_cost=0.0
        )
        self.network.add(
            "Link",
            "rectifier",
            bus0="ac_bus",
            bus1="dc_bus",
            carrier="ac",
            p_nom=999,
            efficiency=self.eff_rectifier,
            marginal_cost=0.0
        )
        self.network.add(
            "Link",
            "electrolyser",
            bus0="ac_bus",
            bus1="h2_bus",
            carrier="ac",
            p_nom=self.electrolyser_p_nom,
            efficiency=1.0,
            marginal_cost=0.0
        )

        self.network_original = self.network.copy()

    def read_data_from_csv(self, file_name_pattern, years):
        # Initialize an empty DataFrame to store all data
        df_all = pd.DataFrame()

        # Loop over each year in the provided list
        for year in years:
            # Construct the file name by appending the year to the file name pattern
            file_name = file_name_pattern + "_" + str(year) + ".csv"

            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_name, index_col=0, parse_dates=True)

            # If the DataFrame does not have 8760 rows (which is the expected number for hourly data in a year),
            # print a warning message in red color
            if len(df) != 8760:
                print(f"\033[91m Warning: check data length in {file_name} \033[0m")

            # Concatenate the DataFrame from the current file with the overall DataFrame
            df_all = pd.concat([df_all, df])
            df_all.index = pd.to_datetime(df_all.index)

        # Return the overall DataFrame containing data from all years
        return df_all

    def plot_network(self):
        pypsatopo.generate(self.network, file_output="mynetworkfinal.cvs")
        self.network.plot(line_colors="r", title="network")

    def run_ltp(self):

        # We want cycling of the battery for the long-term planner. This is set in the creation of the network

        if self.remaining_days_in_delivery > 1:  # This is not the last day
            # Set the snapshots
            self.network.set_snapshots(self.snapshots_ltp)
            # Constraint for hydrogen production
            # **********************************
            self.network.stores_t["e_min_pu"]["electrolyser"] = np.zeros(len(self.snapshots_ltp))
            self.network.stores_t["e_min_pu"].loc[self.snapshots_ltp[-1], "electrolyser"] = self.electrolyser_initial_soc + self.remaining_electrolyser_full_load_hours_per_delivery/self.electrolyser_e_nom
            # **********************************
            # Optimise the network
            self.network.optimize(solver_name="gurobi")

            # Get the total hydrogen for daily planner
            self.electrolyser_full_load_hours_to_dp = self.network.links_t["p0"].loc[self.snapshots_dp, "electrolyser"].sum()
        else:  # This is the last day
            self.electrolyser_full_load_hours_to_dp = self.remaining_electrolyser_full_load_hours_per_delivery
            # Set the snapshots
            self.snapshots_dp = self.snapshots_dp[0:24]

        print("LTP run successfully!")
        print(self.electrolyser_full_load_hours_to_dp)

        return self.electrolyser_full_load_hours_to_dp

    def run_dp(self):
        # Set the snapshots
        self.network.set_snapshots(self.snapshots_dp)

        # If cyclic is used, the dp does not consider the actual initial state of the battery.
        self.network.stores.loc["battery", "e_cyclic"] = False

        # Enforcing e_cycling while taking into account the initial condition of the battery.
        # if self.remaining_days_in_delivery == 1:
        #     self.network.stores_t["e_min_pu"].loc[self.network.stores_t["e_min_pu"].index[-1], "battery"] = self.battery_initial_soc

        if self.remaining_days_in_delivery == 1:  # Enforcing cycling on the last day of the delivery period
            self.network.stores_t["e_min_pu"]["battery"] = np.zeros(len(self.snapshots_dp))
            self.network.stores_t["e_min_pu"].loc[self.snapshots_dp[-1], "battery"] = self.battery_initial_soc

        # Hydrogen production constraint
        self.network.stores_t["e_min_pu"]["electrolyser"] = np.zeros(len(self.snapshots_dp))
        self.network.stores_t["e_min_pu"].loc[self.snapshots_dp[-1], "electrolyser"] = self.electrolyser_initial_soc + self.electrolyser_full_load_hours_to_dp/self.electrolyser_e_nom

        self.network.optimize(solver_name="gurobi")

        self.electrolyser_full_load_hours_to_hpp = self.network.links_t.p0["electrolyser"]

        print("DP run successfully!")
        print(self.electrolyser_full_load_hours_to_hpp)

        return self.electrolyser_full_load_hours_to_hpp

    def run_hpp(self):

        # Set the snapshots
        self.network.set_snapshots(self.snapshots_hpp)
        # self.network.stores.loc["battery", "e_cyclic"] = False

        # Hydrogen production constraint
        self.network.links_t.p_min_pu["electrolyser"] = self.network.links_t.p0["electrolyser"]
        # Battery power flow constraint
        self.network.stores_t.p_set["battery"] = self.network.stores_t.p["battery"]

        self.network.optimize(solver_name="gurobi")

        self.produced_electrolyser_full_load_hours_per_day = self.network.links_t["p0"].loc[
            self.snapshots_hpp, "electrolyser"].sum()

        print(self.produced_electrolyser_full_load_hours_per_day)
        print("HPP run successfully!")
        return self.produced_electrolyser_full_load_hours_per_day

    def simulate_day_to_day(self):

        os.makedirs(f"{outputs_folder}/day_to_day", exist_ok=True)

        for delivery in range(self.n_delivery):
            print(f"Delivery {(delivery+1)}/{self.n_delivery}")
            self.remaining_days_in_delivery = self.n_days_per_delivery
            self.remaining_electrolyser_full_load_hours_per_delivery = self.electrolyser_full_load_hours_per_delivery
            for day in range(self.n_days_per_delivery):
                print(f"Delivery {(delivery + 1)}/{self.n_delivery}")
                print(f"Day {day+1}/{self.n_days_per_delivery}")

                self.update_snapshots()
                self.run_ltp()

                self.run_dp()

                self.run_hpp()

                time_series_df = pd.concat(
                    [self.network.generators_t.p, self.network.stores_t.p, self.network.stores_t.e,
                     self.network.links_t.p0, self.network.links_t.p1], axis=1)

                self.network_time_series_history_df = pd.concat([self.network_time_series_history_df, time_series_df])

                self.remaining_days_in_delivery -= 1
                print(self.remaining_days_in_delivery)
                self.remaining_electrolyser_full_load_hours_per_delivery -= self.produced_electrolyser_full_load_hours_per_day
                print(self.remaining_electrolyser_full_load_hours_per_delivery)

                # ******************************************************************************************************
                # Update section
                # ******************************************************************************************************
                # Updating self.now_idx is critical for the calculation. self.now_idx is used instead of self.now to
                # avoid the leap day issues
                self.now_idx += 24
                # Updating self.now to avoid confusion. This update is not critical for the calculation.
                self.now = self.all_snapshots[self.now_idx]

                # Get battery soc
                self.battery_initial_soc = self.network.stores_t.e["battery"].iloc[-1]/self.battery_e_nom
                self.electrolyser_initial_soc = self.network.stores_t.e["electrolyser"].iloc[-1]/self.electrolyser_e_nom
                self.grid_initial_soc = self.network.stores_t.e["grid"].iloc[-1]/self.grid_e_nom

                # Get original network
                self.network = self.network_original.copy()

                # Update battery state of charge
                self.network.stores.loc["battery", "e_initial"] = self.battery_initial_soc*self.battery_e_nom
                self.network.stores.loc["electrolyser", "e_initial"] = self.electrolyser_initial_soc*self.electrolyser_e_nom
                self.network.stores.loc["grid", "e_initial"] = self.grid_initial_soc*self.grid_e_nom
                # ******************************************************************************************************
                # ******************************************************************************************************
        # Save time series to csv
        # ***********************
        suffixes = 5 * ["_MW"] + 3 * ["_MWh"] + 5 * ["_p0_MW"] + 5 * ["_p1_MW"]
        self.network_time_series_history_df.columns = [f"{col}{suffix}" for col, suffix in zip(self.network_time_series_history_df.columns, suffixes)]
        self.network_time_series_history_df.to_csv(f"{outputs_folder}/'day_to_day'/{self.network_time_series_output_file_name}.csv")
        print("Simulation run successfully!")

    def simulate_benchmark(self):

        os.makedirs(f"{outputs_folder}/benchmark", exist_ok=True)

        self.snapshots_benchmark = self.all_snapshots[self.all_snapshots.year == self.year]
        self.snapshots_benchmark = self.snapshots_benchmark[0:self.n_total_hours]

        # Set the snapshots
        self.network.set_snapshots(self.snapshots_benchmark)


        # Hydrogen production constraint
        e_min_pu_electrolyser = create_array_value_every_nth_element(length=len(self.snapshots_benchmark), x=self.electrolyser_full_load_hours_per_delivery, n=24*self.n_days_per_delivery)
        self.network.stores.loc["electrolyser", "e_nom"] = e_min_pu_electrolyser[-1]
        self.network.stores_t["e_min_pu"]["electrolyser"] = e_min_pu_electrolyser/e_min_pu_electrolyser[-1]
        self.network.optimize(solver_name="gurobi")

        self.network_time_series_history_df = pd.concat(
            [self.network.generators_t.p, self.network.stores_t.p, self.network.stores_t.e, self.network.links_t.p0,
             self.network.links_t.p1], axis=1)

        # Save time series to csv
        # ***********************
        suffixes = 5 * ["_MW"] + 3 * ["_MWh"] + 5 * ["_p0_MW"] + 5 * ["_p1_MW"]
        self.network_time_series_history_df.columns = [f"{col}{suffix}" for col, suffix in
                                                       zip(self.network_time_series_history_df.columns, suffixes)]
        self.network_time_series_history_df.to_csv(
            f"{outputs_folder}/benchmark/{self.network_time_series_output_file_name}.csv")
        print("Simulation run successfully!")


if __name__ == '__main__':
    hps = HydrogenProductionSystem(year=2018, delivery_period="week", h2_kg_annual_target=108000, alpha_co2=0.5,
                                   wind_capacity=1.0, solar_capacity=1.0, electrolyser_capacity=1.0,
                                   battery_capacity=1.0, battery_initial_soc=0.5)
    hps.create_network()

    # hps.simulate_benchmark()

    hps.simulate_day_to_day()
