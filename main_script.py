from class_network import *
from convert_functions import *
import random
import matplotlib
from matplotlib import pyplot as plt
from csv_file_function import *
from function_AC import *
import pandas as pd
import os, sys
import time


def run_system_simulation(year, alpha, time_period):

    with open('data/solar/solar_capacity_factor_2017.csv', 'r') as PF_solar_data:
        PF_solar_2017 = csv_reader_function(PF_solar_data)

    with open('data/wind/wind_capacity_factor_2017.csv', 'r') as PF_wind_data:
        PF_wind_2017 = csv_reader_function(PF_wind_data)

    with open('data/electricity_price/electricity_price_2017.csv', 'r') as price_data:
        price_2017 = np.multiply(csv_reader_function(price_data), 0.001).tolist()

    with open('data/co2_intensity/co2_intensity_2017.csv', 'r') as CO2int_data:
        CO2int_2017 = np.multiply(csv_reader_function(CO2int_data), 0.001).tolist()

    with open('data/solar/solar_capacity_factor_2018.csv', 'r') as PF_solar_data:
        PF_solar_2018 = csv_reader_function(PF_solar_data)

    with open('data/wind/wind_capacity_factor_2018.csv', 'r') as PF_wind_data:
        PF_wind_2018 = csv_reader_function(PF_wind_data)

    with open('data/electricity_price/electricity_price_2018.csv', 'r') as price_data:
        price_2018 = np.multiply(csv_reader_function(price_data), 0.001).tolist()

    with open('data/co2_intensity/co2_intensity_2018.csv', 'r') as CO2int_data:
        CO2int_2018 = np.multiply(csv_reader_function(CO2int_data), 0.001).tolist()

    with open('data/solar/solar_capacity_factor_2019.csv', 'r') as PF_solar_data:
        PF_solar_2019 = csv_reader_function(PF_solar_data)

    with open('data/wind/wind_capacity_factor_2019.csv', 'r') as PF_wind_data:
        PF_wind_2019 = csv_reader_function(PF_wind_data)

    with open('data/electricity_price/electricity_price_2019.csv', 'r') as price_data:
        price_2019 = np.multiply(csv_reader_function(price_data), 0.001).tolist()

    with open('data/co2_intensity/co2_intensity_2019.csv', 'r') as CO2int_data:
        CO2int_2019 = np.multiply(csv_reader_function(CO2int_data), 0.001).tolist()

    with open('data/solar/solar_capacity_factor_2020.csv', 'r') as PF_solar_data:
        PF_solar_2020 = csv_reader_function(PF_solar_data)

    with open('data/wind/wind_capacity_factor_2020.csv', 'r') as PF_wind_data:
        PF_wind_2020 = csv_reader_function(PF_wind_data)

    with open('data/electricity_price/electricity_price_2020.csv', 'r') as price_data:
        price_2020 = np.multiply(csv_reader_function(price_data), 0.001).tolist()

    with open('data/co2_intensity/co2_intensity_2020.csv', 'r') as CO2int_data:
        CO2int_2020 = np.multiply(csv_reader_function(CO2int_data), 0.001).tolist()

    with open('data/solar/solar_capacity_factor_2021.csv', 'r') as PF_solar_data:
        PF_solar_2021 = csv_reader_function(PF_solar_data)

    with open('data/wind/wind_capacity_factor_2021.csv', 'r') as PF_wind_data:
        PF_wind_2021 = csv_reader_function(PF_wind_data)

    with open('data/electricity_price/electricity_price_2021.csv', 'r') as price_data:
        price_2021 = np.multiply(csv_reader_function(price_data), 0.001).tolist()

    with open('data/co2_intensity/co2_intensity_2021.csv', 'r') as CO2int_data:
        CO2int_2021 = np.multiply(csv_reader_function(CO2int_data), 0.001).tolist()

    with open('data/solar/solar_capacity_factor_2022.csv', 'r') as PF_solar_data:
        PF_solar_2022 = csv_reader_function(PF_solar_data)

    with open('data/wind/wind_capacity_factor_2022.csv', 'r') as PF_wind_data:
        PF_wind_2022 = csv_reader_function(PF_wind_data)

    with open('data/electricity_price/electricity_price_2022.csv', 'r') as price_data:
        price_2022 = np.multiply(csv_reader_function(price_data), 0.001).tolist()

    with open('data/co2_intensity/co2_intensity_2022.csv', 'r') as CO2int_data:
        CO2int_2022 = np.multiply(csv_reader_function(CO2int_data), 0.001).tolist()

    # Initialization of variables

    battery_on = True
    # alpha = 0.0001

    # for alpha in alphas:
    print(alpha)
    # [0.0001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.999]:

    # year = 2018

    if time_period == "day":
        number_of_days = 1  # number of days in a delivery period
        simulation_period = 365  # number of delivery periods in the simulation
    elif time_period == "week":
        number_of_days = 7
        simulation_period = 52
    elif time_period == "month":
        number_of_days = 30
        simulation_period = 12
    elif time_period == "year":
        number_of_days = 365
        simulation_period = 1
    elif time_period == "test":
        number_of_days = 1
        simulation_period = 1


    delivery_period = 24*number_of_days
    # delivery_mass = number_of_days*0.7*24*P_to_H2(1000, 33.3)
    # total_mass = (108000/8760)*delivery_period*simulation_period
    total_mass = (108000 / 8760) * delivery_period * simulation_period
    # total_mass = 0
    delivery_mass = total_mass/simulation_period
    initial_battery = 0
    initial_hour = 0


    capital_cost_electrolyzer = Annualized_cost(0.07, 10, 700000, 14000)
    # capital_cost_electrolyzer = 0
    opportunity_cost_nb_0 = 199187.789
    opportunity_cost_nb_1 = 198002.809
    # capital_cost_solar = Annualized_cost(0, 40, 310000, 9500)
    # capital_cost_wind = Annualized_cost(0, 30, 990000, 12600)
    capital_cost_battery = Annualized_cost(0, 40, 700000, 540)
    # capital_cost_electrolyzer = Annualized_cost(0, 10, 700000, 14000)
    # capital_cost_converter = Annualized_cost(0, 15, 40000, 0)
    # capital_cost_grid = Annualized_cost(0, 40, 140000, 2800)

    # capital_cost = (capital_cost_solar+capital_cost_wind+capital_cost_battery+capital_cost_electrolyzer+capital_cost_converter+capital_cost_grid)*delivery_period*simulation_period/8760
    # capital_cost = (capital_cost_solar+capital_cost_wind+capital_cost_electrolyzer+capital_cost_converter+capital_cost_grid)*delivery_period*simulation_period/8760
    # capital_cost = 0


    PF_wind = np.array((PF_wind_2017+PF_wind_2018+PF_wind_2019+PF_wind_2020+PF_wind_2021+PF_wind_2022)[initial_hour+8760*(year-2017)-delivery_period:])
    PF_solar = np.array((PF_solar_2017+PF_solar_2018+PF_solar_2019+PF_solar_2020+PF_solar_2021+PF_solar_2022)[initial_hour+8760*(year-2017)-delivery_period:])
    price = np.array((price_2017+price_2018+price_2019+price_2020+price_2021+price_2022)[initial_hour+8760*(year-2017)-delivery_period:])
    CO2int = np.array((CO2int_2017+CO2int_2018+CO2int_2019+CO2int_2020+CO2int_2021+CO2int_2022)[initial_hour+8760*(year-2017)-delivery_period:])

    ltp_pf_series = pd.DataFrame()
    ltp_H2_series = pd.DataFrame()
    ltp_H2_target = pd.DataFrame()
    dp_pf_series = pd.DataFrame()
    dp_Hydro_plan = pd.DataFrame()
    dr_pf_series = pd.DataFrame()
    dr_H2_prod = pd.DataFrame()
    CO2_h2h = pd.DataFrame()
    solar_pf_series = pd.DataFrame()
    wind_pf_series = pd.DataFrame()
    fromgrid_pf_series = pd.DataFrame()
    togrid_pf_series = pd.DataFrame()
    H2_pf_series = pd.DataFrame()


    i = 0
    j = 0
    time_left = delivery_period
    H2_mass_remaining = delivery_mass
    total_production = 0
    total_emissions = 0
    total_electricity_buying_cost = 0
    total_H2_min_cost = 0
    total_opportunity_cost = opportunity_cost_nb_0
    # total_electricity_net_cost = (Annualized_cost(0, 50, 1000000, 10000)/8760)*delivery_period*simulation_period

    hydrogen_plant = HydrogenProductionSystem(battery_on)


    # Loop for the whole delivery period

    while j < simulation_period:
        while i < number_of_days:
            print(j+1, i+1)
            s_idx = 10+i*24*2+j*delivery_period
            e_idx = delivery_period+34+i*24+j*delivery_period

            operative_PF_wind = PF_wind[s_idx:e_idx]
            operative_PF_solar = PF_solar[s_idx:e_idx]
            operative_price = price[s_idx:e_idx]
            operative_CO2int = CO2int[s_idx:e_idx]

            hydrogen_plant.input(time_left, delivery_period, H2_mass_remaining, initial_battery, alpha, operative_PF_wind, operative_PF_solar, operative_price, operative_CO2int)

            hydrogen_plant.lt_planner()
            hydrogen_plant.daily_planner()
            hydrogen_plant.realisation()
            # hydrogen_plant.realisation_opp()

            ltp_pf_series = pd.concat([ltp_pf_series, hydrogen_plant.ltp_pf_series])
            ltp_H2_series = pd.concat([ltp_H2_series, hydrogen_plant.ltp_H2_series])
            ltp_H2_target = pd.concat([ltp_H2_target, hydrogen_plant.ltp_H2_target])
            dp_pf_series = pd.concat([dp_pf_series, hydrogen_plant.dp_pf_series])
            dp_Hydro_plan = pd.concat([dp_Hydro_plan, hydrogen_plant.dp_Hydro_plan])
            dr_pf_series = pd.concat([dr_pf_series, hydrogen_plant.dr_pf_series])
            dr_H2_prod = pd.concat([dr_H2_prod, hydrogen_plant.dr_H2_prod])

            total_emissions = total_emissions + hydrogen_plant.CO2_emissions
            total_electricity_buying_cost = total_electricity_buying_cost + hydrogen_plant.electricity_cost
            net_cost = hydrogen_plant.electricity_net_cost


            total_production = total_production + hydrogen_plant.total_production
            H2_mass_remaining = delivery_mass*(j+1) - total_production

            # hydrogen_plant = HydrogenProductionSystem(time_left, delivery_period, 0, initial_battery, battery_on, alpha, operative_PF_wind, operative_PF_solar, operative_price, operative_CO2int)
            #
            #
            # hydrogen_plant.lt_planner()
            # hydrogen_plant.daily_planner()
            # hydrogen_plant.realisation()
            # opportunity_cost = -hydrogen_plant.electricity_net_cost
            # total_opportunity_cost = total_opportunity_cost + opportunity_cost
            CO2_d2d = CO2_h2h + hydrogen_plant.dr_h2h_CO2

            # total_H2_min_cost = (total_H2_min_cost + net_cost + opportunity_cost)
            total_H2_min_cost = (total_H2_min_cost + net_cost)

            solar_pf_series = pd.concat([solar_pf_series,hydrogen_plant.dr_pf_series[["Solar"]]])
            wind_pf_series = pd.concat([wind_pf_series,hydrogen_plant.dr_pf_series[["Wind"]]])
            fromgrid_pf_series = pd.concat([fromgrid_pf_series,hydrogen_plant.dr_pf_series[["NetworkImport"]]])
            togrid_pf_series = pd.concat([togrid_pf_series,hydrogen_plant.dr_pf_series[["NetworkExport"]]])
            H2_pf_series = pd.concat([H2_pf_series,hydrogen_plant.dr_pf_series[["H2gen"]]])

            i = i + 1
            time_left = delivery_period - i*24
            initial_battery = hydrogen_plant.battery_left
            # pf_series = pf_series + hydrogen_plant.dr_pf_series




        i = 0
        time_left = delivery_period
        H2_mass_remaining = delivery_mass
        j = j+1

    if battery_on:
        capital_cost = capital_cost_electrolyzer
    else:
        capital_cost = capital_cost_electrolyzer

    # print(f"le coût d'opportunité est {total_H2_min_cost}")
    total_H2_min_cost = total_H2_min_cost + total_opportunity_cost + capital_cost


    # print(ltp_pf_series)
    # print(ltp_H2_target)
    # # print(dp_Hydro_plan)
    # print(dp_pf_series)
    print(dr_H2_prod)
    print(dp_pf_series)
    print(dr_pf_series)
    print(total_production)
    print(delivery_mass*j)
    # print(dr_pf_series[0][1].sum())
    # print(total_electricity_buying_cost/total_production)
    # print(total_opportunity_cost/total_production)
    print(total_H2_min_cost/total_production)
    print(total_emissions/total_production)
    # print(total_opportunity_cost)
    # print(total_H2_min_cost)
    # print(total_emissions)

    # figure, axis = plt.subplots(1, 3)
    # axis[0].plot(dr_H2_prod[0])
    # axis[1].plot(dr_pf_series[0][0])
    # axis[2].plot(CO2_d2d)
    # plt.show()

    snapshots = ["Time (h)"] + list(range(delivery_period * simulation_period))

    solar_pf = ["Solar power (kW)"] + list(solar_pf_series["Solar"])

    wind_pf = ["Wind power (kW)"] + list(wind_pf_series["Wind"])

    powerfromgrid_pf = ["Power from grid (kW)"] + list(fromgrid_pf_series["NetworkImport"])

    power2grid_pf = ["Power to grid (kW)"] + list(togrid_pf_series["NetworkExport"])

    power2hydrogen = ["Power to hydrogen (kW)"] + list(H2_pf_series["H2gen"])

    curtailed_solar = ["Curtailed solar power (kW)"] + (
                np.multiply(PF_solar[delivery_period:delivery_period*(simulation_period+1)], 1000) - np.array(solar_pf_series["Solar"])).tolist()

    curtailed_wind = ["Curtailed wind power (kW)"] + (
                np.multiply(PF_wind[delivery_period:delivery_period*(simulation_period+1)], 1000) - np.array(wind_pf_series["Wind"])).tolist()

    # CO2_prod = ["CO2 emitted (kg)"]+list(hydrogen_plant.benchmark_h2h_CO2.T)
    CO2int = ["CO2 intensity"] + CO2int[delivery_period:].tolist()
    price = ["Electricity price"] + price[delivery_period:].tolist()

    array = [snapshots, solar_pf, wind_pf, powerfromgrid_pf, power2grid_pf, power2hydrogen, curtailed_solar, curtailed_wind,
             CO2int, price, ["Minimum price", total_H2_min_cost/total_production], ["H2 CO2 intensity", total_emissions/total_production]]
    df = pd.DataFrame(array).T

    if time_period == "test":
        df.to_excel(excel_writer=f"test.xlsx")
    else:
        if battery_on:
            df.to_excel(excel_writer=f"d2d_{year}_battery_on_period_{number_of_days}d_alpha{alpha}.xlsx")
        else:
            df.to_excel(excel_writer=f"d2d_{year}_battery_off_period{number_of_days}d_alpha{alpha}.xlsx")


if __name__ == "__main__":

    start_time = time.time()
    year = 2019
    delivery_period = "week"
    # alphas = [0.0001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.999]
    alphas = [0.5]

    for alpha in alphas:
        run_system_simulation(year, alpha, delivery_period)

    end_time = time.time()
    print(f"Execution time {end_time - start_time} seconds")


