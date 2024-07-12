import numpy as np

from class_network import *
from convert_functions import *
import random
from csv_file_function import *
from function_AC import *
from matplotlib import pyplot as plt
import pandas as pd

with open('solar_capacity_factor 2018.csv', 'r') as PF_solar_data:
    PF_solar_2018 = csv_reader_function(PF_solar_data)

with open('wind_capacity_factor 2018.csv', 'r') as PF_wind_data:
    PF_wind_2018 = csv_reader_function(PF_wind_data)

with open('electricity_price 2018.csv', 'r') as price_data:
    price_2018 = np.multiply(csv_reader_function(price_data), 0.001).tolist()

with open('co2_intensity 2018.csv', 'r') as CO2int_data:
    CO2int_2018 = np.multiply(csv_reader_function(CO2int_data), 0.001).tolist()

with open('solar_capacity_factor 2019.csv', 'r') as PF_solar_data:
    PF_solar_2019 = csv_reader_function(PF_solar_data)

with open('wind_capacity_factor 2019.csv', 'r') as PF_wind_data:
    PF_wind_2019 = csv_reader_function(PF_wind_data)

with open('electricity_price 2019.csv', 'r') as price_data:
    price_2019 = np.multiply(csv_reader_function(price_data), 0.001).tolist()

with open('co2_intensity 2019.csv', 'r') as CO2int_data:
    CO2int_2019 = np.multiply(csv_reader_function(CO2int_data), 0.001).tolist()

# capital_cost_solar = Annualized_cost(0.07, 40, 310000, 9500)
# capital_cost_wind = Annualized_cost(0.07, 30, 990000, 12600)
capital_cost_battery = Annualized_cost(0.07, 25,142000, 540)
capital_cost_electrolyzer = Annualized_cost(0.07, 10, 700000, 14000)

# capital_cost_converter = Annualized_cost(0.07, 15, 20000, 0)
# capital_cost_grid = Annualized_cost(0.07, 40, 140000, 2800)

opportunity_cost_nb_0 = 199187.789
# opportunity_cost_nb_1 = 198002.809

# Initialization of variables

battery_on = False
# alpha = 0.001

for alpha in [0.0001]:
    number_of_days = 7  # number of days in a delivery period
    simulation_period = 52  # number of delivery periods in the simulation
    delivery_period = 24*number_of_days
    # delivery_mass = number_of_days*0.7*24*P_to_H2(1000, 33.3)
    total_mass = (108000/8760)*delivery_period*simulation_period
    # total_mass = 0
    delivery_mass = total_mass/simulation_period
    initial_battery = 0
    initial_hour = 0

    # capital_cost = (capital_cost_solar+capital_cost_wind+capital_cost_battery+capital_cost_electrolyzer+capital_cost_converter+capital_cost_grid)*delivery_period*simulation_period/8760
    # capital_cost = (capital_cost_solar+capital_cost_wind+capital_cost_electrolyzer+capital_cost_converter+capital_cost_grid)*delivery_period*simulation_period/8760
    # capital_cost = 0




    # benchmark_H2_prod = []
    benchmark_pf_series = []


    PF_wind = (PF_wind_2018+PF_wind_2019)[initial_hour+8760:initial_hour+8760+delivery_period*simulation_period]
    PF_solar = (PF_solar_2018+PF_solar_2019)[initial_hour+8760:initial_hour+8760+delivery_period*simulation_period]
    price = (price_2018+price_2019)[initial_hour+8760:initial_hour+8760+delivery_period*simulation_period]
    CO2int = (CO2int_2018+CO2int_2019)[initial_hour+8760:initial_hour+8760+delivery_period*simulation_period]

    hydrogen_plant = HydrogenProductionSystem(0,0,delivery_mass,0,battery_on, alpha, PF_wind, PF_solar, price, CO2int)
    hydrogen_plant.benchmark(delivery_period*simulation_period, delivery_period, delivery_mass, PF_wind, PF_solar, price, CO2int)
    net_cost = hydrogen_plant.benchmark_electricity_balance
    production = hydrogen_plant.benchmark_total_production
    emissions = hydrogen_plant.benchmark_CO2_emissions
    benchmark_H2_prod = hydrogen_plant.benchmark_H2_prod
    benchmark_pf_series = hydrogen_plant.benchmark_pf_series


    # hydrogen_plant = HydrogenProductionSystem(0,0,0,0,battery_on, alpha, PF_wind, PF_solar, price, CO2int)
    # hydrogen_plant.benchmark(delivery_period*simulation_period, delivery_period, 0, PF_wind, PF_solar, price, CO2int)
    # opportunity_cost = -hydrogen_plant.benchmark_electricity_balance

    opportunity_cost = opportunity_cost_nb_0


    capital_cost = capital_cost_electrolyzer
    H2_min_price = net_cost + opportunity_cost + capital_cost
    # H2_min_price = net_cost

    print(benchmark_pf_series)
    print(total_mass)
    print(production)
    print(net_cost)
    # print(hydrogen_plant.benchmark_d2d_CO2)
    print(H2_min_price/production)
    print(emissions/production)
    # print(price)
    # print(CO2int)
    # print(H2_min_price)
    # print(emissions)

    # figure, axis = plt.subplots(1, 2)
    # axis[0].plot(benchmark_H2_prod[0])
    # axis[1].plot(benchmark_pf_series[0][0])
    # plt.show()

    snapshots = ["Time (h)"]+list(range(delivery_period*simulation_period))

    solar_pf = ["Solar power (kW)"]+list(benchmark_pf_series[0][0]["Solar"])

    wind_pf = ["Wind power (kW)"]+list(benchmark_pf_series[0][0]["Wind"])

    powerfromgrid_pf = ["Power from grid (kW)"]+list(benchmark_pf_series[0][0]["NetworkImport"])

    power2grid_pf = ["Power to grid (kW)"]+list(benchmark_pf_series[0][0]["NetworkExport"])

    power2hydrogen = ["Power to hydrogen (kW)"]+list(benchmark_pf_series[0][1]["H2gen"])

    curtailed_solar = ["Curtailed solar power (kW)"]+(np.multiply(PF_solar, 1000) - np.array(benchmark_pf_series[0][0]["Solar"])).tolist()

    curtailed_wind = ["Curtailed wind power (kW)"]+(np.multiply(PF_wind, 1000) - np.array(benchmark_pf_series[0][0]["Wind"])).tolist()

    # CO2_prod = ["CO2 emitted (kg)"]+list(hydrogen_plant.benchmark_h2h_CO2.T)
    CO2int = ["CO2 intensity"]+CO2int
    price = ["Electricity price"]+price

    array = [snapshots, solar_pf, wind_pf, powerfromgrid_pf, power2grid_pf, power2hydrogen, curtailed_solar, curtailed_wind, CO2int, price, ["Minimum price", H2_min_price/production], ["H2 CO2 intensity", emissions/production]]
    df = pd.DataFrame(array).T
    if battery_on:
        df.to_excel(excel_writer = "C:/Users/sacha/PycharmProjects/projet_SIRD/résultats/benchmark/battery_on_period{nd}d_alpha{alpha}.xlsx".format(alpha=alpha, nd=number_of_days))
    else:
        df.to_excel(excel_writer = "C:/Users/sacha/PycharmProjects/projet_SIRD/résultats/opportunity_battery_off_period{nd}d_alpha{alpha}.xlsx".format(alpha=alpha, nd=number_of_days))