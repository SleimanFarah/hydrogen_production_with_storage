from class_network import *
from convert_functions import *
import random
from matplotlib import pyplot as plt
from csv_file_function import *
from function_AC import *

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

# Initialization of variables

battery_on = False
alpha = 0

number_of_days = 1  # number of days in a delivery period
simulation_period = 5  # number of delivery periods in the simulation
delivery_period = 24*number_of_days
# delivery_mass = number_of_days*0.7*24*P_to_H2(1000, 33.3)
total_mass = (108000/8760)*delivery_period*simulation_period
# total_mass = 0
delivery_mass = total_mass/simulation_period
initial_battery = 0
initial_hour = 0

capital_cost_solar = Annualized_cost(0, 40, 310000, 9500)
capital_cost_wind = Annualized_cost(0, 30, 990000, 12600)
capital_cost_battery = Annualized_cost(0, 40, 700000, 540)
capital_cost_electrolyzer = Annualized_cost(0, 10, 700000, 14000)
capital_cost_converter = Annualized_cost(0, 15, 40000, 0)
capital_cost_grid = Annualized_cost(0, 40, 140000, 2800)

capital_cost = (capital_cost_solar+capital_cost_wind+capital_cost_battery+capital_cost_electrolyzer+capital_cost_converter+capital_cost_grid)*delivery_period*simulation_period/8760



PF_wind = (PF_wind_2018+PF_wind_2019)[initial_hour+8736-delivery_period:]
PF_solar = (PF_solar_2018+PF_solar_2019)[initial_hour+8736-delivery_period:]
price = (price_2018+price_2019)[initial_hour+8736-delivery_period:]
CO2int = (CO2int_2018+CO2int_2019)[initial_hour+8736-delivery_period:]

ltp_pf_series = []
ltp_H2_series = []
ltp_H2_target = []
dp_pf_series = []
dp_Hydro_plan = []
dr_pf_series = []
dr_H2_prod = []
CO2_d2d = []

i = 0
j = 0
time_left = delivery_period
H2_mass_remaining = delivery_mass
total_production = 0
total_emissions = 0
total_electricity_buying_cost = 0
total_H2_min_cost = capital_cost
total_opportunity_cost = 0
# total_electricity_net_cost = (Annualized_cost(0, 50, 1000000, 10000)/8760)*delivery_period*simulation_period

# Benchmark optimization initialization

benchmark_PF_wind = PF_wind[delivery_period:delivery_period*(simulation_period+1)]
benchmark_PF_solar = PF_solar[delivery_period:delivery_period*(simulation_period+1)]
benchmark_price = price[delivery_period:delivery_period*(simulation_period+1)]
benchmark_CO2int = CO2int[delivery_period:delivery_period*(simulation_period+1)]



# Loop for the whole delivery period

while j < simulation_period:
    while i < number_of_days:
        print(j+1, i+1)
        operative_PF_wind = PF_wind[10+i*24*2+j*delivery_period:delivery_period+34+i*24+j*delivery_period]
        operative_PF_solar = PF_solar[10+i*24*2+j*delivery_period:delivery_period+34+i*24+j*delivery_period]
        operative_price = price[10+i*24*2+j*delivery_period:delivery_period+34+i*24+j*delivery_period]
        operative_CO2int = CO2int[10+i*24*2+j*delivery_period:delivery_period+34+i*24+j*delivery_period]

        hydrogen_plant = HydrogenProductionSystem(time_left, delivery_period, H2_mass_remaining, initial_battery, battery_on, alpha, operative_PF_wind, operative_PF_solar, operative_price, operative_CO2int)

        hydrogen_plant.lt_planner()
        hydrogen_plant.daily_planner()
        hydrogen_plant.realisation()

        ltp_pf_series = ltp_pf_series + hydrogen_plant.ltp_pf_series
        ltp_H2_series = ltp_H2_series + hydrogen_plant.ltp_H2_series
        ltp_H2_target = ltp_H2_target + hydrogen_plant.ltp_H2_target
        dp_pf_series = dp_pf_series + hydrogen_plant.dp_pf_series
        dp_Hydro_plan = dp_Hydro_plan + hydrogen_plant.dp_Hydro_plan
        dr_pf_series = dr_pf_series + hydrogen_plant.dr_pf_series
        dr_H2_prod = dr_H2_prod + hydrogen_plant.dr_H2_prod

        total_emissions = total_emissions + hydrogen_plant.CO2_emissions
        total_electricity_buying_cost = total_electricity_buying_cost + hydrogen_plant.electricity_cost
        net_cost = hydrogen_plant.electricity_net_cost


        total_production = total_production + hydrogen_plant.total_production
        H2_mass_remaining = delivery_mass*(j+1) - total_production

        hydrogen_plant = HydrogenProductionSystem(time_left, delivery_period, 0, initial_battery, battery_on, alpha, operative_PF_wind, operative_PF_solar, operative_price, operative_CO2int)


        hydrogen_plant.lt_planner()
        hydrogen_plant.daily_planner()
        hydrogen_plant.realisation()
        opportunity_cost = -hydrogen_plant.electricity_net_cost
        total_opportunity_cost = total_opportunity_cost + opportunity_cost
        CO2_d2d = CO2_d2d + hydrogen_plant.dr_d2d_CO2

        total_H2_min_cost = (total_H2_min_cost + net_cost + opportunity_cost)

        i = i + 1
        time_left = delivery_period - i*24
        initial_battery = hydrogen_plant.battery_left
    i = 0
    time_left = delivery_period
    H2_mass_remaining = delivery_mass
    j = j+1




# print(ltp_pf_series)
# print(ltp_H2_target)
# # print(dp_Hydro_plan)
# print(dp_pf_series)
print(dr_H2_prod)
print(dr_pf_series)
print(total_production)
print(delivery_mass*j)
# print(dr_pf_series[0][1].sum())
# print(total_electricity_buying_cost/total_production)
print(total_opportunity_cost/total_production)
print(total_H2_min_cost/total_production)
print(total_emissions/total_production)
# print(total_H2_min_cost)
# print(total_emissions)

figure, axis = plt.subplots(1, 3)
axis[0].plot(dr_H2_prod[0])
axis[1].plot(dr_pf_series[0][0])
axis[2].plot(CO2_d2d)
plt.show()

