from class_network import *
from convert_functions import *
import random
from csv_file_function import *
from function_AC import *

year = 2018

with open(f'data\{year}\solar_capacity_factor.csv', 'r') as PF_solar_data:
    PF_solar = csv_reader_function(PF_solar_data)

with open(f'data\{year}\wind_capacity_factor.csv', 'r') as PF_wind_data:
    PF_wind = csv_reader_function(PF_wind_data)

with open(f'data\{year}\electricity_price.csv', 'r') as price_data:
    price = np.multiply(csv_reader_function(price_data), 0.001)

with open(f'data\{year}\co2_intensity.csv', 'r') as CO2int_data:
    CO2int = np.multiply(csv_reader_function(CO2int_data), 0.001)

# Initialization of variables

number_of_days = 1
simulation_period = 1  # number of delivery periods in the simulation
delivery_period = 24*number_of_days
# delivery_mass = number_of_days*0.7*24*P_to_H2(1000, 33.3)
delivery_mass = 150
initial_battery = 0
initial_hour = 4000


PF_wind = PF_wind[initial_hour:]
PF_solar = PF_solar[initial_hour:]
price = price[initial_hour:]
CO2int = CO2int[initial_hour:]

# PF_wind = [random.uniform(0.4, 1) for _ in range(delivery_period+48)]
# PF_solar = [random.uniform(0.4, 1) for _ in range(delivery_period+48)]
# price = [random.uniform(1, 10) for _ in range(delivery_period+48)]
# CO2int = [random.uniform(1, 10) for _ in range(delivery_period+48)]


ltp_pf_series = []
ltp_H2_series = []
ltp_H2_target = []
dp_pf_series = []
dp_Hydro_plan = []
dr_pf_series = []
dr_H2_prod = []

i = 0
j = 0
time_left = delivery_period
H2_mass_remaining = delivery_mass
total_production = 0
total_emissions = 0
total_electricity_cost = 0
total_electricity_balance = (Annualized_cost(0, 50, 1000000, 10000)/8760)*delivery_period
# Loop for the whole delivery period

while j < simulation_period:
    while i < number_of_days:
        print(j, i)
        ltp_PF_wind = PF_wind[10+i*24+j*delivery_period:delivery_period+34+i*24+j*delivery_period]
        ltp_PF_solar = PF_solar[10+i*24+j*delivery_period:delivery_period+34+i*24+j*delivery_period]
        ltp_price = price[10+i*24+j*delivery_period:delivery_period+34+i*24+j*delivery_period]
        ltp_CO2int = CO2int[10+i*24+j*delivery_period:delivery_period+34+i*24+j*delivery_period]

        hydrogen_plant = HydrogenProductionSystem(time_left, delivery_period, H2_mass_remaining, initial_battery, ltp_PF_wind, ltp_PF_solar, ltp_price, ltp_CO2int)

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
        total_electricity_cost = total_electricity_cost + hydrogen_plant.electricity_cost
        total_electricity_balance = total_electricity_balance + hydrogen_plant.electricity_balance

        i = i + 1
        time_left = delivery_period - i*24
        total_production = total_production + hydrogen_plant.total_production
        H2_mass_remaining = delivery_mass*(j+1) - total_production
        initial_battery = hydrogen_plant.battery_left
    i = 0
    time_left = delivery_period
    H2_mass_remaining = delivery_mass
    j = j+1


print(ltp_pf_series)
print(ltp_H2_target)
# print(dp_Hydro_plan)
print(dp_pf_series)
print(dr_pf_series)
print(total_production)
print(delivery_mass*j)
# print(dr_pf_series[0][1].sum())
print(total_electricity_balance/total_production)
print(total_emissions/total_production)
