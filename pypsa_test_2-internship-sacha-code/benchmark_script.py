from class_network import *
from convert_functions import *
import random
from csv_file_function import *
from function_AC import *

with open('solar_capacity_factor.csv', 'r') as PF_solar_data:
    PF_solar = csv_reader_function(PF_solar_data)

with open('wind_capacity_factor.csv', 'r') as PF_wind_data:
    PF_wind = csv_reader_function(PF_wind_data)

with open('electricity_price.csv', 'r') as price_data:
    price = np.multiply(csv_reader_function(price_data), 0.001)

with open('co2_intensity.csv', 'r') as CO2int_data:
    CO2int = np.multiply(csv_reader_function(CO2int_data), 0.001)

number_of_days = 1  # number of days in a delivery period
simulation_period = 1  # number of delivery periods in the simulation
delivery_period = 24*number_of_days
# delivery_mass = number_of_days*0.7*24*P_to_H2(1000, 33.3)
total_mass = 296
delivery_mass = total_mass/simulation_period
initial_battery = 0
initial_hour = 0

PF_wind = PF_wind[initial_hour:]
PF_solar = PF_solar[initial_hour:]
price = price[initial_hour:]
CO2int = CO2int[initial_hour:]