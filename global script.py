from pypsa import *
from gurobipy import *
from function_AC import Annualized_cost
import pypsatopo

# Recuperation of the time series
wind_PF = [2, 0.5]
solar_PF = [2, 1]

