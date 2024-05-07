from pypsa import *
from gurobipy import *
from function_AC import Annualized_cost
import pypsatopo


n = Network()
n.set_snapshots(range(2))
n.add("Bus", "ACBus")
n.add("Bus", "DCBus")
n.add("Bus", "H2Bus")
n.add("Bus", "BatteryBus")
n.add("Bus", "NetworkBus")


n.add("Generator", "Wind", bus="ACBus", p_nom=100, p_max_pu=[2, 0.5], p_min_pu=[2, 0.5], marginal_cost=0)
n.add("Generator", "Solar", bus="DCBus", p_nom=50, p_max_pu=[1, 0], p_min_pu=[1, 0], marginal_cost=0)
n.add("Generator", "GlobalNetwork", bus="NetworkBus", p_nom_extendable=True, p_max_pu=1, p_min_pu=-1)


n.add("Load", "H2gen", bus="H2Bus", p_set=200)


n.add("Store", "Battery", bus="BatteryBus", e_nom_extendable=True, e_initial=0, e_nom_max=1000)


eff_charge = 1
eff_discharge = 1
eff_electrolysis = 1
eff_converter = 1

sell_price = 9
buy_price = 10

n.add("Link", "H2Link", bus0="ACBus", bus1="H2Bus", efficiency=eff_electrolysis, p_nom_extendable=True)
n.add("Link", "ChargeLink", bus0="DCBus", bus1="BatteryBus", efficiency=eff_charge, marginal_cost=-(sell_price+0.01), p_nom_extendable=True)
n.add("Link", "DischargeLink", bus0="BatteryBus", bus1="DCBus", efficiency=eff_discharge, marginal_cost=(buy_price-0.01), p_nom_extendable=True)
n.add("Link", "ACtoDCLink", bus0="ACBus", bus1="DCBus", efficiency=eff_converter, p_nom_extendable=True)
n.add("Link", "DCtoACLink", bus0="DCBus", bus1="ACBus", efficiency=eff_converter, p_nom_extendable=True)
n.add("Link", "BuyLink", bus0="NetworkBus", bus1="ACBus", marginal_cost=buy_price, p_nom_extendable=True)
n.add("Link", "SellLink", bus0="ACBus", bus1="NetworkBus", marginal_cost=-sell_price, p_nom_extendable=True)

m = n.optimize.create_model()


n.optimize(solver_name="gurobi")

# print(n.lines_t.p0)
# print(n.links_t.p0)

# pypsatopo.generate(n, file_output="mynetwork2.cvs")
# n.plot(
#     line_colors="r",
#     title="network",
# )
# print(n.buses)

print(n.generators_t.p)
print(n.loads_t.p)
print(n.stores_t.p)
print(n.stores_t.e)
# n.generator.loc["GlobalNetwork","p"]
