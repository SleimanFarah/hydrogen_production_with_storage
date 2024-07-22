from pypsa import *
from gurobipy import *
from function_AC import Annualized_cost
import pypsatopo


n = Network()
nsnap=2
n.set_snapshots(range(nsnap))
n.add("Bus", "ACBus")
n.add("Bus", "DCBus")
n.add("Bus", "H2Bus")
n.add("Bus", "BatteryBus")
n.add("Bus", "NetworkBus")


n.add("Generator", "Wind", bus="ACBus", p_nom=100, p_max_pu=[4, 2], p_min_pu=[2, 0.5], marginal_cost=0)
n.add("Generator", "Solar", bus="DCBus", p_nom=50, p_max_pu=[3, 3], p_min_pu=[1, 0], marginal_cost=0)
n.add("Generator", "GlobalNetwork", bus="NetworkBus", p_nom_extendable=True, p_max_pu=1, p_min_pu=-1)


# n.add("Store", "H2gen", bus="H2Bus", e_nom_extendable=True)
# n.stores_t.p_min_pu["H2gen"]=-400

n.add("Load", "H2gen", bus="H2Bus")
n.loads_t.p_set["H2gen"] = 400

n.add("Store", "Battery", bus="BatteryBus", e_nom_extendable=True, e_initial=0, e_nom_max=1000)


eff_charge = 1
eff_discharge = 1
eff_electrolysis = 0.5
eff_converter = 1

sell_price = 9
buy_price = 10
CO2_intensity = 10

n.add("Link", "H2Link", bus0="ACBus", bus1="H2Bus", efficiency=eff_electrolysis, p_nom_extendable=True)
n.add("Link", "ChargeLink", bus0="DCBus", bus1="BatteryBus", efficiency=eff_charge, p_nom_extendable=True)
n.add("Link", "DischargeLink", bus0="BatteryBus", bus1="DCBus", efficiency=eff_discharge, p_nom_extendable=True)
n.add("Link", "ACtoDCLink", bus0="ACBus", bus1="DCBus", efficiency=eff_converter, p_nom_extendable=True)
n.add("Link", "DCtoACLink", bus0="DCBus", bus1="ACBus", efficiency=eff_converter, p_nom_extendable=True)
n.add("Link", "BuyLink", bus0="NetworkBus", bus1="ACBus", marginal_cost=buy_price, p_nom_extendable=True)
n.add("Link", "SellLink", bus0="ACBus", bus1="NetworkBus", marginal_cost=-sell_price, p_nom_extendable=True)

alpha=0.5
timestep=1

m = n.optimize.create_model()

total_energy = m.variables["Link-p"].loc[0,"BuyLink"]
for i in range(1,nsnap):
    total_energy = total_energy + m.variables["Link-p"].loc[i,"BuyLink"]*timestep
expr = total_energy*(alpha*CO2_intensity+(1-alpha)*buy_price)
m.add_objective(expr, overwrite=True, sense="min")

# total_hydrolyzer_energy = m.variables["Store-p"].loc[0, "H2gen"]
# for i in range(1, 2):
#     total_hydrolyzer_energy = total_hydrolyzer_energy + m.variables["Store-p"].loc[i, "H2gen"]
# cstr = total_hydrolyzer_energy == -1
# m.add_constraints(cstr, name="target match")

n.optimize.solve_model(solver_name="gurobi")

# print(n.lines_t.p0)n.optimize.solve_model(solver_name="gurobi")
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
print(n.loads_t.p.sum())
