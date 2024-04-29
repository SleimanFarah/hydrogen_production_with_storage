from pypsa import *
from gurobipy import *
from function_AC import Annualized_cost
import pypsatopo


n = Network()
n.set_snapshots(range(2))
vbus = 230
n.add("Bus", "ACBus", v_nom=vbus)
n.add("Bus", "DCBus", v_nom=vbus)
n.add("Bus", "H2Bus", v_nom=vbus)
n.add("Bus", "BatteryBus", v_nom=vbus)


n.add("Generator", "Wind", bus="ACBus", p_nom=100, p_max_pu=[0.5, 2], p_min_pu=[0.5, 2], control="PV", marginal_cost=0)
n.add("Generator", "Solar", bus="DCBus", p_nom=50, p_max_pu=[0, 2], p_min_pu=[0, 2], control="PV", marginal_cost=0)
n.add("Generator", "GlobalNetwork", bus="ACBus", control="Slack", marginal_cost=1, p_nom_extendable=True)


n.add("Load", "H2gen", bus="H2Bus", p_set=[300, 50])


n.add("Store", "Battery", bus="BatteryBus", e_nom_extendable=True, e_nom_max=100)


eff_charge = 0.9
eff_discharge = 0.9
eff_electrolysis = 1
eff_converter = 0.95

n.add("Link", "H2Link", bus0="ACBus", bus1="H2Bus", efficiency=eff_electrolysis, p_nom_extendable=True)
n.add("Link", "ChargeLink", bus0="DCBus", bus1="BatteryBus", efficiency=eff_charge, p_nom_extendable=True)
n.add("Link", "DischargeLink", bus0="BatteryBus", bus1="DCBus", efficiency=eff_discharge, p_nom_extendable=True)
n.add("Link", "ACtoDCLink", bus0="ACBus", bus1="DCBus", efficiency=eff_converter, p_nom_extendable=True)
n.add("Link", "DCtoACLink", bus0="DCBus", bus1="ACBus", efficiency=eff_converter, p_nom_extendable=True)

n.optimize(solver_name="gurobi")

# print(n.lines_t.p0)
# print(n.links_t.p0)

# pypsatopo.generate(n, file_output="mynetwork.cvs")
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
