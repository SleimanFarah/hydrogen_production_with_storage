import pypsa
import gurobipy as gp
from function_AC import Annualized_cost
n=pypsa.Network(name="mynetwork")



for i in range(1):
    n.add("Bus", "bus{}".format(i), v_nom=220)

#n.add("Line", "line01", bus0="bus0", bus1="bus1", x=0.1, r=0.01)
#n.add("Line", "line02", bus0="bus0", bus1="bus2", x=0.1, r=0.01)

n.add("Generator", "gas_gen", bus="bus0", p_nom_max=100, p_nom_extendable=True, carrier="gas", marginal_cost=70, capital_cost=Annualized_cost(0, 40, 100000000, 100000))
n.add("Generator", "wind_gen", bus="bus0", p_nom_max=100, p_nom_extendable=True, carrier="wind", marginal_cost=0, capital_cost=Annualized_cost(0, 25, 200000000, 150000))

n.add("Load", "myload", bus="bus0", p_set=150)

m=n.optimize.create_model()

gas = m.variables["Generator-p"].loc["now","gas_gen"]
wind = m.variables["Generator-p"].loc["now","wind_gen"]

expr= wind >= gas
m.add_constraints(expr, name="renewable_use")
m.add_objective(wind.to_linexpr(), overwrite=True, sense="max")


n.optimize.solve_model(solver_name="gurobi")
#m.solve(solver_name="gurobi")
print(n.generators_t.p)

