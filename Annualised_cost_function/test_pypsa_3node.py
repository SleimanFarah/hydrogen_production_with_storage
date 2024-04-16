import pypsa
import pypsatopo
import numpy as np

network = pypsa.Network()

n_buses=3
for i in range(n_buses):
    network.add("Bus", "My bus {}".format(i), v_nom=20.0)


for i in range(n_buses):
    network.add(
        "Line",
        "My line {}".format(i),
        bus0="My bus {}".format(i),

        bus1="My bus {}".format((i + 1) % n_buses),
        x=0.1,
        r=0.01,
    )

network.lines

network.add("Generator", "My gen", bus="My bus 0", p_set=100, control="PQ")


network.add("Load", "My load", bus="My bus 1", p_set=100)


network.loads.q_set = 100.0



