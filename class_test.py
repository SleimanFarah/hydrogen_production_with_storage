class HydrogenProductionSystem:

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def sum(self):
        self.sum = self.a+self.b
        print(self.sum)


# if __name__ == '__main__':

    # This could be the long-term planner object

ltp_hydrogen_production_system = HydrogenProductionSystem(a=5, b=6)
ltp_hydrogen_production_system.sum()

    # This could be the daily planner object
    # dp_hydrogen_production_system = HydrogenProductionSystem(2, 55).sum()
