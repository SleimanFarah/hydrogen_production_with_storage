import argparse
import multiprocessing as mp
from class_network import run_system_simulation


def helper(args):
    year, alpha, time_period = args
    return run_system_simulation(year, alpha, time_period)


def parallel_run(years, alphas, time_periods):
    pool = mp.Pool()
    pool.map(helper, [(year, alpha, time_period) for year in years for alpha in alphas for time_period in time_periods])


if __name__ == "__main__":

    years = [2018, 2019, 2020, 2021]
    time_periods = ["day"]

    parser = argparse.ArgumentParser(description="parallel_processing_input")
    parser.add_argument('--alphas', type=float, nargs="+", help="alpha values")
    alphas = parser.parse_args().alphas

    parallel_run(years, alphas, time_periods)
