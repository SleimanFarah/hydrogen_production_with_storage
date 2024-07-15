import argparse
import multiprocessing as mp
from main_script import run_system_simulation


def parallel_run(alphas, num_processes):
    pool = mp.Pool(num_processes)
    pool.map(run_system_simulation, alphas)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="parallel_processing_input")
    parser.add_argument('--alphas', type=float, nargs="+", help="alpha values")
    alphas = parser.parse_args().alphas

    # alphas = [0.0001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.999]

    parallel_run(alphas, 3)
