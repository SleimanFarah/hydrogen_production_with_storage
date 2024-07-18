import argparse
import multiprocessing as mp
from main_script import run_system_simulation


def parallel_run(alphas):
    pool = mp.Pool()
    pool.map(run_system_simulation, alphas)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="parallel_processing_input")
    parser.add_argument('--alphas', type=float, nargs="+", help="alpha values")
    alphas = parser.parse_args().alphas

    # alphas = [0.7, 0.8]

    parallel_run(alphas)
