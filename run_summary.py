import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal")
from plotting import plot_summary
import argparse


def run(mass, lifetime, year):
    #summary plots build their own histograms from the ntuples; no datacard needed
    plot_summary.run(mass, lifetime, year)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Summary plots for year, mass, lifetime")
    parser.add_argument("-m", "--mass", type=str, help="mass of sample")
    parser.add_argument("-ct", "--ctau", type=str, help="lifetime of sample")
    parser.add_argument("-y", "--year", type=str, help="year of MC and data (single year or Run2/Run3/2022/2023 aggregates)")

    parser.add_argument("-process_run2", "--process_run2", dest="process_run2", type=int, help="Summary plots for all Run 2. 1=yes, 0=no. Runs all mass/lifetime points")
    parser.add_argument("-process_run3", "--process_run3", dest="process_run3", type=int, help="Summary plots for all Run 3. 1=yes, 0=no. Runs all mass/lifetime points")

    args = parser.parse_args()
    mass = args.mass
    lifetime = args.ctau
    year = args.year
    process_run2 = args.process_run2
    process_run3 = args.process_run3

    if process_run2:
        for m in ["15", "20", "30", "40", "50", "55"]:
            for ct in ["0", "10", "20", "50", "100", "1000"]:
                for year in ["2017", "2018"]:  # no longer using 2016 because the trigger is too inefficient
                    run(m, ct, year)

    elif process_run3:
        for m in ["15", "20", "30", "40", "50", "55"]:
            for ct in ["0", "10", "20", "50", "100", "1000"]:
                for year in ["2022", "2023", "2024"]:
                    run(m, ct, year)

    else:
        run(mass, lifetime, year)
