import ROOT
import subprocess
from datacard.ggHdatacardmaker import main as make_datacard
from run_datacard import combined_datacard
from ggHparameters import signal_path, bkg_path


def get_median_expected_UL(root_file):
    f = ROOT.TFile.Open(root_file)
    limit_tree = f.Get("limit")
    median_UL = None
    for entry in limit_tree:
        if abs(entry.quantileExpected - 0.5) < 1e-4:
            median_UL = entry.limit
    f.Close()
    return median_UL


def scan_mass_lifetime(masses, lifetimes, years, categories, bins,
                        finalstate="4g", physics="ggH", order_fit=4):
    ROOT.gROOT.SetBatch(True)
    results = []

    for mass in masses:
        for ctau in lifetimes:
            for year in years:
                sig = signal_path(mass, ctau, year)
                bkg = bkg_path(year)

                for cat in categories:
                    make_datacard(paths=[sig, bkg], isMC=[1, 0], trees=["ggH4g", "ggH4g"],
                                   var=f"best_4g_corr_mass_m{mass}", categories=[cat], period=year,
                                   bins=bins, lifetime=ctau, mass=mass, finalstate=finalstate,
                                   physics=physics, order_fit=order_fit)

                combined_datacard(year, categories, mass, ctau, finalstate, physics)

                combined_txt = f"datacard_{physics}_{finalstate}_m{mass}_ct{ctau}_combined_{year}.txt"
                combined_root = f"datacard_{physics}_{finalstate}_m{mass}_ct{ctau}_combined_{year}.root"
                subprocess.run(["text2workspace.py", combined_txt, "-o", combined_root])
                subprocess.run(["combine", "-M", "AsymptoticLimits", combined_root, "-m", "125"])

                result_root = f"higgsCombineTest_m{mass}_ct{ctau}_{year}.AsymptoticLimits.mH125.root"
                subprocess.run(["mv", "higgsCombineTest.AsymptoticLimits.mH125.root", result_root])

                median_UL = get_median_expected_UL(result_root)
                print(f"m={mass} ct={ctau} year={year}: median expected UL on r = {median_UL}")
                results.append((mass, ctau, median_UL, year))

    return results


if __name__ == "__main__":
    masses = [15]
    lifetimes = [10]
    years = ["2018"]
    categories = ["low", "asym", "high"]

    results = scan_mass_lifetime(masses, lifetimes, years, categories, bins=[60, 110, 140])
    print(results)
