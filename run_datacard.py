from ggHdatacardmaker import *

if __name__=="__main__":
    parser = argparse.ArgumentParser("Datacard Processing for year, category")
    parser.add_argument("-s","--signal", type=str, help="Signal file")
    parser.add_argument("-b","--background", type=str, help="Background file")
    parser.add_argument("-m","--mass", type=str, help="mass of sample")
    parser.add_argument("-ct","--ctau", type=str, help="lifetime of sample")
    parser.add_argument("-y","--year", type=str, help="year of MC and data")
    parser.add_argument("-c1","--cat1", dest="c1", type=str, help="one of: low, high, asym, all")
    parser.add_argument("-c2","--cat2", dest="c2", type=str, help="one of: low, high, asym, all")
    parser.add_argument("-c3","--cat3", dest="c3", type=str, help="one of: low, high, asym, all")
    parser.add_argument("-c4","--cat4", dest="c4", type=str, help="one of: low, high, asym, all")

    args = parser.parse_args()
    mass=args.mass
    lifetime=args.ctau
    year=args.year

    categories = []
    for c in ["c1", "c2", "c3", "c4"]:
        argu = getattr(args, c)
        if argu is not None:
            categories.append(argu)
    print(categories)
    main(paths=[args.signal, args.background], isMC=[1,0], trees=["ggH4g","ggH4g"], var=f"best_4g_corr_mass_m{mass}", categories=categories,period=year, bins=[33, 120, 130], lifetime=lifetime, mass=mass)
    for cat in categories:
        subprocess.run(["text2workspace.py", f"datacard_ggH_4g_m{mass}_ct{lifetime}_{cat}_{year}.txt", "-o", f"datacard_ggH_4g_m{mass}_ct{lifetime}_{cat}_{year}.root"])
        subprocess.run(["combine", "-M", "FitDiagnostics", f"datacard_ggH_4g_m{mass}_ct{lifetime}_{cat}_{year}.root", "--saveShapes", "--saveWithUncertainties", "--saveNormalizations", "--saveWorkspace", "-m", "125"])
        subprocess.run(["mv", "fitDiagnosticsTest.root", f"fitDiagnosticsTest_m{mass}_ct{lifetime}_{cat}_{year}.root"])
