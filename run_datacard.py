from ggHdatacardmaker import *

def combine_workflow(cat, year, mass, lifetime, finalstate, physics):
    subprocess.run(["text2workspace.py", f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}.txt", "-o", f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}.root"])
    subprocess.run(["combine", f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}.root", "-M", "MultiDimFit", "--saveWorkspace", "--robustFit", "1","--cminDefaultMinimizerStrategy", "0", "-m", "125"])
    subprocess.run(["combine", f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}.root", "-M", "FitDiagnostics", "--saveShapes", "--saveWorkspace", "--saveWithUncertainties", "--saveNormalizations", "--robustFit", "1", "--cminDefaultMinimizerStrategy", "2", "-m", "125"])
    subprocess.run(["mv", "higgsCombineTest.MultiDimFit.mH125.root", f"higgsCombineTest.MultiDimFit.mH125_m{mass}_ct{lifetime}_{cat}_{year}.root"])
    subprocess.run(["mv", "fitDiagnosticsTest.root", f"fitDiagnosticsTest_m{mass}_ct{lifetime}_{cat}_{year}.root"])

def combined_datacard(year, cats, mass, lifetime, finalstate, physics):
    if len(cats)==1:
        print("ERROR: need to combine one or more categories.")
    if len(cats)>3:
        print("ERROR: only 3 categories exist. ")
    if len(cats)==2:
        subprocess.run(["combineCards.py", f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cats[0]}_{year}.txt",f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cats[1]}_{year}.txt", ">", f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cats[0]}_{cats[1]}_{year}.txt"])
    if len(cats)==3:
        subprocess.run(["combineCards.py", f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cats[0]}_{year}.txt",f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cats[1]}_{year}.txt",f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cats[2]}_{year}.txt", ">", f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cats[0]}_{cats[1]}_{cats[2]}_{cats[3]}_{year}.txt"])

def run(signal, bkg, cat, year, mass, lifetime, finalstate, physics, bins):
    main(paths=[signal, bkg], isMC=[1,0], trees=["ggH4g","ggH4g"], var=f"best_4g_corr_mass_m{mass}", categories=categories,period=year, bins=bins, lifetime=lifetime, mass=mass)
    combine_workflow(cat, year,mass, lifetime, finalstate, physics)
    #plot(f"higgsCombineTest.MultiDimFit.mH125_m{mass}_ct{lifetime}_{cat}_{year}.root", f"fitDiagnosticsTest_m{mass}_ct{lifetime}_{cat}_{year}.root",f"{cat}", f"{year}", bins=bins, finalstate=finalstate, physics=physics, order=3)

if __name__=="__main__":
    parser = argparse.ArgumentParser("Datacard Processing for year, category")
    parser.add_argument("-s","--signal", type=str, help="Signal file")
    parser.add_argument("-b","--background", type=str, help="Background file")
    parser.add_argument("-m","--mass", type=str, help="mass of sample")
    parser.add_argument("-ct","--ctau", type=str, help="lifetime of sample")
    parser.add_argument("-y","--year", type=str, help="year of MC and data")
    parser.add_argument("-c1","--cat1", dest="c1", type=str, help="choose one of: low, high, asym, all")
    parser.add_argument("-c2","--cat2", dest="c2", type=str, help="choose one of: low, high, asym, all")
    parser.add_argument("-c3","--cat3", dest="c3", type=str, help="choose one of: low, high, asym, all")
    parser.add_argument("-c4","--cat4", dest="c4", type=str, help="choose one of: low, high, asym, all")

    args = parser.parse_args()
    mass=args.mass
    lifetime=args.ctau
    year=args.year
    sig=args.signal
    bkg=args.background

    categories = []
    for c in ["c1", "c2", "c3", "c4"]:
        argu = getattr(args, c)
        if argu is not None:
            categories.append(argu)
    for cat in categories:
        run(sig, bkg, cat, year, mass, lifetime,"4g", "ggH", bins=[60,110,140])
    
    if len(categories)>1:
        combined_datacard(year,categories,mass,lifetime,"4g","ggH")
