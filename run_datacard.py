import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal")
from datacard.ggHdatacardmaker import *
from ggHparameters import bins

def combine_workflow(cat, year, mass, lifetime, finalstate, physics):
    card_txt = f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}.txt"
    card_root = f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}.root"

    subprocess.run(["text2workspace.py", card_txt, "-o", card_root], check=True)
    subprocess.run(["combine", card_root,"-M", "MultiDimFit","--saveWorkspace","--robustFit", "1","--cminDefaultMinimizerStrategy", "2","-m", "125"], check=True)
    subprocess.run(["combine", card_root,"-M", "FitDiagnostics","--saveShapes","--saveWorkspace","--saveWithUncertainties","--saveNormalizations","--robustFit", "1","--cminDefaultMinimizerStrategy", "2","-m", "125"], check=True)
    subprocess.run(["mv","higgsCombineTest.MultiDimFit.mH125.root",f"higgsCombineTest.MultiDimFit.mH125_m{mass}_ct{lifetime}_{cat}_{year}.root"], check=True)
    subprocess.run(["mv","fitDiagnosticsTest.root",f"fitDiagnosticsTest_m{mass}_ct{lifetime}_{cat}_{year}.root"], check=True)

def combined_datacard(year, cats, mass, lifetime, finalstate, physics):
    if len(cats) < 2:
        print("Need 2+ categories to combine.")
        return
    if len(cats) > 3:
        print("Only 3 categories exist.")
        return

    inputs = [f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{c}_{year}.txt" for c in cats]
    out = f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_combined_{year}.txt"

    with open(out, "w") as f:
        subprocess.run(["combineCards.py", *inputs], stdout=f, check=True)

    green = "\033[1;92m"
    reset = "\033[0m"
    print(f"{green}Success. The following datacards have been made:{reset}")
    for card in inputs:
        print(f"{green}  - {card}{reset}")
    print(f"{green}  - {out}  (combined){reset}")


def run(signal, bkg, cat, year, mass, lifetime, finalstate, physics, bins):
    #only builds the datacard; run combine manually afterwards
    mustard = "\033[38;5;136m"
    teal = "\033[38;5;44m"
    reset = "\033[0m"
    border = "============================================================"
    dirname = f"m{mass}_ct{lifetime}_{cat}_{year}_{finalstate}_{physics}"

    print(f"{teal}{border}{reset}")
    print(f"{mustard}• searching for pre-existing directory: {dirname}{reset}")
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)
        print(f"{mustard}  ↳ found. overriding previous datacard results.{reset}")
    else:
        print(f"{mustard}  ↳ not found, building directory: {dirname}{reset}")

    main(paths=[signal, bkg], isMC=[1,0], trees=["ggH4g","ggH4g"], var=f"best_4g_corr_mass_m{mass}", categories=[cat],period=year, bins=bins, lifetime=lifetime, mass=mass, lumi_scaling=1)

    print(f"{mustard}• updated directory {dirname}.{reset}")

if __name__=="__main__":
    parser = argparse.ArgumentParser("Datacard Processing for year, category")
    parser.add_argument("-s","--signal", type=str, help="Signal file")
    parser.add_argument("-b","--background", type=str, help="Background file")
    parser.add_argument("-m","--mass", type=str, help="mass of sample")
    parser.add_argument("-ct","--ctau", type=str, help="lifetime of sample")
    parser.add_argument("-y","--year", type=str, help="year of MC and data")
    parser.add_argument("-c1","--cat1", dest="c1", type=str, help="choose one of: prompt, displaced, asym, all")
    parser.add_argument("-c2","--cat2", dest="c2", type=str, help="choose one of: prompt, displaced, asym, all")
    parser.add_argument("-c3","--cat3", dest="c3", type=str, help="choose one of: prompt, displaced, asym, all")
    parser.add_argument("-c4","--cat4", dest="c4", type=str, help="choose one of: prompt, displaced, asym, all")

    parser.add_argument("-process_run2","--process_run2", dest="process_run2", type=int, help="Process all Run 2. 1=yes, 0=no. Will run all mass/lifetime points for all categories and combine the cards")
    parser.add_argument("-process_run3","--process_run3", dest="process_run3", type=int, help="Process all Run 3. 1=yes, 0=no. Will run all mass/lifetime points for all categories and combine the cards")


    args = parser.parse_args()
    mass=args.mass
    lifetime=args.ctau
    year=args.year
    sig=args.signal
    bkg=args.background
    process_run2=args.process_run2
    process_run3=args.process_run3

    if process_run2:
        categories=["prompt", "asym", "displaced"]
        for m in [15,20,30,40,50,55]:
            for ct in [0,10,20,50,100,1000]:
                for year in ["2017", "2018"]: # no longer using 2016 because the trigger is too inefficient
                    sig=f"/eos/uscms/store/user/dacampos/analysis/signal/ggH4g_M{m}_ctau{ct}_{year}_0_ggH4g_M{m}_ctau{ct}_{year}_ggH4g.root"
                    bkg=f"/eos/uscms/store/user/dacampos/analysis/data/EGamma_{year}_updated/EGamma_{year}_all_ggH4g.root"
                    for cat in categories:
                        run(sig, bkg, cat, year, m, ct,"4g", "ggH", bins=bins)
                    combined_datacard(year,categories,m,ct,"4g", "ggH")

    elif process_run3:
        categories=["prompt", "asym", "displaced"]
        for m in [15,20,30,40,50,55]:
            for ct in [0,10,20,50,100,1000]:
                for year in ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]:
                    sig=f"/eos/uscms/store/user/dacampos/analysis/signal/ggH4g_M{m}_ctau{ct}_{year}_0_ggH4g_M{m}_ctau{ct}_{year}_ggH4g.root"
                    bkg=f"/eos/uscms/store/user/dacampos/analysis/data/EGamma_{year}_updated/EGamma_{year}_all_ggH4g.root"
                    for cat in categories:
                        run(sig, bkg, cat, year, m, ct,"4g", "ggH", bins=bins)
                    combined_datacard(year,categories,m,ct,"4g", "ggH")

    else:
        categories = []
        for c in ["c1", "c2", "c3", "c4"]:
            argu = getattr(args, c)
            if argu is not None:
                categories.append(argu)
        for cat in categories:
            run(sig, bkg, cat, year, mass, lifetime,"4g", "ggH", bins=bins)
        
        if len(categories)>1:
            combined_datacard(year,categories,mass,lifetime,"4g","ggH")
