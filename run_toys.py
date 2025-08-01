import ggHdatacardmaker
import ROOT 
import subprocess

def main(year, toys, order_gen, order_fit=[]):
    ROOT.gROOT.SetBatch(True)

    for fit_order in order_fit:
        for t in range(1, toys+1):
            period=f"{year}_toy-{t}_orderfit-{fit_order}_ordergen-{order_gen}"
            ggHdatacardmaker.main(paths=["ggH_M30_ctau0_ggH4g.root", f"EGamma_{year}_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["lowlow"],period=period, bins=[90, 110, 140], order_fit=fit_order, order_gen=order_gen, lxy1=50, lxy2=50, lumi_scaling=10)
            ggHdatacardmaker.main(paths=["ggH_M30_ctau0_ggH4g.root", f"EGamma_{year}_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["lowhigh"],period=period, bins=[90, 110, 140], order_fit=fit_order, order_gen=order_gen, lxy1=50, lxy2=50, lumi_scaling=10)
            ggHdatacardmaker.main(paths=["ggH_M30_ctau0_ggH4g.root", f"EGamma_{year}_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["highlow"],period=period, bins=[90, 110, 140], order_fit=fit_order, order_gen=order_gen, lxy1=50, lxy2=50, lumi_scaling=10)
            ggHdatacardmaker.main(paths=["ggH_M30_ctau0_ggH4g.root", f"EGamma_{year}_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=["highhigh"],period=period, bins=[90, 110, 140], order_fit=fit_order, order_gen=order_gen, lxy1=50, lxy2=50, lumi_scaling=10)

            with open("datacard_"+period+"_combined.txt", "w") as f:
                subprocess.run(["combineCards.py", f"datacard_ggH_4g_lowlow_"+period+".txt", f"datacard_ggH_4g_highlow_"+period+".txt", f"datacard_ggH_4g_lowhigh_"+period+".txt", f"datacard_ggH_4g_highhigh_"+period+".txt"], stdout=f)

            subprocess.run(["text2workspace.py", "datacard_"+period+"_combined.txt", "datacard_"+period+"_combined.root"])
            subprocess.run(["combine", "-M", "FitDiagnostics", "--saveShapes", "datacard_"+period+"_combined.root", "--saveWithUncertainties", "--saveNormalizations", "--saveWorkspace", "--cminDefaultMinimizerStrategy","0","--setParameterRanges", "r=-10,10", "-m", "125"])
            subprocess.run(["mv", "fitDiagnosticsTest.root", "fitDiagnosticsTest_"+period+".root"])
            


def fill_txt(year, toys, order_gen, fit_order):
    r_vals=[]
    r_vals_unc_up=[]
    r_vals_unc_down=[]

    for i in range(1, toys+1):
        period=f"{year}_toy-{i}_orderfit-{fit_order}_ordergen-{order_gen}"
        file=ROOT.TFile.Open(f"fitDiagnosticsTest_"+period+".root")
        fit_s=file.Get("fit_s")
        if fit_s== None:
            continue
        r=fit_s.floatParsFinal().find("r")
        r_val=r.getVal()
        r_unc_up=r.getErrorHi()
        r_unc_down=r.getErrorLo()
        r_vals.append(r_val)
        r_vals_unc_up.append(r_unc_up)
        r_vals_unc_down.append(r_unc_down)
        print(f"extracted r val for toy {i}")

    with open(f"r_vals_{year}_orderfit-{fit_order}_ordergen-{order_gen}.txt", "w") as f:
        for r_val, r_val_unc_up, r_val_unc_down in zip(r_vals, r_vals_unc_up, r_vals_unc_down):
            f.write(f"{r_val}, {r_val_unc_up}, {r_val_unc_down} \n")


if __name__ == "__main__":
    main("2018", 100, 3, [2, 3])
    fill_txt("2018",100 ,3, 2)
    fill_txt("2018",100 ,3, 3)



