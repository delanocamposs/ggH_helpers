import ggHtools 
import ggHcardhelper 
import ggHfitter
import ggHsmear
from DataCardMaker import DataCardMaker
import json 
import os 
import ROOT
import argparse, subprocess
import numpy as np
import shutil
import glob

def cleanup(year, finalstate, physics, cat):
    subprocess.run(["mkdir", f"{cat}_{year}_{finalstate}_{physics}"])
    subprocess.run(["mv", f"fit_bkg_{cat}_{year}.root", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"fit_sig_{cat}_{year}.root", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"bkg_parameters_{cat}_{year}.json", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"sig_parameters_{cat}_{year}.json", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"data_obs_{cat}_{year}.root", f"{cat}_{year}_{finalstate}_{physics}/"])
#    subprocess.run(["mv", f"datacardInputs_{physics}_{finalstate}_{cat}_{year}.root", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"rate_histos_{cat}_{year}.root", f"{cat}_{year}_{finalstate}_{physics}/"])
#    subprocess.run(["mv", f"datacard_{physics}_{finalstate}_{cat}_{year}.txt", f"{cat}_{year}_{finalstate}_{physics}/"])
    


def main(paths, isMC, trees, var, categories, period, bins, finalstate="4g", physics="ggH", bkg_weight=True,order=3, lxy1=50, lxy2=50):
    ROOT.gROOT.SetBatch(True)
    year=period
    cat_sum=0
    for cat in categories:
        cat_sum+=len(cat)
    extra_str="="*cat_sum

    print("============================================================"+extra_str)
    print("processing datacard for year {} in categories: {}".format(year, categories))
    print("============================================================"+extra_str)

    ##luminosities in pb^-1 from: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2
    lumis = {"2016":[36310],
             "2017":[41480], 
             "2018":[59830], 
             "Run-2":[101310]} ##this is the sum of 2017 and 2018 only because we do not use 2016


    cat_dict = {"highhigh" : {"cut" : f"(best_4g_phi1_dxy_m30>{lxy1})&&(best_4g_phi2_dxy_m30>{lxy2})", "file" : ""},
                    "highlow" : {"cut" : f"(best_4g_phi1_dxy_m30>{lxy1})&&(best_4g_phi2_dxy_m30<{lxy2})", "file" : ""}, 
                    "lowhigh" : {"cut" : f"(best_4g_phi1_dxy_m30<{lxy1})&&(best_4g_phi2_dxy_m30>{lxy2})", "file" : ""}, 
                    "lowlow" : {"cut" : f"(best_4g_phi1_dxy_m30<{lxy1})&&(best_4g_phi2_dxy_m30<{lxy2})", "file" : ""}, 
                    "none" : {"cut" : "(HLT_passed == 1)", "file" : ""}}

#    cat_dict = {"highhigh" : {"cut" : f"(best_4g_phi1_dxy_m30>{lxy1})&&(best_4g_phi2_dxy_m30>{lxy2}) &&HLT_passed==1&&best_4g_ID_m30==1&&best_4g_phi1_mass_m30>14&&best_4g_phi2_mass_m30>14 ", "file" : ""},
#                    "highlow" : {"cut" : f"(best_4g_phi1_dxy_m30>{lxy1})&&(best_4g_phi2_dxy_m30<{lxy2}) &&HLT_passed==1&&best_4g_ID_m30==1&&best_4g_phi1_mass_m30>14&&best_4g_phi2_mass_m30>14 ", "file" : ""}, 
#                    "lowhigh" : {"cut" : f"(best_4g_phi1_dxy_m30<{lxy1})&&(best_4g_phi2_dxy_m30>{lxy2}) &&HLT_passed==1&&best_4g_ID_m30==1&&best_4g_phi1_mass_m30>14&&best_4g_phi2_mass_m30>14 ", "file" : ""}, 
#                    "lowlow" : {"cut" : f"(best_4g_phi1_dxy_m30<{lxy1})&&(best_4g_phi2_dxy_m30<{lxy2}) &&HLT_passed==1&&best_4g_ID_m30==1&&best_4g_phi1_mass_m30>14&&best_4g_phi2_mass_m30>14 ", "file" : ""}, 
#                    "none" : {"cut" : "(HLT_passed == 1)", "file" : ""}}


    output_names = ["rate_histos_{}_{}.root".format(cat, year) for cat in categories]
    histo_names = []

    for i in range(len(categories)):
        pair = ["hist{}".format(2*i+1), "hist{}".format(2*i+2)]
        histo_names.append(pair)

    #######IMPORTANT#######
    #since the cross section we use to scale the MC is both ggH+VBF combined, these xsec and lumi unc are derived from adding the unc from 
    #each process in quadrature    

    #obtained xsec and PDF from:https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#gluon_gluon_Fusion_Process
    xsec_unc={"ggH":[0.046, -0.067],
              "VBF":[0.004, -0.003]}

    PDF_alphas_unc={"ggH":[0.032, 0.032], 
                    "VBF":[0.021, 0.021]}

    #lumi unc obtained from: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2
    lumi_unc={"2016":[0.012],
              "2017":[0.023],
              "2018":[0.025], 
              "Run-2": [0.03397]} #Run-2 lumi unc obtained from adding in quadrature 2017 and 2018 unc

    xsec_quad_up = np.sqrt((xsec_unc["ggH"][0])**2+(xsec_unc["VBF"][0])**2)
    xsec_quad_down = np.sqrt((xsec_unc["ggH"][1])**2+(xsec_unc["VBF"][1])**2)
    PDF_alphas_unc = np.sqrt((PDF_alphas_unc["ggH"][0])**2 + (PDF_alphas_unc["VBF"][0])**2)
        
    signal_samples={}
    for j in isMC:
        if isMC[j]==1:
            signal_samples[paths[j]]=trees[j]

    selections = []

    for cat in cat_dict:    
        if cat in categories:
            selections.append(cat_dict[cat]["cut"])

    selections.reverse()
    print(selections)

    th1d_files, th1d_filenames, th1d_histos, th1d_histo_obj = ggHtools.sig_bkg_histos(paths, isMC, trees, selections, var, output_names, bins, year, histo_names, bkg_weight)

    i=0
    for cat in categories:
        bkg_rate=th1d_histo_obj[i][1].Integral()
        dcm_cat_year = DataCardMaker(finalstate, cat, period, lumis["Run-2"][0], physics)

        ggHfitter.fitBKG("{}".format(th1d_filenames[i]), "{}".format(th1d_histos[i][1]), "fit_bkg_{}_{}.root".format(cat, year), order=order)
        ggHtools.extract_JSON("fit_bkg_{}_{}.root".format(cat, year), "w", "bkg_parameters_{}_{}.json".format(cat, year))
        
        ggHfitter.fitSIG("{}".format(th1d_filenames[i]), "{}".format(th1d_histos[i][0]), "fit_sig_{}_{}.root".format(cat, year))
        ggHtools.extract_JSON("fit_sig_{}_{}.root".format(cat, year), "w", "sig_parameters_{}_{}.json".format(cat, year))
        dcm_cat_year.addFixedYieldFromFile(name="signal", ID=0, filename=th1d_filenames[i], histoName=th1d_histos[i][0], lumi=True)
        print(f"signal yield for {cat} using file {th1d_filenames[i]} with histogram {th1d_histos[i][0]}")

        ggHcardhelper.addDCB(dcm_cat_year, "signal", "mass", "sig_parameters_{}_{}.json".format(cat, year), resolution={"nuisance_smear_{}_{}".format(cat, year):"0.264"})
        ggHcardhelper.addBernstein(dcm_cat_year, "background", "mass", "bkg_parameters_{}_{}.json".format(cat, year))

        dcm_cat_year.addSystematic(name="nuisance_smear_{}_{}".format(cat, year), kind = "param", values=[0.0, 1.0])

        dcm_cat_year.addSystematic(name="xsec_unc_{}_{}".format(cat, year), kind = "lnN", values={"signal":"{}/{}".format(1+xsec_quad_up,1-xsec_quad_down)})
        
        dcm_cat_year.addSystematic(name="lumi_unc_{}_{}".format(cat, year), kind = "lnN", values={"signal":"{}".format(1+lumi_unc["Run-2"][0])})

        dcm_cat_year.addSystematic(name="PDF_alphas_unc_{}_{}".format(cat, year), kind = "lnN", values={"signal":"{}".format(1+PDF_alphas_unc)})
        dcm_cat_year.addSystematic(name=f"bkg_rate_{cat}_{year}", kind = "rateParam", values=[f"{physics}_{finalstate}_{cat}_{year}", "background", "1", "[0,10]"])

        dcm_cat_year.addFixedYieldFromFile(name="background", ID=1, filename=th1d_filenames[i], histoName=th1d_histos[i][1], lumi=False)

        workspace_file_cat_year = "datacardInputs_"+dcm_cat_year.tag+".root"

        data_cat_year_name = ggHtools.generate_data_hist("fit_bkg_{}_{}.root".format(cat, year), bins_num=bins[0], norm=th1d_histo_obj[i][1].Integral(), output_name="data_obs_{}_{}.root".format(cat, year))
        print(th1d_histo_obj[i][1])

        dcm_cat_year.importBinnedData(data_cat_year_name, "h_pdf__mass", ["mass"])
            
        dcm_cat_year.makeCard()

        print("-------------------------------------------------")
        print("datacard successfully made for {}".format(dcm_cat_year.tag))
        print("-------------------------------------------------")

        cleanup(year, finalstate, physics, cat)

        i+=1


if __name__=="__main__":
    parser = argparse.ArgumentParser("Datacard Processing for year, category")
    parser.add_argument("-y","--year", type=str)
    parser.add_argument("-c1","--cat1", dest="c1", type=str)
    parser.add_argument("-c2","--cat2", dest="c2", type=str)
    parser.add_argument("-c3","--cat3", dest="c3", type=str)
    parser.add_argument("-c4","--cat4", dest="c4", type=str)

    args = parser.parse_args()

    categories = []
    for c in ["c1", "c2", "c3", "c4"]:
        argu = getattr(args, c)
        if argu is not None:
            categories.append(argu)
    print(categories)

    main(paths=["ggH_M30_ctau0_ggH4g.root", "EGamma_2018_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=categories,period=args.year, bins=[90, 110, 140], order=3, lxy1=50, lxy2=50)


        
