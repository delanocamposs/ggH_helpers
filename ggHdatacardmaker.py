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
    subprocess.run(["mv", f"fit_bkg_{cat}_{year}_fit.root", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"fit_bkg_{cat}_{year}_gen.root", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"fit_sig_{cat}_{year}.root", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"bkg_parameters_{cat}_{year}_fit.json", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"bkg_parameters_{cat}_{year}_gen.json", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"sig_parameters_{cat}_{year}.json", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"data_obs_{cat}_{year}.root", f"{cat}_{year}_{finalstate}_{physics}/"])
#    subprocess.run(["mv", f"datacardInputs_{physics}_{finalstate}_{cat}_{year}.root", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"rate_histos_{cat}_{year}.root", f"{cat}_{year}_{finalstate}_{physics}/"])
#    subprocess.run(["mv", f"datacard_{physics}_{finalstate}_{cat}_{year}.txt", f"{cat}_{year}_{finalstate}_{physics}/"])
    


def main(paths, isMC, trees, var, categories, period, bins, finalstate="4g", physics="ggH", bkg_weight=True,order_fit=3, order_gen=3,lxy1=50, lxy2=50, lumi_scaling=1):
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
             "Run-2":[137620],
             "2022":[34748], 
             "2023":[27245], 
             "2024":[108920], 
             "Run-3":[170857]}


    #combining HL and LH into one category now and renaming to simplify
    cat_dict = {"high" : {"cut" : f"(best_4g_phi1_dxy_m30>{lxy1})&&(best_4g_phi2_dxy_m30>{lxy2})", "file" : ""},
                    "asym" : {"cut" : f"(best_4g_phi1_dxy_m30>{lxy1})&&(best_4g_phi2_dxy_m30<{lxy2})||(best_4g_phi1_dxy_m30<{lxy1})&&(best_4g_phi2_dxy_m30>{lxy2})", "file" : ""}, 
                    "low" : {"cut" : f"(best_4g_phi1_dxy_m30<{lxy1})&&(best_4g_phi2_dxy_m30<{lxy2})", "file" : ""}, 
                    "all" : {"cut" : f"(best_4g_phi1_dxy_m30<{lxy1})&&(best_4g_phi2_dxy_m30<{lxy2})||(best_4g_phi1_dxy_m30<{lxy1})&&(best_4g_phi2_dxy_m30>{lxy2})||(best_4g_phi1_dxy_m30>{lxy1})&&(best_4g_phi2_dxy_m30<{lxy2})||(best_4g_phi1_dxy_m30>{lxy1})&&(best_4g_phi2_dxy_m30>{lxy2})", "file" : ""}}


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
    #not sure about 2022, 2023, 2024, Run-3
    lumi_unc={"2016":[0.012],
              "2017":[0.023],
              "2018":[0.025], 
              "Run-2": [0.016], 
              "2022":[0.012],
              "2023":[0.023],
              "2024":[0.025],
              "Run-3":[0.02]} 

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
    #print(selections)
    
    th1d_files, th1d_filenames, th1d_histos, th1d_histo_obj = ggHtools.sig_bkg_histos(paths, isMC, trees, selections, var, output_names, bins, year, histo_names, bkg_weight, lumi_scaling=lumi_scaling)

    i=0
    for cat in categories:
        bkg_rate=th1d_histo_obj[i][1].Integral()
        dcm_cat_year = DataCardMaker(finalstate, cat, period, lumis[year][0], physics)

        #we need to fit bkg twice because one fit is used to generate the data and one fit is used for combine fit
        #two jsons with parameters: fit and gen. gen=used to generate data. fit=used in combine.
        ggHfitter.fitBKG(f"{th1d_filenames[i]}", f"{th1d_histos[i][1]}", f"fit_bkg_{cat}_{year}_fit.root")
        ggHtools.extract_JSON(f"fit_bkg_{cat}_{year}_fit.root", "w", f"bkg_parameters_{cat}_{year}_fit.json")

        ggHfitter.fitBKG(f"{th1d_filenames[i]}", f"{th1d_histos[i][1]}", f"fit_bkg_{cat}_{year}_gen.root")
        ggHtools.extract_JSON(f"fit_bkg_{cat}_{year}_gen.root", "w", f"bkg_parameters_{cat}_{year}_gen.json")

        ggHfitter.fitSIG(f"{th1d_filenames[i]}", f"{th1d_histos[i][0]}", f"fit_sig_{cat}_{year}.root")
        ggHtools.extract_JSON(f"fit_sig_{cat}_{year}.root", "w", f"sig_parameters_{cat}_{year}.json")

        ggHcardhelper.addDCB(dcm_cat_year, "signal", "mass", f"sig_parameters_{cat}_{year}.json", resolution={f"nuisance_smear_{cat}_{year}":"0.264"})

        #the bernstein used in the combine workspace, hence the order=order_fit
        ggHcardhelper.addBernstein(dcm_cat_year, "background", "mass", f"bkg_parameters_{cat}_{year}_fit.json")

        dcm_cat_year.addSystematic(name=f"nuisance_smear_{cat}_{year}", kind = "param", values=[0.0, 1.0])
        dcm_cat_year.addSystematic(name=f"xsec_unc_{cat}_{year}", kind = "lnN", values={"signal":f"{1+xsec_quad_up}/{1-xsec_quad_down}"})
        dcm_cat_year.addSystematic(name=f"lumi_unc_{cat}_{year}", kind = "lnN", values={"signal":"{}".format(1+lumi_unc[year][0])})
        dcm_cat_year.addSystematic(name=f"PDF_alphas_unc_{cat}_{year}", kind = "lnN", values={"signal":"{}".format(1+PDF_alphas_unc)})
        dcm_cat_year.addSystematic(name=f"bkg_rate_{cat}_{year}", kind = "rateParam", values=[f"{physics}_{finalstate}_{cat}_{year}", "background", "1", "[0,10]"])

        dcm_cat_year.addFixedYieldFromFile(name="background", ID=1, filename=th1d_filenames[i], histoName=th1d_histos[i][1], lumi=False)
        dcm_cat_year.addFixedYieldFromFile(name="signal", ID=0, filename=th1d_filenames[i], histoName=th1d_histos[i][0], lumi=True)

        workspace_file_cat_year = "datacardInputs_"+dcm_cat_year.tag+".root"

        data_cat_year_name = ggHtools.generate_data_hist("fit_bkg_{}_{}_gen.root".format(cat, year), bins_num=bins[0], norm=th1d_histo_obj[i][1].Integral(), output_name="data_obs_{}_{}.root".format(cat, year))

        dcm_cat_year.importBinnedData(data_cat_year_name, "h_pdf__mass", ["mass"])
            
        dcm_cat_year.makeCard()

        print("======================================")
        print("datacard successfully made for {}".format(dcm_cat_year.tag))
        print("=====================================")

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

    main(paths=["ggH_M30_ctau0_ggH4g.root", "EGamma_2018_all_ggH4g.root"], isMC=[1,0], trees=["ggH4g","ggH4g"], var="best_4g_corr_mass_m30", categories=categories,period=args.year, bins=[33, 120, 130])


        
