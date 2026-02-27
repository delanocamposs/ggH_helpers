import ROOT 
import uproot
import awkward as ak
from ROOT import RooFit, RooArgList
import sys
import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.special import comb
from scipy.integrate import simps
from scipy.stats import beta
from ggHfitter import fitBKG
ROOT.gROOT.SetBatch(False)

def define_weightMC(file, tree, BR, xsec):
    '''returns a dataframe that has the event weight defined by the BR, xsec, genWeight, PU weight and sumw'''
    df=ROOT.RDataFrame(tree, file)
    sumw=0.0
    weight_formula = f"(genWeight / {sumw}) * {xsec} * {BR} * Pileup_weight"

    with ROOT.TFile.Open(file) as f:
        runs_tree = f.Get("Runs")
        if runs_tree:
            for entry in runs_tree:
                sumw += entry.genEventSumw
    if sumw == 0:
        print("sum of weights is 0")
        sumw = 1.0
    df=df.Define("event_weight", weight_formula)
    return df



def sig_bkg_histos(files, isMC, trees, mass, lifetime, selections, var, output_names, bins, year, histo_names=[], bkg_weight=True, lumi_scaling=1):
    '''
    files = path to the root data/MC files to construct the histos from
    isMC = is the file MC or not (matters because it gets weighted and scaled correctly if it is)
    trees = the trees inside each file that contain the event data 
    selections = cuts to be applied to the files provided 
    var = the independetn variable in the histograms 
    output_names = the names of the new files which hold the new histos 
    histo_names = new names of the new histos made inside the new root output file 
    bkg_weight = the weight we scale the background histogram by (derived SF from looking at blinded data compared to background)
    bins = the binning of the new histograms

    EXAMPLE:

    sig_bkg_histos(["path/to/MC", "path/to/bkg"], [1, 0], ["ggH4g", "ggH4g"], ["best_4g_phi1_dxy_m{mass}>50", "best_4g_phi2_dxy_m{mass}<50"], "best_4g_cor_mass_m{mass}", ["output1.root", "output2.root"], [["hist1", "hist2"], ["hist3", "hist4"]])
    '''

    sample_sigma=52.143
    BR=1e-4
    outputs=[]
    sumw_dict={}
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #this is for the purposes of testing the effects of various triggers (pT cuts on leading/subleading photons) on the limit
    #im simulating 33GeV double photon trigger here
    double_photon=False
    PT_MIN=33
   # ROOT.gInterpreter.Declare(f"""
   # bool passPtCut(float a,float b,float c,float d) {{
   #     float lead=std::max({{a,b,c,d}});
   #     float sub=-1e9f;
   #     if(a!=lead) sub=a;
   #     if(b!=lead && b>sub) sub=b;
   #     if(c!=lead && c>sub) sub=c;
   #     if(d!=lead && d>sub) sub=d;
   #     return lead>={PT_MIN} && sub>={PT_MIN};
   # }}
   # struct Pair {{ float l; float s; }};
   # Pair getLeadSub(float a,float b,float c,float d) {{
   #     float pts[4]={{a,b,c,d}};
   #     std::sort(pts,pts+4,std::greater<float>());
   #     return {{pts[0], pts[1]}};
   # }}
   # """)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #now this is for testing the effects of the different ID'ing methods on the limit: loose EGM ID versus the custom loose ID we use. 
    EGM=False
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

    

    if len(output_names) != len(selections):
        print(f"ERROR: the number of selections must match the number of desired saved output file names. one saved file for each selection")
        return

    selection_num = len(selections)
    file_num = len(files)
    idx=0

    if len(histo_names) == 0:
        histo_names = [[] for _ in range(selection_num)]
        for i in range(selection_num):
            for j in range(file_num):
                histo_names[i].append(f"hist_{i}_{j}")

    histo_obj = [[] for _ in range(selection_num)]


    #adding this raw background histo so i can generate toys from it. need it prior to any category cuts. dont need sig histo 
    output_raw_bkg = ROOT.TFile(f"rate_histos_raw_bkg_{year}.root", "RECREATE")
    raw_rdf_bkg=ROOT.RDataFrame(trees[1], files[1]) 

    raw_rdf_bkg=raw_rdf_bkg.Filter(f"HLT_passed==1 && Photon_preselection[best_4g_idx1_m{mass}]==1 && Photon_preselection[best_4g_idx2_m{mass}]==1 && Photon_preselection[best_4g_idx3_m{mass}]==1 && Photon_preselection[best_4g_idx4_m{mass}]==1 && best_4g_phi1_dxy_m{mass}>-20 && best_4g_phi2_dxy_m{mass}>-20")

    if double_photon: 
        raw_rdf_bkg=raw_rdf_bkg.Define(f"top2_pT_over_{PT_MIN}", f"passPtCut(best_4g_phi1_gamma1_pt_m{mass},best_4g_phi1_gamma2_pt_m{mass},best_4g_phi2_gamma1_pt_m{mass},best_4g_phi2_gamma2_pt_m{mass})").Define("pTpair", f"getLeadSub(best_4g_phi1_gamma1_pt_m{mass},best_4g_phi1_gamma2_pt_m{mass},best_4g_phi2_gamma1_pt_m{mass},best_4g_phi2_gamma2_pt_m{mass})").Define("LeadpT", "pTpair.l").Define("subLeadpT", "pTpair.s")
        raw_rdf_bkg = raw_rdf_bkg.Filter(f"top2_pT_over_{PT_MIN}==1")

    raw_hist_bkg=raw_rdf_bkg.Histo1D(("bkg_hist", f"bkg_hist;{var};Events", bins[0], bins[1], bins[2]),f"{var}")
    raw_hist_bkg.Write()
    output_raw_bkg.Close()
    
    for i in range(len(files)):
        filepath_i = files[i]
        file_i = ROOT.TFile.Open(files[i])
        tree_i = trees[i]

        if isMC[i]:
            sumw=0.0
            with ROOT.TFile.Open(filepath_i) as f:
                runs_tree = f.Get("Runs")
                if runs_tree:
                    for entry in runs_tree:
                        sumw += entry.genEventSumw
            if sumw == 0:
                print("sum of weights is 0")
                sumw = 1.0

            sumw_dict[filepath_i]=sumw


    #these dictionary values are scales for the background distribution according to what year is being processed. 
    #this needs to be done since hte 0.05124.... factor was derived from comparing sidebands of 2018 preselected vs ID data. 
    #if we want to interpolate to any other year having derived the initial scale factor from 2018, all years need their own additional weight normalized by 2018 lumi
    bkg_factor={"2016preVFP":[1/2*(36310/59830)], 
                "2016postVFP":[1/2*(36210/59830)],
                "2017":[41480/59830], 
                "2018":[1], 
                "Run2":[137620/59830], 
                "2022":[34748/59830], 
                "2023":[27245/59830], 
                "2024":[108920/59830], 
                "Run3":[170857/59830]}

    for j in range(len(selections)):
        selection_j = selections[j]
        output_file_j = ROOT.TFile(f"{output_names[j]}", "RECREATE")
        outputs.append(output_file_j)

        for k in range(len(files)):
            rdf_j_k=ROOT.RDataFrame(trees[k], files[k])
    
            #common filters for both signal and backgroound
            rdf_j_k = rdf_j_k.Filter(f"HLT_passed==1 && best_4g_phi1_dxy_m{mass}>-20 && best_4g_phi2_dxy_m{mass}>-20").Filter(selection_j)

            #if we want to test the double photon trigger
            if double_photon:
                rdf_j_k=rdf_j_k.Define(f"top2_pT_over_{PT_MIN}", f"passPtCut(best_4g_phi1_gamma1_pt_m{mass},best_4g_phi1_gamma2_pt_m{mass},best_4g_phi2_gamma1_pt_m{mass},best_4g_phi2_gamma2_pt_m{mass})")
                rdf_j_k=rdf_j_k.Define("pTpair", f"getLeadSub(best_4g_phi1_gamma1_pt_m{mass},best_4g_phi1_gamma2_pt_m{mass},best_4g_phi2_gamma1_pt_m{mass},best_4g_phi2_gamma2_pt_m{mass})").Define("LeadpT", "pTpair.l").Define("subLeadpT", "pTpair.s")
                rdf_j_k = rdf_j_k.Filter(f"top2_pT_over_{PT_MIN}==1")
    

            #this if/else statement then filters beyond the common filters depending on if its signal or background, then writes and saves the histos
            if isMC[k]:
                weight_formula_k = f"(genWeight / {sumw_dict[files[k]]}) * {sample_sigma} * {BR} * Pileup_weight"
                rdf_j_k = rdf_j_k.Filter(f"best_4g_ID_m{mass}==1&&best_4g_passBitMap_loose_iso_m{mass}==1").Define("event_weight", weight_formula_k)
                #rdf_j_k = rdf_j_k.Filter("best_4g_ID_m{mass}==1&&((Photon_isScEtaEB[best_4g_idx1_m{mass}]==1 && Photon_corrIso_m{mass}[best_4g_idx1_m{mass}]<0.1)||(Photon_isScEtaEE[best_4g_idx1_m{mass}]==1 && Photon_corrIso_m{mass}[best_4g_idx1_m{mass}]<0.1)) && ((Photon_isScEtaEB[best_4g_idx2_m{mass}]==1 && Photon_corrIso_m{mass}[best_4g_idx2_m{mass}]<0.1)||(Photon_isScEtaEE[best_4g_idx2_m{mass}]==1 && Photon_corrIso_m{mass}[best_4g_idx2_m{mass}]<0.1)) && ((Photon_isScEtaEB[best_4g_idx3_m{mass}]==1 && Photon_corrIso_m{mass}[best_4g_idx3_m{mass}]<0.1)||(Photon_isScEtaEE[best_4g_idx3_m{mass}]==1 && Photon_corrIso_m{mass}[best_4g_idx3_m{mass}]<0.1)) && ((Photon_isScEtaEB[best_4g_idx4_m{mass}]==1 && Photon_corrIso_m{mass}[best_4g_idx4_m{mass}]<0.1)||(Photon_isScEtaEE[best_4g_idx4_m{mass}]==1 && Photon_corrIso_m{mass}[best_4g_idx4_m{mass}]<0.1))").Define("event_weight", weight_formula_k)

                if EGM:
                    EGM_ID = " && ".join(f"((Photon_isScEtaEB[best_4g_idx{i}_m{mass}]==1 && Photon_hoe[best_4g_idx{i}_m{mass}]<0.04596 && Photon_sieie[best_4g_idx{i}_m{mass}]<0.0106)||(Photon_isScEtaEE[best_4g_idx{i}_m{mass}]==1 && Photon_hoe[best_4g_idx{i}_m{mass}]<0.0590 && Photon_sieie[best_4g_idx{i}_m{mass}]<0.0272))" for i in range(1, 5))
                    rdf_j_k = rdf_j_k.Filter(f"{EGM_ID}")

                hist_j_k = rdf_j_k.Histo1D((f"{histo_names[j][k]}", f"{j}_{k};{var};Events", bins[0], bins[1], bins[2]), f"{var}", "event_weight")

                if double_photon:
                    hist_test1 = rdf_j_k.Histo1D(("leadsig", f"leadsig;LeadpT;Events", 100, 0, 160), "LeadpT")
                    hist_test2 = rdf_j_k.Histo1D(("subleadsig", f"subleadsig;subLeadpT;Events", 100, 0, 160), "subLeadpT")
                    hist_test1.Write()
                    hist_test2.Write()

                print("event weight: ", weight_formula_k)
                hist_j_k.Scale(lumi_scaling)
                hist_j_k.Write()
                histo_obj[j].append(hist_j_k)

            else:
                rdf_j_k = rdf_j_k.Filter(f"Photon_preselection[best_4g_idx1_m{mass}]==1 && Photon_preselection[best_4g_idx2_m{mass}]==1 && Photon_preselection[best_4g_idx3_m{mass}]==1 && Photon_preselection[best_4g_idx4_m{mass}]==1")

                #if EGM:
                #    EGM_ID = " && ".join(f"((Photon_isScEtaEB[best_4g_idx{i}_m{mass}]==1 && Photon_hoe[best_4g_idx{i}_m{mass}]<0.04596 && Photon_sieie[best_4g_idx{i}_m{mass}]<0.0106)||(Photon_isScEtaEE[best_4g_idx{i}_m{mass}]==1 && Photon_hoe[best_4g_idx{i}_m{mass}]<0.0590 && Photon_sieie[best_4g_idx{i}_m{mass}]<0.0272))" for i in range(1, 5))
                #    rdf_j_k = rdf_j_k.Filter(f"{EGM_ID}")


                if double_photon:
                    hist_test3 = rdf_j_k.Histo1D(("leadbkg", f"leadbkg;LeadpT;Events", 100, 0, 160), "LeadpT")
                    hist_test4 = rdf_j_k.Histo1D(("subleadbkg", f"subleadbkg;subLeadpT;Events", 100, 0, 160), "subLeadpT")
                    hist_test3.Write()
                    hist_test4.Write()

                hist_j_k = rdf_j_k.Histo1D((f"{histo_names[j][k]}", f"{j}_{k};{var};Events", bins[0], bins[1], bins[2]), f"{var}")

                if bkg_weight:
                    ##this weight comes from scaling PRESELECTED 4g to the sidebands of FULLY ID'd 4g data in 2018 EGamma data.
                    if EGM:    
                        bkg_w = (1/7099)*(bkg_factor[year][0])
                    else:
                        bkg_w = (15/7099)*(bkg_factor[year][0])
                            
                    hist_j_k.Scale(lumi_scaling*bkg_w)
        
                else:
                    hist_j_k.Scale(lumi_scaling*1*bkg_factor[year][0])

                hist_j_k.Write()
                histo_obj[j].append(hist_j_k)
         
        output_file_j.Close()

    return outputs, output_names, histo_names, histo_obj


def extract_JSON(root_filename, workspace_name, json_filename):

    '''
    creates a JSON file in the format the datcard needs from a root file, the name of the 
    RooWorkspace in the file.
    '''

    f = ROOT.TFile(root_filename)
    ws = f.Get(workspace_name)

    if not ws:
        print(f"Error: Workspace '{workspace_name}' not found in {root_filename}")
        return

    params = {}
    all_vars = ws.allVars()
    it = all_vars.createIterator()
    var = it.Next()

    while var:
        if var.InheritsFrom("RooRealVar"):
            params[var.GetName()] = {"value": var.getVal(), "error": var.getError()}
        var = it.Next()

    with open(json_filename, "w") as json_file:
        json.dump(params, json_file, indent=4)

    f.Close()


def generate_data_hist(file, bins_num, norm, output_name):

    f=ROOT.TFile.Open(file)
    w=f.Get("w")
    pdf =w.pdf("model")
    x=w.var("mass")
    
    h_pdf =pdf.createHistogram("h_pdf", x, ROOT.RooFit.Binning(bins_num))
    h_pdf1 =pdf.createHistogram("h_pdf1", x, ROOT.RooFit.Binning(bins_num))

    h_pdf.Scale(norm)
    h_pdf1.Scale(norm)

    total_expected = float(h_pdf.Integral())
    prob_all_zero = float(np.exp(-total_expected))
    print(f"[generate_data_hist] input={file}, output={output_name}, "f"expected_total={total_expected:.6g}, P(all bins are 0)={prob_all_zero:.6g}")
    if total_expected < 1:
        print("[generate_data_hist] WARNING: expected total events is < 1. A fully empty toy histogram is likely and can be statistically consistent.")

    output_file = ROOT.TFile(output_name, "RECREATE")
    h_pdf1.Write()

    for i in range(1, h_pdf.GetNbinsX()+1):
        density=h_pdf.GetBinContent(i)
        bw=h_pdf.GetXaxis().GetBinWidth(i)
        h_pdf.SetBinContent(i, np.random.poisson(density))

    print(f"[generate_data_hist] toy_total={h_pdf.Integral():.6g}")

    h_pdf.Write()
    output_file.Close()

    return output_name


def generate_bkg_toy(file, bins_num, norm, output_name, toys=1):
    f=ROOT.TFile.Open(file)
    h=f.Get("hist")
    output_file = ROOT.TFile(output_name, "RECREATE")

    for i in range(1, h.GetNbinsX()+1):
        density=h.GetBinContent(i)
        bw=h.GetXaxis().GetBinWidth(i)
        h.SetBinContent(i, (1/100)*np.random.poisson(100*density))
    h.Write()
    output_file.Close()

    return output_name

def clopper_pearson(X, n, alpha=0.05):
    "by default i am doing 95% CL, change alpha to change this: CL=1-alpha"
    lower=beta.ppf(alpha/2, X,n-X+1)
    upper=beta.ppf(1-(alpha/2),X+1,n-X)
    return lower,upper
    

