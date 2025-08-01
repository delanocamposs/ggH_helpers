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
from ggHfitter import fitBKG


def sig_bkg_histos(files, isMC, trees, selections, var, output_names, bins, year, histo_names=[], bkg_weight=True, lumi_scaling=1):
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

    sig_bkg_histos(["path/to/MC", "path/to/bkg"], [1, 0], ["ggH4g", "ggH4g"], ["best_4g_phi1_dxy_m30>50", "best_4g_phi2_dxy_m30<50"], "best_4g_cor_mass_m30", ["output1.root", "output2.root"], [["hist1", "hist2"], ["hist3", "hist4"]])
    '''
    scales_dict={"2016":[36310/59830],"2017":[41480/59830],"2018":[1],"Run-2":[101310/59830]}


    sample_sigma=52.143
    BR=1e-4
    outputs=[]
    sumw_dict={}

    if len(output_names) != len(selections):
        print("ERROR: the number of selections must match the number of desired saved output file names. one saved file for each selection")
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
    raw_rdf_bkg=raw_rdf_bkg.Filter("HLT_passed==1 && Photon_preselection[best_4g_idx1_m30]==1 && Photon_preselection[best_4g_idx2_m30]==1 && Photon_preselection[best_4g_idx3_m30]==1 && Photon_preselection[best_4g_idx4_m30]==1 && best_4g_phi1_dxy_m30>-20 && best_4g_phi2_dxy_m30>-20")
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

    #print(sumw_dict)

    for j in range(len(selections)):
        selection_j = selections[j]
        output_file_j = ROOT.TFile(f"{output_names[j]}", "RECREATE")
        outputs.append(output_file_j)

        for k in range(len(files)):
            rdf_j_k = ROOT.RDataFrame(trees[k], files[k])

            if isMC[k]:
                weight_formula_k = f"(genWeight / {sumw_dict[files[k]]}) * {sample_sigma} * {BR} * Pileup_weight"
                rdf_j_k = rdf_j_k.Filter("HLT_passed==1 && best_4g_ID_m30==1 && best_4g_phi1_dxy_m30>-20 && best_4g_phi2_dxy_m30>-20").Filter(selection_j).Define("event_weight", weight_formula_k)

                hist_j_k = rdf_j_k.Histo1D((f"{histo_names[j][k]}", f"{j}_{k};{var};Events", bins[0], bins[1], bins[2]), f"{var}", "event_weight")
                print("event weight: ", weight_formula_k)
                hist_j_k.Scale(lumi_scaling)
                hist_j_k.Write()
                histo_obj[j].append(hist_j_k)

            else:
                rdf_j_k = rdf_j_k.Filter("HLT_passed==1 && Photon_preselection[best_4g_idx1_m30]==1 && Photon_preselection[best_4g_idx2_m30]==1 && Photon_preselection[best_4g_idx3_m30]==1 && Photon_preselection[best_4g_idx4_m30]==1 && best_4g_phi1_dxy_m30>-20 && best_4g_phi2_dxy_m30>-20").Filter(selection_j)
                #rdf_j_k = rdf_j_k.Filter(selection_j)
                hist_j_k = rdf_j_k.Histo1D((f"{histo_names[j][k]}", f"{j}_{k};{var};Events", bins[0], bins[1], bins[2]), f"{var}")

                if bkg_weight:
                    ##this weight comes from scaling PRESELECTED 4g to the sidebands of FULLY ID'd 4g data in 2018 EGamma data.
                    #bkg_w = 0.9307568438
                    bkg_w = 0.051249577845322525*(101310/59830)
                    hist_j_k.Scale(lumi_scaling*bkg_w)
        
                else:
                    hist_j_k.Scale(lumi_scaling*1)

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


def generate_data_hist(file, bins_num, norm, output_name, lumi_scaling=1):

    f=ROOT.TFile.Open(file)
    w=f.Get("w")
    pdf =w.pdf("model")
    x=w.var("mass")
    
    h_pdf =pdf.createHistogram("h_pdf", x, ROOT.RooFit.Binning(bins_num))
    h_pdf.Scale(lumi_scaling*norm)

    output_file = ROOT.TFile(output_name, "RECREATE")

    for i in range(1, h_pdf.GetNbinsX()+1):
        density=h_pdf.GetBinContent(i)
        bw=h_pdf.GetXaxis().GetBinWidth(i)
        h_pdf.SetBinContent(i, np.random.poisson(density))

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
























#if __name__=="__main__":
#    gen_data(bkg_file="EGamma_2018_all_1bad.root", tree="ggH4g_1bad", var="best_4g_corr_mass_m30", cut="(HLT_passed == 1)&&(best_4g_ID_m30 == 1)&&(best_4g_phi1_mass_m30 > 14)&&(best_4g_phi2_mass_m30>14)&&(best_4g_phi1_dxy_m30>-20)&&(best_4g_phi2_dxy_m30>-20)", output_names=[], categories=["lowlow"])
#    gen_data(bkg_file="EGamma_2018_all_1bad.root", tree="ggH4g_1bad", var="best_4g_corr_mass_m30", cut="(HLT_passed == 1)&&(best_4g_ID_m30 == 1)&&(best_4g_phi1_mass_m30 > 14)&&(best_4g_phi2_mass_m30>14)", output_names=[], categories=["lowlow"])
