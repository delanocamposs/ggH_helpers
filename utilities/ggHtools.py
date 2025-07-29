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


#def scale_MC(source_file, hist_file, hist, year, BR=1e-4,sigma=56):
#    lumis = {"2018":"59830"}
#    rdf=ROOT.RDataFrame("Events", file)

#    f = ROOT.TFile(file)
#    t=f.Get("Runs")
#    sumw=0.0

#    for event in t:
#        sumw+=event.genEventSumw     

#    if sumw==0.0:
#        print("sumw is 0")
#        sumw=1.0

#    weight_formula=f"(genWeight/sumw)*{sigma}*{BR}*Pileup_weight"
#    rdf=rdf.Define("event_weight", weight_formula)
#    f.Close()
    


def sig_bkg_histos(files, isMC, trees, selections, var, output_names, bins, year, histo_names=[], bkg_weight=True):
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

#    if idx != selection_num * file_num:
#        print('ERROR: the histo_names should be a nested array with as many nested arrays as there are selections. eg: [["", ""], ["", ""]]')
#        return

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
                #rdf_j_k = rdf_j_k.Filter(selection_j).Define("event_weight", weight_formula_k)
                hist_j_k = rdf_j_k.Histo1D((f"{histo_names[j][k]}", f"{j}_{k};{var};Events", bins[0], bins[1], bins[2]), f"{var}", "event_weight")
                print("event weight: ", weight_formula_k)
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
                    hist_j_k.Scale(bkg_w)
        
                else:
                    hist_j_k.Scale(1)

                hist_j_k.Write()
                histo_obj[j].append(hist_j_k)
        
         
        output_file_j.Close()

    return outputs, output_names, histo_names, histo_obj



def gen_data(bkg_file, tree, var, cut, output_names,categories=[], lxy1=0, lxy2=0):
    """
    cuts=cuts prior to any category definition. ID, preselection, etc. this applies to all categories!!!!!!
    categories=[] is array of category names for which we want to generate data. eg: ["lowlow", "lowhigh"]

    """
    ROOT.gRandom.SetSeed(0)

    if len(categories)==0:
        print("")
        return

    rdf=ROOT.RDataFrame(tree, bkg_file)
    rdf=rdf.Filter(cut) 
    rdf.Snapshot("baseline", "baseline.root")

    with uproot.open("baseline.root") as f:
        baseline_data=f["baseline"].arrays(library="pd")

    N_exp=len(baseline_data)
    N_toy=ROOT.gRandom.Poisson(N_exp)
    toy_data=baseline_data.sample(n=N_toy, replace=True).reset_index(drop=True)
    print(toy_data)
    print(N_toy)
    print(N_exp)
    
    def assign_category(row, lxy1=lxy1, lxy2=lxy2):
        if row["best_4g_phi1_dxy_m30"]>lxy1 and row["best_4g_phi2_dxy_m30"]>lxy2:
            return 0
        elif row["best_4g_phi1_dxy_m30"]>lxy1 and row["best_4g_phi2_dxy_m30"]<lxy2:
            return 1
        elif row["best_4g_phi1_dxy_m30"]<lxy1 and row["best_4g_phi2_dxy_m30"]>lxy2:
            return 2 
        elif row["best_4g_phi1_dxy_m30"]<lxy1 and row["best_4g_phi2_dxy_m30"]<lxy2:
            return 3

    toy_data["category"]=toy_data.apply(assign_category, axis=1)
    hists=[ROOT.TH1D(f"h_cat{i+1}", f"category {i+1}", 30, 110, 140) for i in range(4)]
    for _,row in toy_data.iterrows():
        cat=row["category"]
        mass=row["best_4g_corr_mass_m30"]    
        #print("row: ", row)
        #print("_: ", _)
        hists[cat].Fill(mass)
        print("mass: ", mass)

    for i in range(len(hists)):
        output_i=ROOT.TFile.Open(f"data_obs_{i}.root", "RECREATE")
        hists[i].Write()
        output_i.Close()
        
    
    
    #colors=[2, 4, 6, 8]
    #c=ROOT.TCanvas("", "", 800, 800)
    #for i,h in enumerate(hists):
    #    h.SetLineColor(colors[i])
    #    if i==0:
    #        h.Draw("HIST")
    #        h.SetMaximum(30)
    #    else:
    #        h.Draw("HIST, SAME")
    #c.SaveAs("cat_data_histos.png")
    
    
    
    
    
    
    


def makePoissonEvents(file, hist, visual=False, cat=None, year=None, new_histo_name=None, data_obs=True):

    '''
    given a file and histogram as arguments, it makes a similar histogram 
    except with poisson distributed event count per bin as a randomizer.

    returnsss: 
    the output file TFile ROOT object, and the name of the output file as a str
    '''

    file = ROOT.TFile(file)
    histo = file.Get(hist)
    if new_histo_name is None:
        clone = histo.Clone("histo")
        print("here 1")
    else: 
        clone = histo.Clone(new_histo_name)
        print("here 2")
   
    random = ROOT.TRandom3() 
    for i in range(1, clone.GetNbinsX()+1):
        height = histo.GetBinContent(i)
        print(height)
        if height > 0:
            clone.SetBinContent(i, random.Poisson(height))
        else:
            clone.SetBinContent(i, 0)

    if cat and year is not None:
        if data_obs:
            output = ROOT.TFile("data_obs_{}_{}.root".format(cat, year), "RECREATE")    
            name = "data_obs_{}_{}.root".format(cat, year)
            output.cd()
            clone.Write()
        else:
            output = ROOT.TFile("poisson_gen_{}_{}.root".format(cat, year), "RECREATE")    
            name = "poisson_gen_{}_{}.root".format(cat, year)
            output.cd()
            clone.Write()
    else:
        if data_obs:
            output = ROOT.TFile("data_obs.root", "RECREATE") 
            name = "data_obs.root"
            output.cd()
            clone.Write()
        else:
            output = ROOT.TFile("poisson_gen.root", "RECREATE") 
            name = "poisson_gen.root"
            output.cd()
            clone.Write()

    if visual:
        c = ROOT.TCanvas("c", "", 800, 600)
        clone.SetLineColor(ROOT.kRed)
        histo.SetLineColor(ROOT.kBlue)
        clone.Draw("HIST")
        histo.Draw("HIST, SAME")
        c.Update()
        input("press enter tyo close canvas")
    

    return output, name

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


def make_hist_from_shape(file_name, workspace_name, var_name, pdf_name, N, output_name, bins=[30, 110, 140], histo_name = "data_obs", see=False):
    """
    takes a RooWorkspace root file with a shape pdf in it and generates a histogram from the shape without allowing any bins to have 
    0 event content. The histogram gets saved to a new file named by the user.

    workspace_file = RooWorkspace root file
    workspace_name = name of RooWorkspace inside file 
    var = observable of plotter in the workspace 
    pdf_name = name of the probability density function inside the RooWorkspace that you are trying to generate the histo from 
    N = number of expected events in the new histogram 
    output_name = the file to save the new histogram to 
    histo_name = name of the histo in the new file 
    n_bins = number of bins of the new histo

    returns:
    output TFile object, histogram object generated from pdf, name of the output file
    
    """

    #var = w.var(var_name)
    f=ROOT.TFile(file_name)
    w=f.Get(workspace_name)
    var=w.var(var_name)
    pdf=w.pdf(pdf_name)
    print(pdf)
   
    new_hist = pdf.createHistogram(histo_name, var, ROOT.RooFit.Binning(bins[0], bins[1], bins[2]))
    new_hist.Scale(N/new_hist.Integral())

    output_obj = ROOT.TFile(output_name, "RECREATE")
    new_hist.Write()

    if see:
        frame = var.frame(RooFit.Title(""))
        dataHist = RooDataHist("Histdata", "binned pdf", RooArgList(var), RooFit.Import(new_hist))
        dataHist.plotOn(frame, RooFit.MarkerStyle(20), RooFit.MarkerSize(0.8), DataError=None)
        pdf.plotOn(frame,RooFit.LineColor(ROOT.kBlue),RooFit.LineWidth(2))
        c2 = ROOT.TCanvas("c2","",800,600)
        frame.GetXaxis().SetTitle("Mass")
        frame.GetYaxis().SetTitle("Events")
        frame.Draw()
        c2.Draw()
    
    output_obj.Close()
    f.Close()
    
    return new_hist, output_name


def histograms_from_json(json_file, bins, total_events, output_name, pdf_name="bernPdf", var_name="x"):
    with open(json_file) as f:
        params = json.load(f)

    coeffs = [params[f"c_{i}"]["value"] for i in range(len(params) - 1)]
    print(coeffs)
    degree = len(coeffs) - 1



    w = ROOT.RooWorkspace("w", "workspace")
    w.factory(f"{var_name}[{xmin},{xmax}]")


    for i, val in enumerate(coeffs):
        w.factory(f"c_{i}[{val},0,100]")
    cList = ROOT.RooArgList()
    for i in range(len(coeffs)):
        cList.add(w.var(f"c_{i}"))

    bern_pdf = ROOT.RooBernsteinFast(degree)(pdf_name, pdf_name, w.var(var_name), cList)
    getattr(w, "import")(bern_pdf, ROOT.RooFit.Rename(pdf_name))

    hist = bern_pdf.createHistogram("shape_hist",w.var(var_name),ROOT.RooFit.Binning(bins[0],bins[1], bins[2]),ROOT.RooFit.Extended(True))

    data_hist = hist.Clone("data_hist")
    print("data hist events: ", data_hist.Integral())
    #norm_factor = total_events / data_hist.Integral()
    #print("norm factor before: ", norm_factor)

    #data_hist.Scale(norm_factor)

    for i in range(1, data_hist.GetNbinsX() + 1):
        mean = data_hist.GetBinContent(i)
        print("mean:  ",mean)
        event_num = ROOT.gRandom.Poisson(mean)
        print("event #: ", event_num)
        data_hist.SetBinContent(i, event_num)

    #norm_factor = total_events/data_hist.Integral()
    #data_hist.Scale(norm_factor)
    #print("norm factor after: ", norm_factor)

    output_file=ROOT.TFile(output_name, "RECREATE")
    hist.Write()
    data_hist.Write()
    output_file.Close()

    return data_hist, output_name


def generate_data_hist(file, bins_num, norm, output_name):

    f=ROOT.TFile.Open(file)
    w=f.Get("w")
    pdf =w.pdf("model")
    x=w.var("mass")
    
    h_pdf =pdf.createHistogram("h_pdf", x, ROOT.RooFit.Binning(bins_num))
    h_pdf.Scale(norm)

    output_file = ROOT.TFile(output_name, "RECREATE")

    #print("H_pdf before: ", h_pdf.Integral())
    for i in range(1, h_pdf.GetNbinsX()+1):
        density=h_pdf.GetBinContent(i)
        bw=h_pdf.GetXaxis().GetBinWidth(i)
        h_pdf.SetBinContent(i, np.random.poisson(density))
        #h_pdf.SetBinContent(i, density*bw)
        #h_pdf.Scale(norm)
        #print("sample: ", np.random.poisson(density*norm))
        #print("norm: ", norm)
        #print("bw*density*norm:  ", density*norm)
    
    #if h_pdf.Integral()!=0:
        #h_pdf.Scale(norm/h_pdf.Integral())
    #    h_pdf.Write()
        #print("data_obs norm: ", h_pdf.Integral())
    #    output_file.Close()
    #else:
    h_pdf.Write()
    output_file.Close()

    #h_pdf.Write()
    #output_file.Close()

    return output_name
























if __name__=="__main__":
    gen_data(bkg_file="EGamma_2018_all_1bad.root", tree="ggH4g_1bad", var="best_4g_corr_mass_m30", cut="(HLT_passed == 1)&&(best_4g_ID_m30 == 1)&&(best_4g_phi1_mass_m30 > 14)&&(best_4g_phi2_mass_m30>14)&&(best_4g_phi1_dxy_m30>-20)&&(best_4g_phi2_dxy_m30>-20)", output_names=[], categories=["lowlow"])
#    gen_data(bkg_file="EGamma_2018_all_1bad.root", tree="ggH4g_1bad", var="best_4g_corr_mass_m30", cut="(HLT_passed == 1)&&(best_4g_ID_m30 == 1)&&(best_4g_phi1_mass_m30 > 14)&&(best_4g_phi2_mass_m30>14)", output_names=[], categories=["lowlow"])
