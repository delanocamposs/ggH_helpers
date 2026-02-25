import ROOT 
import os 
import sys 
import argparse, subprocess
from ggHtools import define_weightMC
from ggHcmsstyle import CMSstyle
ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()

lumi_dict={
        "2016preVFP":36310,
        "2016postVFP":36310,
        "2017":41480,
        "2018":59830,
        "Run 2":137620,
        "2022":34748,
        "2023":27245,
        "2024":59830,
        "Run 3":170857
        }


def run(mass, ctau, year):
    BR=1e-4
    xsec=52.143
    lumi=59830
    cut_string=f"HLT_passed==1&&best_4g_phi1_dxy_m{mass}>-20&&best_4g_phi2_dxy_m{mass}>-20"
    preselection=f"(Photon_preselection[best_4g_idx1_m{mass}]==1)&&(Photon_preselection[best_4g_idx2_m{mass}]==1)&&(Photon_preselection[best_4g_idx3_m{mass}]==1)&&(Photon_preselection[best_4g_idx4_m{mass}]==1)"
    #blind=f"(best_4g_corr_mass_m{mass}<110)"
    blind=f"((best_4g_corr_mass_m{mass}<110)||(best_4g_corr_mass_m{mass}>140))"
    
    
    c=ROOT.TCanvas("", "", 800, 800)
    right=0.08
    left=0.14
    up=0.08
    down=0
    l=ROOT.TLegend(0.5+right-left,0.67+up-down, 0.95+right-left, 0.8+up-down)
    bins=[55, 60, 180]
    var=f"best_4g_corr_mass_m{mass}"
    
    sig=f"ggH4g_M{mass}_ctau{ctau}_{year}_0_ggH4g_M{mass}_ctau{ctau}_{year}_ggH4g.root"
    bkg=f"EGamma_{year}_all_ggH4g.root"
    
    sig_open=ROOT.TFile.Open(sig)
    sumw=0.0
    signal_df=ROOT.RDataFrame("ggH4g", sig)
    with sig_open as f:
        runs_tree = f.Get("Runs")
        for entry in runs_tree:
            sumw += entry.genEventSumw
    weight_formula = f"(genWeight / {sumw}) * {xsec} * {BR} * Pileup_weight"
    signal_df=signal_df.Define("event_weight", weight_formula)
    
    data_ID_df=ROOT.RDataFrame("ggH4g", bkg)
    data_pre_df=ROOT.RDataFrame("ggH4g", bkg)
    
    signal_df=signal_df.Filter(f"{cut_string}&&best_4g_ID_m{mass}==1")
    data_ID_df=data_ID_df.Filter(f"{preselection}&&{cut_string}&&best_4g_ID_m{mass}==1&&{blind}")
    data_pre_df=data_pre_df.Filter(f"{cut_string}&&{preselection}&&{blind}")
    
    signal_histo=signal_df.Histo1D(("hist1_1", f"hist1_1;{var};Events", bins[0], bins[1], bins[2]), f"{var}", "event_weight")
    signal_histo.Scale(lumi)
    
    data_ID_histo=data_ID_df.Histo1D(("hist2", f"hist2;4#gamma mass;Events", bins[0], bins[1], bins[2]), f"{var}")
    data_pre_histo=data_pre_df.Histo1D(("hist3", f"hist3;4#gamma mass;Events", bins[0], bins[1], bins[2]), f"{var}")
    N_data_ID=data_ID_histo.Integral()
    N_data_pre=data_pre_histo.Integral()
    N_sig=signal_histo.Integral()
    print("N_ID", N_data_ID)
    print("N_pre", N_data_pre)
    print("N_sig", N_sig)
    
    data_ID_color=ROOT.kOrange+8
    data_pre_color=ROOT.kPink+1
    signal_color=ROOT.kAzure-4
    data_ID_histo.SetFillStyle(1001)
    data_ID_histo.SetLineColor(ROOT.kBlack)
    data_ID_histo.SetLineStyle(1)
    data_ID_histo.SetFillColor(data_ID_color)
    data_pre_histo.SetFillStyle(1001)
    data_pre_histo.SetLineColor(ROOT.kBlack)
    data_pre_histo.SetLineStyle(1)
    data_pre_histo.SetFillColor(data_pre_color)
    signal_histo.SetFillStyle(1001)
    signal_histo.SetFillColor(signal_color)
    signal_histo.SetLineColor(ROOT.kBlack)
    signal_histo.SetLineStyle(1)
    per_bin=round((bins[2]-bins[1])/bins[0], 2)
    
    data_ID_clone=data_ID_histo.Clone("data_ID_clone")
    data_pre_clone=data_pre_histo.Clone("data_pre_clone")
    signal_histo_clone=signal_histo.Clone("signal_histo_clone")
    signal_histo_clone.SetFillColor(signal_color)
    data_ID_clone.SetFillColor(data_ID_color)
    data_pre_clone.SetFillColor(data_pre_color)
    data_ID_clone.SetFillStyle(1001)
    data_pre_clone.SetFillStyle(1001)
    signal_histo_clone.SetFillStyle(1001)
    data_ID_clone.SetLineColor(ROOT.kBlack)
    data_pre_clone.SetLineColor(ROOT.kBlack)
    signal_histo_clone.SetLineColor(ROOT.kBlack)
    
    data_ID_histo.SetLineColor(ROOT.kBlack)
    data_pre_histo.SetLineColor(ROOT.kBlack)
    signal_histo.SetLineColor(ROOT.kBlack)
    data_ID_histo.SetLineWidth(2)
    data_pre_histo.SetLineWidth(2)
    signal_histo.SetLineWidth(2)
    ymax=100*max(data_ID_histo.GetMaximum(), data_pre_histo.GetMaximum(), signal_histo.GetMaximum())
    data_pre_histo.SetMaximum(ymax)
    data_pre_histo.SetMinimum(0.01)
    
    c.SetLogy()
    data_pre_histo.Draw("HIST")
    data_ID_histo.Draw("HIST,SAME")
    signal_histo.Draw("HIST,SAME")
    
    ROOT.gPad.RedrawAxis()  
    ROOT.gPad.RedrawAxis("g") 
    
    c.Update()
    l.SetTextSize(0.05)
    CMSstyle(c, l,year, lumi_dict[year],[f"H #rightarrow #phi#phi #rightarrow 4#gamma", f"c#tau = {ctau} mm", f"m_{{#phi}} = {mass} #font[42]{{GeV}}"])
    l.AddEntry(signal_histo_clone, "Signal MC", "f")
    l.AddEntry(data_ID_clone, "ID Data (blinded)", "f")
    l.AddEntry(data_pre_clone, "pre Data (blinded)", "f")
    l.Draw("SAME")
    
    c.Update()
    c.SaveAs(f"ID_pre_sideband_{year}_ct{ctau}_m{mass}.png")


if __name__=="__main__":
    parser = argparse.ArgumentParser("Running sideband plot")
    parser.add_argument("-m","--mass", type=str, help="mass of sample")
    parser.add_argument("-ct","--ctau", type=str, help="lifetime of sample")
    parser.add_argument("-y","--year", type=str, help="year of MC and data")
    args = parser.parse_args()
    mass=args.mass
    lifetime=args.ctau
    year=args.year
    run(mass, lifetime, year)
    
