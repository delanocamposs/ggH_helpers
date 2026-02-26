import ROOT 
import numpy as np
import ggHfitter
import os 
import sys 
import argparse, subprocess
from ggHtools import define_weightMC
from ggHcmsstyle import CMSstyle
ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()

def fetchError(q, n):
    l=0
    if n!=0:
        l = ROOT.Math.chisquared_quantile_c(1.0-q, 2.0*n)/2.0
    u=ROOT.Math.chisquared_quantile_c(q, 2.0*n+2)/2.0
    return [l,u]

def getPoissonGraph(h):
    q = (1 - 0.6827) / 2.0
    g = ROOT.TGraphAsymmErrors()
    n = 0
    for i in range(1, h.GetNbinsX() + 1):
        x  = h.GetBinCenter(i)
        if 111<=x<=138:
            continue
        bw = h.GetBinWidth(i)
        N  = h.GetBinContent(i)
        lo, hi = fetchError(q, N)

        g.SetPoint(n, x, N)
        g.SetPointError(n, bw/2.0, bw/2.0, (N - lo), (hi - N))
        n += 1
    return g


lumi_dict={
        #recommendation for pre/postVFP values here: https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDatasetsUL2016
        "2016preVFP":19500,
        "2016postVFP":16800,
        "2016":19500+16800,
        "2017":41480,
        "2018":59830,
        "Run2":137620,
        "2022":34748,
        "2023":27245,
        "2024":59830,
        "Run3":170857
        }

sideband_SFs={
    "2016preVFP":[0],
    "2016postVFP":[0],
    "2017":[55/33242],
    "2018":[2069/30869]
    }

def save_histos(mass, year, ctau, hist1, hist2):
    outfile=ROOT.TFile(f"sig_bkg_summary_histos_m{mass}_ct{ctau}_year{year}.root", "RECREATE")
    hist1.GetValue().Write()
    hist2.GetValue().Write()
    outfile.Close()
    return 


def run(mass, ctau, year):
    BR=1e-4
    scale_factor=1
    xsec=52.143
    cut_string=f"HLT_passed==1&&best_4g_phi1_dxy_m{mass}>-20&&best_4g_phi2_dxy_m{mass}>-20"
    preselection=f"(Photon_preselection[best_4g_idx1_m{mass}]==1)&&(Photon_preselection[best_4g_idx2_m{mass}]==1)&&(Photon_preselection[best_4g_idx3_m{mass}]==1)&&(Photon_preselection[best_4g_idx4_m{mass}]==1)"
    blind=f"((best_4g_corr_mass_m{mass}<110)||(best_4g_corr_mass_m{mass}>140))"
    #blind=f"(best_4g_corr_mass_m{mass}<110)"
    
    
    c=ROOT.TCanvas("", "", 800, 800)
    right=0.1
    left=0.14
    up=0.08
    down=0
    l=ROOT.TLegend(0.5+right-left,0.67+up-down, 0.95+right-left, 0.8+up-down)
    bins=[55, 70, 180]
    var=f"best_4g_corr_mass_m{mass}"
    
    #sig=f"ggH4g_M{mass}_ctau{ctau}_{year}_0_ggH4g_M{mass}_ctau{ctau}_{year}_ggH4g.root"
    sig=f"ggH4g_M{mass}_ctau{ctau}_2018_0_ggH4g_M{mass}_ctau{ctau}_2018_ggH4g.root"
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
    
    signal_df=signal_df.Filter(f"{cut_string}&&best_4g_ID_m{mass}==1&&{preselection}")
    data_ID_df=data_ID_df.Filter(f"{preselection}&&{cut_string}&&best_4g_ID_m{mass}==1&&{blind}")
    
    h_sig=signal_df.Histo1D(("signal_hist", f"signal_hist;{var};Events", bins[0], bins[1], bins[2]), f"best_4g_corr_mass_m{mass}", "event_weight")
    h_sig.Scale(lumi_dict[year])
    
    h_data=data_ID_df.Histo1D(("data_hist", f"data_hist;4#gamma mass;Events", bins[0], bins[1], bins[2]), f"{var}")
    N_data_ID=h_data.Integral()
    N_sig=h_sig.Integral()
    
    data_ID_color=ROOT.kOrange-3
    signal_color=ROOT.kAzure-4
    b_fit_clone=h_data.Clone("b_fit_clone")
    b_fit_clone.SetLineWidth(4)
    h_data.SetFillStyle(1001)
    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetLineStyle(1)
    h_data.SetFillColor(data_ID_color)
    h_sig.SetFillStyle(1001)
    h_sig.SetFillColor(signal_color)
    h_sig.SetLineColor(ROOT.kBlack)
    h_sig.SetLineStyle(1)
    per_bin=round((bins[2]-bins[1])/bins[0], 2)

    ROOT.gStyle.SetEndErrorSize(6)
    data_ID_err = getPoissonGraph(h_data)
    data_ID_err.SetLineColor(ROOT.kBlack)
    data_ID_err.SetMarkerColor(ROOT.kBlack)
    data_ID_err.SetMarkerStyle(20)
    data_ID_err.SetMarkerSize(0.9)

    
    data_ID_clone=h_data.Clone("data_ID_clone")
    h_sig_clone=h_sig.Clone("h_sig_clone")
    h_sig_clone.SetFillColor(signal_color)
    data_ID_clone.SetFillColor(data_ID_color)
    data_ID_clone.SetFillStyle(1001)
    h_sig_clone.SetFillStyle(1001)
    data_ID_clone.SetLineColor(ROOT.kBlack)
    b_fit_clone.SetLineColor(ROOT.kRed)
    h_sig_clone.SetLineColor(ROOT.kBlack)
    
    h_data.SetLineColor(ROOT.kBlack)
    h_sig.SetLineColor(ROOT.kBlack)
    h_data.SetLineWidth(2)
    h_sig.SetLineWidth(2)
    ymax=100*max(h_data.GetMaximum(),h_sig.GetMaximum())
    h_data.SetMaximum(20)
    h_data.SetMinimum(0.01)

    data_ID_points = h_data.Clone("data_ID_points")
    data_ID_points.SetFillStyle(0)
    data_ID_points.SetLineColor(ROOT.kBlack)
    data_ID_points.SetMarkerColor(ROOT.kBlack)
    data_ID_points.SetMarkerStyle(20)
    data_ID_points.SetMarkerSize(0.9)

    ax= h_data.GetXaxis()
    blind_low=110.0
    blind_high=140.0
    first_bin = ax.FindFixBin(blind_low+1e-6)
    last_bin=ax.FindFixBin(blind_high-1e-6)
    x1=ax.GetBinLowEdge(first_bin)
    x2=ax.GetBinLowEdge(last_bin)
    y1=h_data.GetMinimum()
    y2=h_data.GetMaximum()
    line1 = ROOT.TLine(x1, y1, x1, y2)
    line2 = ROOT.TLine(x2, y1, x2, y2)
    for ln in (line1, line2):
        ln.SetLineColor(ROOT.kBlack)
        ln.SetLineWidth(2)
        ln.SetLineStyle(2)

    #c.SetLogy()
    
    save_histos(mass,year, ctau,h_data,h_sig)
    c.cd()
    ggHfitter.fitBKG(f"sig_bkg_summary_histos_m{mass}_ct{ctau}_year{year}.root", "data_hist", f"bkg_fit_result_m{mass}_ct{ctau}_year{year}.root", order=4)
    ggHfitter.fitSIG(f"sig_bkg_summary_histos_m{mass}_ct{ctau}_year{year}.root", "signal_hist", f"sig_fit_result_m{mass}_ct{ctau}_year{year}.root")
    sig_fit_file=ROOT.TFile.Open(f"sig_fit_result_m{mass}_ct{ctau}_year{year}.root")
    wfitsig = sig_fit_file.Get("w")
    xsig = wfitsig.var("mass")
    pdfsig = wfitsig.pdf("model")
    xmin = h_sig.GetXaxis().GetXmin()
    xmax = h_sig.GetXaxis().GetXmax()
    N_sig = h_sig.Integral()
    binw = h_sig.GetXaxis().GetBinWidth(1)
    def pdf_as_events(xx, pp):
        xsig.setVal(xx[0])
        return N_sig*binw*pdfsig.getVal(ROOT.RooArgSet(xsig))
    sig_curve = ROOT.TF1("sig_curve_roofit", pdf_as_events, xmin, xmax, 0)
    sig_curve.SetLineColor(ROOT.kRed)
    sig_curve.SetLineWidth(2)

    bkg_fit_file=ROOT.TFile.Open(f"bkg_fit_result_m{mass}_ct{ctau}_year{year}.root")
    wfitbkg = bkg_fit_file.Get("w")
    xbkg = wfitbkg.var("mass")
    pdfbkg = wfitbkg.pdf("model")
    N_bkg=h_data.Integral()
    def bkg_as_events(xx, pp):
        xbkg.setVal(xx[0])
        return N_bkg * binw * pdfbkg.getVal(ROOT.RooArgSet(xbkg))
    xmin = h_data.GetXaxis().GetXmin()
    xmax = h_data.GetXaxis().GetXmax()
    binw = h_data.GetXaxis().GetBinWidth(1)
    def pdf_as_events(xx, pp):
        xbkg.setVal(xx[0])
        return N_bkg*binw*pdfbkg.getVal(ROOT.RooArgSet(xbkg))
    bkg_curve = ROOT.TF1("bkg_curve_roofit", pdf_as_events, xmin, xmax, 0)
    bkg_curve.SetLineColor(ROOT.kRed)
    bkg_curve.SetLineWidth(2)
    def sig_plus_bkg(xx, pp):
        xbkg.setVal(xx[0])
        xsig.setVal(xx[0])
        b = N_bkg * binw * pdfbkg.getVal(ROOT.RooArgSet(xbkg))
        s = scale_factor*N_sig * binw * pdfsig.getVal(ROOT.RooArgSet(xsig))
        return b + s        
    splusb_curve = ROOT.TF1("splusb_curve", sig_plus_bkg, 110, 138, 0)
    splusb_curve.SetLineColor(signal_color)
    splusb_curve.SetLineWidth(3)
    splusb_curve.SetLineStyle(2)
    c.cd()


    n_points = 200
    x_vals = np.linspace(110, 140, n_points)
    shade = ROOT.TGraph(2 * n_points)
    for i, x in enumerate(x_vals):
        xbkg.setVal(x)
        xsig.setVal(x)
        b = N_bkg*binw*pdfbkg.getVal(ROOT.RooArgSet(xbkg))
        s = scale_factor*N_sig*binw*pdfsig.getVal(ROOT.RooArgSet(xsig))
        shade.SetPoint(i, x, b+s)
    for i, x in enumerate(reversed(x_vals)):
        xbkg.setVal(x)
        b = N_bkg*binw*pdfbkg.getVal(ROOT.RooArgSet(xbkg))
        shade.SetPoint(n_points+i, x, b)
    shade.SetFillColor(ROOT.kAzure-4)
    shade.SetFillStyle(1001)
    shade.SetFillColorAlpha(ROOT.kAzure-4, 0.3)
    shade.SetLineColor(0)

    h_data.Draw("HIST")
    data_ID_points.Draw("P, SAME")
    #h_sig.Draw("HIST,SAME")
    #sig_curve.Draw("SAME")
    bkg_curve.Draw("SAME")
    shade.Draw("F SAME")
    splusb_curve.Draw("SAME")
    line1.Draw("SAME")
    line2.Draw("SAME")
    data_ID_err.Draw("P E1 SAME") 
    
    ROOT.gPad.RedrawAxis()  
    ROOT.gPad.RedrawAxis("g") 
    
    c.Update()
    l.SetTextSize(0.05)
    #CMSstyle(c, l,year, lumi_dict[year],[f"H #rightarrow #phi#phi #rightarrow 4#gamma", f"c#tau = {ctau} mm", f"m_{{#phi}} = {mass} #font[42]{{GeV}}", f"#font[42]{{BR}}(H #rightarrow #phi#phi)(#phi #rightarrow #gamma#gamma) = {BR}", f"#sigma = {xsec} #font[42]{{pb}}"])
    CMSstyle(c, l,year, lumi_dict[year],[f"H #rightarrow #phi#phi #rightarrow 4#gamma", f"c#tau = {ctau} mm", f"m_{{#phi}} = {mass} #font[42]{{GeV}}"])
    l.AddEntry(data_ID_clone, f"Blind {year} Data", "f")
    if scale_factor!=1:
        l.AddEntry(h_sig_clone, f"S+B fit (x{scale_factor})", "f")
    else:
        l.AddEntry(h_sig_clone, f"S+B fit", "f")
    l.AddEntry(b_fit_clone, f"B only fit", "f")
    l.Draw("SAME")
    
    c.Update()
    c.SaveAs(f"summary_fit_plot{year}_ct{ctau}_m{mass}.png")


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

