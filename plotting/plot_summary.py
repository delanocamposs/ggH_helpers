import ROOT
from datacard import ggHfitter
import argparse
import subprocess
from plotting import tdrstyle
from plotting import CMS_lumi
from datacard.ggHdatacardmaker import main
from ggHparameters import lumi
import ggHcuts as cuts
from plotting.plottingtools import fetchError, getPoisson, getPoisson2, save_histos
ROOT.gROOT.SetBatch(True) 

def run(mass, ctau, year):
    tdrstyle.setTDRStyle()
    CMS_lumi.writeExtraText = True
    CMS_lumi.extraText="Preliminary"
    CMS_lumi.lumi_13TeV=f"{year}, {lumi[year]/1000} fb^{-1}"
    ROOT.gStyle.SetEndErrorSize(2)

    
    right=0.1
    left=0.14
    up=0.08
    down=0
    l=ROOT.TLegend(0.5+right-left,0.67+up-down, 0.95+right-left, 0.8+up-down)
    bins_sig = [45, 80, 170]
    bins_data = [45, 80, 170]
    bins=bins_data
    var=f"best_4g_corr_mass_m{mass}"

    run2_years = ["2017", "2018"]
    run3_years = ["2022preEE", "2022postEE", "2023preBPix","2023postBPix", "2024"]
    years_2022 = ["2022preEE", "2022postEE"]
    years_2023 = ["2023preBPix", "2023postBPix"]
    if year == "Run2":
        years_to_process = run2_years
    elif year == "Run3":
        years_to_process = run3_years
    elif year == "2022":
        years_to_process = years_2022
    elif year == "2023":
        years_to_process = years_2023
    else:
        years_to_process = [year]

    bkg=f"/eos/uscms/store/user/dacampos/analysis/data/EGamma_{year}_updated/EGamma_{year}_all_ggH4g.root"
    data_ID_df=ROOT.RDataFrame("ggH4g", bkg)
    data_ID_df=data_ID_df.Filter(cuts.combine(cuts.preselection(mass), cuts.trigger(), cuts.dxy_valid(mass), cuts.full_id(mass), cuts.blind(mass)))
    h_data=data_ID_df.Histo1D(("data_hist", f"data_hist;4#gamma mass;Events", bins_data[0], bins_data[1], bins_data[2]), f"{var}")

    #signal histogram built separately for each year.
    h_sig_total = None
    for y in years_to_process:
        sig_y = f"/eos/uscms/store/user/dacampos/analysis/signal/ggH4g_M{mass}_ctau{ctau}_{y}_0_ggH4g_M{mass}_ctau{ctau}_{y}_ggH4g.root"
        sig_open_y = ROOT.TFile.Open(sig_y)
        sumw_y = 0.0
        with sig_open_y as f:
            runs_tree = f.Get("Runs")
            for entry in runs_tree:
                sumw_y += entry.genEventSumw
        signal_df_y = ROOT.RDataFrame("ggH4g", sig_y)
        weight_y = cuts.mc_weight(sumw_y)
        signal_df_y = signal_df_y.Define("event_weight", weight_y)
        signal_df_y = signal_df_y.Filter(cuts.combine(cuts.trigger(), cuts.dxy_valid(mass), cuts.full_id(mass), cuts.preselection(mass), cuts.pileup()))
        core = signal_df_y.Filter(f"abs({var}-125.0)<0.25")
        print(f"[{y}] max event_weight: all={signal_df_y.Max('event_weight').GetValue():.4f} "
              f"core={core.Max('event_weight').GetValue():.4f}")
        print(f"[{y}] max Pileup_weight core={core.Max('Pileup_weight').GetValue():.3f} "
              f"max genWeight core={core.Max('genWeight').GetValue():.3f}")
        core.Display(["event_weight", "genWeight", "Pileup_weight", var], 8).Print()
        print(f"[{y}] max Pileup_weight ALL = {signal_df_y.Max('Pileup_weight').GetValue():.1f}")
        signal_df_y.Filter("event_weight > 1e-4").Display(["event_weight", "Pileup_weight", var], 10).Print()
        h_y = signal_df_y.Histo1D((f"sig_tmp_{y}", f"sig_tmp_{y}", bins_sig[0], bins_sig[1], bins_sig[2]), var, "event_weight")
        h_y_clone = h_y.GetValue().Clone(f"sig_scaled_{y}")
        h_y_clone.Scale(lumi[y])
        if h_sig_total is None:
            h_sig_total = h_y_clone.Clone("signal_hist")
        else:
            h_sig_total.Add(h_y_clone)
    nonzero = sum(1 for i in range(1, h_sig_total.GetNbinsX()+1) if h_sig_total.GetBinContent(i) > 0)
    print(f"signal nonzero bins={nonzero}, integral={h_sig_total.Integral():.3f}")
    imax = max(range(1, h_sig_total.GetNbinsX()+1), key=lambda i: h_sig_total.GetBinContent(i))
    print(f"peak bin center={h_sig_total.GetBinCenter(imax):.2f} content={h_sig_total.GetBinContent(imax):.3f}")
    lo = h_sig_total.FindBin(120); hi = h_sig_total.FindBin(130)
    print(f"integral in 120-130={h_sig_total.Integral(lo, hi):.3f} of total {h_sig_total.Integral():.3f}")

    class HistProxy:
        def __init__(self, h):
            self._h = h
        def GetValue(self):
            return self._h
        def Scale(self, s):
            self._h.Scale(s)
        def Integral(self):
            return self._h.Integral()
    h_sig = HistProxy(h_sig_total)

    print(f"signal effective entries={h_sig_total.GetEffectiveEntries():.2f} raw entries={h_sig_total.GetEntries():.0f}")
    save_histos(mass, year, ctau, h_data, h_sig)
    sresult, bresult=ggHfitter.fitSIGBKG(f"sig_bkg_summary_histos_m{mass}_ct{ctau}_year{year}.root", "signal_hist", "data_hist", f"SB_fit_result_m{mass}_ct{ctau}_year{year}.root", order=4)
    SB_file=ROOT.TFile.Open(f"SB_fit_result_m{mass}_ct{ctau}_year{year}.root")
    print("got here")
    w = SB_file.Get("w")
    x= w.var("mass")
    s_model = w.pdf("model_s")
    b_model = w.pdf("model_b")

    # Remove fitrange attributes so plotOn draws over the full range
    s_model.removeStringAttribute("fitrange")
    b_model.removeStringAttribute("fitrange")

    # Get normalizations from data
    bkg_norm = h_data.GetValue().Integral()
    sig_norm = h_sig.GetValue().Integral()

    # Build the combined S+B model
    n_sig = ROOT.RooRealVar("n_sig", "n_sig", sig_norm)
    n_bkg = ROOT.RooRealVar("n_bkg", "n_bkg", bkg_norm)
    sb_model = ROOT.RooAddPdf("model_sb", "S+B model",
        ROOT.RooArgList(s_model, b_model),
        ROOT.RooArgList(n_sig, n_bkg))

    #set the proportions and make it beautiful
    can = ROOT.TCanvas("c")
    pad1 = ROOT.TPad("pad1", "pad1", 0,   0.3, 1, 1.0)
    pad1.SetTopMargin(0.08506945)
    pad1.SetBottomMargin(0.00)  
    pad1.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.Draw()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.00, 1, 0.25)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.15)
    pad2.SetRightMargin(0.05)
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    pad2.Draw()
    pad1.cd()
    x.setBins(bins[0])
    plot=x.frame()

    #create data_obs as a TH1 obj so i can look through the contents and throw away the error bars (roofit fucks them up. i compute them manually later for the case of N=0 bins):
    data_obs = w.data("data_b")
    data_obs_binned = ROOT.RooDataHist("data_obs_binned", "binned data", ROOT.RooArgSet(x), data_obs)
    data_obs_TH1 = data_obs_binned.createHistogram("mass")

    cs=[]
    cloned_data = data_obs_TH1.Clone("int_hist")
    for i in range(1, cloned_data.GetNbinsX()+1):
        count = cloned_data.GetBinContent(i)
        cs.append(cloned_data.GetBinContent(i))
        cloned_data.SetBinError(i, 0)
    cloned_data_binned = ROOT.RooDataHist("data_obs_binned", "binned hist", ROOT.RooArgList(x), cloned_data)

    #create points at the locations of the data bins with correct error bars from fetcError earlier in the code:
    gres1=ROOT.TGraphAsymmErrors()
    q = (1-0.6827)/2.0
    bin_centers = [cloned_data.GetBinCenter(i) for i in range(1, cloned_data.GetNbinsX()+1)]
    n_point = 0
    for n in range(cloned_data.GetNbinsX()):
        b_n = bin_centers[n]
        # Skip bins in the blinded signal region
        if 110 <= b_n <= 140:
            continue
        c_n = cs[n]
        error=fetchError(q, c_n)
        gres1.SetPoint(n_point, b_n, c_n)
        gres1.SetPointError(n_point, 0.0, 0.0, (c_n-error[0]), (error[1]-c_n))
        n_point += 1

    gres1.SetLineColor(ROOT.kBlack)
    gres1.SetMarkerColor(ROOT.kBlack)
    gres1.SetMarkerStyle(20)

    #plotting. the order here matters a lot in order to get the brazil plot colors to show up correctly and to get the data points on top of everything
    # Plot data invisibly to set frame normalization and create "data_points" for residHist
    cloned_data_binned.plotOn(plot,ROOT.RooFit.Binning(bins[0], bins[1], bins[2]),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("data_points"),ROOT.RooFit.XErrorSize(0),ROOT.RooFit.Invisible())

    if data_obs_TH1.Integral()!=0:
        n_data = data_obs_TH1.Integral()
        print("data obs integral: ", n_data)
        b_model.plotOn(plot,ROOT.RooFit.VisualizeError(bresult, 2, ROOT.kFALSE),ROOT.RooFit.FillColor(ROOT.kYellow),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("bkg_2sigma"),ROOT.RooFit.DrawOption("F"))
        b_model.plotOn(plot,ROOT.RooFit.VisualizeError(bresult, 1, ROOT.kFALSE),ROOT.RooFit.FillColor(ROOT.kGreen),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("bkg_1sigma"),ROOT.RooFit.DrawOption("F"))
    else:
        print("data is 0. ignoring uncertainty bands because uncertainties on fit parameters are unstable")
        n_data = 0

    norm_to_use = data_obs_TH1.Integral() if data_obs_TH1.Integral() != 0 else 1
    b_model.plotOn(plot,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.LineStyle(2),ROOT.RooFit.Name("bkg_curve"),ROOT.RooFit.Normalization(norm_to_use, ROOT.RooAbsReal.NumEvent))
    sb_model.plotOn(plot,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("sb_curve"),ROOT.RooFit.Normalization(norm_to_use, ROOT.RooAbsReal.NumEvent))

    integral_sb = sb_model.getVal(ROOT.RooArgSet(x))
    print("norm after scaling sb: ", integral_sb)
    integral_b = b_model.getVal(ROOT.RooArgSet(x))
    print("norm after scaling b: ", integral_b)

    y_ax_val=round((bins[2]-bins[1])/(bins[0]),2)
    plot.GetYaxis().SetTitle(f"Events/{y_ax_val} (GeV)")
    plot.GetXaxis().SetLabelSize(0)    
    plot.GetXaxis().SetTitleSize(0)
    plot.GetYaxis().SetTitleOffset(0.65)

    if cs:
        max_y = max(cs[n]+(fetchError(q,cs[n])[1]-cs[n]) for n in range(len(cs)))
        plot.SetMaximum(2.0*max_y)
    else:
        plot.SetMaximum(10) 
    plot.Draw()

    #legend and contents to draw on the canvas
    leg = ROOT.TLegend(0.2, 0.5, 0.68, 0.88)
    cate = ROOT.TLatex()
    cate.SetTextSize(0.06)
    cate.DrawLatexNDC(0.55, 0.84, rf"c#tau = {ctau} mm, m_{{#phi}} = {mass} GeV")
    cate.DrawLatexNDC(0.55, 0.78, f"{year}")
    leg.AddEntry("data_points","Blinded Data", "p")
    leg.AddEntry("sb_curve","S+B fit sum", "L")
    leg.AddEntry("bkg_curve", "B component", "L")
    if data_obs_TH1.Integral()!=0:
        leg.AddEntry("bkg_1sigma",r"\pm 1 \sigma", "f")
        leg.AddEntry("bkg_2sigma",r"\pm 2 \sigma", "f")
    else: 
        print("")
    leg.SetHeader("H #rightarrow #phi#phi #rightarrow 4#gamma")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw("SAME")
    gres1.Draw("p, same")

    #build the signal-only curve for the bottom panel by subtracting the B curve from the S+B curve on the upper panel
    #this guarantees the bottom panel signal shape matches what's shown on top
    sb_curve_obj = plot.getCurve("sb_curve")
    b_curve_obj = plot.getCurve("bkg_curve")

    n_curve_points = sb_curve_obj.GetN()
    g_sig_only = ROOT.TGraph(n_curve_points)
    for i in range(n_curve_points):
        x_i = sb_curve_obj.GetX()[i]
        y_sb = sb_curve_obj.GetY()[i]
        y_b = b_curve_obj.interpolate(x_i)
        g_sig_only.SetPoint(i, x_i, y_sb - y_b)
    g_sig_only.SetLineColor(ROOT.kRed)
    g_sig_only.SetLineWidth(2)

    #cd into bottom panel and begin populating it with s-b curve and datapoints
    pad2.cd()
    resid_hist = plot.residHist("data_points", "bkg_curve")
    nres  = resid_hist.GetN()
    xs = resid_hist.GetX()
    ys = resid_hist.GetY()
    x_vals = [xs[i] for i in range(nres)]
    y_vals = [ys[i] for i in range(nres)]

    #define the residual data points and set errors/locations manually
    g_res=ROOT.TGraphAsymmErrors()
    q = (1-0.6827)/2.0
    n_res_point = 0
    for n in range(nres):
        x_n=xs[n]
        y_n=ys[n]
        # Skip bins in the blinded signal region
        if 110 <= x_n <= 140:
            continue
        c_n = cs[n]
        g_res.SetPoint(n_res_point, x_n, y_n)
        error=fetchError(q, c_n)
        g_res.SetPointError(n_res_point, 0.0, 0.0, (c_n-error[0]), (error[1]-c_n))
        n_res_point += 1
    g_res.SetMarkerStyle(20)
    g_res.SetMarkerColor(ROOT.kBlack)
    g_res.SetLineColor(ROOT.kBlack)

    #fill axis titles on the bottom panel and write on the canvas
    x.setRange(bins[1], bins[2])
    x.setBins(bins[0])
    lower_plot = x.frame(ROOT.RooFit.Range(bins[1], bins[2]))
    lower_plot.GetXaxis().SetTitle("mass(4#gamma) [GeV]")
    lower_plot.GetYaxis().SetTitle("")
    lower_plot.GetXaxis().SetLabelSize(0.09)
    lower_plot.GetYaxis().SetLabelSize(0.07)
    lower_plot.GetXaxis().SetTitleSize(0.12)
    lower_plot.GetYaxis().SetTitleSize(0.10)
    lower_plot.GetXaxis().SetTitleOffset(1.1)
    lower_plot.GetYaxis().SetTitleOffset(0.5)
    lower_plot.GetYaxis().SetNdivisions(505)

    #dynamic y axis scaling
    n_points = g_res.GetN()
    if n_points > 0:
        ymin = float("inf")
        ymax = float("-inf")
        for i in range(n_points):
            y = g_res.GetY()[i]
            y_err_low = g_res.GetErrorYlow(i)
            y_err_high = g_res.GetErrorYhigh(i)
            ymin = min(ymin, y-y_err_low)
            ymax = max(ymax, y+y_err_high)

        y_range = ymax-ymin
        margin = 0.3*y_range if y_range>0 else 1.0
        lower_plot.SetMinimum(ymin-margin)
        lower_plot.SetMaximum(ymax+margin)
    else:
        lower_plot.SetMinimum(-2)
        lower_plot.SetMaximum(2)

    lower_plot.Draw("axis")
    lower_plot.GetYaxis().SetTitle("")
    lower_plot.SetTitle("")

    n_bins =bins[0]
    x_min=bins[1]
    x_max=bins[2]

    dx=(x_max-x_min)/n_bins
    g_bzb = ROOT.TGraph(n_bins)

    for i in range(n_bins):
        x_center=x_min+(i+0.5)*dx
        g_bzb.SetPoint(i, x_center, 0.0)

    g_bzb.SetLineColor(ROOT.kRed)
    g_bzb.SetLineStyle(2)
    g_bzb.SetLineWidth(2)
    g_bzb.Draw("L")
    g_sig_only.Draw("L same")
    g_res.Draw("P")

    leg2 = ROOT.TLegend(0.15, 0.8, 0.68, 0.95)
    leg2.SetHeader("B component subtracted")
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.Draw("Same")
    pad1.cd()
    can.Update()
    can.cd()
    can.Update()
    can.SaveAs(f"summary_plot_m{mass}_ct{ctau}_{year}.png")
    can.SaveAs(f"summary_plot_m{mass}_ct{ctau}_{year}.pdf")
    
    can.Close()
    SB_file.Close()
    del can, pad1, pad2, plot, lower_plot
    del w, x, s_model, b_model, sb_model
    del bresult, sresult
    ROOT.gROOT.GetListOfCanvases().Clear()
    ROOT.gROOT.GetListOfFiles().Clear()

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
