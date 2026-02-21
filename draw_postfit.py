import ROOT
import subprocess
import tdrstyle
import CMS_lumi
from ggHdatacardmaker import main

def fetchError(q, n):
    l=0
    if n!=0:
        l = ROOT.Math.chisquared_quantile_c(1.0-q, 2.0*n)/2.0
    u=ROOT.Math.chisquared_quantile_c(q, 2.0*n+2)/2.0
    return [l,u]

def getPoisson(h):
    q = (1-0.6827)/2.0
    gRate = ROOT.TGraphAsymmErrors()
    n = 0
    for i in range(1, h.GetNbinsX()+1):
        thresh = h.GetBinCenter(i)
        N = h.GetBinContent(i)
        gRate.SetPoint(n, thresh, N)
        error = fetchError(q, N)
        gRate.SetPointError(n, 0, 0, (N-error[0]), (error[1]-N))
        n+=1
    return gRate


def getPoisson2(h, scale):
    q = (1-0.6827)/2.0
    gRate = ROOT.TGraphAsymmErrors()
    n = 0
    for i in range(1, h.GetNbinsX()+1):
        thresh = h.GetBinCenter(i)
        N = h.GetBinContent(i)
        #if N == 0:
        #    continue
        gRate.SetPoint(n, thresh, scale*N)
        error = fetchError(q, N)
        gRate.SetPointError(n, h.GetBinWidth(i)/2.0, h.GetBinWidth(i)/2.0, scale*(N-error[0]), scale*(error[1]-N))
        n+=1
    return gRate

def plot(MultiDimFit, fitDiagnosticsTest, cat, year, bins, finalstate, physics, order=3):
    loop=0
    for j in range(2):
        tdrstyle.setTDRStyle()
        CMS_lumi.writeExtraText = True
        CMS_lumi.extraText    = "Preliminary"

        if year=="2018":
            CMS_lumi.lumi_13TeV   = "2018, 59 fb^{-1}"
        if year=="Run-2":
            CMS_lumi.lumi_13TeV   = "Run-2, 101 fb^{-1}"
            
        ROOT.gStyle.SetEndErrorSize(2)
        f1 = ROOT.TFile(MultiDimFit)
        f2 = ROOT.TFile(fitDiagnosticsTest)
        r1 = f2.Get("fit_s")
        r2 = f2.Get("fit_b")
        w = f1.Get("w")
        
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
        x = w.var("mass")
        x.setBins(bins[0])
        plot = x.frame()

        #create data_obs as a TH1 obj so i can look through the contents and throw away the error bars (roofit fucks them up. i compute them manually later for the case of N=0 bins):
        data_obs = w.data("data_obs")
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
        for n in range(cloned_data.GetNbinsX()):
            b_n = bin_centers[n]
            c_n = cs[n]
            error=fetchError(q, c_n)
            gres1.SetPoint(n, b_n, c_n)
            gres1.SetPointError(n, 0.0, 0.0, (c_n-error[0]), (error[1]-c_n))

        gres1.SetLineColor(ROOT.kBlack)
        gres1.SetMarkerColor(ROOT.kBlack)
        gres1.SetMarkerStyle(20)

        #grab the s+b and b models from the RooWorkspace and get the normalizations of the models after the fits
        sb_model  = w.pdf("model_s").getPdf(f"{physics}_{finalstate}_{cat}_2018")
        b_model   = w.pdf("model_b").getPdf(f"{physics}_{finalstate}_{cat}_2018")
        if loop==0:
            w.loadSnapshot("MultiDimFit")
        norm_fit = f2.Get("norm_fit_s")
        bkg_norm = norm_fit.find(f"{physics}_{finalstate}_{cat}_2018/background").getVal()
        sig_norm = norm_fit.find(f"{physics}_{finalstate}_{cat}_2018/signal").getVal()
        print("bkg norm: ", bkg_norm)
        print("sig norm: ", sig_norm)

        #RooFormulaVar's that scale the pdf's from the RooWorkspace based on the norms of the fits
        b_func = ROOT.RooFormulaVar("bkg_func","N_bkg * B(x)",f"1*@0",ROOT.RooArgList(b_model))
        sb_func = ROOT.RooFormulaVar("sb_func","N_tot*sb_model",f"1*@0",ROOT.RooArgList(sb_model))

        #plotting. the oreder here matters alot in order to get the brazil plot colors to show up correctly and to get the data points on top of everything
        cloned_data_binned.plotOn(plot,ROOT.RooFit.Binning(bins[0], bins[1], bins[2]),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("data_points"), XErrorSize=0, DataError=None)

        if data_obs_TH1.Integral()!=0:
            print("data obs integral: ", data_obs_TH1.Integral())
            b_model.plotOn(plot,ROOT.RooFit.VisualizeError(r2, 2, ROOT.kFALSE),ROOT.RooFit.FillColor(ROOT.kYellow),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("bkg_2sigma"),ROOT.RooFit.DrawOption("F"))
            b_model.plotOn(plot,ROOT.RooFit.VisualizeError(r2, 1, ROOT.kFALSE),ROOT.RooFit.FillColor(ROOT.kGreen),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("bkg_1sigma"),ROOT.RooFit.DrawOption("F"))
        else:
            print("data is 0. ignoring uncertainty bands because uncertainties on fit parameters are unstable")

        b_model.plotOn(plot,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.LineStyle(2),ROOT.RooFit.Name("bkg_curve"))
        sb_model.plotOn(plot,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("sb_curve"))
        cloned_data_binned.plotOn(plot, ROOT.RooFit.Binning(bins[0], bins[1], bins[2]),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("data_points"), XErrorSize=0, DataError=None)

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
        cate.DrawLatexNDC(0.55, 0.84, r"c#tau = 100 mm, m_{#phi} = 30 GeV")
        cate.DrawLatexNDC(0.55, 0.78, f"category: {cat}")
        #cate.DrawLatexNDC(0.55, 0.72, f"order={order}")
        leg.AddEntry("data_points","Background Data", "p")
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
        cate.Draw("SAME")
        gres1.Draw("p, same")

        #same process as before to define the s-b curve for the bottom panel
        s_minus_b = ROOT.RooFormulaVar("s_minus_b", "S-B", f"({sig_norm}+{bkg_norm})*(@0 - @1)", ROOT.RooArgList(sb_model, b_model))

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
        for n in range(nres):
            x_n=xs[n]
            y_n=ys[n]
            c_n = cs[n]
            g_res.SetPoint(n, x_n, y_n)
            error=fetchError(q, c_n)
            g_res.SetPointError(n, 0.0, 0.0, (c_n-error[0]), (error[1]-c_n))
        g_res.SetMarkerStyle(20)
        g_res.SetMarkerColor(ROOT.kBlack)
        g_res.SetLineColor(ROOT.kBlack)

        #fill axis titles on the bottom panel and write on the canvas
        x.setRange(bins[1], bins[2])
        x.setBins(bins[0])
        lower_plot = x.frame(ROOT.RooFit.Range(bins[1], bins[2]))
        lower_plot.GetXaxis().SetTitle("mass(4#gamma) [GeV]")
        lower_plot.GetYaxis().SetTitle("")
        lower_plot.GetXaxis().SetLabelSize(0.11)
        lower_plot.GetYaxis().SetLabelSize(0.11)
        lower_plot.GetXaxis().SetTitleSize(0.14)
        lower_plot.GetYaxis().SetTitleSize(0.14)
        lower_plot.GetXaxis().SetTitleOffset(1.1)
        lower_plot.GetYaxis().SetTitleOffset(0)

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

        s_minus_b.plotOn(lower_plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("sb") )

        n_bins =bins[0]
        x_min=bins[1]
        x_max=bins[2]

        dx=(x_max-x_min)/n_bins
        g_bzb = ROOT.TGraph(n_bins)

        for i in range(n_bins):
            x_center=x_min+(i+0.5)*dx
            g_bzb.SetPoint(i, x_center, 0.0)

        lower_plot.Draw("axis") 
        g_bzb.SetLineColor(ROOT.kRed)
        g_bzb.SetLineStyle(2)
        g_bzb.SetLineWidth(2)
        g_res.Draw("P") 
        g_bzb.Draw("L")

        leg2 = ROOT.TLegend(0.15, 0.8, 0.68, 0.95)
        leg2.SetHeader("B component subtracted")
        leg2.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.Draw("Same")
        pad1.cd()
        #CMS_lumi.CMS_lumi(can, 4,0, relPosX=0.077, lumi_13TeV="59")
        can.Update()
        can.cd()
        can.Update()
        if loop==0:
            can.SaveAs(f"postfit_{cat}_{year}_postfit.png")
            can.SaveAs(f"postfit_{cat}_{year}_postfit.pdf")
        else:
            can.SaveAs(f"postfit_{cat}_{year}_prefit.png")
            can.SaveAs(f"postfit_{cat}_{year}_prefit.pdf")
        loop+=1

