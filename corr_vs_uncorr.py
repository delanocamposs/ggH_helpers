import ROOT 
from ggHcmsstyle import CMSstyle


c=ROOT.TCanvas("", "", 800, 800)
right=0.05
left=0.04
up=0.05
down=0.02
l=ROOT.TLegend(0.5+right-left,0.67+up-down, 0.95+right-left, 0.8+up-down)
file1=ROOT.TFile.Open("ggH_M30_ctau0_ggH4g_double.root")
file2=ROOT.TFile.Open("ggH_M30_ctau0_ggH4g_EGM_double.root")

#file2=ROOT.TFile.Open("none_2018_4g_ggH/rate_histos_none_2018.root")
#corr=file2.Get("hist2")
#file3=ROOT.TFile.Open("none_2018_4g_ggH/fit_bkg_none_2018_fit.root")
#w=file3.Get("w")
#pdf=w.pdf("model")
#x=w.var("mass")
#frame=x.frame()


tree1=file1.Get("ggH4g")
tree2=file2.Get("ggH4g")
uncorr=ROOT.TH1F("", "",100, 110, 140)
corr=ROOT.TH1F("", "", 100, 110, 140)
uncorr_color=ROOT.kOrange+8
#uncorr_color=ROOT.kRed+2
corr_color=ROOT.kAzure-4
#corr_color=ROOT.kOrange-4
uncorr.SetFillStyle(1001)
uncorr.SetFillColor(uncorr_color)
uncorr.SetLineColor(ROOT.kBlack)
uncorr.SetLineStyle(1)
corr.SetFillStyle(1001)
corr.SetFillColor(corr_color)
corr.SetLineColor(ROOT.kBlack)
corr.SetLineStyle(1)
corr.GetXaxis().SetTitle("4#gamma Mass [GeV]")
corr.GetYaxis().SetTitle("Events / 0.33 GeV")

uncorr_clone=uncorr.Clone("uncorr_clone")
corr_clone=corr.Clone("corr_clone")
uncorr_clone.SetFillColor(uncorr_color)
corr_clone.SetFillColor(corr_color)
uncorr_clone.SetFillStyle(1001)
corr_clone.SetFillStyle(1001)
uncorr_clone.SetLineColor(ROOT.kBlack)
corr_clone.SetLineColor(ROOT.kBlack)

runs1=file1.Get("Runs")
sumw1=0.0
for r1 in runs1: 
    sumw1+=getattr(r1, "genEventSumw")

runs2=file2.Get("Runs")
sumw2=0.0
for r2 in runs2: 
    sumw2+=getattr(r2, "genEventSumw")



BR=1e-4
lumi=59830
xsec=52.143

for ev in tree1:
    if not (ev.HLT_passed==1 and ev.best_4g_phi1_dxy_m30 >-20 and ev.best_4g_phi2_dxy_m30>-20 and ev.best_4g_ID_m30==1):
        continue

    event_weight=(BR)*(xsec)*(1/sumw1)*(ev.genWeight)*(ev.Pileup_weight)*(lumi)
    corr.Fill(ev.best_4g_corr_mass_m30, event_weight)
    
for ev in tree2:
    if not (ev.HLT_passed==1 and ev.best_4g_phi1_dxy_m30 >-20 and ev.best_4g_phi2_dxy_m30>-20 and ev.best_4g_ID_m30==1):
        continue

    event_weight=(BR)*(xsec)*(1/sumw2)*(ev.genWeight)*(ev.Pileup_weight)*(lumi)
    uncorr.Fill(ev.best_4g_corr_mass_m30, event_weight)
        
    



uncorr.SetLineColor(ROOT.kBlack)
corr.SetLineColor(ROOT.kBlack)
uncorr.SetLineWidth(2)
corr.SetLineWidth(2)
corr.Draw("HIST")
uncorr.Draw("HIST, SAME")
ymax=max(uncorr.GetMaximum(), corr.GetMaximum())
print(ymax*1.35)
corr.SetMaximum(0.7382053703069688)

#dataHist = ROOT.RooDataHist("dataHist", "dataHist", ROOT.RooArgList(x), corr)
#dataHist.plotOn(frame, ROOT.RooFit.DrawOption("HIST"),ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.FillColor(corr_color), ROOT.RooFit.FillStyle(3001))
#frame.addTH1(corr, "HIST")
#frame.SetYTitle("Events / 0.33 GeV")
#frame.SetXTitle("4#gamma Mass [GeV]")
#pdf.plotOn(frame, ROOT.RooFit.Normalization(5.9449, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.LineColor(ROOT.kRed+2), ROOT.RooFit.Name("model_plot"))
#frame.Draw()

c.Update()
l.AddEntry(corr_clone, "custom ID", "f")
l.AddEntry(uncorr_clone, "EGM ID", "f")
#l.AddEntry(frame.findObject("model_plot"), "fit", "l")
l.SetTextSize(0.05)
CMSstyle(c, l, ["DoublePhoton33_CaloIdL","H #rightarrow #phi#phi #rightarrow 4#gamma", "c#tau = 100 mm", "m_{#phi} = 30 GeV"])
#CMSstyle(c, l, ["TriplePhoton_20_20_20_CaloIdLV2","H #rightarrow #phi#phi #rightarrow 4#gamma", "c#tau = 100 mm", "m_{#phi} = 30 GeV"])
#CMSstyle(c, l, ["Background"])
l.Draw("SAME")
c.Update()
c.SaveAs("corr_vs_uncorr.png")
c.SaveAs("corr_vs_uncorr.pdf")

