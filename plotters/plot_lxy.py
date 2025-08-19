import ROOT
from ggHcmsstyle import CMSstyle
import numpy


ROOT.gStyle.SetPalette(103)
c=ROOT.TCanvas("", "", 800, 800)
#pad1=ROOT.TPad("", "", 0, 0.3, 1, 1)
#pad2=ROOT.TPad("", "", 0, 0, 1, 0.3)
right=0.05
left=0
up=0.05
down=0
l=ROOT.TLegend(0.5+right-left,0.67+up-down, 0.95+right-left, 0.8+up-down)
temps=[]
xmax=150
ymax=150
xmin=-20
ymin=-20
nbins=int((xmax-xmin)/2)

thres=50

signal="ggH_M30_ctau0_ggH4g.root"
sig_tree="ggH4g"
bkg="EGamma_2018_all_ggH4g.root"
bkg_tree="ggH4g"

sig_temps=[]
bkg_temps=[]

histo=ROOT.TH2F("hist", "", nbins, xmin, xmax, nbins, ymin, ymax)
signal_file=ROOT.TFile.Open(signal)

runs=signal_file.Get("Runs")
sumw=0.0
for r in runs:
    sumw+=getattr(r, "genEventSumw")
BR=1e-4
lumi=59830
xsec=52.143
#df=ROOT.RDataFrame(bkg_tree, bkg)
df=ROOT.RDataFrame(sig_tree, signal)
weight_formula=f"(genWeight / {sumw})*{xsec}*{BR}*Pileup_weight"
#df=df.Filter("HLT_passed==1 && Photon_preselection[best_4g_idx1_m30]==1 && Photon_preselection[best_4g_idx2_m30]==1 && Photon_preselection[best_4g_idx3_m30]==1 && Photon_preselection[best_4g_idx4_m30]==1 && best_4g_corr_mass_m30>110 && best_4g_corr_mass_m30<140 && best_4g_phi1_dxy_m30>-20 && best_4g_phi2_dxy_m30>-20")
df=df.Filter("HLT_passed==1 && best_4g_ID_m30==1 && best_4g_corr_mass_m30>110 && best_4g_corr_mass_m30<140 && best_4g_phi1_dxy_m30>-20 && best_4g_phi2_dxy_m30>-20").Define("event_weight", weight_formula)
histo=df.Histo2D(("", "", nbins, xmin, xmax, nbins, ymin, ymax), "best_4g_phi1_dxy_m30", "best_4g_phi2_dxy_m30", "event_weight")
#histo=df.Histo2D(("", "", nbins, xmin, xmax, nbins, ymin, ymax), "best_4g_phi1_dxy_m30", "best_4g_phi2_dxy_m30")

#histo.Scale(0.051249577845322525)
histo.Scale(lumi)
print(histo.Integral())
histo.GetXaxis().SetTitle("#phi_{1} L_{xy} (cm)        ")
histo.GetYaxis().SetTitle("#phi_{2} L_{xy} (cm)        ")
#histo.GetYaxis().SetTitle("   #phi_{2} L_{xy} (cm)")
histo.GetZaxis().SetTitle("Events")


c.Update()
ROOT.gStyle.SetPalette(103)
c.cd()
histo.Draw("COLZ")
#hline=ROOT.TLine(-20, -20, -20, ymax)
#hline.SetLineColor(ROOT.kRed-3)
#hline.SetLineWidth(3)
#hline.Draw("SAME")

#vline=ROOT.TLine(xmax, -20, -20, -20)
#vline.SetLineColor(ROOT.kRed-3)
#vline.SetLineWidth(3)
#vline.Draw("SAME")


hline2=ROOT.TLine(-20, thres, xmax, thres)
vline2=ROOT.TLine(thres, -20, thres, ymax)
hline2.SetLineColor(ROOT.kRed-7)
vline2.SetLineColor(ROOT.kRed-7)
hline2.SetLineWidth(3)
vline2.SetLineWidth(3)
hline2.Draw("SAME")
vline2.Draw("SAME")

#box1=ROOT.TBox(xmin, ymin, -20, ymax)
#box2=ROOT.TBox(-20, ymin, xmax, -20)
#box1.SetFillColor(ROOT.kRed-7)
#box2.SetFillColor(ROOT.kRed-7)
#box1.SetFillStyle(3004)
#box2.SetFillStyle(3004)
#box1.Draw("SAME")
#box2.Draw("SAME")

c.SetLeftMargin(0.08)
c.SetRightMargin(0.2)
CMSstyle(c, l, ["signal", "H #rightarrow #phi#phi #rightarrow 4#gamma", "c#tau = 100 mm", "m_{#phi} = 30 GeV"])
#CMSstyle(c, l, ["background"])
histo.SetStats(0)
#histo.GetXaxis().SetTitleOffset(2)
#histo.GetXaxis().SetLabelSize(0.02)
#histo.GetYaxis().SetLabelSize(0.02)
#histo.GetZaxis().SetLabelSize(0.02)
histo.GetZaxis().SetTitleSize(0.05)
histo.GetYaxis().SetTitleSize(0.04)
histo.GetXaxis().SetTitleSize(0.05)
histo.GetZaxis().SetTitleOffset(1.3)
histo.GetYaxis().SetTitleOffset(1.1)
histo.GetXaxis().SetTitleOffset(1.1)
c.Update()
c.SaveAs("lxy_heatmap.png")
c.SaveAs("lxy_heatmap.pdf")
input()

