import ROOT 
from ggHcmsstyle import CMSstyle


c=ROOT.TCanvas("", "", 800, 800)
right=0.05
left=0
up=0.05
down=0
l=ROOT.TLegend(0.5+right-left,0.67+up-down, 0.95+right-left, 0.8+up-down)
file=ROOT.TFile.Open("ggH_M30_ctau0_ggH4g.root")
tree=file.Get("ggH4g")
uncorr=ROOT.TH1F("", "", 100, 100, 160)
corr=ROOT.TH1F("", "", 100, 100, 160)
#uncorr_color=ROOT.kRed-3
#corr_color=ROOT.kGreen-6
uncorr_color=ROOT.kOrange+8
corr_color=ROOT.kAzure-4
uncorr.SetFillStyle(3001)
uncorr.SetFillColor(uncorr_color)
corr.SetFillStyle(3001)
corr.SetFillColor(corr_color)
uncorr.GetXaxis().SetTitle("4#gamma Mass [GeV]")
uncorr.GetYaxis().SetTitle("Events / 0.6 GeV")

uncorr_clone=uncorr.Clone("uncorr_clone")
corr_clone=corr.Clone("corr_clone")
uncorr_clone.SetFillColor(uncorr_color)
corr_clone.SetFillColor(corr_color)
uncorr_clone.SetFillStyle(1001)
corr_clone.SetFillStyle(1001)
uncorr_clone.SetLineColor(0)
corr_clone.SetLineColor(0)

runs=file.Get("Runs")
sumw=0.0
for r in runs: 
    sumw+=getattr(r, "genEventSumw")
BR=1e-4
lumi=59830
xsec=52.143
for ev in tree:
    if not (ev.HLT_passed==1 and ev.best_4g_phi1_mass_m30 >14 and ev.best_4g_phi2_mass_m30>14 and ev.best_4g_ID_m30==1):
        continue
    event_weight=(BR)*(xsec)*(1/sumw)*(ev.genWeight)*(ev.Pileup_weight)*(lumi)
    uncorr.Fill(ev.best_4g_uncorr_mass_m30, event_weight)
    corr.Fill(ev.best_4g_corr_mass_m30, event_weight)
    
        
    



uncorr.SetLineColor(uncorr_color)
corr.SetLineColor(corr_color)
uncorr.SetLineWidth(2)
corr.SetLineWidth(2)
uncorr.Draw("HIST")
corr.Draw("HIST, SAME")
ymax=max(uncorr.GetMaximum(), corr.GetMaximum())
uncorr.SetMaximum(ymax*1.35)

c.Update()
l.AddEntry(uncorr_clone, "uncorrected", "f")
l.AddEntry(corr_clone, "corrected", "f")
l.SetTextSize(0.05)
CMSstyle(c, l, ["H #rightarrow #phi#phi #rightarrow 4#gamma", "c#tau = 100 mm", "m_{#phi} = 30 GeV"])
l.Draw("SAME")
c.Update()
c.SaveAs("corr_vs_uncorr.png")
c.SaveAs("corr_vs_uncorr.pdf")

