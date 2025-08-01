import math
import ROOT
from ggHcmsstyle import CMSstyle
import CMS_lumi
import tdrstyle
files = [
    ("r_vals_2018_orderfit-2_ordergen-3.txt", "order=2"),
    ("r_vals_2018_orderfit-3_ordergen-3.txt", "order=3")
    #("r_vals_4_100.txt", "order=4"),
    #("r_vals_5_100.txt", "order=5"), 
    #("r_vals_6_100.txt", "order=6")
]
colors = [ROOT.kOrange-2, ROOT.kAzure+1, ROOT.kGreen-3, ROOT.kRed+1] 
histograms = []
for idx, (filename, label) in enumerate(files):
    h = ROOT.TH1F(f"h_ratio_{idx}", r";\frac{r}{\Delta r};# toys", 50, -10, 10)
    with open(filename, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split(",")
            if len(parts) != 3:
                continue
            r_val = float(parts[0])
            r_unc_up = float(parts[1])
            r_unc_down = float(parts[2])
            delta_r = (abs(r_unc_up) + abs(r_unc_down)) / 2.0
            if delta_r > 0:
                ratio = r_val / delta_r
                h.Fill(ratio)
    h.SetLineColor(colors[idx])
    h.SetLineWidth(2)
    histograms.append((h, label))


c = ROOT.TCanvas("c", "r_over_deltaR", 800, 800)
ymax=max(histograms[0][0].GetMaximum(), histograms[1][0].GetMaximum())
#ymax=histograms[0][0].GetMaximum()
histograms[0][0].SetMaximum(1.2*ymax)
histograms[0][0].Draw("HIST")


leg=ROOT.TLegend(0.3, 0.55, 0.8, 0.75)
for h, label in histograms[1:]:
    h.Draw("HIST SAME")

#graphs=[]
#for idx, (h, label) in enumerate(histograms):
#    graph=ROOT.TGraph()
#    n_points=0
#    for i in range(1, h.GetNbinsX()+1):
#        x=h.GetBinCenter(i)
#        y=h.GetBinContent(i)
#        graph.SetPoint(n_points, x, y)
#        n_points+=1
#    graph.SetMarkerStyle(20)
#    graph.SetMarkerColor(colors[idx])
#    graph.Draw("P, SAME")
#    graphs.append(graph)
#    c.Update()


h1clone=histograms[0][0].Clone("h1fill")
h2clone=histograms[1][0].Clone("h2fill")
#h3clone=histograms[2][0].Clone("h3fill")
#h4clone=histograms[3][0].Clone("h3fill")

h1clone.SetFillStyle(1001)
h1clone.SetFillColor(colors[0])
h2clone.SetFillStyle(1001)
h2clone.SetFillColor(colors[1])
#h3clone.SetFillStyle(1001)
#h3clone.SetFillColor(colors[2])
#h4clone.SetFillStyle(1001)
#h4clone.SetFillColor(colors[3])

mean1=histograms[0][0].GetMean()
mean2=histograms[1][0].GetMean()
#mean3=histograms[2][0].GetMean()
#mean4=histograms[3][0].GetMean()

rms1=histograms[0][0].GetRMS()
rms2=histograms[1][0].GetRMS()
#rms3=histograms[2][0].GetRMS()
#rms4=histograms[3][0].GetRMS()
hi1=round(mean1/rms1, 2)
hi2=round(mean2/rms2, 2)
#hi3=round(mean3/rms3, 2)
#hi4=round(mean4/rms4, 2)

leg.AddEntry(h1clone,f"order=2, mean/RMS={hi1}", "f")
leg.AddEntry(h2clone,f"order=3, mean/RMS={hi2}", "f")
#leg.AddEntry(h3clone,f"order=5, mean/RMS={hi3}", "f")
#leg.AddEntry(h4clone,f"order=6, mean/RMS={hi4}", "f")
leg.SetTextSize(0.03)
leg.Draw()
        
histograms[0][0].GetXaxis().SetTitleSize(0.05)
histograms[0][0].GetYaxis().SetTitleSize(0.05)


CMSstyle(c, leg, ["toys=100", "r=0 injected data"])


c.Update()
c.SaveAs("r_over_deltaR_CMS_multi.png")
c.SaveAs("r_over_deltaR_CMS_multi.pdf")
input()
