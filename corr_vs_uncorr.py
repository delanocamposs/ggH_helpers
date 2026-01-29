import ROOT 
from ggHtools import define_weightMC
from ggHcmsstyle import CMSstyle



double_photon=False
PT_MIN=33
ROOT.gInterpreter.Declare(f"""
bool passPtCut(float a,float b,float c,float d) {{
    float lead=std::max({{a,b,c,d}});
    float sub=-1e9f;
    if(a!=lead) sub=a;
    if(b!=lead && b>sub) sub=b;
    if(c!=lead && c>sub) sub=c;
    if(d!=lead && d>sub) sub=d;
    return lead>={PT_MIN} && sub>={PT_MIN};
}}
struct Pair {{ float l; float s; }};
Pair getLeadSub(float a,float b,float c,float d) {{
    float pts[4]={{a,b,c,d}};
    std::sort(pts,pts+4,std::greater<float>());
    return {{pts[0], pts[1]}};
}}
""")



c=ROOT.TCanvas("", "", 800, 800)
right=0.02
left=0.14
up=0.04
down=0
l=ROOT.TLegend(0.5+right-left,0.67+up-down, 0.95+right-left, 0.8+up-down)
bins=[30, 110,140]
var="best_4g_corr_mass_m30"

file1="ggH_M30_ctau0_ggH4g.root"
file2="EGamma_2018_all_ggH4g.root"

file1_open=ROOT.TFile.Open(file1)
file2_open=ROOT.TFile.Open(file2)

tree1 = file1_open.Get("ggH4g")
tree2 = file2_open.Get("ggH4g")

BR=1e-4
xsec=52.143
lumi=59830
EGM_ID = " && ".join(f"((Photon_isScEtaEB[best_4g_idx{i}_m30]==1 && Photon_hoe[best_4g_idx{i}_m30]<0.04596 && Photon_sieie[best_4g_idx{i}_m30]<0.0106)||(Photon_isScEtaEE[best_4g_idx{i}_m30]==1 && Photon_hoe[best_4g_idx{i}_m30]<0.0590 && Photon_sieie[best_4g_idx{i}_m30]<0.0272))" for i in range(1, 5))
cut_string="HLT_passed==1&&best_4g_phi1_dxy_m30>-20&&best_4g_phi2_dxy_m30>-20"
preselection="(Photon_preselection[best_4g_idx1_m30]==1)&&(Photon_preselection[best_4g_idx2_m30]==1)&&(Photon_preselection[best_4g_idx3_m30]==1)&&(Photon_preselection[best_4g_idx4_m30]==1)"
blind="((best_4g_corr_mass_m30<110)||(best_4g_corr_mass_m30>140))"

sumw=0.0
custom_signal_df=ROOT.RDataFrame("ggH4g", file1)
EGM_signal_df=ROOT.RDataFrame("ggH4g", file1)
with file1_open as f:
    runs_tree = f.Get("Runs")
    for entry in runs_tree:
        sumw += entry.genEventSumw
weight_formula = f"(genWeight / {sumw}) * {xsec} * {BR} * Pileup_weight"
custom_signal_df=custom_signal_df.Define("event_weight", weight_formula)
EGM_signal_df=EGM_signal_df.Define("event_weight", weight_formula)

#custom_signal_df=define_weightMC(file1, "ggH4g", BR, xsec)
#EGM_signal_df=define_weightMC(file1, "ggH4g", BR, xsec)
custom_sideband_df=ROOT.RDataFrame("ggH4g", file2)
EGM_sideband_df=ROOT.RDataFrame("ggH4g", file2)

if double_photon:
    custom_signal_df=custom_signal_df.Define(f"top2_pT_over_{PT_MIN}", "passPtCut(best_4g_phi1_gamma1_pt_m30,best_4g_phi1_gamma2_pt_m30,best_4g_phi2_gamma1_pt_m30,best_4g_phi2_gamma2_pt_m30)").Define("pTpair", "getLeadSub(best_4g_phi1_gamma1_pt_m30,best_4g_phi1_gamma2_pt_m30,best_4g_phi2_gamma1_pt_m30,best_4g_phi2_gamma2_pt_m30)").Define("LeadpT", "pTpair.l").Define("subLeadpT", "pTpair.s")
    custom_signal_df = custom_signal_df.Filter(f"top2_pT_over_{PT_MIN}==1")

    EGM_signal_df=EGM_signal_df.Define(f"top2_pT_over_{PT_MIN}", "passPtCut(best_4g_phi1_gamma1_pt_m30,best_4g_phi1_gamma2_pt_m30,best_4g_phi2_gamma1_pt_m30,best_4g_phi2_gamma2_pt_m30)").Define("pTpair", "getLeadSub(best_4g_phi1_gamma1_pt_m30,best_4g_phi1_gamma2_pt_m30,best_4g_phi2_gamma1_pt_m30,best_4g_phi2_gamma2_pt_m30)").Define("LeadpT", "pTpair.l").Define("subLeadpT", "pTpair.s")
    EGM_signal_df = EGM_signal_df.Filter(f"top2_pT_over_{PT_MIN}==1")

    custom_sideband_df=custom_sideband_df.Define(f"top2_pT_over_{PT_MIN}", "passPtCut(best_4g_phi1_gamma1_pt_m30,best_4g_phi1_gamma2_pt_m30,best_4g_phi2_gamma1_pt_m30,best_4g_phi2_gamma2_pt_m30)").Define("pTpair", "getLeadSub(best_4g_phi1_gamma1_pt_m30,best_4g_phi1_gamma2_pt_m30,best_4g_phi2_gamma1_pt_m30,best_4g_phi2_gamma2_pt_m30)").Define("LeadpT", "pTpair.l").Define("subLeadpT", "pTpair.s")
    custom_sideband_df = custom_sideband_df.Filter(f"top2_pT_over_{PT_MIN}==1")

    EGM_sideband_df=EGM_sideband_df.Define(f"top2_pT_over_{PT_MIN}", "passPtCut(best_4g_phi1_gamma1_pt_m30,best_4g_phi1_gamma2_pt_m30,best_4g_phi2_gamma1_pt_m30,best_4g_phi2_gamma2_pt_m30)").Define("pTpair", "getLeadSub(best_4g_phi1_gamma1_pt_m30,best_4g_phi1_gamma2_pt_m30,best_4g_phi2_gamma1_pt_m30,best_4g_phi2_gamma2_pt_m30)").Define("LeadpT", "pTpair.l").Define("subLeadpT", "pTpair.s")
    EGM_sideband_df = EGM_sideband_df.Filter(f"top2_pT_over_{PT_MIN}==1")



custom_signal_df = custom_signal_df.Filter(f"{cut_string}&&best_4g_ID_m30==1&&best_4g_passBitMap_loose_iso_m30==1")
EGM_signal_df = EGM_signal_df.Filter(f"{cut_string}&&best_4g_ID_m30==1&&{EGM_ID}&&best_4g_passBitMap_loose_iso_m30==1")
custom_sideband_df = custom_sideband_df.Filter(f"{cut_string}&&best_4g_ID_m30==1&&{blind}")
#custom_sideband_df = custom_sideband_df.Filter(f"{cut_string}&&{preselection}")
EGM_sideband_df = EGM_sideband_df.Filter(f"{cut_string}&&{preselection}")
#EGM_sideband_df = EGM_sideband_df.Filter(f"{cut_string}&&{preselection}")

custom_signal_histo = custom_signal_df.Histo1D(("hist1_1", f"hist1_1;{var};Events", bins[0], bins[1], bins[2]), f"{var}", "event_weight")
EGM_signal_histo = EGM_signal_df.Histo1D(("hist1_2", f"hist1_2;{var};Events", bins[0], bins[1], bins[2]), f"{var}", "event_weight")
custom_signal_histo.Scale(lumi)
EGM_signal_histo.Scale(lumi)

custom_sideband_histo = custom_sideband_df.Histo1D(("hist2", f"hist2;{var};Events", bins[0], bins[1], bins[2]), f"{var}")
EGM_sideband_histo = EGM_sideband_df.Histo1D(("hist2", f"hist2;{var};Events", bins[0], bins[1], bins[2]), f"{var}")
#custom_sideband_histo.Scale(15/7099)
#EGM_sideband_histo.Scale(1/7099)
N_cus_bkg=custom_sideband_histo.Integral()
N_cus_sig=custom_signal_histo.Integral()
N_EGM_bkg=EGM_sideband_histo.Integral()
N_EGM_sig=EGM_signal_histo.Integral()
print("N custom bkg: ", N_cus_bkg)
print("N custom sig: ", N_cus_sig)
print("N EGM bkg: ", N_EGM_bkg)
print("N EGM sig: ", N_EGM_sig)

sideband_custom_color=ROOT.kOrange+8
#sideband_EGM_color=ROOT.kRed+2
signal_EGM_color=ROOT.kMagenta-10
signal_custom_color=ROOT.kAzure-4
sideband_EGM_color=ROOT.kBlue-8
#signal_EGM_color=ROOT.kOrange-4
custom_sideband_histo.SetFillStyle(1001)
custom_sideband_histo.SetFillColor(signal_EGM_color)
custom_sideband_histo.SetLineColor(ROOT.kBlack)
custom_sideband_histo.SetLineStyle(1)
custom_signal_histo.SetFillStyle(1001)
custom_signal_histo.SetFillColor(signal_custom_color)
custom_signal_histo.SetLineColor(ROOT.kBlack)
custom_signal_histo.SetLineStyle(1)
EGM_sideband_histo.SetFillStyle(1001)
EGM_sideband_histo.SetFillColor(sideband_EGM_color)
EGM_sideband_histo.SetLineColor(ROOT.kBlack)
EGM_sideband_histo.SetLineStyle(1)
EGM_signal_histo.SetFillStyle(1001)
EGM_signal_histo.SetFillColor(sideband_custom_color)
EGM_signal_histo.SetLineColor(ROOT.kBlack)
EGM_signal_histo.SetLineStyle(1)

EGM_sideband_histo.GetXaxis().SetTitle("4#gamma Mass [GeV]")
per_bin=round((bins[2]-bins[1])/bins[0], 2)
EGM_sideband_histo.GetYaxis().SetTitle(f"Events / {per_bin} GeV")

sideband_custom_clone=custom_sideband_histo.Clone("sideband_custom_clone")
sideband_EGM_clone=EGM_sideband_histo.Clone("uncorrbkg_clone")
signal_custom_clone=custom_signal_histo.Clone("signal_custom_clone")
signal_EGM_clone=EGM_signal_histo.Clone("signal_EGM_clone")
sideband_custom_clone.SetFillColor(signal_EGM_color)
signal_custom_clone.SetFillColor(signal_custom_color)
sideband_custom_clone.SetFillStyle(1001)
signal_custom_clone.SetFillStyle(1001)
sideband_custom_clone.SetLineColor(ROOT.kBlack)
signal_custom_clone.SetLineColor(ROOT.kBlack)
sideband_EGM_clone.SetFillColor(sideband_EGM_color)
signal_EGM_clone.SetFillColor(sideband_custom_color)
sideband_EGM_clone.SetFillStyle(1001)
signal_EGM_clone.SetFillStyle(1001)
sideband_EGM_clone.SetLineColor(ROOT.kBlack)
signal_EGM_clone.SetLineColor(ROOT.kBlack)

custom_sideband_histo.SetLineColor(ROOT.kBlack)
EGM_sideband_histo.SetLineColor(ROOT.kBlack)
custom_signal_histo.SetLineColor(ROOT.kBlack)
EGM_signal_histo.SetLineColor(ROOT.kBlack)
custom_sideband_histo.SetLineWidth(2)
EGM_sideband_histo.SetLineWidth(2)
custom_signal_histo.SetLineWidth(2)
EGM_signal_histo.SetLineWidth(2)
#custom_signal_histo.Draw("HIST")
#custom_sideband_histo.Draw("HIST")
#EGM_signal_histo.Draw("HIST, SAME")
EGM_sideband_histo.Scale(15/7099)
EGM_sideband_histo.Draw("HIST")
#custom_sideband_histo.Draw("HIST, SAME")
ymax=max(custom_sideband_histo.GetMaximum(), EGM_sideband_histo.GetMaximum())
ymax=ymax*1.3
print(ymax)
#custom_sideband_histo.SetMaximum(ymax)
EGM_sideband_histo.SetMinimum(0.1)

c.Update()
#l.AddEntry(signal_custom_clone, "signal: custom ID", "f")
#l.AddEntry(signal_EGM_clone, "signal: EGM ID", "f")
#l.AddEntry(sideband_custom_clone, "full ID", "f")
l.AddEntry(sideband_EGM_clone, "preselected (with SF)", "f")
l.SetTextSize(0.05)
#CMSstyle(c, l, ["DoublePhoton33_CaloIdL","H #rightarrow #phi#phi #rightarrow 4#gamma", "c#tau = 100 mm", "m_{#phi} = 30 GeV"])
#CMSstyle(c, l, ["TriplePhoton_20_20_20_CaloIdLV2","H #rightarrow #phi#phi #rightarrow 4#gamma", "c#tau = 100 mm", "m_{#phi} = 30 GeV"])
#CMSstyle(c, l, ["DoublePhoton33_CaloIdL", "Background"])
#CMSstyle(c, l, ["TriplePhoton_20_20_20_CaloIdLV2", "Background"])
CMSstyle(c, l, ["Preselected 2018 EGamma Data"])
#c.SetLogy()
l.Draw("SAME")
c.Update()
c.SaveAs("plot.png")
c.SaveAs("plot.pdf")

