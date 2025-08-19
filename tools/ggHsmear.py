import ROOT, math, random

def smear(file, treename, cat, year, xsec=56, BR=1e-4, bins=100, xmin=100, xmax=160):
#    cuts= {"lowlow":"", "":"", "":"", "":""}
    f = ROOT.TFile.Open(file)
    runs = f.Get("Runs")
    sumw = 0.0

    for r in runs:
        sumw += getattr(r, "genEventSumw")
    if sumw == 0:
        raise RuntimeError("sum(genEventSumw) is zero!")

    #print(sumw)
    norm = BR*xsec * 1 / sumw
    #print(norm)
    h_smear  = ROOT.TH1D("h_smear",  "m_{4γ} smeared; m [GeV]; events", bins, xmin, xmax)
    h_none = ROOT.TH1D("h_none", "m_{4γ} no smear; m [GeV]; events", bins, xmin, xmax)
    tree = f.Get(treename)
    i=0
    for ev in tree:
        if not (ev.HLT_passed==1 and ev.best_4g_ID_m30==1):
            continue

        m1, m2 = ev.best_4g_phi1_mass_m30, ev.best_4g_phi2_mass_m30
        lxy1, lxy2 = ev.best_4g_phi1_dxy_m30, ev.best_4g_phi2_dxy_m30

        if cat=="lowlow":
            passes = (m1>14 and m2>14 and lxy1<50 and lxy2<50)
        elif cat=="highlow":
            passes = (m1>14 and m2>14 and lxy1>50 and lxy2<50)
        elif cat=="lowhigh":
            passes = (m1>14 and m2>14 and lxy1<50 and lxy2>50)
        else:
            passes = (m1>14 and m2>14 and lxy1>50 and lxy2>50)

        if not passes:
            continue

        i+=1
    #    print("got here:  ",i)
        event_weight = ev.genWeight * norm * ev.Pileup_weight
        #print("genWeight: ", ev.genWeight)
        #print("PU: ", ev.Pileup_weight)
        m0 = ev.best_4g_corr_mass_m30
    #    print(m0)
        h_none.Fill(m0, event_weight)
        pts = [ev.best_4g_phi1_gamma1_pt_m30,
               ev.best_4g_phi1_gamma2_pt_m30,
               ev.best_4g_phi2_gamma1_pt_m30,
               ev.best_4g_phi2_gamma2_pt_m30]
        idxs = [int(ev.best_4g_idx1_m30),
                int(ev.best_4g_idx2_m30),
                int(ev.best_4g_idx3_m30),
                int(ev.best_4g_idx4_m30)]
        deltasUp = [ev.Photon_dEsigmaUp[i] for i in idxs]
        smeared_pts = [p * (1 + random.gauss(0,2)*d)
                       for p,d in zip(pts, deltasUp)]
        total = ROOT.TLorentzVector()
        for p, n, phi in zip(smeared_pts,
                          [ev.best_4g_phi1_gamma1_eta_m30,
                           ev.best_4g_phi1_gamma2_eta_m30,
                           ev.best_4g_phi2_gamma1_eta_m30,
                           ev.best_4g_phi2_gamma2_eta_m30],
                          [ev.best_4g_phi1_gamma1_phi_m30,
                           ev.best_4g_phi1_gamma2_phi_m30,
                           ev.best_4g_phi2_gamma1_phi_m30,
                           ev.best_4g_phi2_gamma2_phi_m30]):
            E = p*math.cosh(n)
            v = ROOT.TLorentzVector()
            v.SetPtEtaPhiE(p,n,phi, E)
            total += v

        h_smear.Fill(total.M(), event_weight)

    #c= ROOT.TCanvas("c","c",800,600)
    #h_none.SetLineColor(ROOT.kGreen)
    #h_smear.SetLineColor(ROOT.kRed)
    #dcb_fit.SetLineColor(ROOT.kBlue)
    #h_down.SetLineColor(ROOT.kBlue)
    #ymax = max(h_smear.GetMaximum(), h_none.GetMaximum())*1.1
    #h_none.SetStats(0)
    #h_none.SetMaximum(ymax)
    #h_none.Draw("HIST")
    #h_smear.Draw("HIST, SAME")
    #dcb_fit.Draw("SAME")
    #h_none.SetTitle("4 photon invariant mass distributions")
    #h_none.GetXaxis().SetTitle("4 photon mass")
    #h_none.GetYaxis().SetTitle("Events")
    #h_down.Draw("SAME")

    #leg = ROOT.TLegend(0.60,0.70,0.90,0.90)
    #leg.AddEntry(h_none,  "no smearing",    "l")
    #leg.AddEntry(h_smear,   "smearing up",  "l")
    #leg.AddEntry(h_down, "Smear Down",  "l")
    #leg.Draw()
    #c.Update()
    file=ROOT.TFile("smeared_signal_histograms_{}_{}.root".format(cat, year), "RECREATE")
    filename="smeared_signal_histograms_{}_{}.root".format(cat, year)
    file.cd()
    h_smear.Write("smeared_up")
    h_none.Write("no_smearing")
    #c.SaveAs("smeared_not_smeared_histo_{}_{}.png".format(cat, year))
    file.Write()
    file.Close()
    f.Close()

    return h_smear, h_none, filename



if __name__=="__main__":
    smear("ggH_M30_ctau0_ggH4g.root", "ggH4g","lowlow", "2018")
