import ROOT
import subprocess
import json
import os
from datacard.ggHdatacardmaker import main as make_datacard
from run_datacard import combined_datacard
from ggHparameters import signal_path, bkg_path, lumi
from plotting.style import tdrstyle


# AsymptoticLimits stores one tree entry per quantile. Map the quantileExpected
# value combine writes to the keys we use everywhere downstream.
_QUANTILE_KEYS = {-1.0: "obs", 0.025: "exp-2", 0.16: "exp-1",
                  0.5: "exp0", 0.84: "exp+1", 0.975: "exp+2"}


def get_limits(root_file):
    """Return all 95% CL upper limits on r from an AsymptoticLimits root file.

    Returns a dict like {"obs": .., "exp-2": .., "exp-1": .., "exp0": ..,
    "exp+1": .., "exp+2": ..}. Missing quantiles are simply absent.
    """
    f = ROOT.TFile.Open(root_file)
    tree = f.Get("limit") if f else None
    out = {}
    if not tree:
        print(f"[get_limits] WARNING: no 'limit' tree in {root_file}; skipping point.")
        if f:
            f.Close()
        return out
    for entry in tree:
        for q, key in _QUANTILE_KEYS.items():
            if abs(entry.quantileExpected - q) < 1e-3:
                out[key] = entry.limit
    f.Close()
    return out


def save_results(results, filename):
    with open(filename, "w") as fout:
        json.dump(results, fout, indent=2)
    print(f"[save_results] wrote {filename}")


def load_results(filename):
    with open(filename) as fin:
        return json.load(fin)


def scan_mass_lifetime(masses, lifetimes, years, categories, bins, finalstate="4g",
                       physics="ggH", order_fit=4, results_json="limits_UL_vs_mass.json"):
    """Build/combine datacards and run AsymptoticLimits over the mass/lifetime grid.

    Returns (and saves to JSON) a nested dict:
        results[year][str(ctau)][str(mass)] = {"obs":.., "exp0":.., ...}   # raw r
    Raw r limits are stored; BR scaling is applied at plot time via br_scale.
    """
    ROOT.gROOT.SetBatch(True)
    results = {y: {} for y in years}
    for mass in masses:
        for ctau in lifetimes:
            for year in years:
                sig = signal_path(mass, ctau, year)
                bkg = bkg_path(year)
                for cat in categories:
                    make_datacard(paths=[sig, bkg], isMC=[1, 0], trees=["ggH4g", "ggH4g"],
                                  var=f"best_4g_corr_mass_m{mass}", categories=[cat], period=year,
                                  bins=bins, lifetime=ctau, mass=mass, finalstate=finalstate,
                                  physics=physics, order_fit=order_fit)

                combined_datacard(year, categories, mass, ctau, finalstate, physics)

                cats_combined = [c for c in categories if c != "none"]
                combined_tag = f"combined_{'_'.join(cats_combined)}"
                combined_txt = f"datacard_{physics}_{finalstate}_m{mass}_ct{ctau}_{combined_tag}_{year}.txt"
                combined_root = f"datacard_{physics}_{finalstate}_m{mass}_ct{ctau}_{combined_tag}_{year}.root"
                subprocess.run(["text2workspace.py", combined_txt, "-o", combined_root])
                subprocess.run(["combine", "-M", "AsymptoticLimits", combined_root, "-m", "125"])

                result_root = f"higgsCombineTest_m{mass}_ct{ctau}_{year}.AsymptoticLimits.mH125.root"
                subprocess.run(["mv", "higgsCombineTest.AsymptoticLimits.mH125.root", result_root])

                lims = get_limits(result_root)
                if lims:
                    results[year].setdefault(str(ctau), {})[str(mass)] = lims
                print(f"m={mass} ct={ctau} year={year}: expected UL on r = {lims.get('exp0')}")

    save_results(results, results_json)
    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _band_graph(masses, lo, hi, color):
    """Closed-polygon TGraph for a Brazil band between lo[] and hi[] vs masses[]."""
    n = len(masses)
    g = ROOT.TGraph(2 * n)
    for i in range(n):
        g.SetPoint(i, masses[i], hi[i])
    for i in range(n):
        g.SetPoint(n + i, masses[n - 1 - i], lo[n - 1 - i])
    g.SetFillColor(color)
    g.SetLineColor(color)
    return g


def _line_graph(masses, vals, style, width=2):
    g = ROOT.TGraph(len(masses))
    for i, (x, y) in enumerate(zip(masses, vals)):
        g.SetPoint(i, x, y)
    g.SetLineColor(ROOT.kBlack)
    g.SetLineStyle(style)
    g.SetLineWidth(width)
    g.SetMarkerSize(0)
    return g


def _panel_arrays(panel, blind):
    """From results[year][ctau] (mass->limits dict) build sorted arrays.

    Keeps only masses that have the full expected band. Returns
    (masses, dict_of_arrays). 'obs' present only if available and not blind.
    """
    masses = sorted(int(m) for m in panel)
    needed = ["exp-2", "exp-1", "exp0", "exp+1", "exp+2"]
    masses = [m for m in masses if all(k in panel[str(m)] for k in needed)]
    arrays = {k: [panel[str(m)][k] for m in masses] for k in needed}
    if not blind and masses and all("obs" in panel[str(m)] for m in masses):
        arrays["obs"] = [panel[str(m)]["obs"] for m in masses]
    return masses, arrays


def plot_UL_vs_mass(results, year, ctaus, br_scale=1.0, blind=False,
                    ytitle="95% upper limit on BR(H#rightarrow#Phi#Phi)",
                    extra_labels=("BR(#Phi#rightarrow#gamma#gamma) = 0.5",),
                    yrange=None, outname=None):
    """Multi-panel 95% CL UL vs m_Phi, one panel per lifetime (Brazil bands).

    results   : nested dict from scan_mass_lifetime / load_results
    year      : which year key to plot
    ctaus     : list of lifetimes (panel order, left->right)
    br_scale  : multiply raw r limits by this to get the plotted quantity (BR)
    """
    ROOT.gROOT.SetBatch(True)
    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetOptStat(0)
    # lighter, less distracting grid lines
    ROOT.gStyle.SetGridColor(ROOT.kGray)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridWidth(1)

    year_res = results[year] if year in results else results
    ctaus = [ct for ct in ctaus if str(ct) in year_res]
    n = len(ctaus)
    if n == 0:
        print("[plot_UL_vs_mass] no lifetimes with data to plot.")
        return

    # gather scaled arrays per panel + global ranges
    panels, all_y, all_m = [], [], []
    for ct in ctaus:
        masses, arrays = _panel_arrays(year_res[str(ct)], blind)
        arrays = {k: [v * br_scale for v in vals] for k, vals in arrays.items()}
        panels.append((ct, masses, arrays))
        for vals in arrays.values():
            all_y += vals
        all_m += masses

    if not all_y:
        print("[plot_UL_vs_mass] no complete points to plot.")
        return

    ymin,ymax =(min(all_y) * 0.5, max(all_y) * 2.0) if yrange is None else yrange
    xlo,xhi= min(all_m), max(all_m)
    lm,rm,tm,bm = 0.17, 0.035, 0.07, 0.14
    denom = (1.0 / (1 - lm)) + max(n - 2, 0) + (1.0 / (1 - rm) if n > 1 else 0)
    if n == 1:
        denom = 1.0 / (1 - lm - rm)
    fw = 1.0 / denom
    widths = []
    for i in range(n):
        if n == 1:
            widths.append(fw / (1 - lm - rm))
        elif i == 0:
            widths.append(fw / (1 - lm))
        elif i == n - 1:
            widths.append(fw / (1 - rm))
        else:
            widths.append(fw)
    edges = [0.0]
    for w in widths:
        edges.append(edges[-1] + w)
    edges = [e / edges[-1] for e in edges]

    canv = ROOT.TCanvas("UL_vs_mass", "UL_vs_mass", 260 * max(n, 4), 600)
    keep = []
    green, yellow = ROOT.kGreen + 1, ROOT.kOrange

    for i, (ct, masses, arrays) in enumerate(panels):
        canv.cd()
        pad = ROOT.TPad(f"pad{i}", "", edges[i], 0.0, edges[i + 1], 1.0)
        pad.SetLeftMargin(lm if (i == 0 or n == 1) else 0.0)
        pad.SetRightMargin(rm if (i == n - 1 or n == 1) else 0.0)
        pad.SetTopMargin(tm)
        pad.SetBottomMargin(bm)
        pad.SetLogy()
        pad.SetTicks(1, 1)
        pad.SetGridx()
        pad.SetGridy()
        pad.Draw()
        keep.append(pad)
        pad.cd()

        frame = pad.DrawFrame(xlo, ymin, xhi, ymax)
        frame.GetXaxis().SetNdivisions(505)
        frame.GetXaxis().SetTitleSize(0.06)
        frame.GetXaxis().SetLabelSize(0.055)
        frame.GetXaxis().SetLabelOffset(-0.005)
        frame.GetYaxis().SetTitle(ytitle if i == 0 else "")
        frame.GetYaxis().SetTitleSize(0.055)
        frame.GetYaxis().SetTitleOffset(1.5)
        if i != 0 and n > 1:
            frame.GetYaxis().SetLabelSize(0.0)
        else:
            frame.GetYaxis().SetLabelSize(0.045)
        keep.append(frame)

        if masses:
            g2 = _band_graph(masses, arrays["exp-2"], arrays["exp+2"], yellow)
            g1 = _band_graph(masses, arrays["exp-1"], arrays["exp+1"], green)
            g2.Draw("F same")
            g1.Draw("F same")
            gexp = _line_graph(masses, arrays["exp0"], 2)  # dashed expected
            gexp.Draw("L same")
            keep += [g2, g1, gexp]
            gobs = None
            if "obs" in arrays:
                gobs = _line_graph(masses, arrays["obs"], 1)  # solid observed
                gobs.Draw("L same")
                keep.append(gobs)

        pad.RedrawAxis()
        pad.RedrawAxis("g")  # grid lines on top of the bands

        cx = pad.GetLeftMargin() + 0.5 * (1 - pad.GetLeftMargin() - pad.GetRightMargin())
        tl = ROOT.TLatex()
        tl.SetNDC()
        tl.SetTextFont(42)
        tl.SetTextAlign(22)
        tl.SetTextSize(0.07)
        tl.DrawLatex(cx, bm + 0.04, f"c#tau = {ct} mm")
        keep.append(tl)

        # legend + extra labels only in the first panel
        if i == 0 and masses:
            leg = ROOT.TLegend(lm + 0.04, 0.65, lm + 0.92, 0.90)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetTextSize(0.06)
            leg.AddEntry(g1, "#pm1#sigma", "f")
            leg.AddEntry(g2, "#pm2#sigma", "f")
            leg.AddEntry(gexp, "Expected", "l")
            if gobs is not None:
                leg.AddEntry(gobs, "Observed", "l")
            leg.Draw()
            keep.append(leg)

            lab = ROOT.TLatex()
            lab.SetNDC()
            lab.SetTextFont(42)
            lab.SetTextSize(0.05)
            for j, txt in enumerate(extra_labels):
                lab.DrawLatex(lm + 0.04, 0.60 - 0.06 * j, txt)
            keep.append(lab)

    # header drawn on the canvas, spanning full width
    canv.cd()
    head = ROOT.TLatex()
    head.SetNDC()
    head.SetTextAlign(13)
    cms_size = 0.05
    head.SetTextFont(61)
    head.SetTextSize(cms_size)
    head.DrawLatex(0.04, 0.985, "CMS")
    # place "Preliminary" right after "CMS": its NDC width scales with the (wide)
    # canvas aspect ratio, so a fixed NDC offset would leave a large pixel gap.
    cms_w = cms_size * (canv.GetWh() / canv.GetWw()) * 1.8
    head.SetTextFont(52)
    head.SetTextSize(0.038)
    head.DrawLatex(0.04 + cms_w + 0.006, 0.983, "Preliminary")
    head.SetTextFont(42)
    head.SetTextAlign(33)
    head.SetTextSize(0.038)
    head.DrawLatex(1 - rm, 0.985, f"{lumi[year] / 1000:.2f} fb^{{-1}} (13 TeV)")
    # global x-axis title near the bottom-right
    head.SetTextAlign(31)
    head.SetTextSize(0.045)
    head.DrawLatex(1 - rm, 0.035, "m_{#Phi} (GeV)")
    keep.append(head)

    canv.Update()
    if outname is None:
        outname = f"UL_vs_mass_{year}"
    canv.SaveAs(f"{outname}.png")
    canv.SaveAs(f"{outname}.pdf")
    return canv


if __name__ == "__main__":
    masses = [15, 20, 30, 40, 50, 55]
    lifetimes = [0, 10, 20, 50, 100, 1000]
    years = ["2018"]
    categories = ["prompt", "asym", "displaced"]

    results_json = "limits_UL_vs_mass.json"

    # Re-use existing limits if present; otherwise run the (slow) combine scan.
    if os.path.exists(results_json):
        print(f"loading cached limits from {results_json} (delete it to re-scan)")
        results = load_results(results_json)
    else:
        results = scan_mass_lifetime(masses, lifetimes, years, categories,
                                     bins=[30, 110, 140], results_json=results_json)

    for year in years:
        # br_scale: r->BR conversion applied to every y value (= reference BR=1e-4).
        plot_UL_vs_mass(results, year, lifetimes, br_scale=1e-4)
