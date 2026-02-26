from Fitter import *
ROOT.gROOT.SetBatch(False)


def fitBKG(file, hist, output_name, order=3, POI="mass", verbose=False):
    f = ROOT.TFile(file)
    bkg = f.Get(hist)
    #fit the background with a Berntsein Polynomial of order 3 
    fitter = Fitter([POI])
    fitter.setRange("lower", "mass", 70,110)
    fitter.setRange("upper", "mass", 140,180)
    fitter.bernstein('model',POI,order=order)
    fitter.importBinnedData(bkg,[POI])
    fitter.fit("model","data",fitRange="lower,upper")
    chi2 = fitter.projection("model","data",POI,filename=output_name)
    output = ROOT.TFile(output_name, "UPDATE")
    output.cd()
    fitter.w.Write("w")
    output.Close()
    print
    if verbose: 
        print("bkg chi-squared={}".format(chi2))

def fitSIG(file, hist, output_name, POI="mass", verbose=False):
    f = ROOT.TFile(file)
    signal = f.Get(hist)
    #fit the signal with a Double Crystalball
    fitterS = Fitter([POI])
    fitterS.doubleCB('model',POI)
    fitterS.importBinnedData(signal,[POI],name = "data")
    #fit twice
    fitterS.fit("model","data")
    chi2 = fitterS.projection("model","data",POI,filename=output_name)
    if verbose:
        print("signal chi-squared={}".format(chi2))
    output = ROOT.TFile(output_name, "UPDATE")
    output.cd()
    fitterS.w.Write("w")
    output.Close()


