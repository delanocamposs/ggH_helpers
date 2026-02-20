from Fitter import *
ROOT.gROOT.SetBatch(False)


def fitBKG(file, hist, output_name, order=3, POI="mass", verbose=False):
    f = ROOT.TFile(file)
    bkg = f.Get(hist)
    #fit the background with a Berntsein Polynomial of order 3 
    fitter = Fitter([POI])
    fitter.bernstein('model',POI,order=order)
    fitter.importBinnedData(bkg,[POI])
    fitter.fit("model","data")
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


if __name__=="__main__":
    fitBKG("rate_histos_lowlow_2018.root", "hist2", "bkgFit_scaled_lowlow.root")
    #fitSIG("smeared_signal_histograms.root", "smeared_up", "smeared_up_fit.root", POI="mass", verbose=True)
#    fitSIG("smeared_signal_histograms.root", "no_smearing", "no_smearing_fit.root", POI="mass", verbose=True)
#    fitSIG("output_lowlow.root", "hist1", "sigFit_scaled_lowlow.root")

    
 
