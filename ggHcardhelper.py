import ROOT
import json
from DataCardMaker import DataCardMaker




#this function adds a double crystal ball to the RooWorkspace for the ggH shape analysis
def addDCB(self, name, variable, jsonFile, scale={}, resolution={}):
    #set the mass of the higgs to be a constant at 125 GeV in the workspace
    self.w.factory("MH[125.0]")
    self.w.var("MH").setConstant(1)
    scaleStr = '0'
    resolutionStr = '0'
    scaleSysts = []
    resolutionSysts = []

    #here we define nuisance parameters that are defined from the scale and resolutuion dictionaries
    #the key is the nuisance parameter name and the value is the actual value to the parameter 
    for syst, factor in scale.items():
        self.w.factory(f"{syst}[0,-0.1,0.1]")
        scaleStr += f"+{factor}*{syst}"
        scaleSysts.append(syst)
    for syst, factor in resolution.items():
        self.w.factory(f"{syst}[0,-5,5]")
        resolutionStr += f"+{factor}*{syst}"
        resolutionSysts.append(syst)
    self.w.factory(f"{variable}[100,200]")

    #load the optimized parameters from the JSON input file which are derived from the fitter in RooFit
    #model parameters derived from maximum likelihood analysis
    with open(jsonFile) as f:
        info = json.load(f)

        SCALEVar="_".join(["MEAN",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',{vv_systs})".format(
        name=SCALEVar,
        param=info['mean']['value'],
        vv_syst=scaleStr,
        vv_systs=','.join(scaleSysts)))

        SIGMAVar="_".join(["SIGMA",name,self.tag])
        self.w.factory("expr::{name}('({param})*(1+{vv_syst})',{vv_systs})".format(
        name=SIGMAVar,
        param=info['sigma']['value'],
        vv_syst=resolutionStr,
        vv_systs=','.join(resolutionSysts)))

        ALPHA1Var="_".join(["ALPHA1",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(
        name=ALPHA1Var,
        param=info['alpha1']['value']))

        ALPHA2Var="_".join(["ALPHA2",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(
        name=ALPHA2Var,
        param=info['alpha2']['value']))

        N1Var="_".join(["N1",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(
        name=N1Var,
        param=info['n1']['value']))

        N2Var="_".join(["N2",name,self.tag])
        self.w.factory("expr::{name}('MH*0+{param}',MH)".format(
        name=N2Var,
        param=info['n2']['value']))        

        pdfName="_".join([name,self.tag])


    fourPhotonMass = ROOT.RooDoubleCB(pdfName, pdfName,
                                      self.w.var(variable),
                                      self.w.function(SCALEVar),
                                      self.w.function(SIGMAVar),
                                      self.w.function(ALPHA1Var),
                                      self.w.function(N1Var),
                                      self.w.function(ALPHA2Var),
                                      self.w.function(N2Var))
    self.w.Import(fourPhotonMass)
    #getattr(self.w, 'import')(fourPhotonMass)





#this function adds a bernstein polynomial of degree 3 to the Rooworkspace for the ggH shape analysis
def addBernstein(self, name, variable, jsonFile, order):
#    self.w.factory(f"{variable}[100,200]")

    with open(jsonFile) as f:
        params = json.load(f)    

    coeffNames=[]
    keys = []    
    for i in range(order):
        coeffNames.append(f"c_{i}_bkg")
        keys.append(f"c_{i}")    

    clist=ROOT.RooArgList()
    for coeffName, key in zip(coeffNames, keys):
        self.w.factory("{name}[{val}, 0, 50]".format(name=coeffName, val=params[key]['value']))
        #self.w.var(coeffName).setConstant(True)    
        clist.add(self.w.var(coeffName))
    
    pdfName = "_".join([name, self.tag])
    bernstein_pdf = ROOT.RooBernsteinFast(order)(pdfName,pdfName,self.w.var(variable),clist)
    self.w.Import(bernstein_pdf)
    #getattr(self.w, 'import')(bernstein_pdf)





#finalstate = "4gamma"
#category = "DDP"
#period = "13TeV"
#luminosity = 59830
#physics = "ggH"


#dcm = DataCardMaker(finalstate, category, period, luminosity, physics)
#dcm.importBinnedData("data_obs_lowlowLXY.root", "histo", ["mass"])
#addDCB(dcm, "signal", "mass", "signal_fit_params.json", resolution={"res_uncertainty_lowlowLXY":"0.1"})
#addBernstein(dcm, "background", "mass", "background_fit_params.json")

#dcm.addSystematic(name="QCD_ggH_highhighLXY", kind="lnN", values={"signal":1.0458}) 
#dcm.addSystematic(name="PDFalphaS_ggH_highhighLXY", kind="lnN", values={"signal":1.0383}) 
#dcm.addSystematic(name="res_uncertainty_lowlowLXY", kind="param", values=[0.0,1.0]) 

#dcm.addFixedYieldFromFile(name="signal", ID="0", filename="TH1D_4g_3g1b_scaled_lowlowLXY.root", histoName="hist_4good")
#dcm.addFixedYieldFromFile(name="background", ID="1", filename="TH1D_4g_3g1b_scaled_lowlowLXY.root", histoName="hist_3good1bad", lumi=False)
#dcm.makeCard()
#print("Datacard and workspace written for tag:", dcm.tag)



