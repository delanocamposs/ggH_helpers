import ROOT
import json

class Fitter(object):
    def __init__(self,poi = ['x']):
        self.cache=ROOT.TFile("cache.root","RECREATE")
        self.cache.cd()

        self.w=ROOT.RooWorkspace("w","w")
        self.dimensions = len(poi)
        self.poi=poi
        for v in poi:
            self.w.factory(v+"[100,160]")

    def setRange(self,name,poi,low,high):
        self.w.var(poi).setRange(name,low,high)

    def factory(self,expr):
        self.w.factory(expr)

    def bernstein(self,name = 'model',poi='x',order=1):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

        cList = ROOT.RooArgList()
        for i in range(0,order):
            self.w.factory("c_"+str(i)+"[0,100]")
            cList.add(self.w.var("c_"+str(i)))
        bernsteinPDF = ROOT.RooBernsteinFast(order)(name,name,self.w.var(poi),cList)
        getattr(self.w,'import')(bernsteinPDF)


    def doubleCB(self,name = 'model',poi='x'):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        self.w.factory("mean[125,120,130]")
        self.w.factory("sigma[2,0.1,20]")
        self.w.factory("alpha1[2,1,20]")
        self.w.factory("n1[2,1,50]")
        self.w.factory("alpha2[2,1,20]")
        self.w.factory("n2[2,1,50]")
        doubleCB = ROOT.RooDoubleCB(name,name,self.w.var(poi),self.w.var("mean"),self.w.var("sigma"),self.w.var("alpha1"),self.w.var("n1"),self.w.var("alpha2"),self.w.var("n2"))
        getattr(self.w,'import')(doubleCB)

    def importBinnedData(self,histogram,poi = ["x"],name = "data"):
        cList = ROOT.RooArgList()
        for i,p in enumerate(poi):
            cList.add(self.w.var(p))
            if i==0:
                axis=histogram.GetXaxis()
            elif i==1:
                axis=histogram.GetYaxis()
            elif i==2:
                axis=histogram.GetZaxis()
            else:
                print( 'Asking for more than 3 D . ROOT doesnt support that, use unbinned data instead')
                return
            mini=axis.GetXmin()
            maxi=axis.GetXmax()
            bins=axis.GetNbins()
            self.w.var(p).setMin(mini)
            self.w.var(p).setMax(maxi)
            self.w.var(p).setBins(bins)
        dataHist=ROOT.RooDataHist(name,name,cList,histogram)
        getattr(self.w,'import')(dataHist,ROOT.RooFit.Rename(name))



    def fit(self, model="model", data="data", options=[], fitRange=None, sumW2=False):
        opts = list(options)
        opts.append(ROOT.RooFit.Save(True))
        opts.append(ROOT.RooFit.PrintLevel(-1))
        opts.append(ROOT.RooFit.Verbose(False))
        if fitRange:
            opts.append(ROOT.RooFit.Range(fitRange))
        if sumW2:
            opts.append(ROOT.RooFit.SumW2Error(True))
        return self.w.pdf(model).fitTo(self.w.data(data), *opts)


    def projection(self,model = "model",data="data",poi="x",filename="fit.root"):
        
        self.frame=self.w.var(poi).frame()
        self.w.data(data).plotOn(self.frame)
        self.w.pdf(model).plotOn(self.frame)
        self.c=ROOT.TCanvas("c","c")
        self.c.cd()
        self.frame.Draw()
        self.c.SaveAs(filename)
        return self.frame.chiSquare()

    def DCBandBernstein(self, order, dcbname="s_model", bname="b_model", poi="x", model_name="sb_model"):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        self.w.factory("mean[125,120,130]")
        self.w.factory("sigma[2,0.1,20]")
        self.w.factory("alpha1[2,1,20]")
        self.w.factory("n1[2,1,50]")
        self.w.factory("alpha2[2,1,20]")
        self.w.factory("n2[2,1,50]")
        cList = ROOT.RooArgList()
        for i in range(0,order):
            self.w.factory("c_"+str(i)+"[0,100]")
            cList.add(self.w.var("c_"+str(i)))

        bernsteinPDF = ROOT.RooBernsteinFast(order)(bname,bname,self.w.var(poi),cList)
        doubleCB = ROOT.RooDoubleCB(dcbname,dcbname,self.w.var(poi),self.w.var("mean"),self.w.var("sigma"),self.w.var("alpha1"),self.w.var("n1"),self.w.var("alpha2"),self.w.var("n2"))
        getattr(self.w,'import')(bernsteinPDF)
        getattr(self.w,'import')(doubleCB)
        self.w.factory(f"SUM::{model_name}(nsig[0,0,10000]*{dcbname}, nbkg[100,0,100000]*{bname})")
        return





