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

    def function(self,name,function,dependents):
        self.w.factory("expr::"+name+"('"+function+"',"+','.join(dependents)+")")


    def bernstein(self,name = 'model',poi='x',order=1):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

        cList = ROOT.RooArgList()
        for i in range(0,order):
            self.w.factory("c_"+str(i)+"[0,100]")
            cList.add(self.w.var("c_"+str(i)))
        bernsteinPDF = ROOT.RooBernsteinFast(order)(name,name,self.w.var(poi),cList)
        getattr(self.w,'import')(bernsteinPDF,ROOT.RooFit.Rename(name))


    def bernsteinPlusGaus(self,name = 'model',poi='x',order=1):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        self.w.factory("RooGaussian::"+name+"G(x,c_0[-10,-15,15],c_1[3,0,10])")


        cList = ROOT.RooArgList()
        for i in range(3,order):
            self.w.factory("c_"+str(i)+"[0.1,0,100]")
            cList.add(self.w.var("c_"+str(i)))
        bernsteinPDF = ROOT.RooBernsteinFast(order)(name+"B",name,self.w.var(poi),cList)
        getattr(self.w,'import')(bernsteinPDF,ROOT.RooFit.Rename(name+"B"))

        self.w.factory("SUM::"+name+"(c_2[0.5,0,1]*"+name+"G,"+name+"B)")

    def expo(self,name = 'model',poi='x'):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        self.w.factory("RooExponential::"+name+"("+poi+",c_0[-1,-1000,0])")

    def gaus(self,name = 'model',poi='x'):
        self.w.factory("RooGaussian::"+name+"("+poi+",c_0[50,0,10000],c_1[30,0,10000])")

    def doubleCB(self,name = 'model',poi='x'):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

        self.w.factory("mean[125,120,130]")
        self.w.factory("sigma[2,0,20]")
        self.w.factory("alpha1[2,1,20]")
        self.w.factory("n1[2,0,20]")
        self.w.factory("alpha2[2,1,20]")
        self.w.factory("n2[2,0,20]")

        doubleCB = ROOT.RooDoubleCB(name,name,self.w.var(poi),self.w.var("mean"),self.w.var("sigma"),self.w.var("alpha1"),self.w.var("n1"),self.w.var("alpha2"),self.w.var("n2"))
        getattr(self.w,'import')(doubleCB,ROOT.RooFit.Rename(name))

    def backgroundFast(self):
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

        self.w.factory("M[1000,20000]")
        self.w.factory("m[25,175]")



        self.w.factory("alpha_0[7,0,100]")
        self.w.factory("alpha_1[0,-0.1,0.1]")
        self.w.factory("expr::alpha('alpha_0+alpha_1*M',alpha_0,alpha_1,M)")

        self.w.factory("beta_0[6,0,100]")
        self.w.factory("beta_1[0,-0.1,0.1]")
        self.w.factory("expr::beta('beta_0+beta_1*M',beta_0,beta_1,M)")       

        self.w.factory("gamma_0[1,0,100]")
        self.w.factory("gamma_1[0,-0.1,0.1]")
        self.w.factory("expr::gamma('gamma_0+gamma_1*M',gamma_0,gamma_1,M)")       
        self.w.factory("delta[1.4,0,1000]")


        cList = ROOT.RooArgList()
        cList.add(self.w.function("alpha"))
        cList.add(self.w.function("beta"))
        cList.add(self.w.function("gamma"))
        cList.add(self.w.var("delta"))


        softDrop = ROOT.RooBernsteinFast(4)("modelJJ","modelJJ",self.w.var('m'),cList)
        getattr(self.w,'import')(softDrop,ROOT.RooFit.Rename("modelJJ"))
        self.w.factory("p0[20,0,100]")
        self.w.factory("p1[0.5,0,100]")
        self.w.factory("p2[0.0001,0,10]")

        qcd = ROOT.RooQCDPdf("modelQ","",self.w.var("M"),self.w.var("p0"),self.w.var("p1"),self.w.var("p2"))
        getattr(self.w,'import')(qcd,ROOT.RooFit.Rename("modelQ"))

        self.w.factory("PROD::model(modelJJ|M,modelQ)")

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



    def importUnbinnedData(self,tree,name,variables,poi,accept):
        cList = ROOT.RooArgSet()
        for i,p in enumerate(poi):
            cList.add(self.w.var(p))
        
        data=ROOT.RooDataSet(name,name,cList)

        for event in tree:
            if not accept(event):
                continue
            for i,p in enumerate(poi):
                val =  getattr(event,variables[i])
                cList.find(p).setVal(val[0])

            data.add(cList)


        getattr(self.w,'import')(data,ROOT.RooFit.Rename(name))
            


    def fit(self, model="model", data="data", options=[], fitRange=None):
        opts = list(options)
        if fitRange:
            opts.append(ROOT.RooFit.Range(fitRange))
        self.w.pdf(model).fitTo(self.w.data(data), *opts)


    def fetch(self,var):
        return (self.w.var(var).getVal(), self.w.var(var).getError())

    def projection(self,model = "model",data="data",poi="x",filename="fit.root"):
        
        self.frame=self.w.var(poi).frame()
        self.w.data(data).plotOn(self.frame)
        self.w.pdf(model).plotOn(self.frame)
        self.c=ROOT.TCanvas("c","c")
        self.c.cd()
        self.frame.Draw()
        self.c.SaveAs(filename)
        return self.frame.chiSquare()

