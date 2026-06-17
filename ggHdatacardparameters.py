order_fit = 4
order_gen = 4

signal_window = (110, 140)
n_bins = 60
bins = [n_bins, signal_window[0], signal_window[1]]

lxy1 = 50
lxy2 = 50

smear_resolution = 0.264

##this weight comes from scaling PRESELECTED 4g to the sidebands of FULLY ID'd 4g data in 2018 EGamma data.
bkg_scale_factor = 15 / 7099
bkg_scale_factor_egm = 1 / 7099

##luminosities in pb^-1 from: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2
lumis = {"2016preVFP": [1/2*(36310)],
         "2016postVFP": [1/2*(36310)],
         "2017": [41480],
         "2018": [59830],
         "Run2": [137620],
         "2022": [34748],
         "2023": [27245],
         "2024": [108920],
         "Run3": [170857]}

#these dictionary values are scales for the background distribution according to what year is being processed.
#this needs to be done since hte 0.05124.... factor was derived from comparing sidebands of 2018 preselected vs ID data.
#if we want to interpolate to any other year having derived the initial scale factor from 2018, all years need their own additional weight normalized by 2018 lumi
bkg_factor = {"2016preVFP": [1/2*(36310/59830)],
              "2016postVFP": [1/2*(36210/59830)],
              "2017": [41480/59830],
              "2018": [1],
              "Run2": [137620/59830],
              "2022": [34748/59830],
              "2023": [27245/59830],
              "2024": [108920/59830],
              "Run3": [170857/59830]}

#######IMPORTANT#######
#since the cross section we use to scale the MC is both ggH+VBF combined, these xsec and lumi unc are derived from adding the unc from
#each process in quadrature

#obtained xsec and PDF from:https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#gluon_gluon_Fusion_Process
xsec_unc = {"ggH": [0.046, -0.067],
            "VBF": [0.004, -0.003]}

pdf_alphas_unc = {"ggH": [0.032, 0.032],
                  "VBF": [0.021, 0.021]}

#lumi unc obtained from: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2
lumi_unc = {"2016preVFP": [0.012],
            "2016postVFP": [0.012],
            "2017": [0.023],
            "2018": [0.025],
            "Run2": [0.016]}
#          "2022":[0.012],
#          "2023":[0.023],
#         "2024":[0.025],
#          "Run3":[0.02]}

dcb_mean = (125, 120, 130)
dcb_sigma = (2, 0.1, 20)
dcb_alpha1 = (2, 1, 20)
dcb_n1 = (2, 1, 50)
dcb_alpha2 = (2, 1, 20)
dcb_n2 = (2, 1, 50)

bernstein_coeff = (0, 100)
bernstein_coeff_card = (0, 50)
