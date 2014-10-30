#!/usr/bin/env python

import ROOT
import numpy as np
import os
import re
import sys
import math


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

def add_histograms_all(analysis, params, smear=0.):
    # Quadratic slope = lastbinwidth/firstbinwidth
    # ROOT.LinearHistogram(file_name, hist_name, n_bins, obs_min, obs_max)
    # ROOT.QuadraticHistogram(file_name, hist_name, n_bins, obs_min, obs_max, slope)
    # ROOT.SmearedLinearHistogram(file_name, hist_name, n_bins, obs_min, obs_max, smear_fac in [0, 1])
    # ROOT.SmearedQuadraticHistogram(file_name, hist_name, n_bins, obs_min, obs_max, slope, smear_fac in [0, 1])

    # List histogram takes a vector of bin edges
    # bin_edges = ROOT.std.vector('double')()
    # map(bin_edges.push_back, [1, 2, 3, 4, 5, 6])  # 5 bins in [1, 6)
    # ROOT.ListHistogram(file_name, hist_name, bin_edges)

    name = "l_s%g" % (smear,)
    filename = (params.output % name) + '.hist'

    def Histogram(*args):
        print "Histogram", args
        if smear == 0.:
            return ROOT.LinearHistogram(*args)
        else:
            # add smearing parameter
            newargs = args + (smear,)
            return ROOT.SmearedLinearHistogram(*newargs)

    histdefs = [
        ["higgs_pt"             , ("higgs_pt_60", 60, 0, 300)],
        ["higgs_pt"             , ("higgs_pt_200", 200, 0, 500)],
        ["higgs_pt"             , ("higgs_pt_100", 100, 0, 500)],
        ["higgs_pt"             , ("higgs_pt_50", 50, 0, 500)],
        ["higgs_eta"            , ("higgs_eta_lh", 20, -5, 5)],
        ["higgs_eta"            , ("higgs_eta_20", 20, -4.4, 4.4)],
        ["higgs_eta"            , ("higgs_eta_40", 40, -4.4, 4.4)],
        ["higgs_y"              , ("higgs_y_lh", 20, -5, 5)],
        ["higgs_y"              , ("higgs_y_20", 20, -4.4, 4.4)],
        ["higgs_y"              , ("higgs_y_40", 40, -4.4, 4.4)],
        ["jjj_ystar"            , ("jjj_ystar_20", 20, -5, 5)],
        ["jet_pt_n"             , ("jet_pT_%d_60", 60, 0, 300)],
        ["jet_pt_n"             , ("jet_pT_%d_200", 200, 0, 500)],
        ["jet_pt_n"             , ("jet_pT_%d_100", 100, 0, 500)],
        ["jet_pt_n"             , ("jet_pT_%d_50", 50, 0, 500)],
        ["jet_eta_n"            , ("jet_eta_%d_lh", 20, -5, 5)],
        ["jet_eta_n"            , ("jet_eta_%d_20", 20, -4.4, 4.4)],
        ["jet_eta_n"            , ("jet_eta_%d_40", 40, -4.4, 4.4)],
        ["jet_y_n"              , ("jet_y_%d_lh", 20, -5, 5)],
        ["jet_y_n"              , ("jet_y_%d_20", 20, -4.4, 4.4)],
        ["jet_y_n"              , ("jet_y_%d_40", 40, -4.4, 4.4)],
        ["jet_jet_mass_ij"      , ("jet_jet_mass_%d%d_25", 25, 0, 1000)],
        ["jet_jet_mass_ij"      , ("jet_jet_mass_%d%d_100", 100, 0, 1000)],
        ["jet_jet_dy_ij"        , ("jet_jet_dy_%d%d_20", 20, 0, 10)],
        ["jet_jet_dphi_ij"      , ("jet_jet_dphi_%d%d_20", 20, 0, math.pi+1e-10)],
        ["jet_jet_dR_ij"        , ("jet_jet_dR_%d%d_10", 10, 0, 5)],
        ["jet_jet_dR_ij"        , ("jet_jet_dR_%d%d_50", 50, 0, 5)],
        ["higgs_dijet_pt_ij"    , ("higgs_dijet_pt_%d%d_60", 60, 0, 300)],
        ["higgs_dijet_dphi_ij"  , ("higgs_dijet_dphi_%d%d_20", 20, 0, math.pi+1e-10)],
        ["higgs_dijet_ystar_ij" , ("higgs_dijet_ystar_%d%d_20", 20, 0, 10)],
    ]

    for hname, hparam in sorted(histdefs):
        histparam = hparam[1:]
        assert len(histparam) >= 3
        if '%d%d' in hparam[0]:
            for j in range(1, params.njet+1):
                for i in range(j):
                    histname = hparam[0] % (i+1, j+1)
                    n = (j-1)*j/2 + i
                    getattr(analysis, hname)[n].push_back(Histogram(filename, histname, *histparam))
        elif '%d' in hparam[0]:
            for i in range(params.njet+1):
                histname = hparam[0] % (i+1)
                getattr(analysis, hname)[i].push_back(Histogram(filename, histname, *histparam))
        else:
            histname = hparam[0]
            getattr(analysis, hname).push_back(Histogram(filename, histname, *histparam))

    return

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# main function which is called from hammer.py
def initialize(params, selector):
    print "Using higgs analysis %s" % __file__

    analysis = ROOT.HiggsJetsAnalysis.create()
    analysis.runname = params.runname
    analysis.setJetNumber(params.njet)
    analysis.setAntiKt(0.4)
    analysis.jet_ptmin = 30
    analysis.jet_etamax = 4.4
    #analysis.min_dijet_m = 400
    #analysis.min_dijet_y = 2.8
    selector.opt_alphas_ignore = 2  # two powers of alphaS are not reweighted
    selector.opt_ignore_scale = 125.  # two powers of alphaS are not reweighted

    # Extract smear value from the output pattern
    smearpat = r'-smear(\d+\.?\d*)-'
    m = re.match(r".*?%s.*?" % smearpat, params.output)
    if m:
        s = float(m.group(1))
        params.output = re.sub(smearpat, '-', params.output)
    else:
        s = 0.

    # Add analysis histograms (see the function above)
    add_histograms_all(analysis, params, smear=s)

    # assign to selector
    selector.analysis = analysis


if __name__ == '__main__':
    print "Analysis module can't be run alone"
    print "pass --analysis=%s to hammer.py" % __file__
