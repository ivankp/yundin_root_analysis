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

# Grid helper
def create_grid(name, edges, params):
    if not params.grids:
        return None
    if params.warmup and os.path.exists(name):
        bakname = name + '.bak'
        print "WARNING: warmup file exists, moving '%s' to '%s'" % (name, bakname)
        os.rename(name, bakname)
    if not os.path.exists(name):
        if not params.warmup:
            print "ERROR: No --warmup option, but grid '%s' is not found" % name
            sys.exit(2)
        obs = ROOT.std.vector('double')()
        for ed in edges:
            obs.push_back(ed)
        grid = ROOT.Grid.capture(ROOT.Grid(name, obs))
    else:
        grid = ROOT.Grid.capture(ROOT.Grid(name))
    if not grid.isWarmup():
        assert not params.warmup
        grid.setFilename(name.replace('.root', '-o.root'))
    return grid

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

def add_histograms_all(analysis, params, smear=[0.]):
    # Quadratic slope = lastbinwidth/firstbinwidth
    # ROOT.LinearHistogram(file_name, hist_name, n_bins, obs_min, obs_max)
    # ROOT.QuadraticHistogram(file_name, hist_name, n_bins, obs_min, obs_max, slope)
    # ROOT.SmearedLinearHistogram(file_name, hist_name, n_bins, obs_min, obs_max, smear_fac in [0, 1])
    # ROOT.SmearedQuadraticHistogram(file_name, hist_name, n_bins, obs_min, obs_max, slope, smear_fac in [0, 1])

    # List histogram takes a vector of bin edges
    # bin_edges = ROOT.std.vector('double')()
    # map(bin_edges.push_back, [1, 2, 3, 4, 5, 6])  # 5 bins in [1, 6)
    # ROOT.ListHistogram(file_name, hist_name, bin_edges)

    def get_filename(s):
        name = "l_s%g" % (s,)
        filename = (params.output % name) + '.hist'
        return filename

    def Histogram(s, *args):
        if s == 0.:
            print "LinearHistogram ", args
            return ROOT.LinearHistogram(*args)
        else:
            args += (s,)  # add smearing parameter
            print "SmearedHistogram", args
            return ROOT.SmearedLinearHistogram(*args)

    def AddHistogram(hvec, *args):
        for s in smear:
            filename = get_filename(s)
            hvec.push_back(Histogram(s, filename, *args))

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
                    AddHistogram(getattr(analysis, hname)[n], histname, *histparam)
        elif '%d' in hparam[0]:
            for i in range(params.njet+1):
                histname = hparam[0] % (i+1)
                AddHistogram(getattr(analysis, hname)[i], histname, *histparam)
        else:
            histname = hparam[0]
            AddHistogram(getattr(analysis, hname), histname, *histparam)

    return

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

def add_grids_all(analysis, params):
    # total xs grid
    if True:
        obs = [-0.5, 0.5]  # only one bin to save memory
        filename = (params.output % 'totxs') + '.root'  # has to end with '.root'
        analysis.g_jet_inclusive = create_grid(filename, obs, params)

    # inclusive jets grid
    if False:
        obs = (lambda n: np.linspace(-0.5, n+0.5, n+2))(params.njet+1)
        filename = (params.output % 'incl') + '.root'  # has to end with '.root'
        analysis.g_jet_inclusive = create_grid(filename, obs, params)

    if False:
        obs = np.linspace(0, 500, 15+1)
        filename = (params.output % 'phmass') + '.root'  # has to end with '.root'
        analysis.g_photon_mass = create_grid(filename, obs, params)

    return

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# main function which is called from hammer.py
def initialize(params, selector):
    print "Using higgs analysis %s" % __file__

    jet_R = 0.4
    jetrpat = r'-R(\d+\.?\d*)-'
    m = re.match(r".*?%s.*?" % jetrpat, params.output)
    if m:
        jet_R = float(m.group(1))

    analysis = ROOT.HiggsJetsAnalysis.create()
    analysis.runname = params.runname
    analysis.setJetNumber(params.njet)
    analysis.setAntiKt(jet_R)
    analysis.jet_ptmin = 30
    analysis.jet_etamax = 4.4
    if 'vbf' in params.output:
        analysis.min_dijet_m = 400
        analysis.min_dijet_y = 2.8
    selector.opt_extra_scale = 125
    selector.opt_extra_alphas = 2
    selector.opt_extra_factor = 1
    if params.rescaler in ['mult_higgs', 'fixed_higgs', 'hthat_higgs']:
        selector.opt_extra_alphas = 0
    if params.rescaler == 'multiplicative':
        selector.opt_extra_factor = selector.opt_rescale_factor

    # Extract smear value from the output pattern, e.g. -smear0.3- or -smear0,0.1,0.3-
    smearpat = r'-smear(\d+\.?\d*|[\d.,]+)-'
    m = re.match(r".*?%s.*?" % smearpat, params.output)
    if m:
        params.output = re.sub(smearpat, '-', params.output)
        val = m.group(1)
        try:
            smear = [float(val)]
        except ValueError:
            smear = eval(val)
    else:
        smear = [0.]

    if not params.grids:
        # Add analysis histograms (see the function above)
        add_histograms_all(analysis, params, smear=smear)
    else:
        # Setting up grids
        ROOT.Grid.nloops = 1                     # number of loops, 0 - LO, 1 - NLO
        ROOT.Grid.pdf_function = "ntuplejets"    # 'ntuplephjets' for photons, 'ntuplejets' for jets
        ROOT.Grid.aparam = 5.
        ROOT.Grid.born_alphaspower = selector.opt_born_alphaspower
        # set the limits on x1, x2 and Q2
        fac = selector.opt_rescale_factor
        if 'vbf' in params.output:
            ROOT.Grid.def_opts = ROOT.GridOpts(100, (90*fac)**2, (4000*fac)**2, 5,
                                               100, 0.0044, 1., 5)
        else:
            ROOT.Grid.def_opts = ROOT.GridOpts(100, (90*fac)**2, (4000*fac)**2, 5,
                                               100, 0.00054, 1., 5)
        add_grids_all(analysis, params)

    # assign to selector
    selector.analysis = analysis


if __name__ == '__main__':
    print "Analysis module can't be run alone"
    print "pass --analysis=%s to hammer.py" % __file__
