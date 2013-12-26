#!/usr/bin/env python

import ROOT
import numpy as np
import os
import sys
import math

def get(i, arr):
    if i < len(arr):
        return arr[i]
    else:
        return arr[-1]

# helper
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

def add_histograms_all(analysis, params):
    # Quadratic slope = lastbinwidth/firstbinwidth
    # ROOT.LinearHistogram(file_name, hist_name, n_bins, obs_min, obs_max)
    # ROOT.QuadraticHistogram(file_name, hist_name, n_bins, obs_min, obs_max, slope)
    # ROOT.SmearedLinearHistogram(file_name, hist_name, n_bins, obs_min, obs_max, smear_fac in [0, 1])
    # ROOT.SmearedQuadraticHistogram(file_name, hist_name, n_bins, obs_min, obs_max, slope, smear_fac in [0, 1])

    # List histogram takes a vector of bin edges
    # bin_edges = ROOT.std.vector('double')()
    # map(bin_edges.push_back, [1, 2, 3, 4, 5, 6])  # 5 bins in [1, 6)
    # ROOT.ListHistogram(file_name, hist_name, bin_edges)

    filename = (params.output % "l64_l20") + '.hist'

    # if there is not enough limits, last one is taken for excess elements
    minpt = [analysis.jet_ptmin]
    maxpt = [1420, 1400, 800]

    # jet-pT histograms
    for i in range(params.njet+1):
        pmin = get(i, minpt)
        pmax = get(i, maxpt)
        analysis.jet_pt_n[i].push_back(
            ROOT.LinearHistogram(filename, "jet_pT_%d" % (i+1),
                                 64, get(i, minpt), get(i, maxpt))
        )

    # jet-eta histograms
    for i in range(params.njet+1):
        maxeta = analysis.jet_etamax
        analysis.jet_eta_n[i].push_back(
            ROOT.LinearHistogram(filename, "jet_eta_%d" % (i+1),
                                 20, -maxeta, maxeta)
        )

    # vboson histograms
    analysis.vboson_pt.push_back(
        ROOT.LinearHistogram(filename, "vboson_pt", 15, 0, 500)
    )
    analysis.vboson_eta.push_back(
        ROOT.LinearHistogram(filename, "vboson_eta", 15, -4, 4)
    )

    return

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

def add_grids_all(analysis, params):
    # inclusive jets grid
    if True:
        obs = (lambda n: np.linspace(-0.5, n+0.5, n+2))(params.njet+1)
        filename = (params.output % 'incl') + '.root'  # has to end with '.root'
        analysis.g_jet_inclusive = create_grid(filename, obs, params)

    return

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# main function which is called from hammer.py
def initialize(params, selector):
    print "Using vector boson analysis %s" % __file__

    analysis = ROOT.VJetAnalysis.create()
    analysis.runname = params.runname
    analysis.setJetNumber(params.njet)  # resets all histograms in analysis
    analysis.setAntiKt(0.5)
    analysis.jet_ptmin = 25
    analysis.jet_etamax = 3.0
    analysis.lepton_ptmin = 20
    analysis.lepton_etamax = 2.5
    analysis.vboson_onshell_mass = 80.419
    analysis.vboson_mass_max = 7e3
    analysis.vboson_mass_min = 0.0
    analysis.etmiss_min = 20.0
    analysis.lepton_jet_Rsep = 0.0;

    # Add analysis histograms (see the function above)
    add_histograms_all(analysis, params)

    # Setting up grids
    ROOT.Grid.nloops = 1                     # number of loops, 0 - LO, 1 - NLO
    ROOT.Grid.pdf_function = "ntuplevjets"  # 'ntuplephjets' for photons, 'ntuplejets' for jets
    ROOT.Grid.aparam = 5.
    ROOT.Grid.born_alphapower = selector.born_alphapower
    # set the limits on x1, x2 and Q2
    fac = selector.rescale_factor
    ROOT.Grid.def_opts = ROOT.GridOpts(50, (22*fac)**2, (3900*fac)**2, 5,
                                       50, 0.0008, 1., 5)
    # Add grids
    if params.grids:
        add_grids_all(analysis, params)

    # assign to selector
    selector.analysis = analysis

if __name__ == '__main__':
    print "Analysis module can't be run alone"
    print "pass --analysis=%s to hammer.py" % __file__
