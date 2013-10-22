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

    # photon histograms
    analysis.photon_mass.push_back(
        ROOT.LinearHistogram(filename, "photon_mass", 15, 0, 500)
    )
    analysis.photon_pt.push_back(
        ROOT.LinearHistogram(filename, "photon_pt", 15, 0, 500)
    )
    maxeta = analysis.photon_etamax
    analysis.photon_eta.push_back(
        ROOT.LinearHistogram(filename, "photon_eta", 15, -maxeta, maxeta)
    )
    analysis.jet_jet_phi12.push_back(
        ROOT.LinearHistogram(filename, "jet_jet_phi12", 31, 0, math.pi)
    )
    analysis.photon_jet_R11.push_back(
        ROOT.LinearHistogram(filename, "photon_jet_R11", 30, 0, 5)
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
    print "Using diphoton analysis %s" % __file__

    analysis = ROOT.DiPhotonAnalysis.create()
    analysis.runname = params.runname
    analysis.setJetNumber(params.njet)  # resets all histograms in analysis
    analysis.setAntiKt(0.5)
    analysis.jet_ptmin = 30
    analysis.jet_etamax = 4.7
    analysis.photon_R = 0.4
    analysis.photon_pt1min = 40
    analysis.photon_pt2min = 25
    analysis.photon_etamax = 2.5
    analysis.photon_jet_Rsep = 0.5
    analysis.photon_photon_Rsep = 0.45

    # Add analysis histograms (see the function above)
    add_histograms_all(analysis, params)

    # Setting up grids
    ROOT.Grid.nloops = 1                     # number of loops, 0 - LO, 1 - NLO
    ROOT.Grid.pdf_function = "ntuplephjets"  # 'ntuplephjets' for photons, 'ntuplejets' for jets
    ROOT.Grid.aparam = 5.
    ROOT.Grid.born_alphapower = selector.born_alphapower
    # set the limits on x1, x2 and Q2
    fac = selector.rescale_factor
    ROOT.Grid.def_opts = ROOT.GridOpts(50, (22*fac)**2, (3900*fac)**2, 5,
                                       50, 0.0008, 1., 5)
    # Add grids
    add_grids_all(analysis, params)

    # assign to selector
    selector.analysis = analysis


if __name__ == '__main__':
    print "Analysis module can't be run alone"
    print "pass --analysis=%s to hammer.py" % __file__
