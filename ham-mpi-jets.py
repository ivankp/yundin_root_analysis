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

    filename = (params.output % "l32_l15") + '.hist'

    # if there is not enough limits, last one is taken for excess elements
    minpt = [analysis.jet_pt1min, analysis.jet_ptmin]
    maxpt = [1420, 1400, 800]

    # jet-pT histograms
    for i in range(params.njet+1):
        pmin = get(i, minpt)
        pmax = get(i, maxpt)
        analysis.jet_pt_n[i].push_back(
            ROOT.SmearedLinearHistogram(filename, "jet_pT_%d" % (i+1),
                                        32, get(i, minpt), get(i, maxpt), 0.5)
        )

    # jet-eta histograms
    for i in range(params.njet+1):
        maxeta = analysis.jet_etamax
        analysis.jet_eta_n[i].push_back(
            ROOT.SmearedLinearHistogram(filename, "jet_eta_%d" % (i+1),
                                        15, -maxeta, maxeta, 0.5)
        )

    ### effective double diff. type histogram ###
    analysis.jets_d12_d34.push_back(
        ROOT.LinearHistogram2D(filename, "jets_d12_d34",
          15, 0, 1e6, 15, 0, 1e6)
        )

    return

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# main function which is called from hammer.py
def initialize(params, selector):
    print "Using jets analysis %s" % __file__

    analysis = ROOT.FourJetMPIAnalysis.create()
    analysis.runname = params.runname
    analysis.setJetNumber(params.njet)
    analysis.setAntiKt(0.4)
    analysis.jet_ptmin = 60
    analysis.jet_etamax = 2.8
    analysis.jet_pt1min = 80

    # Add analysis histograms (see the function above)
    add_histograms_all(analysis, params)

    # assign to selector
    selector.analysis = analysis


if __name__ == '__main__':
    print "Analysis module can't be run alone"
    print "pass --analysis=%s to hammer.py" % __file__
