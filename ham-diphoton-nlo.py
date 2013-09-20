#!/usr/bin/env python

import ROOT
import numpy as np
import os

# helper
def create_grid(name, edges):
    if not os.path.exists(name):
        obs = ROOT.std.vector('double')()
        for ed in edges:
            obs.push_back(ed)
        return ROOT.Grid.capture(ROOT.Grid(name, obs))
    else:
        return ROOT.Grid.capture(ROOT.Grid(name))

# main function which is called from hammer.py
def initialize(params, selector):
    print "Using diphoton analysis %s" % __file__

    analysis = ROOT.DiPhotonAnalysis.create()
    analysis.runname = params.runname
    analysis.setJetNumber(params.njet)
    analysis.setAntiKt(0.5)
    analysis.jet_ptmin = 30
    analysis.jet_etamax = 4.7
    analysis.photon_R = 0.4
    analysis.photon_pt1min = 40
    analysis.photon_pt2min = 25
    analysis.photon_etamax = 2.5
    analysis.photon_jet_Rsep = 0.5
    analysis.photon_photon_Rsep = 0.45

    # grids
    ROOT.Grid.aparam = 5.
    fac = selector.rescale_factor
    ROOT.Grid.def_opts = ROOT.GridOpts(75, (20*fac)**2, (4000*fac)**2, 5,
                                       75, 1e-8, 1., 5)
    ROOT.Grid.born_alphapower = selector.born_alphapower
    ROOT.Grid.nloops = 1
    ROOT.Grid.pdf_function = "ntuplephjets"

    obs = (lambda n: np.linspace(-0.5, n+0.5, n+2))(params.njet+1)
    filename = (params.output % 'incl') + '.root'
    analysis.g_jet_inclusive = create_grid(filename, obs)
    if not analysis.g_jet_inclusive.isWarmup():
        analysis.g_jet_inclusive.setFilename(filename.replace('.root', '-o.root'))

    # if more than 3 jets, extra jets use the last defined value (here 800)
    maxpt = ROOT.std.vector('double')()
    for pt in [1420, 1400, 800]:
        maxpt.push_back(pt)

    # add histograms
    # minpt is default to jet_ptmin
    analysis.addPtLinearHistograms(params.output % "l64_l20", 64, None, maxpt)
    # eta limits are default to +-jet_etamax
    analysis.addEtaLinearHistograms(params.output % "l64_l20", 20)

    # assign to selector
    selector.analysis = analysis


if __name__ == '__main__':
    print "Analysis module can't be run alone"
    print "pass --analysis=%s to hammer.py" % __file__
