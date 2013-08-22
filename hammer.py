#!/usr/bin/env python

import getopt
import sys
import re

def process(params):
    import ROOT

    # we want to handle Ctrl+C
    sh = ROOT.TSignalHandler(ROOT.kSigInterrupt, False)
    sh.Add()
    sh.Connect("Notified()", "TROOT", ROOT.gROOT, "SetInterrupt()")

    # load libraries
    ROOT.gSystem.Load("libfastjet.so")
    ROOT.gSystem.Load("libLHAPDF.so")

    # load macros
    ROOT.gROOT.LoadMacro("LHAPDF.h+")
    ROOT.gROOT.LoadMacro("SherpaAlphaS.C+")
    ROOT.gROOT.LoadMacro("SelectorCommon.C+")

    # create a chain
    chain = ROOT.TChain("t3")
    for name in params.inputs:
        chain.Add(name)

    selector = ROOT.SelectorCommon()

    maxpt = ROOT.std.vector('double')()
    for pt in [1420, 1400, 800, 800, 800]:
        maxpt.push_back(pt)

    # Initialize analysis
    if False:
        analysis = ROOT.JetAnalysis.create()
        analysis.setAntiKt(0.4)
        analysis.jet_ptmin = 60
        analysis.jet_etamax = 2.8
        analysis.jet_pt1min = 80
    else:
        analysis = ROOT.DiPhotonAnalysis.create()
        analysis.setAntiKt(0.5)
        analysis.jet_ptmin = 30
        analysis.jet_etamax = 4.7
        analysis.photon_pt1min = 40
        analysis.photon_pt2min = 25
        analysis.photon_etamax = 2.5
        analysis.photon_jet_Rsep = 0.5
        analysis.photon_photon_Rsep = 0.45

    selector.analysis = analysis
    selector.analysis.setJetNumber(params.njet)
    selector.analysis.runname = params.runname

    # Initialize reweighting
    selector.FROMPDF = 0
    selector.TOPDF = 0
    selector.rescale_factor = 1
    selector.alphapower = 0
    if params.qfilter:
        selector.filter_inq = params.qfilter[0]
        selector.filter_nq = params.qfilter[1]
    if params.frompdf is not None:
        ROOT.LHAPDF.setVerbosity(0)  # comment out for more output

        selector.alphapower = params.power
        selector.rescale_factor = params.scale

        # scale change
        selector.rescale_n = params.rescale_n
        if params.rescaler == 'simple':
            if selector.rescale_factor != 1:
                selector.setrescaler_multiplicative()
        elif params.rescaler == 'ht':
            selector.setrescaler_ht()
        elif params.rescaler == 'hthat':
            selector.setrescaler_hthat()
        elif params.rescaler == 'htn':
            selector.setrescaler_htn()
        elif params.rescaler == 'htnhat':
            selector.setrescaler_htnhat()
        elif params.rescaler == 'sumpt2':
            selector.setrescaler_sumpt2()
        elif params.rescaler == 'sumpt2hat':
            selector.setrescaler_sumpt2hat()

        # FROMPDF is always initialized
        if True:
            selector.FROMPDF = 1
            pdf, m = params.frompdf, 0
            if params.frompdf.find(':') >= 0:
                pdf, m = params.frompdf.split(':')
            if m != 0:
                print "Selected FROMPDF member = %s" % m
            ROOT.LHAPDF.initPDFSet(selector.FROMPDF, pdf, ROOT.LHAPDF.LHGRID, int(m))

        # TOPDF is initialized is it is different
        if params.topdf == params.frompdf:
            selector.TOPDF = selector.FROMPDF
        else:
            selector.TOPDF = selector.FROMPDF + 1
            pdf, m = params.topdf, 0
            if params.topdf.find(':') >= 0:
                pdf, m = params.topdf.split(':')
            if m != 0:
                print "Selected TOPDF member = %s" % m
            ROOT.LHAPDF.initPDFSet(selector.TOPDF, pdf,  ROOT.LHAPDF.LHGRID, int(m))

        print "Scale: '%s' x %f, AlphaPow %d" % (params.rescaler, selector.rescale_factor, selector.alphapower)
        print "------------- FROMPDF %d - %s (Nf=%d) ---------------" % (selector.FROMPDF, params.frompdf, ROOT.LHAPDF.getNf(selector.FROMPDF))
        print "QMASS %s" % repr([ROOT.LHAPDF.getQMass(selector.FROMPDF, qn) for qn in [1,2,3,4,5,6]])
        print "QTHRE %s" % repr([ROOT.LHAPDF.getThreshold(selector.FROMPDF, qn) for qn in [1,2,3,4,5,6]])
        ROOT.LHAPDF.getDescription(selector.FROMPDF)
        print "------------- TOPDF %d - %s (Nf=%d) -----------------" % (selector.TOPDF, params.topdf, ROOT.LHAPDF.getNf(selector.TOPDF))
        print "QMASS %s" % repr([ROOT.LHAPDF.getQMass(selector.TOPDF, qn) for qn in [1,2,3,4,5,6]])
        print "QTHRE %s" % repr([ROOT.LHAPDF.getThreshold(selector.TOPDF, qn) for qn in [1,2,3,4,5,6]])
        ROOT.LHAPDF.getDescription(selector.TOPDF)
        print "--------------------------------------------------"

        if params.debug:
            print "WARNING! using SHERPA running AlphaS!"
            order = ROOT.LHAPDF.getOrderAlphaS(selector.FROMPDF)
            # MUST set mZ, asMZ and qmass exactly as in Sherpa
            mZ = 91.188
            asMZ = 0.12018  # MSTW2008
            asMZ = 0.117982
            qmass = ROOT.std.vector('double')()
            for m in [0., 0.01, 0.005, 0.2, 1.42, 4.8, 1e10]:
                qmass.push_back(m)
            print "ALPHA_S(mZ) = %f, mZ = %f, order = %d, qM = %s" % (asMZ,
                mZ, order, repr(["%e" % x for x in qmass]))
            selector.initAS(order, asMZ, mZ*mZ, qmass)
            selector.use_sherpa_alphas = True

        if params.beta0fix:
            selector.beta0fix = 99
            print "WARNING! nonsense beta0 fix enabled with for m_oqcd = %d" % selector.beta0fix

    # add histograms
    selector.analysis.addPtLinearHistograms(params.output % "l64_l20", 64, maxpt)
    selector.analysis.addEtaLinearHistograms(params.output % "l64_l20", 20)

    selector.stat_step = params.stat
    chain.Process(selector)

    if selector.stat_step:
        import matplotlib.pyplot as plt
        import numpy as np

        yval = []
        for x in selector.xsvals:
            yval.append(x)
        yerr = []
        for x in selector.xserrs:
            yerr.append(x)

        xval = [1+i*selector.stat_step for i in range(len(yval))]

        fig = plt.figure()
        axs = plt.subplot(111)
        axs.grid(True)
        axs.plot(xval, yval, lw=2, c='k')
        axs.plot(xval, np.array(yval)+np.array(yerr), c='r')
        axs.plot(xval, np.array(yval)-np.array(yerr), c='b')
        plt.show()

def usage():
    print """\
Usage: hammer [OPTION...] [FILE]
Reweight events
  -n, --njet                Number of jets
  -r, --runname='Test'      Run name
  -o, --output='%s.hist'    Output name pattern

  -s, --scale               Multiplicative scale factor (1. does nothing)
  -p, --power               Alpha_s power (0 does nothing)
  -f, --frompdf             From PDF set "name.LHgrid:member"
  -t, --topdf               To PDF set "name.LHgrid:member"

  -b, --beta0fix            Fix comix beta0 weight
  -d, --debug               Use sherpa alphas

  --qfilter=inq:nq           Filter input by 'incoming quarks':'total quarks'

  --stat=<N>                Eventoscope with step N
  --rescaler=<name>         Use rescaler 'simple', 'ht', 'hthat', 'htn', 'htnhat'

Other options:
  -h, --help                show this help message
"""


class Params:
    def __init__(self):
        try:
            opts, args = getopt.getopt(sys.argv[1:], "n:s:p:o:r:f:t:bdh",
                                 ["njet=", "scale=", "power=", "output=", "runname=",
                                  "frompdf=", "topdf=", "beta0fix", "debug", "help",
                                  "stat=", "rescaler=", "qfilter="])
        except getopt.GetoptError, err:
            print str(err)
            usage()
            sys.exit(2)

        self.njet = None
        self.runname = 'Test'
        self.output = '%s.hist'
        self.scale = None
        self.power = None
        self.frompdf = None
        self.topdf = None
        self.beta0fix = False
        self.debug = False
        self.qfilter = None
        self.stat = 0
        self.rescaler = 'simple'
        self.rescale_n = None

        for op, oparg in opts:
            if op in ("-h", "--help"):
                usage()
                sys.exit()
            elif op in ("-n", "--njet"):
                self.njet = int(oparg)
            elif op in ("-s", "--scale"):
                self.scale = float(oparg)
            elif op in ("-p", "--power"):
                self.power = int(oparg)
            elif op in ("-o", "--output"):
                self.output = oparg
            elif op in ("-r", "--runname"):
                self.runname = re.sub(r'\s+', r'_', oparg)
            elif op in ("-f", "--frompdf"):
                self.frompdf = oparg.replace('.LHgrid', '')
            elif op in ("-t", "--topdf"):
                self.topdf = oparg.replace('.LHgrid', '')
            elif op in ("-b", "--beta0fix"):
                self.beta0fix = True
            elif op in ("-d", "--debug"):
                self.debug = True
            elif op in ("--qfilter"):
                self.qfilter = oparg
            elif op in ("--stat"):
                self.stat = int(oparg)
            elif op in ("--rescaler"):
                self.rescaler = oparg.lower()
                if self.rescaler.find(':') >= 0:
                    self.rescaler, self.rescale_n = self.rescaler.split(':')
                    self.rescale_n = int(self.rescale_n)
            else:
                assert False, "unhandled option"

        if not self.njet:
            print "Error: njet is not set"
            usage()
            sys.exit(2)

        if self.output.count(r'%s') != 1:
            print "Error: wrong output pattern"
            usage()
            sys.exit(2)

        if self.rescaler not in ['simple', 'ht', 'hthat', 'htn', 'htnhat']:
            print "Unknown value for rescaler: %s" % self.rescaler
            usage()
            sys.exit(2)

        rewparam = [self.scale, self.power, self.frompdf, self.topdf]
        print rewparam
        if any(v is not None for v in rewparam) and any(v is None for v in rewparam):
            print "Error: must set --scale, --power, --frompdf and --topdf"
            usage()
            sys.exit(2)

        if self.rescale_n is None:
            self.rescale_n = self.njet

        if self.rescale_n > self.njet:
            print "Error: 'rescale_n' %d cannont be larger than 'njet' %d" % (self.rescale_n, self.njet)
            usage()
            sys.exit(2)

        try:
            if self.qfilter is not None:
                self.qfilter = [int(x) for x in self.qfilter.split(':')]
                assert len(self.qfilter) == 2
                assert self.qfilter[0] <= 2
                assert 0 <= self.qfilter[0] <= self.qfilter[1]
                assert self.qfilter[1] % 2 == 0
                assert 0 <= self.qfilter[1] <= self.njet+2
        except:
            print "Error: wrong qfilter pattern"
            usage()
            sys.exit(2)

        if len(args) > 0:
            self.inputs = args[:]
        else:
            print "Error: missing input files"
            usage()
            sys.exit(2)


def main():
    params = Params()
    process(params)


if __name__ == '__main__':
    main()
