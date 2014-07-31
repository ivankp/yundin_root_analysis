#!/usr/bin/env python

import getopt
import sys
import os
import re
import imp
import glob


# python 2.4 does not have any and all
try: any, all
except NameError:
    any = lambda x: reduce(lambda a,b: a or b, x)
    all = lambda x: reduce(lambda a,b: a and b, x)

# all custom rescalers
def set_rescaler(selector, params):
    selector.opt_rescale_n = params.rescale_n
    if params.rescaler == 'simple':
        if selector.opt_rescale_factor != 1:
            selector.setrescaler_multiplicative()
    elif params.rescaler == 'ht':
        selector.setrescaler_ht()
    elif params.rescaler == 'hthat':
        selector.setrescaler_hthat()
    elif params.rescaler == 'sumpt2':
        selector.setrescaler_sumpt2()
    elif params.rescaler == 'sumpt2hat':
        selector.setrescaler_sumpt2hat()
    elif params.rescaler == 'maaht':
        selector.setrescaler_maaht()
    elif params.rescaler == 'maahthat':
        selector.setrescaler_maahthat()
    elif params.rescaler == 'maa2sumpt2':
        selector.setrescaler_maa2sumpt2()
    elif params.rescaler == 'maa2sumpt2hat':
        selector.setrescaler_maa2sumpt2hat()
    else:
        name = 'setrescaler_%s' % params.rescaler
        try:
            func = getattr(selector, name)
            func()
        except:
            print "Rescaler %s not found" % (name)
            sys.exit(2)


def get_pdfname(pdfopt, name=''):
    pdf, m = pdfopt, 0
    if pdfopt.find(':') >= 0:
        pdf, m = pdfopt.split(':')
    if m != 0:
        print "Selected %s: %s (%s)" % (name, pdf, m)
    return (pdf, m)


def process(params):
    import ROOT

    # we want to handle Ctrl+C
    sh = ROOT.TSignalHandler(ROOT.kSigInterrupt, False)
    sh.Add()
    sh.Connect("Notified()", "TROOT", ROOT.gROOT, "SetInterrupt()")

    # add hammer.py directory to macro path
    hammer_path = os.path.dirname(__file__)
    ROOT.gROOT.SetMacroPath(ROOT.gROOT.GetMacroPath().rstrip(':') + ':' + hammer_path)

    if params.noapplgrid:
        ROOT.gSystem.AddIncludePath("-DDISABLE_APPLGRID")

    if params.noloopsim:
        ROOT.gSystem.AddIncludePath("-DDISABLE_LOOPSIM")

    # load libraries
    ROOT.gSystem.Load("libfastjet.so")
    ROOT.gSystem.Load("libLHAPDF.so")
    if not params.noapplgrid:
        ROOT.gSystem.Load("libAPPLgrid.so")
    if not params.noloopsim:
        ROOT.gSystem.Load("libLoopSim.so")

    # load macros
    ROOT.gSystem.AddIncludePath("-I%s" % hammer_path)
    ROOT.gROOT.LoadMacro("LHAPDF.h+")
    ROOT.gROOT.LoadMacro("SherpaAlphaS.C+")
    ROOT.gROOT.LoadMacro("FlavourKT.cpp+")
    ROOT.gROOT.LoadMacro("SelectorCommon.C+")
    if not params.noapplgrid:
        ROOT.gROOT.LoadMacro("appl_grid.h+")

    # create a chain
    chain = ROOT.TChain("t3")
    for name in params.inputs:
        chain.Add(name)

    selector = ROOT.SelectorCommon()

    # Initialize reweighting
    selector.opt_frompdf = 0
    selector.opt_topdf = 0
    selector.opt_rescale_factor = 1
    if params.power is not None:
        selector.opt_born_alphaspower = params.power
    if params.nborn is not None:
        selector.opt_loopsim_nborn = params.nborn

    if params.qfilter:
        selector.opt_filter_inq = params.qfilter[0]
        selector.opt_filter_nq = params.qfilter[1]
    if params.frompdf is not None:
        ROOT.LHAPDF.setVerbosity(0)  # comment out for more output

        selector.opt_rescale_factor = params.scale

        # FROMPDF is always initialized
        if True:
            selector.opt_frompdf = 1
            pdf, m = get_pdfname(params.frompdf, 'FROMPDF')
            ROOT.LHAPDF.initPDFSet(selector.opt_frompdf, pdf, ROOT.LHAPDF.LHGRID, int(m))

        # TOPDF is initialized and it is different
        if params.topdf == params.frompdf:
            selector.opt_topdf = selector.opt_frompdf
        else:
            selector.opt_topdf = selector.opt_frompdf + 1
            pdf, m = get_pdfname(params.topdf, 'TOPDF')
            ROOT.LHAPDF.initPDFSet(selector.opt_topdf, pdf,  ROOT.LHAPDF.LHGRID, int(m))

        # set rescaler after PDFs
        set_rescaler(selector, params)

        print "Scale: '%s (%d)' x %f" % (
            params.rescaler, selector.opt_rescale_n, selector.opt_rescale_factor)

        print "------------- FROMPDF %d - %s (Nf=%d) ---------------" % (
            selector.opt_frompdf, params.frompdf, ROOT.LHAPDF.getNf(selector.opt_frompdf))
        print "QMASS %s" % repr([ROOT.LHAPDF.getQMass(selector.opt_frompdf, qn) for qn in [1,2,3,4,5,6]])
        print "QTHRE %s" % repr([ROOT.LHAPDF.getThreshold(selector.opt_frompdf, qn) for qn in [1,2,3,4,5,6]])
        ROOT.LHAPDF.getDescription(selector.opt_frompdf)

        print "------------- TOPDF %d - %s (Nf=%d) -----------------" % (
            selector.opt_topdf, params.topdf, ROOT.LHAPDF.getNf(selector.opt_topdf))
        print "QMASS %s" % repr([ROOT.LHAPDF.getQMass(selector.opt_topdf, qn) for qn in [1,2,3,4,5,6]])
        print "QTHRE %s" % repr([ROOT.LHAPDF.getThreshold(selector.opt_topdf, qn) for qn in [1,2,3,4,5,6]])
        ROOT.LHAPDF.getDescription(selector.opt_topdf)
        print "--------------------------------------------------"

        if params.debug:
            print "WARNING! using SHERPA running AlphaS!"
            order = ROOT.LHAPDF.getOrderAlphaS(selector.opt_frompdf)
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

        if params.beta0fix is True:
            params.beta0fix = 99
        if params.beta0fix:
            selector.opt_beta0fix = params.beta0fix
            if params.beta0fix > 0:
                print "WARNING! nonsense beta0 fix enabled with for m_oqcd = %d" % params.beta0fix
            elif params.beta0fix < 0:
                print "WARNING! scheme dependence is disabled"

        if params.pi2o12fix:
            selector.opt_pi2o12fix = params.pi2o12fix
            print "WARNING! pi^2/12 factor included"

        if params.cdr2fdhfix is not None:
            selector.opt_cdr2fdhfix = params.cdr2fdhfix
            print "WARNING! CDR-DRED conversion factor included = %d" % params.cdr2fdhfix

    # call analysis module to initialize all settings
    params.analysis_mod.initialize(params, selector)

    selector.opt_stat_step = params.stat
    chain.Process(selector)
    selector.stat_report()

    if selector.opt_stat_step:
        import matplotlib.pyplot as plt
        import numpy as np

        yval = []
        for x in selector.xsvals:
            yval.append(x)
        yerr = []
        for x in selector.xserrs:
            yerr.append(x)

        xval = [1+i*selector.opt_stat_step for i in range(len(yval))]

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
Basic options:
  -a, --analysis            Analysis module
  -n, --njet                Number of jets
  -r, --runname='Test'      Run name for AidaPath
  -o, --output='%s.hist'    Output name pattern

  -s, --scale               Multiplicative scale factor
  --rescaler=<name>         Use rescaler 'simple', 'ht', 'hthat', 'sumpt2', 'sumpt2hat,
                            'maaht', 'maahthat', 'maa2sumpt2', 'maa2sumpt2hat'

  -p, --power               Born Alpha_s power (needed only for APPLgrid)
  -f, --frompdf             From PDF set "name.LHgrid:member"
  -t, --topdf               To PDF set "name.LHgrid:member"

  --nborn=n                 Set nborn for LoopSim

  --grids                   Specify to activate APPLgrid (for both warmup and fill)
  --warmup                  Select "warmup" mode for APPLgrid (otherwise "fill" mode)

Expert options:
  --beta0fix=true/[n]       Fix beta0 weight (0 turns off scheme dep)
  --debug                   Use sherpa alphas

  --qfilter=inq:nq          Filter input by 'incoming quarks':'total quarks'
  --stat=<N>                Eventoscope with step N

  --cdr2fdhfix=n            Convert CDR to FDH
  --pi2o12fix               Change prefactor conventions

Other options:
  -h, --help                show this help message
"""


class Params:
    def __init__(self):
        try:
            opts, args = getopt.getopt(sys.argv[1:], "a:n:s:p:o:r:f:t:h",
                                 ["analysis=", "njet=", "scale=", "power=", "output=", "runname=",
                                  "frompdf=", "topdf=", "beta0fix=", "cdr2fdhfix=", "pi2o12fix", "debug", "help",
                                  "stat=", "rescaler=", "qfilter=", "grids", "warmup", "nborn=",
                                  "noapplgrid", "noloopsim"])
        except getopt.GetoptError, err:
            print str(err)
            usage()
            sys.exit(2)

        self.analysis_mod = None
        self.njet = None
        self.runname = 'Test'
        self.output = '%s.hist'
        self.scale = None
        self.power = None
        self.frompdf = None
        self.topdf = None
        self.beta0fix = False
        self.cdr2fdhfix = None
        self.pi2o12fix = None
        self.debug = False
        self.grids = False
        self.warmup = False
        self.noapplgrid = False
        self.noloopsim = False
        self.qfilter = None
        self.stat = 0
        self.rescaler = 'simple'
        self.rescale_n = None
        self.nborn = None

        for op, oparg in opts:
            if op in ("-h", "--help"):
                usage()
                sys.exit()
            elif op in ("-a", "--analysis"):
                self.analysis_mod = oparg
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
            elif op in ("--beta0fix"):
                try:
                    self.beta0fix = int(oparg)
                except:
                    self.beta0fix = oparg.lower() in ['true', 'yes', 'on']
            elif op in ("--cdr2fdhfix"):
                try:
                    self.cdr2fdhfix = int(oparg)
                except:
                    if oparg.lower() in ['true', 'yes', 'on']:
                        self.cdr2fdhfix = 0
            elif op in ("--pi2o12fix"):
                self.pi2o12fix = True
            elif op in ("--debug"):
                self.debug = True
            elif op in ("--grids"):
                self.grids = True
            elif op in ("--warmup"):
                self.warmup = True
            elif op in ("--noapplgrid"):
                self.noapplgrid = True
            elif op in ("--noloopsim"):
                self.noloopsim = True
            elif op in ("--qfilter"):
                self.qfilter = oparg
            elif op in ("--nborn"):
                self.nborn = int(oparg)
            elif op in ("--stat"):
                self.stat = int(oparg)
            elif op in ("--rescaler"):
                self.rescaler = oparg.lower()
                if self.rescaler.find(':') >= 0:
                    self.rescaler, self.rescale_n = self.rescaler.split(':')
                    self.rescale_n = int(self.rescale_n)
            else:
                assert False, "unhandled option"

        if self.output.endswith('.hist'):
            self.output = self.output[:-5]

        if not self.analysis_mod:
            print "Error: --analysis option is mandatory"
            usage()
            sys.exit(2)
        try:
            # python 2.4 doesn't have sys.dont_write_bytecode
            if 'dont_write_bytecode' not in sys.__dict__:
                sys.dont_write_bytecode = False
            dwbc = sys.dont_write_bytecode
            sys.dont_write_bytecode = True
            if not self.analysis_mod.startswith('/'):
                self.analysis_mod = os.path.join(os.path.dirname(__file__), self.analysis_mod)
            if not os.path.exists(self.analysis_mod) and not self.analysis_mod.endswith('.py'):
                self.analysis_mod += '.py'
            self.analysis_mod = imp.load_source('analysis_mod', self.analysis_mod)
            sys.dont_write_bytecode = dwbc
        except Exception, e:
            print "Error: must give analysis module (%s)" % str(e)
            usage()
            sys.exit(2)

        if self.njet is None:
            print "Error: njet is not set"
            usage()
            sys.exit(2)

        if self.output.count(r'%s') != 1:
            print "Error: wrong output pattern"
            usage()
            sys.exit(2)

        if self.rescaler not in ['simple', 'ht', 'hthat', 'sumpt2', 'sumpt2hat',
                                 'maaht', 'maahthat', 'maa2sumpt2', 'maa2sumpt2hat',
                                 'minlo']:
            print "Unknown value for rescaler: %s" % self.rescaler
            usage()
            sys.exit(2)

        rewparam = [self.scale, self.frompdf, self.topdf]
        print rewparam
        if any(v is not None for v in rewparam) and any(v is None for v in rewparam):
            print "Error: must set --scale, --frompdf and --topdf together"
            usage()
            sys.exit(2)

        if self.grids and not self.scale:
            print "Error: has to specify --scale with --grids"
            usage()
            sys.exit(2)

        if self.grids and self.noapplgrid:
            print "Error: cannot use --noapplgrid with --grids"
            usage()
            sys.exit(2)

        if self.power is not None and not self.grids:
            print "Warning: --power option without --grids is deprecated (taken from event data)"

        if self.rescale_n is None:
            self.rescale_n = -1

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
            self.inputs = []
            for a in args:
                self.inputs.extend(glob.glob(a))
            self.inputs.sort()
        else:
            print "Error: missing input files"
            usage()
            sys.exit(2)


def main():
    params = Params()
    process(params)


if __name__ == '__main__':
    main()
