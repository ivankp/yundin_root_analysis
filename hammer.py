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
    ROOT.gROOT.LoadMacro("T3selector.C+")
    ROOT.gROOT.LoadMacro("LHAPDF.h+")

    # create a chain
    chain = ROOT.TChain("t3")
    for name in params.inputs:
        chain.Add(name)

    selector = ROOT.T3selector()

    maxpt = ROOT.std.vector(int)()
    for pt in [1420, 1400, 800, 800, 800]:
        maxpt.push_back(pt)

    # Initialize analysis
    selector.analysis.setN(params.njet)
    selector.analysis.init("%s! algo=antikt R=0.4 Pt=60 Pt1=80 Eta=2.8" % params.runname)

    # Initialize reweighting
    selector.FROMPDF = 1
    selector.TOPDF = 1
    selector.scalefactor = 1
    selector.alphapower = 0
    if params.frompdf is not None:
        ROOT.LHAPDF.setVerbosity(0)  # comment out for more output

        selector.scalefactor = params.scale
        selector.alphapower = params.power

        # FROMPDF is always initialized
        if True:
            selector.FROMPDF = 1
            pdf, m = params.frompdf, 0
            if params.frompdf.find(':') >= 0:
                pdf, m = params.frompdf.split(':')
            ROOT.LHAPDF.initPDFByName(selector.FROMPDF, pdf, int(m));

        # TOPDF is initialized is it is different
        if params.topdf == params.frompdf:
            selector.TOPDF = 1
        else:
            selector.TOPDF = 2
            pdf, m = params.topdf, 0
            if params.topdf.find(':') >= 0:
                pdf, m = params.topdf.split(':')
            ROOT.LHAPDF.initPDFByName(selector.TOPDF, pdf, int(m));

        print "Scale: %f, AlphaPow %d, PDF1 %d PDF2 %d" % (selector.scalefactor,
                        selector.alphapower, selector.FROMPDF, selector.TOPDF)
        ROOT.LHAPDF.getDescription(selector.TOPDF)
        ROOT.LHAPDF.getDescription(selector.FROMPDF)

    # add histograms
    selector.analysis.addPtLinearHistograms(params.output % "l64_l20", 64, maxpt)
    selector.analysis.addEtaLinearHistograms(params.output % "l64_l20", 20)

    chain.Process(selector)


def usage():
    print """\
Usage: hammer [OPTION...] [FILE]
Reweight events
  -n, --njet                Number of jets
  -r, --runname='Test'      Run name
  -o, --output='%s.hist'    Output name pattern

  -s, --scale               Multiplicative scale factor (1. does nothing)
  -p, --power               Alpha_s power (0 does nothing)
  -f, --frompdf             From PDF set
  -t, --topdf               To PDF set (equal to "frompdf" does nothing)

Other options:
  -h, --help                show this help message
"""


class Params:
    def __init__(self):
        try:
            opts, args = getopt.getopt(sys.argv[1:], "n:s:p:o:r:f:t:h",
                                 ["njet", "scale", "power", "output", "runname", "frompdf", "topdf", "help"])
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
                self.frompdf = oparg
            elif op in ("-t", "--topdf"):
                self.topdf = oparg
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

        rewparam = [self.scale, self.power, self.frompdf, self.topdf]
        print rewparam
        if any(v is not None for v in rewparam) and any(v is None for v in rewparam):
            print "Error: must set --scale, --power, --frompdf and --topdf"
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
