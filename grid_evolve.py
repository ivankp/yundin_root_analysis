#!/usr/bin/env python

import getopt
import sys
import re
from array import array
import os
import math
import ROOT
import numpy as np


def get_pdfname(pdfopt, name=''):
    pdf, m = pdfopt, 0
    if pdfopt.find(':') >= 0:
        pdf, m = pdfopt.split(':')
    if m != 0:
        print "Selected %s: %s (%s)" % (name, pdf, m)
    return (pdf, m)

def process(params):
    # we want to handle Ctrl+C
    #sh = ROOT.TSignalHandler(ROOT.kSigInterrupt, False)
    #sh.Add()
    #sh.Connect("Notified()", "TROOT", ROOT.gROOT, "SetInterrupt()")
    # load libraries
    ROOT.gSystem.Load("libfastjet.so")
    ROOT.gSystem.Load("libLHAPDF.so")
    ROOT.gSystem.Load("libAPPLgrid.so")
    #AddDynamicPath

    # load macros
    # gSystem->AddIncludePath(" -I$HOME/mypackage/include ")
    # gSystem->AddLinkedLibs("-L/my/path -lanylib");
    ROOT.gROOT.LoadMacro("LHAPDF.h+")
    ROOT.gROOT.LoadMacro("SherpaAlphaS.C+")
    ROOT.gROOT.LoadMacro("SelectorCommon.C+")

    ROOT.gROOT.LoadMacro("appl_grid.h+")
    ROOT.gROOT.LoadMacro("appl_helper.h+")
    ROOT.gROOT.LoadMacro("ntuple_pdf.h+")

    pdf_functions = ["ntupleall", "ntuplephjets", "ntuplejets"]
    pdf_objects = [getattr(ROOT, "%s_pdf" % p)() for p in pdf_functions]

    grid = ROOT.appl.grid(params.inputs[0])
    grid.trim()

    pdf, m = get_pdfname('CT10') # params.frompdf)
    ROOT.LHAPDF.initPDFSet(pdf, ROOT.LHAPDF.LHGRID, int(m))

    xsec = ROOT.vconvolute(grid, 0, 1, 1)
    print [(x, x/2.379820259160597e+00) for x in xsec[2:]]
    xsec = ROOT.vconvolute(grid, 0, 2, 2)
    print [(x, x/1.898278301459246e+00) for x in xsec[2:]]
    xsec = ROOT.vconvolute(grid, 0, 0.5, 0.5)
    print [(x, x/3.039590253832915e+00) for x in xsec[2:]]

    if False:
        xval = np.logspace(-1.5, 2.1, 1000)
        yval = []
        for f in xval:
            xsec = ROOT.vconvolute(grid, 0, f, f)
            yval.append(xsec[2])
        yval = np.array(yval)
        xyval = np.array([xval, yval])
        np.savetxt('diphoton-lo.txt', xyval)
    else:
        xval, yval = np.loadtxt('diphoton-lo.txt')
        #numpy.interp(x, xval, yval, left=None, right=None)

    if True:
        import matplotlib.pyplot as plt

        fig = plt.figure()
        axs = plt.subplot(111)
        axs.grid(True)
        def vline(x, *args, **kwargs):
            axs.plot([x,x], [min(yval), max(yval)], *args, **kwargs)
        vline(x=1, color='k', ls='--', lw=1.5)
        vline(x=0.5, color='gray', ls='--', lw=1.5)
        vline(x=2, color='gray', ls='--', lw=1.5)
        axs.set_xscale('log')
        #axs.set_ylim(0, 6)
        #axs.set_xlim(0.1, 8)
        axs.plot(xval, yval, lw=2, c='r', marker='x')
        #axs.plot(xval, np.array(yval)+np.array(yerr), c='r')
        #axs.plot(xval, np.array(yval)-np.array(yerr), c='b')
        #plt.xticks([0.2,1.,4.], ['0.2', '1', '4'])
        plt.show()


def usage():
    print """\
Usage: analyze [OPTION...] [FILE]
Reweight events
  -o, --output='name.root'  Output name
  -d, --debug               debug

Other options:
  -h, --help                show this help message
"""


class Params:
    def __init__(self):
        try:
            opts, args = getopt.getopt(sys.argv[1:], "o:dh",
                                 ["output=", "debug", "help"])
        except getopt.GetoptError, err:
            print str(err)
            usage()
            sys.exit(2)

        self.output = 'test.root'
        self.debug = False

        for op, oparg in opts:
            if op in ("-h", "--help"):
                usage()
                sys.exit()
            elif op in ("-o", "--output"):
                self.output = oparg
            elif op in ("-d", "--debug"):
                self.debug = True
            else:
                assert False, "unhandled option"

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
