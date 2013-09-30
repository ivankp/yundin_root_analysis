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

    pdfabbr, pdfname = ('ct10', 'CT10')
    #pdfabbr, pdfname = ('nnpdf', 'NNPDF23_nlo_FFN_NF5_as_0118')
    #pdfabbr, pdfname = ('mstw08', 'MSTW2008nlo68cl')
    pdf, m = get_pdfname(pdfname) # 'CT10') # params.frompdf)
    ROOT.LHAPDF.initPDFSet(pdf, ROOT.LHAPDF.LHGRID, int(m))

    grids = [ROOT.appl.grid(name) for name in params.inputs]
    for g in grids:
      g.trim()

    order = g.nloops()
    firstname = params.inputs[0].replace('.root', '')
    firstname += pdfabbr
    print "Order ", order, " Name ", firstname

    if False:
        for g,n in zip(grids, params.inputs):
            print n
            print ROOT.vconvolute(g, 0, 1, 1)[2]
            if order > 0:
                print ROOT.vconvolute(g, 1, 1, 1)[2]

    # sum grids
    for g in grids[1:]:
      print g
      grids[0] += g
    g = grids[0]

    if True:
        xsec = ROOT.vconvolute(g, order, 1, 1)
        print 1, [x for x in xsec]
        xsec = ROOT.vconvolute(g, order, 2, 2)
        print 2, [x for x in xsec]
        xsec = ROOT.vconvolute(g, order, 0.5, 0.5)
        print 0.5, [x for x in xsec]

    try:
        xval_lo, yval_lo = np.loadtxt(firstname + '-lo.txt')
        if order > 0:
            xval_nlo, yval_nlo = np.loadtxt(firstname + '-nlo.txt')
    except IOError:
        xval_lo = np.logspace(-1.5, 1, 100)
        yval_lo = []
        for f in xval_lo:
            xsec_lo = ROOT.vconvolute(g, 0, f, f)
            yval_lo.append(xsec_lo[2])
        yval_lo = np.array(yval_lo)
        np.savetxt(firstname + '-lo.txt', np.array([xval_lo, yval_lo]))
        if order > 0:
            xval_nlo = xval_lo
            yval_nlo = []
            for f in xval_nlo:
                print f
                xsec_nlo = ROOT.vconvolute(g, 1, f, f)
                yval_nlo.append(xsec_nlo[2])
            yval_nlo = np.array(yval_nlo)
            np.savetxt(firstname + '-nlo.txt', np.array([xval_nlo, yval_nlo]))

    if True:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        axs = plt.subplot(111)
        axs.grid(True)
        def vline(x, *args, **kwargs):
            axs.plot([x, x], [min(yval_lo), max(yval_lo)], *args, **kwargs)
        vline(x=1, color='k', ls='--', lw=1.5)
        vline(x=0.5, color='gray', ls='--', lw=1.5)
        vline(x=2, color='gray', ls='--', lw=1.5)
        #axs.set_xscale('log')
        #axs.set_ylim(0, 6)
        #axs.set_xlim(0.1, 8)
        axs.plot(xval_lo, yval_lo, lw=2, c='b') #, marker='x')
        if order > 0:
            axs.plot(xval_nlo, yval_nlo, lw=2, c='r') #, marker='x')
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
