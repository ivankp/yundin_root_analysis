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

def evolve(name, g, order, xval_lo, xval_nlo, fac=1):
    yval_lo = []
    for f in xval_lo:
        xsec_lo = ROOT.vconvolute(g, 0, f*fac, f*fac)
        yval_lo.append(list(xsec_lo))
        print 'LO', f, yval_lo[-1]
    yval_lo = np.array(yval_lo)
    np.savetxt(name % 'lo', np.column_stack((xval_lo,yval_lo)))
    if order > 0:
        yval_nlo = []
        for f in xval_nlo:
            xsec_nlo = ROOT.vconvolute(g, 1, f*fac, f*fac)
            yval_nlo.append(list(xsec_nlo))
            print 'NLO', f, yval_nlo[-1]
        yval_nlo = np.array(yval_nlo)
        np.savetxt(name % 'nlo', np.column_stack((xval_nlo,yval_nlo)))
    return yval_lo, yval_nlo


def process(params):
    # we want to handle Ctrl+C
    #sh = ROOT.TSignalHandler(ROOT.kSigInterrupt, False)
    #sh.Add()
    #sh.Connect("Notified()", "TROOT", ROOT.gROOT, "SetInterrupt()")

    # add hammer.py directory to macro path
    hammer_path = os.path.dirname(__file__)
    ROOT.gROOT.SetMacroPath(ROOT.gROOT.GetMacroPath().rstrip(':') + ':' + hammer_path)

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

    if params.pdf == 'CT10':
        pdfabbr, pdfname = ('ct10', 'CT10')
    elif params.pdf == 'NNPDF23':
        pdfabbr, pdfname = ('nnpdf', 'NNPDF23_nlo_FFN_NF5_as_0118')
    elif params.pdf == 'MSTW08':
        pdfabbr, pdfname = ('mstw08', 'MSTW2008nlo68cl')
    else:
        pdfabbr, pdfname = params.pdf, params.pdf
    pdf, m = get_pdfname(pdfname)
    ROOT.LHAPDF.initPDFSet(pdf, ROOT.LHAPDF.LHGRID, int(m))

    grids = [ROOT.appl.grid(name) for name in params.inputs]
    for g in grids:
      g.trim()

    order = g.nloops()
    firstname = params.inputs[0].replace('.root', '')
    firstname += pdfabbr
    print "Order ", order, " Name ", firstname

    if True:
        for g,n in zip(grids, params.inputs):
            print n
            print list(ROOT.vconvolute(g, 0, 1, 1))
            if order > 0:
                print list(ROOT.vconvolute(g, 1, 1, 1))

    # sum grids
    g = grids[0]
    for i in range(1, len(grids)):
      g += grids[i]
      grids[i] = None

    if True:
        xsec = ROOT.vconvolute(g, order, 1, 1)
        print 1, list(xsec)
        xsec = ROOT.vconvolute(g, order, 2, 2)
        print 2, list(xsec)
        xsec = ROOT.vconvolute(g, order, 0.5, 0.5)
        print 0.5, list(xsec)

    try:
        tmp = np.loadtxt(firstname + '-lo.txt')
        xval_lo = tmp[...,0]
        yval_lo = tmp[...,1:]
        if order > 0:
            tmp = np.loadtxt(firstname + '-nlo.txt')
            xval_nlo = tmp[...,0]
            yval_nlo = tmp[...,1:]
    except IOError:
        xval_lo = np.sort(np.append(np.logspace(-1.5, 1, 50), [0.5, 1, 2]))
        fac = 1
        if firstname.find('0.025') > 0:
            print "0.025 detected"
            fac = 1./0.025
        if firstname.find('0.1') > 0:
            print "0.1 detected"
            fac = 1./0.1
        xval_nlo = xval_lo
        yval_lo, yval_nlo = evolve(firstname + '-%s.txt', g, order, xval_lo, xval_nlo, fac)
    # pick total XS
    yval_lo = yval_lo[...,0]
    yval_nlo = yval_nlo[...,0]

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

        #axs.plot(xval_lo2, yval_lo2, lw=1, c='c') #, marker='x')
        #if order > 0:
        #    axs.plot(xval_nlo2, yval_nlo2, lw=1, c='k') #, marker='x')
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
            opts, args = getopt.getopt(sys.argv[1:], "o:dhp:",
                                 ["output=", "debug", "help", "pdf="])
        except getopt.GetoptError, err:
            print str(err)
            usage()
            sys.exit(2)

        self.output = 'test.root'
        self.debug = False
        self.pdf = None

        for op, oparg in opts:
            if op in ("-h", "--help"):
                usage()
                sys.exit()
            elif op in ("-o", "--output"):
                self.output = oparg
            elif op in ("-d", "--debug"):
                self.debug = True
            elif op in ("-p", "--pdf"):
                self.pdf = oparg
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
