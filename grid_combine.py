#!/usr/bin/env python

import getopt
import sys
import re
from array import array
import os
import math
import ROOT
import numpy as np
import gc

def process(params):
    # add hammer.py directory to macro path
    hammer_path = os.path.dirname(__file__)
    ROOT.gROOT.SetMacroPath(ROOT.gROOT.GetMacroPath().rstrip(':') + ':' + hammer_path)

    # load libraries
    ROOT.gSystem.Load("libAPPLgrid.so")

    ROOT.gROOT.LoadMacro("appl_grid.h+")
    ROOT.gROOT.LoadMacro("ntuple_pdf.h+")

    pdf_functions = ["ntupleall", "ntuplephjets", "ntuplejets"]
    pdf_objects = [getattr(ROOT, "%s_pdf" % p)() for p in pdf_functions]

    grid = ROOT.appl.grid(params.inputs[0])
    for name in params.inputs[1:]:
        g = ROOT.appl.grid(name)
        grid += g
        g = None
        gc.collect()

    print "Writing out %s" % params.output
    grid.Write(params.output)


def usage():
    print """\
Usage: grid_combine [OPTION...] [FILE]
Add together several grid files

  -o, --output='name.root'  Output file name
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

        self.output = 'combined.root'
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
