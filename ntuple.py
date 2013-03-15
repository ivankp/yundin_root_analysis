#!/usr/bin/env python

import ROOT

ROOT.gSystem.Load('libfastjet.so')

chain = ROOT.TChain("t3")
chain.Add('RS*.root')
#f = ROOT.TFile('Bevent1.root')
#t3 = f.Get('t3')
#t3.SetBranchStatus("*", 0)
#t3.SetBranchStatus("id", 1);
#t3.SetBranchStatus("weight", 1);
#t3.SetBranchStatus("nuwgt", 1);

chain.Process("T3selector.C+", "N=3 algo=antikt Pt=30 Eta=2.8")

#seen = set()

#n = 0
#w = 0
#for e in chain:
    #n, w, nw = e.id, e.weight, e.nuwgt
    #if nw not in seen:
        #print nw
        #seen.add(nw)
    ##else:
        ##print n
