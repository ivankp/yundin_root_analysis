#!/usr/bin/env python

import ROOT

# load fastjet
ROOT.gSystem.Load("libfastjet.so")

# we want to handle Ctrl+C
sh = ROOT.TSignalHandler(ROOT.kSigInterrupt, False)
sh.Add()
sh.Connect("Notified()", "TROOT", ROOT.gROOT, "SetInterrupt()")

chain = ROOT.TChain("t3")
#chain.Add("B*.root")
chain.Add("B_HT0.5_mstw08_j80.0_j60.0_e2.8/*.root")

#f = ROOT.TFile('Bevent1.root')
#t3 = f.Get('t3')
#t3.SetBranchStatus("*", 0)
#t3.SetBranchStatus("id", 1);
#t3.SetBranchStatus("weight", 1);
#t3.SetBranchStatus("nuwgt", 1);

ROOT.gROOT.LoadMacro("T3selector.C+")
selector = ROOT.T3selector()

maxpt = ROOT.std.vector(int)()
for pt in [1420, 1400, 800, 800, 800]:
  maxpt.push_back(pt)

selector.analysis.setN(2) # IMPORTANT!!! Set it first!
selector.analysis.init("Test_analysis! algo=antikt R=0.4 Pt=60 Pt1=80 Eta=2.8")
#selector.analysis.addPtLinearHistograms("hist_linear_20.txt", 20, maxpt)
#selector.analysis.addEtaLinearHistograms("hist_linear_20.txt", 20)
#selector.analysis.addPtLinearHistograms("hist_linear_40.txt", 40, maxpt)
#selector.analysis.addEtaLinearHistograms("hist_linear_40.txt", 40)
selector.analysis.addPtLinearHistograms("hist_nail_cmp.txt", 64, maxpt)
selector.analysis.addEtaLinearHistograms("hist_nail_cmp.txt", 20)

chain.Process(selector) #, "Test_analysis! algo=antikt R=0.4 Pt=60 Pt1=80 Eta=2.8") # call init above


