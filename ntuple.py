#!/usr/bin/env python

import ROOT

# load fastjet
ROOT.gSystem.Load("libfastjet.so")

# we want to handle Ctrl+C
sh = ROOT.TSignalHandler(ROOT.kSigInterrupt, False)
sh.Add()
sh.Connect("Notified()", "TROOT", ROOT.gROOT, "SetInterrupt()")

chain = ROOT.TChain("t3")
chain.Add("B*.root")
#f = ROOT.TFile('Bevent1.root')
#t3 = f.Get('t3')
#t3.SetBranchStatus("*", 0)
#t3.SetBranchStatus("id", 1);
#t3.SetBranchStatus("weight", 1);
#t3.SetBranchStatus("nuwgt", 1);

ROOT.gROOT.LoadMacro("T3selector.C+")
selector = ROOT.T3selector()

maxpt = ROOT.std.vector(int)()
maxpt.push_back(1000)
maxpt.push_back(700)

selector.analysis.setN(3)
selector.analysis.addPtLinearHistograms("hist_linear_20.txt", 20, maxpt)
selector.analysis.addEtaLinearHistograms("hist_linear_20.txt", 20)
selector.analysis.addPtLinearHistograms("hist_linear_40.txt", 40, maxpt)
selector.analysis.addEtaLinearHistograms("hist_linear_40.txt", 40)

chain.Process(selector, "Test analysis! algo=antikt Pt=30 Eta=2.8")


