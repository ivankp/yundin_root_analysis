import ROOT

# main function which is called from hammer.py
def initialize(params, selector):
    print "Using jets analysis %s" % __file__

    analysis = ROOT.JetAnalysis.create()
    analysis.runname = params.runname
    analysis.setJetNumber(params.njet)
    analysis.setAntiKt(0.4)
    analysis.jet_ptmin = 60
    analysis.jet_etamax = 2.8
    analysis.jet_pt1min = 80

    maxpt = ROOT.std.vector('double')()
    for pt in [1420, 1400, 800, 800, 800]:
        maxpt.push_back(pt)

    # add histograms
    analysis.addPtLinearHistograms(params.output % "l64_l20", 64, maxpt)
    analysis.addEtaLinearHistograms(params.output % "l64_l20", 20)

    # assign to selector
    selector.analysis = analysis


if __name__ == '__main__':
    print "Analysis module can't be run alone"
    print "pass --analysis=%s to hammer.py" % __file__
