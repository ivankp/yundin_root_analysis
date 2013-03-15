#define T3selector_cxx
// The class definition in T3selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("T3selector.C")
// Root > T->Process("T3selector.C","some options")
// Root > T->Process("T3selector.C+")
//

#include "T3selector.h"
#include <TH2.h>
#include <TStyle.h>

#ifndef NDEBUG
  #include <iostream>
#endif

// --------------------------------------------------------------------------- //
// Histograms
// --------------------------------------------------------------------------- //

class Histogram
{
  protected:
    Histogram(string name_, int nbin_, double x1_, double x2_)
      : name(name_), nbin(nbin_), prevevt(-1), lastidx(0),
      x1(x1_), x2(x2_), x12(x2-x1),
      curidx(nbin), events(nbin), curwgt(nbin),
      wgt(nbin), wgt2(nbin), bwidth(nbin), edge(nbin+1)
    { }

  public:
    void print(string runname, double count, bool unweight=true) {
      flush();
      cout << "# BEGIN HISTOGRAM /" << runname << "/" << name << endl;
      cout << "AidaPath=/" << runname << "/" << name << endl;
      cout << "Title=" << endl;
      double area = 0.;
      for (int i=0; i<nbin; i++) {
        area += wgt[i];
      }
      area /= count;
      cout << "## Area: " << area << endl;
      cout << "## Num bins: " << nbin << endl;
      cout << "## xlow  \txhigh   \tval    \terrminus\terrplus" << endl;

      for (int i=0; i<nbin; i++) {
        const double value = wgt[i]/count;
        double error = 0;
        if (unweight) {
          double variance = wgt2[i]; //events[i] > 1 ? (wgt2[i]-wgt[i]*wgt[i]/events[i])/(events[i]-1) : wgt[i]*wgt[i];
          error = variance > 0. ? sqrt(variance)/count : 0.;
        } else {
          double variance = (wgt2[i]-wgt[i]*wgt[i]/count)/(count-1);
          error = variance > 0. ? sqrt(variance/count) : 0.;
        }

        cout << edge[i] << "\t"
             << edge[i+1] << "\t"
             << value/bwidth[i] << "\t"
             << error/bwidth[i] << "\t"
             << error/bwidth[i] << endl;
      }
      cout << "# END HISTOGRAM" << endl << endl << endl;
    }

  protected:
    void flush() {
      for (int i=0; i<lastidx; i++) {
        int idx = curidx[i];
        wgt[idx] += curwgt[i];
        wgt2[idx] += curwgt[i]*curwgt[i];
        events[idx]++;
        curwgt[i] = 0.;
      }
      lastidx = 0;
    }

    void setedges() {
      edge[0] = x1;
      for (int i=0; i<nbin; i++) {
        edge[i+1] = edge[i] + bwidth[i];
      }
    }

    inline
    void fill(int evt, int n, double w) {
      assert(0 <= n && n < nbin);

      if (evt != prevevt) {
        flush();
        prevevt = evt;
      }

      int i;
      for (i=0; i<lastidx; i++) {
        if (curidx[i] == n) {
          break;
        }
      }
      curidx[i] = n;
      curwgt[i] += w;
      if (i == lastidx) {
        lastidx++;
      }
    }

    string name;
    double x1, x2, x12;
    const int nbin;
    int prevevt, lastidx;
    vector<int> curidx;
    vector<int> events;
    vector<double> curwgt;
    vector<double> wgt;
    vector<double> wgt2;
    vector<double> bwidth;
    vector<double> edge;
};

class LinearHistogram : public Histogram
{
  public:
    LinearHistogram(string name_, int nbin_, double x1_, double x2_)
      : Histogram(name_, nbin_, x1_, x2_)
    {
      const double step = x12/nbin;
      for (int i=0; i<nbin; i++) {
        bwidth[i] = step;
      }
      setedges();
    }

    void bin(int nextevt, double x, double w) {
//       cout << name << ": E(" << evt << ") LE (" << lastevt << ") LI(" << lastidx << ")" << endl; cout.flush();
      if (x < x1 or x > x2) return;
      int n = static_cast<int>(nbin*(x-x1)/x12);
      assert(0 <= n && n < nbin);
      fill(nextevt, n, w);
    }
};


class QuadraticHistogram : public Histogram
{
  public:
    QuadraticHistogram(string name_, int nbin_, double x1_, double x2_, double f)
      : Histogram(name_, nbin_, x1_, x2_),
        step(2*x12/((1 + f)*nbin)),
        slope((f - 1)/(nbin - 1))
    {
      for (int i=0; i<nbin; i++) {
        bwidth[i] = step*(1 + i*slope);
      }
      setedges();
    }

    void bin(int nextevt, double x, double w) {
//       cout << name << ": E(" << evt << ") LE (" << lastevt << ") LI(" << lastidx << ")" << endl; cout.flush();
      if (x < x1 or x > x2) return;
      const double dn = (slope-2 + sqrt((slope-2)*(slope-2) + (8*slope*(x - x1))/step))/(2*slope);
      const int n = static_cast<int>(dn);
      assert(0 <= n && n < nbin);
      fill(nextevt, n, w);
    }
  private:
    const double step;
    const double slope;
};

// --------------------------------------------------------------------------- //
// Parameters
// --------------------------------------------------------------------------- //

double T3params::OptF(const TString s)
{
  return TString(s(s.Index("=")+1, s.Length())).Atof();
}

int T3params::OptI(const TString s)
{
  return TString(s(s.Index("=")+1, s.Length())).Atoi();
}

void T3params::init(const TString* opt_)
{
  njet = 2;
  ptmin = 30.;
  jet1ptmin = ptmin;
  etacut = 5.;
  jetR = 0.4;
  jetalgo = fastjet::antikt_algorithm;

  if (opt_) {
    const TString& opt = *opt_;

    TSubString subopt = opt(TRegexp("N=[0-9]+"));
    if (subopt.Length()) {
      njet = OptI(subopt);
      std::cout << "njet = " << njet << std::endl;
    }

    subopt = opt(TRegexp("Pt=[0-9.eE]+"));
    if (subopt.Length()) {
      ptmin = OptF(subopt);
      std::cout << "ptmin = " << ptmin << std::endl;
    }

    subopt = opt(TRegexp("Pt1=[0-9.eE]+"));
    if (subopt.Length()) {
      jet1ptmin = OptF(subopt);
      std::cout << "jet1ptmin = " << jet1ptmin << std::endl;
    }

    subopt = opt(TRegexp("Eta=[0-9.eE]+"));
    if (subopt.Length()) {
      etacut = OptF(subopt);
      std::cout << "etacut = " << etacut << std::endl;
    }

    subopt = opt(TRegexp("R=[0-9.eE]+"));
    if (subopt.Length()) {
      jetR = OptF(subopt);
      std::cout << "jetR = " << jetR << std::endl;
    }

    if (opt.Contains("algo=kt")) {
      jetalgo = fastjet::kt_algorithm;
    }
    if (opt.Contains("algo=antikt")) {
      jetalgo = fastjet::antikt_algorithm;
    }
  }
  jetdef = fastjet::JetDefinition(jetalgo, jetR);
}

// --------------------------------------------------------------------------- //
// Selector
// --------------------------------------------------------------------------- //

void T3selector::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}

void T3selector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  fastjet::ClusterSequence::set_fastjet_banner_stream(0); // silence fastjet
  TString option = GetOption();
  params.init(&option);

}

Bool_t T3selector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either T3selector::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  GetEntry(entry);

  std::vector<fastjet::PseudoJet> input;
  for (int i=0; i<nparticle; i++) {
    input.push_back(fastjet::PseudoJet(E[i], px[i], py[i], pz[i]));
  }
  fastjet::ClusterSequence cs(input, params.jetdef);
  vector<fastjet::PseudoJet> fjets = fastjet::sorted_by_pt(cs.inclusive_jets(params.ptmin));
  vector<fastjet::PseudoJet> jets;
  for (unsigned i=0; i<fjets.size(); i++) {
    if (std::fabs(fjets[i].eta()) <= params.etacut) {
      jets.push_back(fjets[i]);
    }
  }
  if (jets.size() > 1 && jets[0].pt() >= params.jet1ptmin) {
    std::cout << jets.size() << std::endl;
//     a.jet_exclusive->bin(e.event, jets.size(), e.weight);
//     for (int i=1; i<=jets.size(); i++) {
//       a.jet_inclusive->bin(e.event, i, e.weight);
//     }
//     for (int i=0; i<jets.size(); i++) {
//       assert(i<a.jet_pt_n.size());
//       a.jet_pt_n[i]->bin(e.event, jets[i].pt(), e.weight);
//       a.jet_eta_n[i]->bin(e.event, jets[i].eta(), e.weight);
//     }
  }

  return kTRUE;
}

void T3selector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void T3selector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
