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

#include <LHAPDF.h>

// --------------------------------------------------------------------------- //
// Histograms
// --------------------------------------------------------------------------- //

Histogram::Histogram(const TString& filename_, const TString& name_,
                     int nbin_, double x1_, double x2_)
  : filename(filename_), name(name_),
    x1(x1_), x2(x2_), x12(x2-x1),
    nbin(nbin_), prevevt(-1), lastidx(0),
    curidx(nbin), events(nbin), curwgt(nbin),
    wgt(nbin), wgt2(nbin), bwidth(nbin), edge(nbin+1)
{ }

TString Histogram::getFile() const
{
  return filename;
}

void Histogram::print(std::ostream& stream, const TString& runname, double count, bool unweight)
{
  flush();
  stream << "# BEGIN HISTOGRAM /" << runname << "/" << name << std::endl;
  stream << "AidaPath=/" << runname << "/" << name << std::endl;
  stream << "Title=" << std::endl;

  double area = 0.;
  for (int i=0; i<nbin; i++) {
    area += wgt[i];
  }
  area /= count;
  stream << "## Area: " << area << std::endl;
  stream << "## Num bins: " << nbin << std::endl;
  stream << "## xlow  \txhigh   \tval    \terrminus\terrplus" << std::endl;

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

    stream << edge[i] << "\t"
          << edge[i+1] << "\t"
          << value/bwidth[i] << "\t"
          << error/bwidth[i] << "\t"
          << error/bwidth[i] << std::endl;
  }
  stream << "# END HISTOGRAM" << std::endl << std::endl << std::endl;
}

void Histogram::flush()
{
  for (int i=0; i<lastidx; i++) {
    int idx = curidx[i];
    wgt[idx] += curwgt[i];
    wgt2[idx] += curwgt[i]*curwgt[i];
    events[idx]++;
    curwgt[i] = 0.;
  }
  lastidx = 0;
}

void Histogram::setedges() {
  edge[0] = x1;
  for (int i=0; i<nbin; i++) {
    edge[i+1] = edge[i] + bwidth[i];
  }
}

inline
void Histogram::fill(int evt, int n, double w)
{
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

LinearHistogram::LinearHistogram(const TString& filename_, const TString& name_,
                                 int nbin_, double x1_, double x2_)
  : Histogram(filename_, name_, nbin_, x1_, x2_)
{
  const double step = x12/nbin;
  for (int i=0; i<nbin; i++) {
    bwidth[i] = step;
  }
  setedges();
}

void LinearHistogram::bin(int nextevt, double x, double w)
{
//       std::cout << name << ": E(" << evt << ") LE (" << lastevt << ") LI(" << lastidx << ")" << std::endl; std::cout.flush();
  if (x < x1 or x > x2) return;
  int n = static_cast<int>(nbin*(x-x1)/x12);
  assert(0 <= n && n < nbin);
  fill(nextevt, n, w);
}

QuadraticHistogram::QuadraticHistogram(const TString& filename_, const TString& name_,
                                       int nbin_, double x1_, double x2_, double f)
  : Histogram(filename_, name_, nbin_, x1_, x2_),
    step(2*x12/((1 + f)*nbin)),
    slope((f - 1)/(nbin - 1))
{
  for (int i=0; i<nbin; i++) {
    bwidth[i] = step*(1 + i*slope);
  }
  setedges();
}

void QuadraticHistogram::bin(int nextevt, double x, double w)
{
//       std::cout << name << ": E(" << evt << ") LE (" << lastevt << ") LI(" << lastidx << ")" << std::endl; std::cout.flush();
  if (x < x1 or x > x2) return;
  const double dn = (slope-2 + sqrt((slope-2)*(slope-2) + (8*slope*(x - x1))/step))/(2*slope);
  const int n = static_cast<int>(dn);
  assert(0 <= n && n < nbin);
  fill(nextevt, n, w);
}

// --------------------------------------------------------------------------- //
// Parameters
// --------------------------------------------------------------------------- //

T3analysis::T3analysis()
  : jet_exclusive(0), jet_inclusive(0)
{
  init(TString(""));
}

void T3analysis::clear()
{
  if (jet_exclusive) {
    delete jet_exclusive;
    jet_exclusive = 0;
  }
  if (jet_inclusive) {
    delete jet_inclusive;
    jet_inclusive = 0;
  }

  for (unsigned i=0; i<jet_pt_n.size(); i++) {
    for (unsigned j=0; j<jet_pt_n[i].size(); j++) {
      if (jet_pt_n[i][j]) {
        delete jet_pt_n[i][j];
      }
    }
    jet_pt_n[i].clear();
  }
  jet_pt_n.clear();

  for (unsigned i=0; i<jet_eta_n.size(); i++) {
    for (unsigned j=0; j<jet_eta_n[i].size(); j++) {
      if (jet_eta_n[i][j]) {
        delete jet_eta_n[i][j];
      }
    }
    jet_eta_n[i].clear();
  }
  jet_eta_n.clear();
}

void T3analysis::reset()
{
  clear();
  jet_exclusive = new LinearHistogram("!", "jet_multi_exclusive", njet+3, -0.5, njet+3-0.5);
  jet_inclusive = new LinearHistogram("!", "jet_multi_inclusive", njet+3, -0.5, njet+3-0.5);
  jet_pt_n.resize(njet+1);
  jet_eta_n.resize(njet+1);
  event_count = 0.;
}

void T3analysis::addPtLinearHistograms(TString filename, int nbins, std::vector<int> ptlimits)
{
  double maxpt = 2000;
  for (unsigned i=0; i<=njet; i++) {
    if (i < ptlimits.size()) {
      maxpt = ptlimits[i];
    }
    std::stringstream ptname;
    ptname << "jet_pT_" << i+1;
    jet_pt_n[i].push_back(new LinearHistogram(filename, ptname.str(), nbins, ptmin, maxpt));
  }
  outputfiles.insert(filename);
}

void T3analysis::addPtQuadraticHistograms(TString filename, int nbins, double f, std::vector<int> ptlimits)
{
  double maxpt = 2000;
  for (unsigned i=0; i<=njet; i++) {
    if (i < ptlimits.size()) {
      maxpt = ptlimits[i];
    }
    std::stringstream ptname;
    ptname << "jet_pT_" << i+1;
    jet_pt_n[i].push_back(new QuadraticHistogram(filename, ptname.str(), nbins, ptmin, maxpt, f));
  }
  outputfiles.insert(filename);
}

void T3analysis::addEtaLinearHistograms(TString filename, int nbins)
{
  for (unsigned i=0; i<=njet; i++) {
    std::stringstream etaname;
    etaname << "jet_eta_" << i+1;
    jet_eta_n[i].push_back(new LinearHistogram(filename, etaname.str(), nbins, -etacut, etacut));
  }
  outputfiles.insert(filename);
}

void T3analysis::addEtaQuadraticHistograms(TString filename, int nbins, double f)
{
  for (unsigned i=0; i<=njet; i++) {
    std::stringstream etaname;
    etaname << "jet_eta_" << i+1;
    jet_eta_n[i].push_back(new QuadraticHistogram(filename, etaname.str(), nbins, -etacut, etacut, f));
  }
  outputfiles.insert(filename);
}

double T3analysis::OptF(const TString s)
{
  return TString(s(s.Index("=")+1, s.Length())).Atof();
}

int T3analysis::OptI(const TString s)
{
  return TString(s(s.Index("=")+1, s.Length())).Atoi();
}

void T3analysis::setN(const unsigned njet_)
{
  njet = njet_;
  reset();
}

void T3analysis::init(const TString& opt)
{
  ptmin = 30.;
  jet1ptmin = ptmin;
  etacut = 5.;
  jetR = 0.4;
  jetalgo = fastjet::antikt_algorithm;

  if (opt.Length()) {

    TSubString subopt = opt(TRegexp("^[^!]+"));
    if (subopt.Length()) {
      runname = TString(subopt);
      std::cout << "Runname: " << runname << std::endl;
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
      std::cout << "algo = kt" << std::endl;
      jetalgo = fastjet::kt_algorithm;
    }
    if (opt.Contains("algo=antikt")) {
      std::cout << "algo = antikt" << std::endl;
      jetalgo = fastjet::antikt_algorithm;
    }
  }

  jetdef = fastjet::JetDefinition(jetalgo, jetR);
}

void T3analysis::analysis_bin(const Int_t id, const Double_t weight,
                              const std::vector<fastjet::PseudoJet>& jets)
{
  jet_exclusive->bin(id, jets.size(), weight);
  for (unsigned i=1; i<=jets.size(); i++) {
    jet_inclusive->bin(id, i, weight);
  }
  for (unsigned i=0; i<jets.size(); i++) {
    assert(i<jet_pt_n.size());
    for (unsigned k=0; k<jet_pt_n[i].size(); k++) {
      jet_pt_n[i][k]->bin(id, jets[i].pt(), weight);
    }
    for (unsigned k=0; k<jet_eta_n[i].size(); k++) {
      jet_eta_n[i][k]->bin(id, jets[i].eta(), weight);
    }
  }
}

void T3analysis::analysis_finalize()
{
  std::cout << "Finalize\n";
  std::set<TString>::iterator it;
  for (it=outputfiles.begin(); it!=outputfiles.end(); ++it) {
    std::cout << "File: " << it->Data() << std::endl;
    std::ofstream outfile = std::ofstream(it->Data());
    outfile.precision(7);
    outfile.setf(std::ios::scientific);
    jet_exclusive->print(outfile, runname, event_count);
    jet_inclusive->print(outfile, runname, event_count);

    for (unsigned i=0; i<=njet; i++) {
      for (unsigned k=0; k<jet_pt_n[i].size(); k++) {
        if (*it == jet_pt_n[i][k]->getFile()) {
          jet_pt_n[i][k]->print(outfile, runname, event_count);
        }
      }
      for (unsigned k=0; k<jet_eta_n[i].size(); k++) {
        if (*it == jet_eta_n[i][k]->getFile()) {
          jet_eta_n[i][k]->print(outfile, runname, event_count);
        }
      }
    }
    outfile.close();
  }
}

// --------------------------------------------------------------------------- //
// Selector
// --------------------------------------------------------------------------- //

void T3selector::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

//   TString option = GetOption();

}

void T3selector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  fastjet::ClusterSequence::set_fastjet_banner_stream(0); // silence fastjet
  TString option = GetOption();
}

static LHAPDF::Flavour pdg2lha(int pdgnum)
{
  switch (pdgnum) {
    case 21:
      return LHAPDF::GLUON;
    case 1:
      return LHAPDF::DOWN;
    case -1:
      return LHAPDF::DBAR;
    case 2:
      return LHAPDF::UP;
    case -2:
      return LHAPDF::UBAR;
    case 3:
      return LHAPDF::STRANGE;
    case -3:
      return LHAPDF::SBAR;
    case 4:
      return LHAPDF::CHARM;
    case -4:
      return LHAPDF::CBAR;
    case 5:
      return LHAPDF::BOTTOM;
    case -5:
      return LHAPDF::BBAR;
    case 6:
      return LHAPDF::TOP;
    case -6:
      return LHAPDF::TBAR;
    default:
      throw;
  }
}

double T3selector::beta0pole2()
{
  static const double CA = 3.;
  static const double CF = 4./3.;
  static const double Nf = 5.;

  static const double b0 = 0.5*(11./3.*CA - 2./3.*Nf);

  double sum = 0.;
  for (int i=0; i<nparticle; i++) {
    sum += kf[i] == 21 ? CA : CF;
  }
  sum += id1 == 21 ? CA : CF;
  sum += id2 == 21 ? CA : CF;

  return b0/(-2.*sum);
}

void T3selector::reweight()
{
//   fac_scale *= 2.;
//   ren_scale *= 2.;
  if (FROMPDF == TOPDF and scalefactor == 1.) {
    return;
  }
  const int flav1 = pdg2lha(id1);
  const int flav2 = pdg2lha(id2);

  const double fx1 = LHAPDF::xfx(FROMPDF, x1, fac_scale, flav1)/x1;
  const double fx2 = LHAPDF::xfx(FROMPDF, x2, fac_scale, flav2)/x2;

  const double calc_weight = me_wgt*(fx1*fx2);
  if (std::fabs((calc_weight - weight)/weight) > 1e-10 and nuwgt != 18) {  // this simple check does not work for I-part
    std::cout << "Check your FROMPDF! (" << id << ") " << calc_weight << " != " << weight << std::endl;
  }

  const double new_fx1 = LHAPDF::xfx(TOPDF, x1, fac_scale*scalefactor, pdg2lha(id1))/x1;
  const double new_fx2 = LHAPDF::xfx(TOPDF, x2, fac_scale*scalefactor, pdg2lha(id2))/x2;

  double alphafactor = 0;
  if (use_sherpa_alphas) {
    alphafactor = sherpa_alphas->AlphaS(ren_scale*ren_scale*scalefactor*scalefactor)/alphas;
  } else {
    alphafactor = LHAPDF::alphasPDF(TOPDF, ren_scale*scalefactor)/alphas;
  }
  if (nuwgt == 0) {
    weight = me_wgt*(new_fx1*new_fx2);
  } else if (nuwgt == 2) {
    const double lr = log(scalefactor*scalefactor);  // log(murnew^2/murold^2)

    if (beta0fix) {
      usr_wgts[0] -= (beta0fix - (alphapower-1))*2.*usr_wgts[1]*beta0pole2();
    }

    weight = (me_wgt  + usr_wgts[0]*lr + 0.5*usr_wgts[1]*lr*lr)*(new_fx1*new_fx2);
  } else if (nuwgt == 18) {
    const double lr = log(scalefactor*scalefactor);  // log(murnew^2/murold^2)
    const double lf = log(scalefactor*scalefactor);  // log(mufnew^2/mufold^2)
    double w[8];
    for (int i=0; i<8; i++) {
      w[i] = usr_wgts[2+i] + usr_wgts[10+i]*lf;
    }
    double pdfx1[13]; LHAPDF::xfx(TOPDF, x1, fac_scale*scalefactor, pdfx1);
    double pdfx2[13]; LHAPDF::xfx(TOPDF, x2, fac_scale*scalefactor, pdfx2);

    double pdfx1p[13]; LHAPDF::xfx(TOPDF, x1/x1p, fac_scale*scalefactor, pdfx1p);
    double pdfx2p[13]; LHAPDF::xfx(TOPDF, x2/x2p, fac_scale*scalefactor, pdfx2p);

    double f1[4];
    double f2[4];

    // no top pdf below
    f1[0] = (flav1 != 0) ? pdfx1[6+flav1] : (  pdfx1[7] + pdfx1[8] + pdfx1[9] + pdfx1[10] + pdfx1[11]
                                             + pdfx1[1] + pdfx1[2] + pdfx1[3] + pdfx1[4] + pdfx1[5]);
    f2[0] = (flav2 != 0) ? pdfx2[6+flav2] : (  pdfx2[7] + pdfx2[8] + pdfx2[9] + pdfx2[10] + pdfx2[11]
                                             + pdfx2[1] + pdfx2[2] + pdfx2[3] + pdfx2[4] + pdfx2[5]);

    f1[1] = (flav1 != 0) ? pdfx1p[6+flav1] : (  pdfx1p[7] + pdfx1p[8] + pdfx1p[9] + pdfx1p[10] + pdfx1p[11]
                                              + pdfx1p[1] + pdfx1p[2] + pdfx1p[3] + pdfx1p[4] + pdfx1p[5]);
    f2[1] = (flav2 != 0) ? pdfx2p[6+flav2] : (  pdfx2p[7] + pdfx2p[8] + pdfx2p[9] + pdfx2p[10] + pdfx2p[11]
                                              + pdfx2p[1] + pdfx2p[2] + pdfx2p[3] + pdfx2p[4] + pdfx2p[5]);

    f1[2] = pdfx1[6];
    f2[2] = pdfx2[6];

    f1[3] = pdfx1p[6];
    f2[3] = pdfx2p[6];

    if (beta0fix) {
      usr_wgts[0] -= (beta0fix - (alphapower-1))*2.*usr_wgts[1]*beta0pole2();
    }

    weight = (me_wgt  + usr_wgts[0]*lr + 0.5*usr_wgts[1]*lr*lr)*(new_fx1*new_fx2)
           + (w[0]*f1[0]/x1 + w[1]*f1[1]/x1 + w[2]*f1[2]/x1 + w[3]*f1[3]/x1)*new_fx2
           + (w[4]*f2[0]/x2 + w[5]*f2[1]/x2 + w[6]*f2[2]/x2 + w[7]*f2[3]/x2)*new_fx1;
  } else {
    std::cout << "Unknown value for nuwgt = " << nuwgt << std::endl;
  }
  weight *= pow(alphafactor, alphapower);
}

void T3selector::initAS(const int order, const double asMZ, const double mZ2,
                        const std::vector<double>& qmasses)
{
  sherpa_alphas = new SHERPA::One_Running_AlphaS(order, asMZ, mZ2, qmasses);
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
  // of thid and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  // check for Ctrl+C
  if (gROOT->IsInterrupted()) {
    Abort("Keyboard interrupt");
  }

  GetEntry(entry);
  analysis.event_count += 1;

  std::vector<fastjet::PseudoJet> input;
  for (int i=0; i<nparticle; i++) {
    input.push_back(fastjet::PseudoJet(px[i], py[i], pz[i], E[i]));
  }
  fastjet::ClusterSequence cs(input, analysis.jetdef);
  std::vector<fastjet::PseudoJet> fjets = fastjet::sorted_by_pt(cs.inclusive_jets(analysis.ptmin));
  std::vector<fastjet::PseudoJet> jets;
  for (unsigned i=0; i<fjets.size(); i++) {
    if (std::fabs(fjets[i].eta()) <= analysis.etacut) {
      jets.push_back(fjets[i]);
    }
  }
  if (jets.size() > 1 && jets[0].pt() >= analysis.jet1ptmin) {
    reweight(); // reweight event in-place
    analysis.analysis_bin(id, weight, jets);
  }

  return kTRUE;
}

void T3selector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  analysis.analysis_finalize();
}

void T3selector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
