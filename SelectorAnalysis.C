
#include <utility>

#include "FlavourKT.h"

static inline int getFlavour(const fastjet::PseudoJet& vec)
{
  return FlavourKTPlugin::getFlavour(vec);
}

// analog of delta_R using pseudorapidity
static double DeltaEtaPhi(const fastjet::PseudoJet& a, const fastjet::PseudoJet& b)
{
  double dphi = abs(a.phi() - b.phi());
  if (dphi > fastjet::pi) {
    dphi = fastjet::twopi - dphi;
  }
  const double deta = a.eta() - b.eta();
  return sqrt(dphi*dphi + deta*deta);
}

// --------------------------------------------------------------------------- //
// Parameters
// --------------------------------------------------------------------------- //

Analysis::Analysis()
  : jet_number(0), jet_ptmin(0), jet_etamax(0), jet_ymax(0),
    jet_exclusive(0), jet_inclusive(0), g_jet_inclusive(0)
{
}

Analysis::~Analysis()
{
  clear();
  clear_var(g_jet_inclusive);
}

template <typename T>
void Analysis::clear_var(T*& var)
{
  if (var) {
    delete var;
    var = 0;
  }
}

void Analysis::clear()
{
  outputfiles.clear();

  clear_var(jet_exclusive);
  clear_var(jet_inclusive);

  clear_histvec(scale_wgt);
  clear_histvec(scale_nowgt);
  clear_histvec(jet_ht);

  clear_histvec(jet_pt_n);
  clear_histvec(jet_eta_n);
  clear_histvec(jet_y_n);

}

void Analysis::reset()
{
  clear();

  call_count = 0.;
  event_count = 0.;
  event_binned = 0.;

  jet_exclusive = new LinearHistogram("!", "jet_multi_exclusive", jet_number+3, -0.5, jet_number+3-0.5);
  jet_inclusive = new LinearHistogram("!", "jet_multi_inclusive", jet_number+3, -0.5, jet_number+3-0.5);

  jet_pt_n.resize(jet_number+1);
  jet_eta_n.resize(jet_number+1);
  jet_y_n.resize(jet_number+1);

  g_jet_pt_n.resize(jet_number+1);
  g_jet_eta_n.resize(jet_number+1);
}

void Analysis::clear_histvec(std::vector<HistogramBase*>& histvec)
{
  for (unsigned j=0; j<histvec.size(); j++) {
    clear_var(histvec[j]);
  }
  histvec.clear();
}

void Analysis::clear_histvec(std::vector<std::vector<HistogramBase*> >& histvecvec)
{
  for (unsigned n=0; n<histvecvec.size(); n++) {
    clear_histvec(histvecvec[n]);
  }
  histvecvec.clear();
}

void Analysis::append_output_filename(const TString& name)
{
  outputfiles.insert(name);
}

void Analysis::output_histvec(const std::vector<HistogramBase*>& histvec,
                              const TString& filename, std::ofstream& stream,
                              bool dryrun)
{
  for (unsigned i=0; i<histvec.size(); i++) {
    if (dryrun) {
      append_output_filename(histvec[i]->getFile());
    } else {
      if (filename == histvec[i]->getFile()) {
        histvec[i]->print(stream, runname, event_count);
      }
    }
  }
}

void Analysis::bin_histvec(const std::vector<HistogramBase*>& histvec,
                          int nextevt, double x, double w)
{
  for (unsigned i=0; i<histvec.size(); i++) {
    histvec[i]->bin(nextevt, x, w);
  }
}

void Analysis::bin_histvec(const std::vector<HistogramBase*>& histvec,
                          int nextevt, double x, double y, double w)
{
  for (unsigned i=0; i<histvec.size(); i++) {
    histvec[i]->bin2d(nextevt, x, y, w);
  }
}

// Setting parameters

void Analysis::setJetNumber(const unsigned n)
{
  jet_number = n;
  reset();
}

void Analysis::setAntiKt(double R)
{
  jet_definition = fastjet::JetDefinition(fastjet::antikt_algorithm, R);
}

void Analysis::setKt(double R)
{
  jet_definition = fastjet::JetDefinition(fastjet::kt_algorithm, R);
}

void Analysis::set_input(PseudoJetVector newinput)
{
  input.clear();
  input.swap(newinput);
}

bool Analysis::check_cuts(const SelectorCommon* /*event*/)
{
  jets.clear();

  PseudoJetVector jetinput;
  for (unsigned i=0; i<input.size(); i++) {
    const int lhid = getFlavour(input[i]);
    // gluons, light quarks, everything above 80
    if (lhid == 21 or abs(lhid) <= 5 or abs(lhid) >= 81) {
      jetinput.push_back(input[i]);
    }
  }
  if (jetinput.size() > 0) {
    fastjet::ClusterSequence cs(jetinput, jet_definition);
    PseudoJetVector alljets = fastjet::sorted_by_pt(cs.inclusive_jets(jet_ptmin));
    for (unsigned i=0; i<alljets.size(); i++) {
      if ((jet_etamax == 0. or abs(alljets[i].eta()) <= jet_etamax) and
          (jet_ymax == 0. or abs(alljets[i].rap()) <= jet_ymax)) {
        jets.push_back(alljets[i]);
      }
    }
  }

  return jets.size() >= jet_number;
}

#ifdef DISABLE_APPLGRID
void Analysis::fill_grid(Grid* /*grid*/, int /*nextevt*/, double /*x*/, double /*w*/, const SelectorCommon* /*event*/)
{
}
#else // ! DISABLE_APPLGRID
void Analysis::fill_grid(Grid* grid, int nextevt, double x, double w, const SelectorCommon* event)
{
  static const int lhaids[] = {5,-5,4,-4,3,-3,2,-2,1,-1};

  if (grid) {
    const int id1 = event->get_lhaid1();
    const int id2 = event->get_lhaid2();
    const double x1 = event->get_x1();
    const double x2 = event->get_x2();
    const double x1r = event->get_x1r();
    const double x2r = event->get_x2r();
    const double fac_scale = event->get_fac_scale();
    const double* pdfx1 = event->pdfx1;
    const double* pdfx2 = event->pdfx2;
    const double* pdfx1p = event->pdfx1p;
    const double* pdfx2p = event->pdfx2p;
    const int order = event->get_event_order();
    grid->fill(nextevt, id1, id2, x1, x2, fac_scale, pdfx1, pdfx2,
               x, event->naked_weight, w, order);
    if (event->coll_weights_count == 8) { // I-part, bin collinear CTs separately
      // f1q
      // f1qx
      if (id1 != 0) {
        grid->fill(nextevt, id1, id2, x1, x2, fac_scale, pdfx1, pdfx2,
                   x, event->coll_weights[0], 0, order);
        grid->fill(nextevt, id1, id2, x1r, x2, fac_scale, pdfx1p, pdfx2,
                   x, event->coll_weights[1], 0, order);
      } else {
        if (event->coll_weights[0] != 0.) {
          for (int i=0; i<10; i++) {
            grid->fill(nextevt, lhaids[i], id2, x1, x2, fac_scale, pdfx1, pdfx2,
                       x, event->coll_weights[0], 0, order);
          }
        }
        if (event->coll_weights[1] != 0.) {
          for (int i=0; i<10; i++) {
            grid->fill(nextevt, lhaids[i], id2, x1r, x2, fac_scale, pdfx1p, pdfx2,
                       x, event->coll_weights[1], 0, order);
          }
        }
      }
      // f1g
      grid->fill(nextevt, 0, id2, x1, x2, fac_scale, pdfx1, pdfx2,
                 x, event->coll_weights[2], 0, order);
      // f1gx
      grid->fill(nextevt, 0, id2, x1r, x2, fac_scale, pdfx1p, pdfx2,
                 x, event->coll_weights[3], 0, order);
      // f2q
      // f2qx
      if (id2 != 0) {
        grid->fill(nextevt, id1, id2, x1, x2, fac_scale, pdfx1, pdfx2,
                   x, event->coll_weights[4+0], 0, order);
        grid->fill(nextevt, id1, id2, x1, x2r, fac_scale, pdfx1, pdfx2p,
                   x, event->coll_weights[4+1], 0, order);
      } else {
        if (event->coll_weights[4+0] != 0.) {
          for (int i=0; i<10; i++) {
            grid->fill(nextevt, id1, lhaids[i], x1, x2, fac_scale, pdfx1, pdfx2,
                       x, event->coll_weights[4+0], 0, order);
          }
        }
        if (event->coll_weights[4+1] != 0.) {
          for (int i=0; i<10; i++) {
            grid->fill(nextevt, id1, lhaids[i], x1, x2r, fac_scale, pdfx1, pdfx2p,
                       x, event->coll_weights[4+1], 0, order);
          }
        }
      }
      // f2g
      grid->fill(nextevt, id1, 0, x1, x2, fac_scale, pdfx1, pdfx2,
                 x, event->coll_weights[4+2], 0, order);
      // f2gx
      grid->fill(nextevt, id1, 0, x1, x2r, fac_scale, pdfx1, pdfx2p,
                 x, event->coll_weights[4+3], 0, order);
    }
  }
}
#endif // DISABLE_APPLGRID


bool Analysis::photonIsolation(const SelectorCommon* event, double photon_R,
                               double photon_n, double photon_eps) const
{
  for (int i=0; i<event->get_nparticle(); i++) {
    if (getFlavour(input[i]) != 22) continue;
    const double Eeps = input[i].Et()*photon_eps/pow(1. - cos(photon_R), photon_n);
    std::vector<std::pair<double,double> > hadronic;
    for (int j=0; j<event->get_nparticle(); j++) {
      if (abs(getFlavour(input[j])) <= 6) {
        const double Rij = input[i].delta_R(input[j]);
        if (Rij <= photon_R) {
          hadronic.push_back(std::make_pair(Rij, input[j].Et()));
        }
      }
    }
    sort(hadronic.begin(), hadronic.end());
    double energy = 0.;
    for (unsigned j=0; j<hadronic.size(); j++) {
      energy += hadronic[j].second;
      if (energy > Eeps*pow(1. - cos(hadronic[j].first), photon_n)) {
        return false;
      }
    }
  }
  return true;
}


void Analysis::analysis_bin(const SelectorCommon* event)
{
  const Int_t id = event->get_event_id();
  const Double_t weight = event->get_event_weight();
  const Double_t scale = event->get_ren_scale();

  event_binned += 1;

  jet_exclusive->bin(id, jets.size(), weight);
  for (unsigned i=0; i<=jets.size(); i++) {
    jet_inclusive->bin(id, i, weight);
    fill_grid(g_jet_inclusive, id, i, weight, event);
  }

  bin_histvec(scale_wgt, id, scale, weight);
  bin_histvec(scale_nowgt, id, scale, Double_t(1.));

  double jetht = 0.;
  for (unsigned i=0; i<jets.size(); i++) {
    const double jetpt = jets[i].pt();
    jetht += jetpt;

    if (i >= jet_pt_n.size()) continue;
    const double jeteta = jets[i].eta();
    const double jety = jets[i].rapidity();
    bin_histvec(jet_pt_n[i], id, jetpt, weight);
    bin_histvec(jet_eta_n[i], id, jeteta, weight);
    bin_histvec(jet_y_n[i], id, jety, weight);
    fill_grid(g_jet_pt_n[i], id, jetpt, weight, event);
    fill_grid(g_jet_eta_n[i], id, jeteta, weight, event);
  }
  bin_histvec(jet_ht, id, jetht, weight);
}

void Analysis::analysis_finalize(const SelectorCommon* event)
{
  finalize_stat(std::cout, event);

  std::set<TString>::iterator it;
  std::ofstream null;
  output_histograms(*it, null, true); // dryrun to get outputfiles
  for (it=outputfiles.begin(); it!=outputfiles.end(); ++it) {
    std::cout << "File: " << it->Data() << std::endl;
    std::ofstream outfile(it->Data());
    outfile.precision(15);
    outfile.setf(std::ios::scientific);

    output_histograms(*it, outfile, false);

    finalize_stat(outfile << "### ", event);
    outfile.close();
  }
  output_grids();
}

void Analysis::finalize_stat(std::ostream& stream, const SelectorCommon* event)
{
  stream << "Finalize: "
         << "groups " << event->event_groups
         << ", events " << long(call_count)
         << " (binned " << long(event_binned) << ")"
         << " [trials " << long(event_count) << "]"
         << std::endl;
}

void Analysis::output_histograms(const TString& filename, std::ofstream& stream, bool dryrun)
{
  jet_exclusive->print(stream, runname, event_count);
  jet_inclusive->print(stream, runname, event_count);

  output_histvec(scale_wgt, filename, stream, dryrun);
  output_histvec(scale_nowgt, filename, stream, dryrun);
  output_histvec(jet_ht, filename, stream, dryrun);

  for (unsigned i=0; i<=jet_number; i++) {
    output_histvec(jet_pt_n[i], filename, stream, dryrun);
    output_histvec(jet_eta_n[i], filename, stream, dryrun);
    output_histvec(jet_y_n[i], filename, stream, dryrun);
  }
}

void Analysis::output_grids()
{
  if (g_jet_inclusive) {
    g_jet_inclusive->write(event_count);
  }
  for (unsigned i=0; i<=jet_number; i++) {
    if (g_jet_pt_n[i]) {
      g_jet_pt_n[i]->write(event_count);
    }
    if (g_jet_eta_n[i]) {
      g_jet_eta_n[i]->write(event_count);
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//   JetAnalysis
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

JetAnalysis::JetAnalysis()
  : jet_pt1min(0.), jet_pt2min(0.), jet_pt3min(0.), jet_ht2min(0)
{
}

void JetAnalysis::clear()
{
  Analysis::clear();

  clear_histvec(jet_pt12ave);
}

bool JetAnalysis::check_cuts(const SelectorCommon* event)
{
  if (not Analysis::check_cuts(event)) {
    return false;
  }

  double ht2 = jets[0].pt()+jets[1].pt();
  if (ht2 < jet_ht2min) {
    return false;
  }

  if (jet_number >= 1 and jet_pt1min > 0. and jets[0].pt() < jet_pt1min) {
    return false;
  }
  if (jet_number >= 2 and jet_pt2min > 0. and jets[1].pt() < jet_pt2min) {
    return false;
  }
  if (jet_number >= 3 and jet_pt3min > 0. and jets[2].pt() < jet_pt3min) {
    return false;
  }

  return true;
}

void JetAnalysis::analysis_bin(const SelectorCommon* event)
{
  const Int_t id = event->get_event_id();
  const Double_t weight = event->get_event_weight();
  Analysis::analysis_bin(event);

  double pt12ave = 0.5*(jets[0].pt()+jets[1].pt());
  bin_histvec(jet_pt12ave, id, pt12ave, weight);
}

void JetAnalysis::output_histograms(const TString& filename, std::ofstream& stream, bool dryrun)
{
  Analysis::output_histograms(filename, stream, dryrun);
  output_histvec(jet_pt12ave, filename, stream, dryrun);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//   Jet3Analysis -- three jet analysis with beta23 observable
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Jet3Analysis::Jet3Analysis()
  : jet_eta1max(10), jet_eta2max(10), jet_eta2min(0),
    jet_jet_DR23min(0), jet_jet_DR23max(10), jet_jet_M12min(0)
{
}

void Jet3Analysis::clear()
{
  BaseClass::clear();

  clear_histvec(jet_jet_eta23);
  clear_histvec(jet_jet_phi23);
  clear_histvec(jet_jet_beta23);
}

bool Jet3Analysis::check_cuts(const SelectorCommon* event)
{
  if (not JetAnalysis::check_cuts(event)) {
    return false;
  }

  if (abs(jets[0].eta()) > jet_eta1max) {
    return false;
  }

  if (abs(jets[1].eta()) > jet_eta2max or abs(jets[1].eta()) < jet_eta2min) {
    return false;
  }

  double M12 = (jets[0]+jets[1]).m();
  if (M12 < jet_jet_M12min) {
    return false;
  }

  double DR23 = jets[1].delta_R(jets[2]);
  if (DR23 < jet_jet_DR23min or DR23 > jet_jet_DR23max) {
    return false;
  }

  return true;
}

void Jet3Analysis::analysis_bin(const SelectorCommon* event)
{
  JetAnalysis::analysis_bin(event);

  const Int_t id = event->get_event_id();
  const Double_t weight = event->get_event_weight();

  if (jets.size() >= 3) {
    const double jet2phi = jets[1].phi();
    const double jet3phi = jets[2].phi();
    double jj_phi23 = jet3phi - jet2phi;
    if (jj_phi23 > M_PI) {
      jj_phi23 = jj_phi23 - 2.*M_PI;
    } else if (jj_phi23 < -M_PI) {
      jj_phi23 = jj_phi23 + 2.*M_PI;
    }
    const double jet2eta = jets[1].eta();
    const double jet3eta = jets[2].eta();
    double jj_eta23 = jet3eta - jet2eta;
    if (jet2eta < 0) {
      jj_eta23 = -jj_eta23;
    }

    double jj_beta23 = atan(abs(jj_phi23)/jj_eta23);
    if (jj_beta23<0.) jj_beta23 += M_PI;

    bin_histvec(jet_jet_eta23, id, jj_eta23, weight);
    bin_histvec(jet_jet_phi23, id, jj_phi23, weight);
    bin_histvec(jet_jet_beta23, id, jj_beta23, weight);
  }
}

void Jet3Analysis::output_histograms(const TString& filename, std::ofstream& stream, bool dryrun)
{
  JetAnalysis::output_histograms(filename, stream, dryrun);
  output_histvec(jet_jet_eta23, filename, stream, dryrun);
  output_histvec(jet_jet_phi23, filename, stream, dryrun);
  output_histvec(jet_jet_beta23, filename, stream, dryrun);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//   JetMAnalysis
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

JetMAnalysis::JetMAnalysis()
  : ystar_min(0), ystar_max(1e10)
{
  g_jet_mass_jjj = 0;
}

void JetMAnalysis::clear()
{
  BaseClass::clear();

  clear_histvec(jet_mass_jjj);
}

bool JetMAnalysis::check_cuts(const SelectorCommon* event)
{
  if (not BaseClass::check_cuts(event)) {
    return false;
  }

  if (jets.size() >= 3) {
    const double y1 = jets[0].rap();
    const double y2 = jets[1].rap();
    const double y3 = jets[2].rap();
    const double ystar = abs(y1 - y2) + abs(y2 - y3) + abs(y1 - y3);
    if (ystar < ystar_min or ystar >= ystar_max) {
    return false;
    }
  }

  return true;
}

void JetMAnalysis::analysis_bin(const SelectorCommon* event)
{
  BaseClass::analysis_bin(event);

  const Int_t id = event->get_event_id();
  const Double_t weight = event->get_event_weight();

  if (jets.size() >= 3) {
    const double mass_jjj = (jets[0] + jets[1] + jets[2]).m();
    bin_histvec(jet_mass_jjj, id, mass_jjj, weight);
    fill_grid(g_jet_mass_jjj, id, mass_jjj, weight, event);
  }
}

void JetMAnalysis::output_histograms(const TString& filename, std::ofstream& stream, bool dryrun)
{
  BaseClass::output_histograms(filename, stream, dryrun);

  output_histvec(jet_mass_jjj, filename, stream, dryrun);
}

void JetMAnalysis::output_grids()
{
  BaseClass::output_grids();

  if (g_jet_mass_jjj) {
    g_jet_mass_jjj->write(event_count);
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//   FourJetMPIAnalysis -- four jet analysis for dble diff. d12 d34 observable
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

FourJetMPIAnalysis::FourJetMPIAnalysis() :
  mpivars_d12_bins(1), mpivars_d12_bin_low(0.), mpivars_d12_bin_high(100.)
{
}

void FourJetMPIAnalysis::clear()
{
  BaseClass::clear();

  clear_histvec(jets_d12_d34);
}

bool FourJetMPIAnalysis::check_cuts(const SelectorCommon* event)
{
  if (not JetAnalysis::check_cuts(event)) {
    return false;
  }

  return true;
}

void FourJetMPIAnalysis::analysis_bin(const SelectorCommon* event)
{
  JetAnalysis::analysis_bin(event);

  const Int_t id = event->get_event_id();
  const Double_t weight = event->get_event_weight();

  if (jets.size() == 4) {
    double dij[6];
    for (int i=0; i<3; i++) {
      for (int j=i+1; j<4; j++) {
        const double jet1phi = jets[i].phi();
        const double jet2phi = jets[j].phi();
        double jj_phi12 = jet1phi - jet2phi;
        if (jj_phi12 > M_PI) {
          jj_phi12 = jj_phi12 - 2.*M_PI;
        } else if (jj_phi12 < -M_PI) {
          jj_phi12 = jj_phi12 + 2.*M_PI;
        }
        const double pt1 = jets[i].pt();
        const double pt2 = jets[j].pt();
        dij[i+j*(j-1)/2] = pt1*pt1 + pt2*pt2 + pt1*pt2*cos(jj_phi12);
      }
    }

    double d12 = dij[0];
    unsigned int pos1 = 0, pos2 = 1;
    for (int i=0; i<3; i++) {
      for (int j=i+1; j<4; j++) {
        int idx = i+j*(j-1)/2;
        if (dij[idx]<d12) {
          d12 = dij[idx];
          pos1 = i;
          pos2 = j;
        }
      }
    }

    double d34 = dij[5-pos1-pos2*(pos2-1)/2];

    bin_histvec(jets_d12_d34, id, d12, d34, weight);
  }
}

void FourJetMPIAnalysis::output_histograms(const TString& filename, std::ofstream& stream, bool dryrun)
{
  JetAnalysis::output_histograms(filename, stream, dryrun);
  output_histvec(jets_d12_d34, filename, stream, dryrun);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//   VJetAnalysis
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

VJetAnalysis::VJetAnalysis()
  : lepton_ptmin(0), lepton_etamax(0),
    lepton_etagap_min(0), lepton_etagap_max(0),
    etmiss_min(0),
    vboson_mass_min(0), vboson_mass_max(0),
    lepton_jet_Rsep(0), lepton_lepton_Rsep(0)
{
  g_vboson_pt = 0;
  g_vboson_eta = 0;
}

void VJetAnalysis::clear()
{
  Analysis::clear();

  clear_histvec(vboson_pt);
  clear_histvec(vboson_eta);
}

void VJetAnalysis::reset()
{
  Analysis::reset();
}

bool VJetAnalysis::check_cuts(const SelectorCommon* event)
{
  if (not Analysis::check_cuts(event)) {
    return false;
  }

  // check for W/Z bosons in first two flavours
  PseudoJetVector leptons;
  const int flavour0 = getFlavour(input[0]);
  const int flavour1 = getFlavour(input[1]);
  if (flavour0 == 11 || flavour0 == -11) {
    leptons.push_back(input[0]);
  }
  if (flavour1 == 11 || flavour1 == -11) {
    leptons.push_back(input[1]);
  }

  PseudoJetVector neutrino;
  if (flavour0 == 12 || flavour0 == -12) {
    neutrino.push_back(input[0]);
  }
  if (flavour1 == 12 || flavour1 == -12) {
    neutrino.push_back(input[1]);
  }
  if (leptons.size() + neutrino.size() == 2) {
    vboson = input[0] + input[1];
  } else {
    bool found = false;
    for (unsigned i = 0; i < input.size(); i++) {
      const int flav_i = getFlavour(input[i]);
      if (flav_i == 23 || flav_i == 24) {
        vboson = input[i];
        found = true;
        break;
      }
    }
    if (not found) {
      return false;
    }
  }

  // missing energy/neutrino cuts
  for (unsigned int j=0; j<neutrino.size(); j++) {
    if (neutrino[j].pt() < etmiss_min) {
      return false;
    }
  }

  // lepton cuts
  for (unsigned int j=0; j<leptons.size(); j++) {
    if (leptons[j].pt() < lepton_ptmin) {
      return false;
    }
    if (abs(leptons[j].eta()) > lepton_etamax) {
      return false;
    }
    if (abs(leptons[j].eta()) > lepton_etagap_min && abs(leptons[j].eta()) < lepton_etagap_max) {
      return false;
    }

    // lepton-jets R-separation
    for (unsigned i=0; i<jets.size(); i++) {
      if (leptons[j].delta_R(jets[i]) < lepton_jet_Rsep) {
        return false;
      }
    }
  }

  if (leptons.size() == 2) {
    if (leptons[0].delta_R(leptons[1]) < lepton_lepton_Rsep) {
      return false;
    }
  }

  // vector boson mass cut
  double vmass = vboson.m();
  if (vmass < vboson_mass_min || vmass > vboson_mass_max) {
    return false;
  }

  return true;
}

void VJetAnalysis::analysis_bin(const SelectorCommon* event)
{
  Analysis::analysis_bin(event);

  const Int_t id = event->get_event_id();
  const Double_t weight = event->get_event_weight();

  const double LLpt = vboson.pt();
  const double LLeta = vboson.eta();

  bin_histvec(vboson_pt, id, LLpt, weight);
  bin_histvec(vboson_eta, id, LLeta, weight);

  fill_grid(g_vboson_pt, id, LLpt, weight, event);
  fill_grid(g_vboson_eta, id, LLeta, weight, event);
}

void VJetAnalysis::output_histograms(const TString& filename, std::ofstream& stream, bool dryrun)
{
  // all jet histograms are already in the base class
  Analysis::output_histograms(filename, stream, dryrun);

  output_histvec(vboson_pt, filename, stream, dryrun);
  output_histvec(vboson_eta, filename, stream, dryrun);
}

void VJetAnalysis::output_grids()
{
  Analysis::output_grids();

  if (g_vboson_pt) {
    g_vboson_pt->write(event_count);
  }
  if (g_vboson_eta) {
    g_vboson_eta->write(event_count);
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//   PhotonJetAnalysis
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

PhotonJetAnalysis::PhotonJetAnalysis()
  : jet_pt1min(0), photon_R(0), photon_n(0), photon_eps(0),
    photon_ptmin(0), photon_etamax(0), photon_jet_Rsep(0)
{
  g_photon_pt = 0;
  g_photon_eta = 0;
  g_photon_jet_R11 = 0;
  g_jet_jet_phi12 = 0;
}

void PhotonJetAnalysis::clear()
{
  Analysis::clear();

  clear_histvec(photon_pt);
  clear_histvec(photon_eta);
  clear_histvec(photon_jet_R11);
  clear_histvec(jet_jet_phi12);
}

void PhotonJetAnalysis::reset()
{
  Analysis::reset();
}

bool PhotonJetAnalysis::check_cuts(const SelectorCommon* event)
{
  if (not Analysis::check_cuts(event)) {
    return false;
  }

  assert(getFlavour(input[0]) == 22 and getFlavour(input[1]) != 22);  // safe-guard against di-photon

  const fastjet::PseudoJet& photon = input[0];

  // photon cuts
  if (photon.pt() < photon_ptmin) {
    return false;
  }
  if (abs(photon.eta()) > photon_etamax) {
    return false;
  }

  // keep only observed jets (some maybe hidden by the photon)
  PseudoJetVector obsjets;
  for (unsigned i=0; i<jets.size(); i++) {
    if (photon.delta_R(jets[i]) > photon_R) {
      obsjets.push_back(jets[i]);
    }
  }
  jets.swap(obsjets);
  if (jets.size() < jet_number) {
    return false;
  }

  // photon-jets R-separation
  for (unsigned i=0; i<jets.size(); i++) {
    if (photon.delta_R(jets[i]) < photon_jet_Rsep) {
      return false;
    }
  }

  // photon isolation
  if (photon_n > 0 and photon_eps > 0) {
    if (not photonIsolation(event, photon_R, photon_n, photon_eps)) {
      return false;
    }
  }

  // leading jet pt-cut
  return jets[0].pt() >= jet_pt1min;
}

void PhotonJetAnalysis::analysis_bin(const SelectorCommon* event)
{
  Analysis::analysis_bin(event);

  const Int_t id = event->get_event_id();
  const Double_t weight = event->get_event_weight();

  const fastjet::PseudoJet& photon = input[0];
  const double AApt = photon.pt();
  const double AAeta = photon.eta();

  bin_histvec(photon_pt, id, AApt, weight);
  bin_histvec(photon_eta, id, AAeta, weight);

  fill_grid(g_photon_pt, id, AApt, weight, event);
  fill_grid(g_photon_eta, id, AAeta, weight, event);

  if (jets.size() >= 1) {
    const double R11 = photon.delta_R(jets[0]);
    bin_histvec(photon_jet_R11, id, R11, weight);
    fill_grid(g_photon_jet_R11, id, R11, weight, event);
  }

  if (jets.size() >= 2) {
    const double jet1phi = jets[0].phi();
    const double jet2phi = jets[1].phi();
    double jj_phi12 = jet1phi > jet2phi ? jet1phi - jet2phi : jet2phi - jet1phi;
    if (jj_phi12 > M_PI) {
      jj_phi12 = 2.*M_PI - jj_phi12;
    }
    bin_histvec(jet_jet_phi12, id, jj_phi12, weight);
    fill_grid(g_jet_jet_phi12, id, jj_phi12, weight, event);
  }
}

void PhotonJetAnalysis::output_histograms(const TString& filename, std::ofstream& stream, bool dryrun)
{
  // all jet histograms are already in the base class
  Analysis::output_histograms(filename, stream, dryrun);

  output_histvec(photon_pt, filename, stream, dryrun);
  output_histvec(photon_eta, filename, stream, dryrun);
  output_histvec(photon_jet_R11, filename, stream, dryrun);
  output_histvec(jet_jet_phi12, filename, stream, dryrun);
}

void PhotonJetAnalysis::output_grids()
{
  Analysis::output_grids();

  if (g_photon_pt) {
    g_photon_pt->write(event_count);
  }
  if (g_photon_eta) {
    g_photon_eta->write(event_count);
  }
  if (g_photon_jet_R11) {
    g_photon_jet_R11->write(event_count);
  }
  if (g_photon_jet_R11) {
    g_photon_jet_R11->write(event_count);
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//   DiPhotonAnalysis
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

DiPhotonAnalysis::DiPhotonAnalysis()
  : jet_pt1min(0),
    photon_R(0), photon_n(0), photon_eps(0),
    photon_pt1min(0), photon_pt2min(0), photon_etamax(0),
    photon_photon_Rsep(0), photon_jet_Rsep(0),
    photon_photon_mass_min(0), photon_photon_mass_max(1e10)
{
  g_photon_mass = 0;
  g_photon_pt = 0;
  g_photon_eta = 0;
  g_photon_jet_R11 = 0;
  g_jet_jet_phi12 = 0;
}

void DiPhotonAnalysis::clear()
{
  Analysis::clear();

  clear_histvec(photon_mass);
  clear_histvec(photon_pt);
  clear_histvec(photon_eta);
  clear_histvec(photon_jet_R11);
  clear_histvec(jet_jet_phi12);
}

void DiPhotonAnalysis::reset()
{
  Analysis::reset();
}

bool DiPhotonAnalysis::check_cuts(const SelectorCommon* event)
{
  if (not Analysis::check_cuts(event)) {
    return false;
  }

  if (not(getFlavour(input[0]) == 22 and getFlavour(input[1]) == 22)) {
    return false;
  }

  double pt1 = input[0].pt();
  double pt2 = input[1].pt();
  if (pt1 < pt2) {
    std::swap(pt1, pt2);
    std::swap(input[0], input[1]);
  }

  if (pt1 < photon_pt1min) {
    return false;
  }
  if (pt2 < photon_pt2min) {
    return false;
  }
  if (abs(input[0].eta()) > photon_etamax or
      abs(input[1].eta()) > photon_etamax) {
    return false;
  }
  if (input[0].delta_R(input[1]) < photon_photon_Rsep) {
    return false;
  }

  PseudoJetVector obsjets;
  for (unsigned i=0; i<jets.size(); i++) {
    if (input[0].delta_R(jets[i]) > photon_R and
        input[1].delta_R(jets[i]) > photon_R) {
      obsjets.push_back(jets[i]);
    }
  }
  jets.swap(obsjets);
  if (jets.size() < jet_number) {
    return false;
  }

  for (unsigned i=0; i<jets.size(); i++) {
    if (input[0].delta_R(jets[i]) < photon_jet_Rsep or
        input[1].delta_R(jets[i]) < photon_jet_Rsep) {
      return false;
    }
  }

  double AAmass = (input[0]+input[1]).m();
  if (AAmass < photon_photon_mass_min or AAmass > photon_photon_mass_max) {
    return false;
  }

  // photon isolation
  if (photon_n > 0 and photon_eps > 0) {
    if (not photonIsolation(event, photon_R, photon_n, photon_eps)) {
      return false;
    }
  }

  // leading jet pt-cut
  if (jets.size() > 0) {
    return jets[0].pt() >= jet_pt1min;
  }
  return true;
}

void DiPhotonAnalysis::analysis_bin(const SelectorCommon* event)
{
  Analysis::analysis_bin(event);

  const Int_t id = event->get_event_id();
  const Double_t weight = event->get_event_weight();

  const fastjet::PseudoJet AAmom = input[0] + input[1];
  const double AAmass = AAmom.m();
  const double AApt = AAmom.pt();
  const double AAeta = AAmom.eta();

  bin_histvec(photon_mass, id, AAmass, weight);
  bin_histvec(photon_pt, id, AApt, weight);
  bin_histvec(photon_eta, id, AAeta, weight);

  fill_grid(g_photon_mass, id, AAmass, weight, event);
  fill_grid(g_photon_pt, id, AApt, weight, event);
  fill_grid(g_photon_eta, id, AAeta, weight, event);

  if (jets.size() >= 1) {
    const double R11 = input[0].delta_R(jets[0]);
    bin_histvec(photon_jet_R11, id, R11, weight);
    fill_grid(g_photon_jet_R11, id, R11, weight, event);
  }

  if (jets.size() >= 2) {
    double jj_phi12 = abs(jets[0].phi() - jets[1].phi());
    jj_phi12 = jj_phi12 > M_PI ? 2.*M_PI - jj_phi12 : jj_phi12;
    bin_histvec(jet_jet_phi12, id, jj_phi12, weight);
    fill_grid(g_jet_jet_phi12, id, jj_phi12, weight, event);

    const fastjet::PseudoJet JJmom = jets[0] + jets[1];
    const double JJmass = JJmom.m();
    bin_histvec(jet_jet_mass, id, JJmass, weight);

    const double jj_eta12 = abs(jets[0].eta() - jets[1].eta());
    bin_histvec(jet_jet_eta12, id, jj_eta12, weight);

    double aa_jj_phi12 = abs(JJmom.phi() - AAmom.phi());
    aa_jj_phi12 = aa_jj_phi12 > M_PI ? 2.*M_PI - aa_jj_phi12 : aa_jj_phi12;
    bin_histvec(diphoton_dijet_phi, id, aa_jj_phi12, weight);

    const double aa_j_j_ystar = AAmom.rapidity() - (jets[0].rapidity() + jets[1].rapidity())/2.;
    bin_histvec(diphoton_dijet_ystar, id, aa_j_j_ystar, weight);
  }
}

void DiPhotonAnalysis::output_histograms(const TString& filename, std::ofstream& stream, bool dryrun)
{
  // all jet histograms are already in the base class
  Analysis::output_histograms(filename, stream, dryrun);

  output_histvec(photon_mass, filename, stream, dryrun);
  output_histvec(photon_pt, filename, stream, dryrun);
  output_histvec(photon_eta, filename, stream, dryrun);
  output_histvec(photon_jet_R11, filename, stream, dryrun);
  output_histvec(jet_jet_phi12, filename, stream, dryrun);
  output_histvec(jet_jet_mass, filename, stream, dryrun);
  output_histvec(jet_jet_eta12, filename, stream, dryrun);
  output_histvec(diphoton_dijet_phi, filename, stream, dryrun);
  output_histvec(diphoton_dijet_ystar, filename, stream, dryrun);
}

void DiPhotonAnalysis::output_grids()
{
  Analysis::output_grids();

  if (g_photon_mass) {
    g_photon_mass->write(event_count);
  }
  if (g_photon_pt) {
    g_photon_pt->write(event_count);
  }
  if (g_photon_eta) {
    g_photon_eta->write(event_count);
  }
  if (g_photon_jet_R11) {
    g_photon_jet_R11->write(event_count);
  }
  if (g_photon_jet_R11) {
    g_photon_jet_R11->write(event_count);
  }
}

bool DiPhotonAnalysisBH::check_cuts(const SelectorCommon* event)
{
  if (not Analysis::check_cuts(event)) {
    return false;
  }

  assert(getFlavour(input[0]) == 22 and getFlavour(input[1]) == 22);

  double pt1 = input[0].pt();
  double pt2 = input[1].pt();
  if (pt1 < pt2) {
    std::swap(pt1, pt2);
    std::swap(input[0], input[1]);
  }

  if (pt1 < photon_pt1min) {
    return false;
  }
  if (pt2 < photon_pt2min) {
    return false;
  }
  if (abs(input[0].eta()) > photon_etamax or
      abs(input[1].eta()) > photon_etamax) {
    return false;
  }
  if (input[0].delta_R(input[1]) < photon_photon_Rsep) {
    return false;
  }

  for (unsigned i=0; i<jet_number; i++) {
    if (DeltaEtaPhi(input[0], jets[i]) < photon_jet_Rsep or
        DeltaEtaPhi(input[1], jets[i]) < photon_jet_Rsep) {
      return false;
    }
  }
  for (unsigned i=jet_number; i<jets.size(); i++) {
    if (DeltaEtaPhi(input[0], jets[i]) < photon_R or
        DeltaEtaPhi(input[1], jets[i]) < photon_R) {
      return false;
    }
  }

  // photon isolation
  if (photon_n > 0 and photon_eps > 0) {
    if (not photonIsolation(event, photon_R, photon_n, photon_eps)) {
      return false;
    }
  }

  // leading jet pt-cut
  return jets[0].pt() >= jet_pt1min;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//   HiggsJetsAnalysis
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

HiggsJetsAnalysis::HiggsJetsAnalysis()
  : min_dijet_m(0), min_dijet_y(0)
{
}

void HiggsJetsAnalysis::clear()
{
  Analysis::clear();

  clear_histvec(jet_jet_mass_ij);
  clear_histvec(jet_jet_dy_ij);
  clear_histvec(jet_jet_dphi_ij);
  clear_histvec(jet_jet_dR_ij);

  clear_histvec(higgs_pt);
  clear_histvec(higgs_eta);
  clear_histvec(higgs_y);
  clear_histvec(jjj_ystar);

  clear_histvec(higgs_dijet_pt_ij);
  clear_histvec(higgs_dijet_dphi_ij);
  clear_histvec(higgs_dijet_ystar_ij);
}

void HiggsJetsAnalysis::reset()
{
  Analysis::reset();

  const int n = (jet_number+1)*jet_number/2;

  jet_jet_mass_ij.resize(n);
  jet_jet_dy_ij.resize(n);
  jet_jet_dphi_ij.resize(n);
  jet_jet_dR_ij.resize(n);

  higgs_dijet_pt_ij.resize(n);
  higgs_dijet_dphi_ij.resize(n);
  higgs_dijet_ystar_ij.resize(n);
}

bool HiggsJetsAnalysis::check_cuts(const SelectorCommon* event)
{
  if (not Analysis::check_cuts(event)) {
    return false;
  }

  assert(getFlavour(input[0]) == 25);

  if (min_dijet_m != 0. and
      (jets[0] + jets[1]).m() < min_dijet_m) {
    return false;
  }
  if (min_dijet_y != 0. and
      abs(jets[0].rapidity() - jets[1].rapidity()) < min_dijet_y) {
    return false;
  }

  return true;
}

void HiggsJetsAnalysis::analysis_bin(const SelectorCommon* event)
{
  Analysis::analysis_bin(event);

  const Int_t id = event->get_event_id();
  const Double_t weight = event->get_event_weight();

  const fastjet::PseudoJet& Hmom = input[0];
  const double Hpt = Hmom.pt();
  const double Heta = Hmom.eta();
  const double Hy = Hmom.rapidity();

  bin_histvec(higgs_pt, id, Hpt, weight);
  bin_histvec(higgs_eta, id, Heta, weight);
  bin_histvec(higgs_y, id, Hy, weight);

  if (jets.size() >= 3) {
    const double ystar = 0.5*(jets[0].rap() + jets[1].rap()) - jets[2].rap();
    bin_histvec(jjj_ystar, id, ystar, weight);
  }

  for (unsigned j=1; j<jets.size(); j++) {
    for (unsigned i=0; i<j; i++) {
      const fastjet::PseudoJet& jet1 = jets[i];
      const fastjet::PseudoJet& jet2 = jets[j];
      const double jjm = (jet1 + jet2).m();
      const double jjdphi = abs(jet1.delta_phi_to(jet2));
      const double jjdy = abs(jet1.rapidity() - jet2.rapidity());
      const double jjdR = jet1.delta_R(jet2);
      const double Hjjpt = (Hmom + jet1 + jet2).pt();
      const double Hjjdphi = abs(Hmom.delta_phi_to(jet1 + jet2));
      const double Hjjystar = abs(Hmom.rap() - 0.5*(jet1.rap() + jet2.rap()));
      const int n = (j-1)*j/2 + i;
      bin_histvec(jet_jet_mass_ij[n], id, jjm, weight);
      bin_histvec(jet_jet_dy_ij[n], id, jjdy, weight);
      bin_histvec(jet_jet_dphi_ij[n], id, jjdphi, weight);
      bin_histvec(jet_jet_dR_ij[n], id, jjdR, weight);
      bin_histvec(higgs_dijet_pt_ij[n], id, Hjjpt, weight);
      bin_histvec(higgs_dijet_dphi_ij[n], id, Hjjdphi, weight);
      bin_histvec(higgs_dijet_ystar_ij[n], id, Hjjystar, weight);
    }
  }
}

void HiggsJetsAnalysis::output_histograms(const TString& filename, std::ofstream& stream, bool dryrun)
{
  // all jet histograms are already in the base class
  Analysis::output_histograms(filename, stream, dryrun);

  output_histvec(higgs_pt, filename, stream, dryrun);
  output_histvec(higgs_eta, filename, stream, dryrun);
  output_histvec(higgs_y, filename, stream, dryrun);
  output_histvec(jjj_ystar, filename, stream, dryrun);

  for (unsigned j=1; j<jet_number+1; j++) {
    for (unsigned i=0; i<j; i++) {
      const int n = (j-1)*j/2 + i;
      output_histvec(jet_jet_mass_ij[n], filename, stream, dryrun);
      output_histvec(jet_jet_dy_ij[n], filename, stream, dryrun);
      output_histvec(jet_jet_dphi_ij[n], filename, stream, dryrun);
      output_histvec(jet_jet_dR_ij[n], filename, stream, dryrun);
      output_histvec(higgs_dijet_pt_ij[n], filename, stream, dryrun);
      output_histvec(higgs_dijet_dphi_ij[n], filename, stream, dryrun);
      output_histvec(higgs_dijet_ystar_ij[n], filename, stream, dryrun);
    }
  }
}
