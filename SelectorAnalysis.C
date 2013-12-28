
// --------------------------------------------------------------------------- //
// Parameters
// --------------------------------------------------------------------------- //

Analysis::Analysis()
  : jet_number(0), jet_ptmin(0), jet_etamax(100),
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

  for (unsigned i=0; i<jet_pt_n.size(); i++) {
    clear_histvec(jet_pt_n[i]);
  }
  jet_pt_n.clear();

  for (unsigned i=0; i<jet_eta_n.size(); i++) {
    clear_histvec(jet_eta_n[i]);
  }
  jet_eta_n.clear();
}

void Analysis::reset()
{
  clear();

  event_count = 0.;
  event_binned = 0.;

  jet_exclusive = new LinearHistogram("!", "jet_multi_exclusive", jet_number+3, -0.5, jet_number+3-0.5);
  jet_inclusive = new LinearHistogram("!", "jet_multi_inclusive", jet_number+3, -0.5, jet_number+3-0.5);

  jet_pt_n.resize(jet_number+1);
  jet_eta_n.resize(jet_number+1);

  g_jet_pt_n.resize(jet_number+1);
  g_jet_eta_n.resize(jet_number+1);
}

void Analysis::clear_histvec(std::vector<Histogram*>& histvec)
{
  for (unsigned j=0; j<histvec.size(); j++) {
    clear_var(histvec[j]);
  }
  histvec.clear();
}

void Analysis::append_output_filename(const TString& name)
{
  outputfiles.insert(name);
}

void Analysis::output_histvec(const std::vector<Histogram*>& histvec,
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

void Analysis::bin_histvec(const std::vector<Histogram*>& histvec,
                          int nextevt, double x, double w)
{
  for (unsigned i=0; i<histvec.size(); i++) {
    histvec[i]->bin(nextevt, x, w);
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


bool Analysis::check_cuts(SelectorCommon* event)
{
  input.clear();
  jets.clear();

  PseudoJetVector jetinput;
  for (int i=0; i<event->nparticle; i++) {
    const fastjet::PseudoJet vec = fastjet::PseudoJet(event->px[i], event->py[i], event->pz[i], event->E[i]);
    input.push_back(vec);
    if (event->kf[i] == 21 or abs(event->kf[i]) <= 6) {
      jetinput.push_back(vec);
    }
  }
  if (jetinput.size() > 0) {
    fastjet::ClusterSequence cs(jetinput, jet_definition);
    PseudoJetVector alljets = fastjet::sorted_by_pt(cs.inclusive_jets(jet_ptmin));
    for (unsigned i=0; i<alljets.size(); i++) {
      if (abs(alljets[i].eta()) <= jet_etamax) {
        jets.push_back(alljets[i]);
      }
    }
  }

  return jets.size() >= jet_number;
}

void Analysis::fill_grid(Grid* grid, int nextevt, double x, double w, SelectorCommon* event)
{
  static const int lhaids[] = {5,-5,4,-4,3,-3,2,-2,1,-1};

  if (grid) {
    grid->fill(nextevt,
               event->lhaid1, event->lhaid2,
               event->x1, event->x2, event->fac_scale,
               event->pdfx1, event->pdfx2,
               x, event->naked_weight, w, event->event_order());
    if (event->coll_weights_count == 8) { // I-part, bin collinear CTs separately
      const int id1 = event->lhaid1;
      const int id2 = event->lhaid2;
      const double x1 = event->x1;
      const double x2 = event->x2;
      const double x1r = event->x1r;
      const double x2r = event->x2r;
      const double fac_scale = event->fac_scale;
      const double* pdfx1 = event->pdfx1;
      const double* pdfx2 = event->pdfx2;
      const double* pdfx1p = event->pdfx1p;
      const double* pdfx2p = event->pdfx2p;
      const int order = event->event_order();
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


void Analysis::analysis_bin(SelectorCommon* event)
{
  const Int_t id = event->id;
  const Double_t weight = event->weight;
  const Double_t scale = event->ren_scale;

  event_binned += 1;

  jet_exclusive->bin(id, jets.size(), weight);
  for (unsigned i=0; i<=jets.size(); i++) {
    jet_inclusive->bin(id, i, weight);
    fill_grid(g_jet_inclusive, id, i, weight, event);
  }

  bin_histvec(scale_wgt, id, scale, weight);
  bin_histvec(scale_nowgt, id, scale, Double_t(1.));

  for (unsigned i=0; i<jets.size(); i++) {
    assert(i<jet_pt_n.size());
    const double jetpt = jets[i].pt();
    const double jeteta = jets[i].eta();
    bin_histvec(jet_pt_n[i], id, jetpt, weight);
    bin_histvec(jet_eta_n[i], id, jeteta, weight);
    fill_grid(g_jet_pt_n[i], id, jetpt, weight, event);
    fill_grid(g_jet_eta_n[i], id, jeteta, weight, event);
  }
}

void Analysis::analysis_finalize()
{
  std::cout << "Finalize " << long(event_count) << " events (binned " << long(event_binned) << ")\n";
  std::set<TString>::iterator it;
  std::ofstream null;
  output_histograms(*it, null, true); // dryrun to get outputfiles
  for (it=outputfiles.begin(); it!=outputfiles.end(); ++it) {
    std::cout << "File: " << it->Data() << std::endl;
    std::ofstream outfile = std::ofstream(it->Data());
    outfile.precision(15);
    outfile.setf(std::ios::scientific);

    output_histograms(*it, outfile, false);

    outfile.close();
  }
  output_grids();
}

void Analysis::output_histograms(const TString& filename, std::ofstream& stream, bool dryrun)
{
  jet_exclusive->print(stream, runname, event_count);
  jet_inclusive->print(stream, runname, event_count);

  output_histvec(scale_wgt, filename, stream, dryrun);
  output_histvec(scale_nowgt, filename, stream, dryrun);

  for (unsigned i=0; i<=jet_number; i++) {
    output_histvec(jet_pt_n[i], filename, stream, dryrun);
    output_histvec(jet_eta_n[i], filename, stream, dryrun);
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
  : jet_pt1min(0)
{
}

bool JetAnalysis::check_cuts(SelectorCommon* event)
{
  if (not Analysis::check_cuts(event)) {
    return false;
  }

  return jets[0].pt() >= jet_pt1min;
}

void JetAnalysis::analysis_bin(SelectorCommon* event)
{
//   const Int_t id = event->id;
//   const Double_t weight = event->weight;
  Analysis::analysis_bin(event);
}

void JetAnalysis::output_histograms(const TString& filename, std::ofstream& stream, bool dryrun)
{
  // all jet histograms are already in the base class
  Analysis::output_histograms(filename, stream, dryrun);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//   PhotonJetAnalysis
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

PhotonJetAnalysis::PhotonJetAnalysis()
  : jet_pt1min(0), photon_R(0), photon_ptmin(0), photon_etamax(0),
    photon_jet_Rsep(0)
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

bool PhotonJetAnalysis::check_cuts(SelectorCommon* event)
{
  if (not Analysis::check_cuts(event)) {
    return false;
  }

  assert(event->kf[0] == 22 and event->kf[1] != 22);  // safe-guard against di-photon

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

  // leading jet pt-cut
  return jets[0].pt() >= jet_pt1min;
}

void PhotonJetAnalysis::analysis_bin(SelectorCommon* event)
{
  Analysis::analysis_bin(event);

  const Int_t id = event->id;
  const Double_t weight = event->weight;

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
    photon_R(0), photon_pt1min(0), photon_pt2min(0), photon_etamax(0),
    photon_photon_Rsep(0), photon_jet_Rsep(0)
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

bool DiPhotonAnalysis::check_cuts(SelectorCommon* event)
{
  if (not Analysis::check_cuts(event)) {
    return false;
  }

  assert(event->kf[0] == 22 and event->kf[1] == 22);

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

  // leading jet pt-cut
  return jets[0].pt() >= jet_pt1min;
}

void DiPhotonAnalysis::analysis_bin(SelectorCommon* event)
{
  Analysis::analysis_bin(event);

  const Int_t id = event->id;
  const Double_t weight = event->weight;

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
