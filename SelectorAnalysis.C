
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
  clear_var(jet_exclusive);
  clear_var(jet_inclusive);

  for (unsigned i=0; i<jet_pt_n.size(); i++) {
    for (unsigned j=0; j<jet_pt_n[i].size(); j++) {
      clear_var(jet_pt_n[i][j]);
    }
    jet_pt_n[i].clear();
  }
  jet_pt_n.clear();

  for (unsigned i=0; i<jet_eta_n.size(); i++) {
    for (unsigned j=0; j<jet_eta_n[i].size(); j++) {
      clear_var(jet_eta_n[i][j]);
    }
    jet_eta_n[i].clear();
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
}

template <typename T>
void Analysis::addPtHistograms(TString filename, int nbins,
                               double param1, double param2, double param3,
                               double low, double high,
                               std::vector<double>* lowlimits,
                               std::vector<double>* highlimits)
{
  for (unsigned i=0; i<=jet_number; i++) {
    if (lowlimits and i < lowlimits->size()) {
      low = (*lowlimits)[i];
    }
    if (highlimits and i < highlimits->size()) {
      high = (*highlimits)[i];
    }
    std::stringstream name;
    name << "jet_pT_" << i+1;
    jet_pt_n[i].push_back(new T(filename, name.str(), nbins, low, high, param1, param2, param3));
  }
  outputfiles.insert(filename);
}

void Analysis::addPtLinearHistograms(TString filename, int nbins,
                                     std::vector<double>* ptlowlimits,
                                     std::vector<double>* pthighlimits)
{
  addPtHistograms<LinearHistogram>(filename, nbins, 0, 0, 0, jet_ptmin, 2000., ptlowlimits, pthighlimits);
}

void Analysis::addPtSmearedLinearHistograms(TString filename, int nbins, double s,
                                            std::vector<double>* ptlowlimits,
                                            std::vector<double>* pthighlimits)
{
  addPtHistograms<SmearedLinearHistogram>(filename, nbins, s, 0, 0, jet_ptmin, 2000., ptlowlimits, pthighlimits);
}

void Analysis::addPtQuadraticHistograms(TString filename, int nbins, double f,
                                        std::vector<double>* ptlowlimits,
                                        std::vector<double>* pthighlimits)
{
  addPtHistograms<QuadraticHistogram>(filename, nbins, f, 0, 0, jet_ptmin, 2000., ptlowlimits, pthighlimits);
}

void Analysis::addPtSmearedQuadraticHistograms(TString filename, int nbins, double f, double s,
                                               std::vector<double>* ptlowlimits,
                                               std::vector<double>* pthighlimits)
{
  addPtHistograms<SmearedQuadraticHistogram>(filename, nbins, f, s, 0, jet_ptmin, 2000., ptlowlimits, pthighlimits);
}

template <typename T>
void Analysis::addEtaHistograms(TString filename, int nbins,
                                double param1, double param2, double param3,
                                double low, double high,
                                std::vector<double>* lowlimits,
                                std::vector<double>* highlimits)
{
  for (unsigned i=0; i<=jet_number; i++) {
    if (lowlimits and i < lowlimits->size()) {
      low = (*lowlimits)[i];
    }
    if (highlimits and i < highlimits->size()) {
      high = (*highlimits)[i];
    }
    std::stringstream name;
    name << "jet_eta_" << i+1;
    jet_eta_n[i].push_back(new T(filename, name.str(), nbins, low, high, param1, param2, param3));
  }
  outputfiles.insert(filename);
}

void Analysis::addEtaLinearHistograms(TString filename, int nbins,
                                        std::vector<double>* etalowlimits,
                                        std::vector<double>* etahighlimits)
{
  addEtaHistograms<LinearHistogram>(filename, nbins, 0, 0, 0, -jet_etamax, jet_etamax, etalowlimits, etahighlimits);
}

void Analysis::addEtaSmearedLinearHistograms(TString filename, int nbins, double s,
                                        std::vector<double>* etalowlimits,
                                        std::vector<double>* etahighlimits)
{
  addEtaHistograms<SmearedLinearHistogram>(filename, nbins, s, 0, 0, -jet_etamax, jet_etamax, etalowlimits, etahighlimits);
}

void Analysis::addEtaQuadraticHistograms(TString filename, int nbins, double f,
                                        std::vector<double>* etalowlimits,
                                        std::vector<double>* etahighlimits)
{
  addEtaHistograms<QuadraticHistogram>(filename, nbins, f, 0, 0, -jet_etamax, jet_etamax, etalowlimits, etahighlimits);
}

void Analysis::addEtaSmearedQuadraticHistograms(TString filename, int nbins, double f, double s,
                                        std::vector<double>* etalowlimits,
                                        std::vector<double>* etahighlimits)
{
  addEtaHistograms<SmearedQuadraticHistogram>(filename, nbins, f, s, 0, -jet_etamax, jet_etamax, etalowlimits, etahighlimits);
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
  fastjet::ClusterSequence cs(jetinput, jet_definition);
  PseudoJetVector alljets = fastjet::sorted_by_pt(cs.inclusive_jets(jet_ptmin));
  for (unsigned i=0; i<alljets.size(); i++) {
    if (abs(alljets[i].eta()) <= jet_etamax) {
      jets.push_back(alljets[i]);
    }
  }

  return jets.size() >= jet_number;
}

void Analysis::analysis_bin(SelectorCommon* event)
{
  const Int_t id = event->id;
  const Double_t weight = event->weight;

  event_binned += 1;

  jet_exclusive->bin(id, jets.size(), weight);
  for (unsigned i=1; i<=jets.size(); i++) {
    jet_inclusive->bin(id, i, weight);
    if (g_jet_inclusive) {
      g_jet_inclusive->fill(id,
                          event->lhaid1, event->lhaid2,
                          event->x1, event->x2, event->fac_scale,
                          event->pdfx1, event->pdfx2,
                          i, event->naked_weight, weight, event->event_order());
    }
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

void Analysis::analysis_finalize()
{
  std::cout << "Finalize " << long(event_count) << " events (binned " << long(event_binned) << ")\n";
  std::set<TString>::iterator it;
  for (it=outputfiles.begin(); it!=outputfiles.end(); ++it) {
    std::cout << "File: " << it->Data() << std::endl;
    std::ofstream outfile = std::ofstream(it->Data());
    outfile.precision(10);
    outfile.setf(std::ios::scientific);

    output_histograms(*it, outfile);

    outfile.close();
  }
}

void Analysis::output_histograms(const TString& filename, std::ofstream& stream)
{
  jet_exclusive->print(stream, runname, event_count);
  jet_inclusive->print(stream, runname, event_count);

  if (g_jet_inclusive) {
    g_jet_inclusive->write(event_count);
  }

  for (unsigned i=0; i<=jet_number; i++) {
    for (unsigned k=0; k<jet_pt_n[i].size(); k++) {
      if (filename == jet_pt_n[i][k]->getFile()) {
        jet_pt_n[i][k]->print(stream, runname, event_count);
      }
    }
    for (unsigned k=0; k<jet_eta_n[i].size(); k++) {
      if (filename == jet_eta_n[i][k]->getFile()) {
        jet_eta_n[i][k]->print(stream, runname, event_count);
      }
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

void JetAnalysis::output_histograms(const TString& filename, std::ofstream& stream)
{
  // all jet histograms are already in the base class
  Analysis::output_histograms(filename, stream);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//   DiPhotonAnalysis
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

DiPhotonAnalysis::DiPhotonAnalysis()
  : photon_pt1min(0), photon_pt2min(0), photon_etamax(0),
    photon_photon_Rsep(0), photon_jet_Rsep(0),
    photon_mass(0), photon_jet_R11(0), jet_jet_phi12(0)
{
}

void DiPhotonAnalysis::clear()
{
  Analysis::clear();

  clear_var(photon_mass);
  clear_var(photon_jet_R11);
  clear_var(jet_jet_phi12);
}

void DiPhotonAnalysis::reset()
{
  Analysis::reset();

  jet_jet_phi12 = new LinearHistogram("!", "jet_jet_phi12", 31, 0, M_PI);
  photon_mass = new LinearHistogram("!", "photon_mass", 15, 0, 500);
  photon_jet_R11 = new LinearHistogram("!", "photon_jet_R11", 30, 0, 5);
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

  return true;
}

void DiPhotonAnalysis::analysis_bin(SelectorCommon* event)
{
  Analysis::analysis_bin(event);

  const Int_t id = event->id;
  const Double_t weight = event->weight;

  double mass = (input[0]+input[1]).m();
  photon_mass->bin(id, mass, weight);

  double R11 = input[0].delta_R(jets[0]);
  photon_jet_R11->bin(id, R11, weight);

  double jet1phi = jets[0].phi();
  double jet2phi = jets[1].phi();
  double jj_phi12 = jet1phi > jet2phi ? jet1phi - jet2phi : jet2phi - jet1phi;
  if (jj_phi12 > M_PI) {
    jj_phi12 = 2.*M_PI - jj_phi12;
  }
  jet_jet_phi12->bin(id, jj_phi12, weight);
}

void DiPhotonAnalysis::output_histograms(const TString& filename, std::ofstream& stream)
{
  // all jet histograms are already in the base class
  Analysis::output_histograms(filename, stream);

  photon_mass->print(stream, runname, event_count);
  photon_jet_R11->print(stream, runname, event_count);
  jet_jet_phi12->print(stream, runname, event_count);
}
