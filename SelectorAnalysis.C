
// --------------------------------------------------------------------------- //
// Parameters
// --------------------------------------------------------------------------- //

Analysis::Analysis()
  : jet_number(0), jet_ptmin(0), jet_etamax(100),
    jet_exclusive(0), jet_inclusive(0)
{
}

Analysis::~Analysis()
{
  clear();
}

void Analysis::clear()
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

template <class T>
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

void Analysis::addPtLinearHistograms(TString filename, int nbins, std::vector<double> ptlimits)
{
  addPtHistograms<LinearHistogram>(filename, nbins, 0, 0, 0, jet_ptmin, 2000., 0, &ptlimits);
}

void Analysis::addPtSmearedLinearHistograms(TString filename, int nbins, double s, std::vector<double> ptlimits)
{
  addPtHistograms<SmearedLinearHistogram>(filename, nbins, s, 0, 0, jet_ptmin, 2000., 0, &ptlimits);
}

void Analysis::addPtQuadraticHistograms(TString filename, int nbins, double f, std::vector<double> ptlimits)
{
  addPtHistograms<QuadraticHistogram>(filename, nbins, f, 0, 0, jet_ptmin, 2000., 0, &ptlimits);
}

void Analysis::addPtSmearedQuadraticHistograms(TString filename, int nbins, double f, double s, std::vector<double> ptlimits)
{
  addPtHistograms<SmearedQuadraticHistogram>(filename, nbins, f, s, 0, jet_ptmin, 2000., 0, &ptlimits);
}

template <class T>
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

void Analysis::addEtaLinearHistograms(TString filename, int nbins)
{
  addEtaHistograms<LinearHistogram>(filename, nbins, 0, 0, 0, -jet_etamax, jet_etamax);
}

void Analysis::addEtaSmearedLinearHistograms(TString filename, int nbins, double s)
{
  addEtaHistograms<SmearedLinearHistogram>(filename, nbins, s, 0, 0, -jet_etamax, jet_etamax);
}

void Analysis::addEtaQuadraticHistograms(TString filename, int nbins, double f)
{
  addEtaHistograms<QuadraticHistogram>(filename, nbins, f, 0, 0, -jet_etamax, jet_etamax);
}

void Analysis::addEtaSmearedQuadraticHistograms(TString filename, int nbins, double f, double s)
{
  addEtaHistograms<SmearedQuadraticHistogram>(filename, nbins, f, s, 0, -jet_etamax, jet_etamax);
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
    outfile.precision(7);
    outfile.setf(std::ios::scientific);

    output_histograms(*it, outfile);

    outfile.close();
  }
}

void Analysis::output_histograms(const TString& filename, std::ofstream& stream)
{
  jet_exclusive->print(stream, runname, event_count);
  jet_inclusive->print(stream, runname, event_count);

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
// JetAnalysis
// ---------------------------------------------------------------------------

JetAnalysis::JetAnalysis()
  : jet_pt1min(0)
{
}

bool JetAnalysis::check_cuts(SelectorCommon* event)
{
  if (Analysis::check_cuts(event)) {
    return jets[0].pt() >= jet_pt1min;
  } else {
    return false;
  }
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
