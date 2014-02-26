
// Root > T->Process("SelectorCommon.C")
// Root > T->Process("SelectorCommon.C","some options")
// Root > T->Process("SelectorCommon.C+")

#include "SelectorCommon.h"
#include <TH2.h>
#include <TStyle.h>

#ifndef NDEBUG
  #include <iostream>
#endif

#include <LHAPDF.h>

// --------------------------------------------------------------------------- //
// Selector
// --------------------------------------------------------------------------- //

SelectorCommon::SelectorCommon(TTree* /*tree*/)
  : fChain(0), alphasPower(-1), // ROOT
    analysis(0),
    stat_Q2_min(1e100), stat_Q2_max(0.),
    stat_x1_min(1e100), stat_x1_max(0.),
    stat_x2_min(1e100), stat_x2_max(0.)
{
  // no rescaler by default
  rescaler = 0;
  rescale_factor = 1.;
  rescale_n = -1;

  // print event number every 1e6 events
  print_event_step = 1e6;
  // set pdf warning threshold to 1e-9
  pdf_warning_thresh = 1e-9;
  // limit pdf warnings to 10000 (to avoid large log files)
  pdf_warning_limit = 10000;
  pdf_warning_count = 0;

  // eventoscope disabled
  stat_step = 0;

  // quark-filter disabled
  filter_inq = -1;
  filter_nq = -1;

  // advanced settings disabled
  use_sherpa_alphas = false;
  sherpa_alphas = 0;
  beta0fix = 0;
  cdr2fdhfix = -1;
  pi2o12fix = 0;
}


SelectorCommon::~SelectorCommon()
{
  if (analysis) {
    delete analysis;
    analysis = 0;
  }
}

void SelectorCommon::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("id", &id, &b_id);
  fChain->SetBranchAddress("nparticle", &nparticle, &b_nparticle);
  fChain->SetBranchAddress("px", px, &b_px);
  fChain->SetBranchAddress("py", py, &b_py);
  fChain->SetBranchAddress("pz", pz, &b_pz);
  fChain->SetBranchAddress("E", E, &b_E);
  fChain->SetBranchAddress("alphas", &alphas, &b_alphas);
  fChain->SetBranchAddress("kf", kf, &b_kf);
  fChain->SetBranchAddress("weight", &weight, &b_weight);
  fChain->SetBranchAddress("weight2", &weight2, &b_weight2);
  fChain->SetBranchAddress("me_wgt", &me_wgt, &b_me_wtg);
  fChain->SetBranchAddress("me_wgt2", &me_wgt2, &b_me_wtg2);
  fChain->SetBranchAddress("x1", &x1, &b_x1);
  fChain->SetBranchAddress("x2", &x2, &b_x2);
  fChain->SetBranchAddress("x1p", &x1p, &b_x1p);
  fChain->SetBranchAddress("x2p", &x2p, &b_x2p);
  fChain->SetBranchAddress("id1", &id1, &b_id1);
  fChain->SetBranchAddress("id2", &id2, &b_id2);
  fChain->SetBranchAddress("fac_scale", &fac_scale, &b_fac_scale);
  fChain->SetBranchAddress("ren_scale", &ren_scale, &b_ren_scale);
  fChain->SetBranchAddress("nuwgt", &nuwgt, &b_nuwgt);
  fChain->SetBranchAddress("usr_wgts", &usr_wgts, &b_usr_wgts);
  fChain->SetBranchAddress("alphasPower", &alphasPower, &b_alphasPower);
  fChain->SetBranchAddress("part", part, &b_part);
}

Bool_t SelectorCommon::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  if (fChain && fChain->GetCurrentFile()) {
    std::cout << "File: " << fChain->GetCurrentFile()->GetName() << std::endl;
  }

  return kTRUE;
}

void SelectorCommon::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

//   TString option = GetOption();

}

void SelectorCommon::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  fastjet::ClusterSequence::set_fastjet_banner_stream(0); // silence fastjet
  TString option = GetOption();
}

//---------------------------------------------------------------------
// static members
//---------------------------------------------------------------------

#if !defined(__MAKECINT__)
const double SelectorCommon::Nf = 5.;
const double SelectorCommon::CA = 3.;
const double SelectorCommon::CF = 4./3.;

const double SelectorCommon::b0 = (33. - 2.*SelectorCommon::Nf)/(12.*M_PI);
const double SelectorCommon::b1 = (153. - 19.*SelectorCommon::Nf)/(24.*M_PI*M_PI);
#endif

int SelectorCommon::pdg2lha(int pdgnum)
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

double SelectorCommon::adim(Int_t flav)
{
  if (flav == 21) {
    return CA;
  } else if (abs(flav) <= 6) {
    return CF;
  }
  return 0;
}

// Misc functions

double SelectorCommon::pole2(int id1_, int id2_, int n_, const Int_t* kf_)
{
  double sum = 0.;

  for (int i=0; i<n_; i++) {
    sum += adim(kf_[i]);
  }
  sum += adim(id1_);
  sum += adim(id2_);

  return -2.*sum;
}

double SelectorCommon::beta0pole2(int id1_, int id2_, int n_, const Int_t* kf_)
{
  return (2.*M_PI*b0)/pole2(id1_, id2_, n_, kf_);
}

double SelectorCommon::cdr2fdh(int id1_, int id2_, int n_, const Int_t* kf_)
{
  int nQ = 0;
  int nG = 0;

  for (int i=0; i<n_; i++) {
    nG += int(kf_[i] == 21);
    nQ += int(abs(kf_[i]) <= 6);
  }
  nG += int(id1_ == 21);
  nQ += int(abs(id1_) <= 6);
  nG += int(id2_ == 21);
  nQ += int(abs(id2_) <= 6);

  return nG*CA/3. + nQ*CF;
}

// ----------------------------------------------------------------------------
// Various scale functions
// ----------------------------------------------------------------------------
double SelectorCommon::rescaler_multiplicative(const double scale,
                                               const PseudoJetVector& /*partons*/,
                                               const PseudoJetVector& /*jets*/)
{
  double newscale = scale*rescale_factor;
  return newscale;
}

double SelectorCommon::rescaler_ht(const double /*scale*/,
                                   const PseudoJetVector& /*partons*/,
                                   const PseudoJetVector& jets)
{
  double newscale = 0;
  const int imax = rescale_n >= 0 ? rescale_n : jets.size();
  for (int i=0; i<imax; i++) {
    newscale += jets[i].pt();
  }
  newscale *= 0.5*rescale_factor;
  return newscale;
}

double SelectorCommon::rescaler_hthat(const double /*scale*/,
                                      const PseudoJetVector& partons,
                                      const PseudoJetVector& /*jets*/)
{
  double newscale = 0;
  const int imax = rescale_n >= 0 ? rescale_n : partons.size();
  for (int i=0; i<imax; i++) {
    newscale += partons[i].pt();
  }
  newscale *= 0.5*rescale_factor;
  return newscale;
}

double SelectorCommon::rescaler_sumpt2(const double /*scale*/,
                                       const PseudoJetVector& /*partons*/,
                                       const PseudoJetVector& jets)
{
  double newscale = 0;
  const int imax = rescale_n >= 0 ? rescale_n : jets.size();
  for (int i=0; i<imax; i++) {
    newscale += jets[i].pt2();
  }
  newscale = sqrt(newscale);
  newscale *= 0.5*rescale_factor;
  return newscale;
}

double SelectorCommon::rescaler_sumpt2hat(const double /*scale*/,
                                          const PseudoJetVector& partons,
                                          const PseudoJetVector& /*jets*/)
{
  double newscale = 0;
  const int imax = rescale_n >= 0 ? rescale_n : partons.size();
  for (int i=0; i<imax; i++) {
    newscale += partons[i].pt2();
  }
  newscale = sqrt(newscale);
  newscale *= 0.5*rescale_factor;
  return newscale;
}

double SelectorCommon::rescaler_maaht(const double /*scale*/,
                                   const PseudoJetVector& input,
                                   const PseudoJetVector& jets)
{
  double newscale = (input[0]+input[1]).m();
  const int imax = rescale_n >= 0 ? rescale_n : jets.size();
  for (int i=0; i<imax; i++) {
    newscale += jets[i].pt();
  }
  newscale *= 0.5*rescale_factor;
  return newscale;
}

double SelectorCommon::rescaler_maahthat(const double /*scale*/,
                                      const PseudoJetVector& input,
                                      const PseudoJetVector& /*jets*/)
{
  double newscale = (input[0]+input[1]).m();
  const int imax = rescale_n >= 0 ? 2+rescale_n : input.size();
  for (int i=2; i<imax; i++) {
    newscale += input[i].pt();
  }
  newscale *= 0.5*rescale_factor;
  return newscale;
}

double SelectorCommon::rescaler_maa2sumpt2(const double /*scale*/,
                                       const PseudoJetVector& input,
                                       const PseudoJetVector& jets)
{
  double newscale = (input[0]+input[1]).m2();
  const int imax = rescale_n >= 0 ? rescale_n : jets.size();
  for (int i=0; i<imax; i++) {
    newscale += jets[i].pt2();
  }
  newscale = sqrt(newscale);
  newscale *= 0.5*rescale_factor;
  return newscale;
}

double SelectorCommon::rescaler_maa2sumpt2hat(const double /*scale*/,
                                          const PseudoJetVector& input,
                                          const PseudoJetVector& /*jets*/)
{
  double newscale = (input[0]+input[1]).m2();
  const int imax = rescale_n >= 0 ? 2+rescale_n : input.size();
  for (int i=2; i<imax; i++) {
    newscale += input[i].pt2();
  }
  newscale = sqrt(newscale);
  newscale *= 0.5*rescale_factor;
  return newscale;
}

double SelectorCommon::rescaler_minlo(const double /*scale*/,
                                      const PseudoJetVector& input,
                                      const PseudoJetVector& /*jets*/)
{
  PseudoJetVector ktinput;
  PseudoJetVector primary;

  // inputs 0,1 are dummy beams
  {
    fastjet::PseudoJet beamA = fastjet::PseudoJet(0., 0., 0., 0.);
    FlavourKTPlugin::addFlavour(beamA, abs(id1) > 6 ? id1 : -id1);
    ktinput.push_back(beamA);
    fastjet::PseudoJet beamB = fastjet::PseudoJet(0., 0., 0., 0.);
    FlavourKTPlugin::addFlavour(beamB, abs(id2) > 6 ? id2 : -id2);
    ktinput.push_back(beamB);
  }

  for (unsigned i=0; i<input.size(); i++) {
    const int flav = kf[i];
    fastjet::PseudoJet parton = input[i];
    FlavourKTPlugin::addFlavour(parton, flav);
    if (/*onlyqcd and */abs(flav) > 6 and flav != 21) { // non-qcd stuff goes into primary system
      primary.push_back(parton);
      continue;
    }
    ktinput.push_back(parton);
  }

  fastjet::ClusterSequence cs(ktinput, clustering_def);
  minlo_scales.clear();

  // determine clust_first index of the first clustering
  int imax = ktinput.size()-1;
  int clust_first = 0;
  if (part[0] == 'R') {
    static long unsigned maxinputlen = 0;
    maxinputlen = std::max(maxinputlen, ktinput.size());
    if (ktinput.size() == maxinputlen) {
      imax -= 1;
      clust_first += 1;
    }
  }

  // determine clust_last index of the first clustering (the rest is primary system)
  int clust_last = clust_first;
  double rdij = 0.;
  for (int i = imax; i >= 0; i--) {
    const double nextrdij = sqrt(cs.exclusive_dmerge(i));
    if (nextrdij > 0.) {
      if (nextrdij < rdij) {
        break;
      }
      clust_last += 1;  // number of clusterings, can be limited to keep more stuff in primary system
      rdij = nextrdij;
      minlo_scales.push_back(rdij);
    }
  }

  const std::vector<fastjet::ClusterSequence::history_element>& cshist = cs.history();
  const std::vector<fastjet::PseudoJet>& csjets = cs.jets();

  // check that fastjet didn't do anything unexpected
  for (int i = 0; cshist[i].jetp_index >= 0; i++) {
    assert(i == cshist[i].jetp_index);
  }

  // add un-clustered non-beam partons to the primary system
  clust_first += ktinput.size();  // everything before clust_first is input
  clust_last += ktinput.size();   // everything from clust_last and up is primary system
  for (int i = clust_last; cshist[i].jetp_index >= 0; i++) {
    assert(cshist[i].parent1 >= 0);
    const fastjet::PseudoJet& parent1 = csjets[cshist[i].parent1];
    if (parent1.E() > 0.) {
      primary.push_back(parent1);
    }
    assert(cshist[i].parent2 >= 2);
    const fastjet::PseudoJet& parent2 = csjets[cshist[i].parent2];
    if (parent2.E() > 0.) {
      primary.push_back(parent2);
    }
  }

  double minlo_Q0 = 2.*lambda;  // NP cutoff = 2*lambda
  minlo_Q0 = std::max(minlo_Q0, minlo_scales.front());
  minlo_scales.front() = minlo_Q0;

  double minlo_Q = 0.;
  for (unsigned i = 0; i < primary.size(); i++) {
    minlo_Q += primary[i].pt();
  }
  minlo_Q *= 0.5;  // primary system scale = sum pt/2
  minlo_Q = std::max(minlo_Q, minlo_scales.back());  // round up to q_n if needed

  // compute alpha factor
  alphafactor = 1.;
  double minlo_mur = 1.;
  double minlo_alpha = 0.;
  // first n powers corresponding to n clusterings
  for (unsigned i = 0; i < minlo_scales.size(); i++) {
    const double minlo_qi = minlo_scales[i];
    minlo_mur *= minlo_qi;
    const double as = getAlphaS(rescale_factor*minlo_qi);
    minlo_alpha += as;
    alphafactor *= as/alphas;
  }
  // then m powers corresponding to primary system
  for (int i = minlo_scales.size(); i < born_alphapower; i++) {
    minlo_mur *= minlo_Q;
    const double as = getAlphaS(rescale_factor*minlo_Q);
    minlo_alpha += as;
    alphafactor *= as/alphas;
  }
  minlo_alpha /= born_alphapower; // a_s(n+m+1) as in (3.2)
  if (event_alphapower > born_alphapower) {  // if NLO multiply a_s(n+m+1)
    alphafactor *= minlo_alpha/alphas;
  }

  // set mF and mR scale factors
  fac_scalefactor = rescale_factor*minlo_Q0/fac_scale;
  ren_scalefactor = rescale_factor*pow(minlo_mur, 1./born_alphapower)/ren_scale;

  // compute sudakov FFs for external legs
  double sudakov1 = 1.;
  double sudakov1sub = 0.;
  for (int i = 0; i < clust_first; i++) {
    const int child_idx = cshist[i].child;
    const int flav = FlavourKTPlugin::getFlavour(csjets[cshist[i].jetp_index]);
    double scale1 = 0;
    if (child_idx < clust_first) {  // NLO radiation
      // ignore
    } else if (child_idx < clust_last) {  // clustering scale
      scale1 = sqrt(cshist[child_idx].dij);
    } else {  // hard scale1
      scale1 = minlo_Q;
    }
    if (scale1 != 0. and scale1 != minlo_Q0) {
      sudakov1 *= Deltaf(minlo_Q0, scale1, flav);
      if (part[0] == 'B') {
        sudakov1sub += Deltaf1(minlo_Q0, scale1, flav);
      }
    }
  }

  // compute sudakov FFs for internal sceleton lines
  double sudakov2 = 1.;
  double sudakov2sub = 0.;
  for (int i = clust_first; i < clust_last; i++) {
    const int child_idx = cshist[i].child;
    const int flav = FlavourKTPlugin::getFlavour(csjets[cshist[i].jetp_index]);
    double scale1 = sqrt(cshist[i].dij);
    double scale2 = 0;
    if (child_idx < clust_last) {  // clustering scale
      scale2 = sqrt(cshist[child_idx].dij);
    } else {  // hard scale
      scale2 = minlo_Q;
    }
    assert(scale2 > scale1 or (scale2 == scale1 and scale1 == minlo_Q));
    if (scale1 == minlo_Q0) {
      sudakov2 *= Deltaf(minlo_Q0, scale2, flav);
      if (part[0] == 'B') {
        sudakov2sub += Deltaf1(minlo_Q0, scale2, flav);
      }
    } else {
      sudakov2 *= Deltaf(minlo_Q0, scale2, flav)/Deltaf(minlo_Q0, scale1, flav);
      if (part[0] == 'B') {
        sudakov2sub += Deltaf1(minlo_Q0, scale2, flav) - Deltaf1(minlo_Q0, scale1, flav);
      }
    }
  }

  // adjust alpha factor with wudakov FFs (Note that a_s(m+n+1) is included here, not in Deltaf1)
  alphafactor *= sudakov1*sudakov2*(1. - minlo_alpha*sudakov2sub - minlo_alpha*sudakov1sub);

  if (alphafactor != alphafactor) {
    std::cout << id1 << " " << id2 << " -> "; for (unsigned i=0; i<input.size(); i++) std::cout << kf[i] << " "; std::cout << "\n";
    std::cout << "MM ";for (unsigned i=0; i<minlo_scales.size(); i++) std::cout << minlo_scales[i] << " "; std::cout << "\n";
    for (unsigned i = 0; i < cshist.size(); i++) {
      std::cout << cshist[i].parent1 << " " << cshist[i].parent2 << " " << cshist[i].child << " " << cshist[i].jetp_index << " " << cshist[i].dij << "\n";
    }
    std::cout << "S " << sudakov1 << " " << sudakov2 << " " << sudakov2sub << " " << sudakov1sub << std::endl;
    std::cout << "alphafactor = " << alphafactor << " fac = " << fac_scalefactor << " ren = " << ren_scalefactor << "\n";
  }
  return 0.; // zero means that we use fac_scalefactor/ren_scalefactor
}


double SelectorCommon::Deltaf(double Q0sq, double Qsq, int flav)
{
  const double C = adim(flav);
  const double B = flav == 21 ? M_PI*b0/CA : 3./4.;

  const double Lsq = lambda*lambda;

  return exp(-C/(M_PI*b0)*
             (log(log(Qsq/Lsq)/log(Q0sq/Lsq))*(0.5*log(Qsq/Lsq)-B)-0.5*log(Qsq/Q0sq))
            );
}

double SelectorCommon::Deltaf1(double Q0sq, double Qsq, int flav)
{
  const double C = adim(flav);
  const double B = flav == 21 ? M_PI*b0/CA : 3./4.;

  const double logQoQ0 = log(Qsq/Q0sq);
  return -C/M_PI*(0.25*logQoQ0*logQoQ0 - logQoQ0*B);
}

// ----------------------------------------------------------------------------
// Reweighting: Scale, PDF and AlphaS change
// ----------------------------------------------------------------------------

double SelectorCommon::LambdaQCD(double muR, double aS) const
{
  if (aS < 0.) {
    aS = LHAPDF::alphasPDF(TOPDF, muR);
  }
  double Lam = 0.2;
  double diff = 1;
  while (diff > 1e-7) {
    double t = log((muR*muR)/(Lam*Lam));
    double f = aS*b0*b0*t + (b1*log(t))/(b0*t) - b0;
    double fp = -(2.*(b1 + aS*b0*b0*b0*t*t - b1*log(t)))/(b0*Lam*t*t);
    double newLam = Lam - f/fp;
    diff = abs(2.*(newLam-Lam)/(newLam+Lam));
    Lam = newLam;
  }
  return Lam;
}

double SelectorCommon::getAlphaS(double mur)
{
  if (use_sherpa_alphas) {
    return sherpa_alphas->AlphaS(mur);
  } else {
    return LHAPDF::alphasPDF(TOPDF, mur);
  }
}

void SelectorCommon::reweight(const PseudoJetVector& input,
                              const PseudoJetVector& jets)
{
  if (FROMPDF == 0 and TOPDF == 0 and rescaler == 0) {
    return;
  }

  lhaid1 = pdg2lha(id1);
  lhaid2 = pdg2lha(id2);

  const double fx1 = LHAPDF::xfx(FROMPDF, x1, fac_scale, lhaid1)/x1;
  const double fx2 = LHAPDF::xfx(FROMPDF, x2, fac_scale, lhaid2)/x2;

  // this simple check does not work for I-part (hence != 18)
  if (pdf_warning_count < pdf_warning_limit and nuwgt != 18) {
    const double pdf_calc_weight = me_wgt*(fx1*fx2);
    const double pdf_rel_diff = (pdf_calc_weight - weight)/weight;
    if (abs(pdf_rel_diff) > pdf_warning_thresh) {
      pdf_warning_count += 1;
      std::cout << "Check your FROMPDF! " << pdf_warning_count << " (" << id << ") "
                << pdf_calc_weight << " != " << weight << " (" << pdf_rel_diff << ")\n";
      if (pdf_warning_count == pdf_warning_limit) {
        std::cout << "Check your FROMPDF! Reached warning limit " << pdf_warning_count
                  << " no further warnings will be printed\n";
        std::cout << "Check your FROMPDF! Fraction of suspicious events "
                  << pdf_warning_count/analysis->event_count << "\n";
      }
    }
  }

  if (alphasPower >= 0) {
    event_alphapower = int(alphasPower);
  } else {
    event_alphapower = born_alphapower;
  }
  if (event_order() < 0 or event_order() > 1) {
    std::cout << "Check your alpha power " << int(alphasPower) << " != " << born_alphapower << std::endl;
  }

  double scalefactor = 1.;
  if (rescaler) {
    scalefactor = (*this.*rescaler)(fac_scale, input, jets)/fac_scale;
  }
  if (scalefactor != 0.) {  // zero means that fac_scalefactor/ren_scalefactor are set in rescaler
    fac_scalefactor = scalefactor;
    ren_scalefactor = scalefactor;
  }

  // below new scales
  fac_scale *= fac_scalefactor;
  ren_scale *= ren_scalefactor;

  const double log_r = log(ren_scalefactor*ren_scalefactor);  // log(murnew^2/murold^2)
  const double log_f = log(fac_scalefactor*fac_scalefactor);  // log(mufnew^2/mufold^2)

  LHAPDF::xfx(TOPDF, x1, fac_scale, pdfx1);
  LHAPDF::xfx(TOPDF, x2, fac_scale, pdfx2);
  const double new_fx1 = pdfx1[lhaid1+6]/x1;
  const double new_fx2 = pdfx2[lhaid2+6]/x2;

  if (scalefactor != 0.) {  // zero means that alphafactor is set in rescaler
    alphafactor = pow(getAlphaS(ren_scale)/alphas, event_alphapower);
  }

  coll_weights_count = 0;
  if (nuwgt == 0) {
    weight = me_wgt*(new_fx1*new_fx2);
  } else if (nuwgt == 2) {

    // pi^2 scheme conversion
    if (pi2o12fix) {
      me_wgt += -M_PI*M_PI/12.*usr_wgts[1];
    }
    // CDR to DRED conversion
    if (cdr2fdhfix >= 0) {
      me_wgt += -(cdr2fdhfix - cdr2fdh(id1, id2, nparticle, kf))*usr_wgts[1]/pole2(id1, id2, nparticle, kf);
    }
    // fix for wrong alpha power in old Sherpa ntuples
    if (beta0fix > 0) {
      usr_wgts[0] -= (beta0fix - (event_alphapower-1))*2.*usr_wgts[1]*beta0pole2(id1, id2, nparticle, kf);
    }

    if (beta0fix >= 0) {
      me_wgt += usr_wgts[0]*log_r + 0.5*usr_wgts[1]*log_r*log_r;
    }

    weight = me_wgt*(new_fx1*new_fx2);
  } else if (nuwgt == 18) {
    coll_weights_count = 8;
    for (int i=0; i<8; i++) {
      coll_weights[i] = usr_wgts[2+i] + usr_wgts[10+i]*log_f;
    }

    coll_weights[1] /= x1p;
    coll_weights[3] /= x1p;
    coll_weights[5] /= x2p;
    coll_weights[7] /= x2p;

    x1r = x1/x1p;
    x2r = x2/x2p;

    LHAPDF::xfx(TOPDF, x1r, fac_scale, pdfx1p);
    LHAPDF::xfx(TOPDF, x2r, fac_scale, pdfx2p);

    // no top pdf below
    pdf_f1[0] = (lhaid1 != 0) ? pdfx1[6+lhaid1]/x1 :
                              ( pdfx1[7] + pdfx1[8] + pdfx1[9] + pdfx1[10] + pdfx1[11]
                              + pdfx1[1] + pdfx1[2] + pdfx1[3] + pdfx1[4]  + pdfx1[5])/x1;
    pdf_f2[0] = (lhaid2 != 0) ? pdfx2[6+lhaid2]/x2 :
                              ( pdfx2[7] + pdfx2[8] + pdfx2[9] + pdfx2[10] + pdfx2[11]
                              + pdfx2[1] + pdfx2[2] + pdfx2[3] + pdfx2[4]  + pdfx2[5])/x2;

    pdf_f1[1] = (lhaid1 != 0) ? pdfx1p[6+lhaid1]/x1r :
                              ( pdfx1p[7] + pdfx1p[8] + pdfx1p[9] + pdfx1p[10] + pdfx1p[11]
                              + pdfx1p[1] + pdfx1p[2] + pdfx1p[3] + pdfx1p[4]  + pdfx1p[5])/x1r;
    pdf_f2[1] = (lhaid2 != 0) ? pdfx2p[6+lhaid2]/x2r :
                              ( pdfx2p[7] + pdfx2p[8] + pdfx2p[9] + pdfx2p[10] + pdfx2p[11]
                              + pdfx2p[1] + pdfx2p[2] + pdfx2p[3] + pdfx2p[4]  + pdfx2p[5])/x2r;

    pdf_f1[2] = pdfx1[6]/x1;
    pdf_f2[2] = pdfx2[6]/x2;

    pdf_f1[3] = pdfx1p[6]/x1r;
    pdf_f2[3] = pdfx2p[6]/x2r;

    // pi^2 scheme conversion
    if (pi2o12fix) {
      me_wgt += -M_PI*M_PI/12.*usr_wgts[1];
    }
    // CDR to DRED conversion
    if (cdr2fdhfix >= 0) {
      me_wgt += -(cdr2fdhfix - cdr2fdh(id1, id2, nparticle, kf))*usr_wgts[1]/pole2(id1, id2, nparticle, kf);
    }
    // fix for wrong alpha power in old Sherpa ntuples
    if (beta0fix) {
      usr_wgts[0] -= (beta0fix - (event_alphapower-1))*2.*usr_wgts[1]*beta0pole2(id1, id2, nparticle, kf);
    }

    if (beta0fix >= 0) {
      me_wgt += usr_wgts[0]*log_r + 0.5*usr_wgts[1]*log_r*log_r;
    }
    weight = me_wgt*(new_fx1*new_fx2)
           + (coll_weights[0]*pdf_f1[0] + coll_weights[1]*pdf_f1[1] +
              coll_weights[2]*pdf_f1[2] + coll_weights[3]*pdf_f1[3])*new_fx2
           + (coll_weights[4]*pdf_f2[0] + coll_weights[5]*pdf_f2[1] +
              coll_weights[6]*pdf_f2[2] + coll_weights[7]*pdf_f2[3])*new_fx1;
  } else {
    std::cout << "Unknown value for nuwgt = " << nuwgt << std::endl;
  }
  const double appl_alphas = pow(alphas/(2.*M_PI), event_alphapower);
  for (int i=0; i<coll_weights_count; i++) {
    coll_weights[i] /= appl_alphas;
  }
  naked_weight = me_wgt/appl_alphas;
  weight *= alphafactor;

  statUpdate();
}

void SelectorCommon::statUpdate()
{
  if (ren_scale < stat_Q2_min) {
    stat_Q2_min = ren_scale;
  } else if (ren_scale > stat_Q2_max) {
    stat_Q2_max = ren_scale;
  }

  if (x1 < stat_x1_min) {
    stat_x1_min = x1;
  } else if (x1 > stat_x1_max) {
    stat_x1_max = x1;
  }

  if (x2 < stat_x2_min) {
    stat_x2_min = x2;
  } else if (x2 > stat_x2_max) {
    stat_x2_max = x2;
  }
}

void SelectorCommon::statReport()
{
  std::cout << "STAT Q2 in [" << stat_Q2_min << ", " << stat_Q2_max << "]\n";
  std::cout << "STAT x1 in [" << stat_x1_min << ", " << stat_x1_max << "]\n";
  std::cout << "STAT x2 in [" << stat_x2_min << ", " << stat_x2_max << "]\n";
}

void SelectorCommon::initAS(const int order, const double asMZ, const double mZ2,
                        const std::vector<double>& qmasses)
{
  sherpa_alphas = new SHERPA::One_Running_AlphaS(order, asMZ, mZ2, qmasses);
}

Bool_t SelectorCommon::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either SelectorCommon::GetEntry() or TBranch::GetEntry()
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
  analysis->event_count += 1;

  // runtime statistics
  if (long(analysis->event_count) % print_event_step == 0) {
    std::cout << "Processing event " << analysis->event_count << std::endl;
  }

  // quark-filter
  if (filter_inq >= 0 and filter_nq >= 0) {
    int inq = int(abs(id1) <= 6) + int(abs(id2) <= 6);
    int nq = inq;
    for (int i=0; i<nparticle; i++) {
      nq += int(abs(kf[i]) <= 6);
    }
    if (inq != filter_inq or nq != filter_nq) {
      return kTRUE;
    }
  }

  if (analysis->check_cuts(this)) {
    reweight(analysis->input, analysis->jets);  // reweight event in-place
    analysis->analysis_bin(this);

    // evento-scope
    if (stat_step) {
      xsval_cur += weight;
      xserr_cur += weight*weight;
      if (int(analysis->event_count) % stat_step == 0) {
        xsvals.push_back(xsval_cur/analysis->event_count);
        xserrs.push_back(sqrt(xserr_cur)/analysis->event_count);
      }
    }
  }

  return kTRUE;
}

void SelectorCommon::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  analysis->analysis_finalize();
}

void SelectorCommon::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}

#include "SelectorHistograms.C"
#include "SelectorAnalysis.C"
