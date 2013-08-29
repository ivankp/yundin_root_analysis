
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

static inline double adim(Int_t flav)
{
  static const double CA = 3.;
  static const double CF = 4./3.;
  if (flav == 21) {
    return CA;
  } else if (abs(flav) <= 6) {
    return CF;
  }
  return 0;
}

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
  static const double CA = 3.;
//   static const double CF = 4./3.;
  static const double Nf = 5.;

  static const double b0 = 0.5*(11./3.*CA - 2./3.*Nf);

  return b0/pole2(id1_, id2_, n_, kf_);
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
  const unsigned imax = rescale_n ? rescale_n : jets.size();
  for (unsigned i=0; i<imax; i++) {
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
  const unsigned imax = rescale_n ? rescale_n : partons.size();
  for (unsigned i=0; i<imax; i++) {
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
  const unsigned imax = rescale_n ? rescale_n : jets.size();
  for (unsigned i=0; i<imax; i++) {
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
  const unsigned imax = rescale_n ? rescale_n : partons.size();
  for (unsigned i=0; i<imax; i++) {
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
  const unsigned imax = rescale_n ? rescale_n : jets.size();
  for (unsigned i=0; i<imax; i++) {
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
  const unsigned imax = rescale_n ? rescale_n : input.size();
  for (unsigned i=0; i<imax; i++) {
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
  const unsigned imax = rescale_n ? rescale_n : jets.size();
  for (unsigned i=0; i<imax; i++) {
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
  const unsigned imax = rescale_n ? rescale_n : input.size();
  for (unsigned i=0; i<imax; i++) {
    newscale += input[i].pt2();
  }
  newscale = sqrt(newscale);
  newscale *= 0.5*rescale_factor;
  return newscale;
}

// ----------------------------------------------------------------------------
// Reweighting: Scale, PDF and AlphaS change
// ----------------------------------------------------------------------------

void SelectorCommon::reweight(const PseudoJetVector& input,
                              const PseudoJetVector& jets)
{
  if (FROMPDF == 0 and TOPDF == 0 and rescaler == 0) {
    return;
  }
  double scalefactor = 1.;
  if (rescaler) {
    scalefactor = (*this.*rescaler)(fac_scale, input, jets)/fac_scale;
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
      usr_wgts[0] -= (beta0fix - (alphapower-1))*2.*usr_wgts[1]*beta0pole2(id1, id2, nparticle, kf);
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
      usr_wgts[0] -= (beta0fix - (alphapower-1))*2.*usr_wgts[1]*beta0pole2(id1, id2, nparticle, kf);
    }

    weight = (me_wgt  + usr_wgts[0]*lr + 0.5*usr_wgts[1]*lr*lr)*(new_fx1*new_fx2)
           + (w[0]*f1[0]/x1 + w[1]*f1[1]/x1 + w[2]*f1[2]/x1 + w[3]*f1[3]/x1)*new_fx2
           + (w[4]*f2[0]/x2 + w[5]*f2[1]/x2 + w[6]*f2[2]/x2 + w[7]*f2[3]/x2)*new_fx1;
  } else {
    std::cout << "Unknown value for nuwgt = " << nuwgt << std::endl;
  }
  weight *= pow(alphafactor, alphapower);
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
