
#ifndef SELECTOR_COMMON_H
#define SELECTOR_COMMON_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Stuff

#include <TRegexp.h>
#include <fastjet/ClusterSequence.hh>
#include <fstream>

// Math
#include <cstdlib>
#include <cmath>
using std::abs;
using std::pow;
using std::sqrt;

#include "SherpaAlphaS.h"
#include "SelectorHistograms.h"
#include "SelectorAnalysis.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

#define MAXNPARTICLE 16
#define MAXNUWEIGHT 128

class SelectorCommon : public TSelector
{
  public :
    // -----------------------------------------------------------------------
    // ROOT stuff BEGIN         ROOT stuff BEGIN        ROOT stuff BEGIN
    // -----------------------------------------------------------------------
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain

    // Declaration of leaf types
    Int_t           id;
    Int_t           nparticle;
    Float_t         px[MAXNPARTICLE];   //[nparticle]
    Float_t         py[MAXNPARTICLE];   //[nparticle]
    Float_t         pz[MAXNPARTICLE];   //[nparticle]
    Float_t         E[MAXNPARTICLE];   //[nparticle]
    Double_t        alphas;
    Int_t           kf[MAXNPARTICLE];   //[nparticle]
    Double_t        weight;
    Double_t        weight2;
    Double_t        me_wgt;
    Double_t        me_wgt2;
    Double_t        x1;
    Double_t        x2;
    Double_t        x1p;
    Double_t        x2p;
    Int_t           id1;
    Int_t           id2;
    Double_t        fac_scale;
    Double_t        ren_scale;
    Int_t           nuwgt;
    Double_t        usr_wgts[MAXNUWEIGHT];   //[nuwgt]
    Char_t          alphasPower;
    Char_t          part[2];

    // List of branches
    TBranch        *b_id;   //!
    TBranch        *b_nparticle;   //!
    TBranch        *b_px;   //!
    TBranch        *b_py;   //!
    TBranch        *b_pz;   //!
    TBranch        *b_E;   //!
    TBranch        *b_alphas;   //!
    TBranch        *b_kf;   //!
    TBranch        *b_weight;   //!
    TBranch        *b_weight2;   //!
    TBranch        *b_me_wtg;   //!
    TBranch        *b_me_wtg2;   //!
    TBranch        *b_x1;   //!
    TBranch        *b_x2;   //!
    TBranch        *b_x1p;   //!
    TBranch        *b_x2p;   //!
    TBranch        *b_id1;   //!
    TBranch        *b_id2;   //!
    TBranch        *b_fac_scale;   //!
    TBranch        *b_ren_scale;   //!
    TBranch        *b_nuwgt;   //!
    TBranch        *b_usr_wgts;   //!
    TBranch        *b_alphasPower;   //!
    TBranch        *b_part;   //!

    virtual Int_t   Version() const {
      return 2;
    }
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) {
      return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
    }
    virtual void    SetOption(const char *option) {
      fOption = option;
    }
    virtual void    SetObject(TObject *obj) {
      fObject = obj;
    }
    virtual void    SetInputList(TList *input) {
      fInput = input;
    }
    virtual TList  *GetOutputList() const {
      return fOutput;
    }
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    ClassDef(SelectorCommon, 0);

    // -----------------------------------------------------------------------
    // ROOT stuff END           ROOT stuff END          ROOT stuff END
    // -----------------------------------------------------------------------

    SelectorCommon(TTree * /*tree*/ =0)
    : fChain(0), // ROOT
      analysis(0),
      rescale_factor(1.), rescale_n(1), rescaler(0),
      filter_inq(-1), filter_nq(-1),
      use_sherpa_alphas(false), sherpa_alphas(0), beta0fix(0),
      stat_step()
    { }

    ~SelectorCommon();

    // analysis
    typedef std::vector<fastjet::PseudoJet> PseudoJetVector;
    Analysis* analysis;

    // scale change
    typedef double (SelectorCommon::*RescalerType)(const double scale,
                                                   const PseudoJetVector& partons,
                                                   const PseudoJetVector& jets);
    double rescaler_multiplicative(const double scale,
                                   const PseudoJetVector& partons,
                                   const PseudoJetVector& jets);
    double rescaler_ht(const double scale,
                       const PseudoJetVector& partons,
                       const PseudoJetVector& jets);
    double rescaler_hthat(const double scale,
                          const PseudoJetVector& partons,
                          const PseudoJetVector& jets);
    double rescaler_htn(const double scale,
                        const PseudoJetVector& partons,
                        const PseudoJetVector& jets);
    double rescaler_htnhat(const double scale,
                           const PseudoJetVector& partons,
                           const PseudoJetVector& jets);
    double rescaler_sumpt2(const double scale,
                           const PseudoJetVector& partons,
                           const PseudoJetVector& jets);
    double rescaler_sumpt2hat(const double scale,
                              const PseudoJetVector& partons,
                              const PseudoJetVector& jets);
    void setrescaler_none() { rescaler = 0; }
    void setrescaler_multiplicative() { rescaler = &SelectorCommon::rescaler_multiplicative; }
    void setrescaler_ht() { rescaler = &SelectorCommon::rescaler_ht; }
    void setrescaler_hthat() { rescaler = &SelectorCommon::rescaler_hthat; }
    void setrescaler_htn() { rescaler = &SelectorCommon::rescaler_htn; }
    void setrescaler_htnhat() { rescaler = &SelectorCommon::rescaler_htnhat; }
    void setrescaler_sumpt2() { rescaler = &SelectorCommon::rescaler_sumpt2; }
    void setrescaler_sumpt2hat() { rescaler = &SelectorCommon::rescaler_sumpt2hat; }

    double rescale_factor;
    unsigned rescale_n;
    RescalerType rescaler;

    int alphapower;

    // qfilter
    int filter_inq;
    int filter_nq;

    // pdf sets
    int FROMPDF, TOPDF;

    // alphas
    bool use_sherpa_alphas;
    void initAS(const int order, const double asMZ, const double mZ2,
                const std::vector<double>& qmasses);
    SHERPA::One_Running_AlphaS* sherpa_alphas;

    // beta0 workaround
    int beta0fix;
    static double beta0pole2(int id1_, int id2_, int n_, const Int_t* kf_);
    static double pole2(int id1_, int id2_, int n_, const Int_t* kf_);

    // eventoscope
    int stat_step;
    double xsval_cur, xserr_cur;
    std::vector<double> xsvals;
    std::vector<double> xserrs;

  protected:
    void reweight(const PseudoJetVector& input,
                  const PseudoJetVector& jets);
};

#if defined(__MAKECINT__)
#pragma link C++ class fastjet::JetDefinition;
#endif

#endif
