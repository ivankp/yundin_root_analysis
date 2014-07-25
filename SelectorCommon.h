
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
#include "FlavourKT.h"
#include <LoopSim.hh>

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
    Int_t           ntuple_id;
    Int_t           ntuple_nparticle;
    Float_t         ntuple_px[MAXNPARTICLE];   //[nparticle]
    Float_t         ntuple_py[MAXNPARTICLE];   //[nparticle]
    Float_t         ntuple_pz[MAXNPARTICLE];   //[nparticle]
    Float_t         ntuple_E[MAXNPARTICLE];   //[nparticle]
    Double_t        ntuple_alphas;
    Int_t           ntuple_kf[MAXNPARTICLE];   //[nparticle]
    Double_t        ntuple_weight;
    Double_t        ntuple_weight2;
    Double_t        ntuple_me_wgt;
    Double_t        ntuple_me_wgt2;
    Double_t        ntuple_x1;
    Double_t        ntuple_x2;
    Double_t        ntuple_x1p;
    Double_t        ntuple_x2p;
    Int_t           ntuple_id1;
    Int_t           ntuple_id2;
    Double_t        ntuple_fac_scale;
    Double_t        ntuple_ren_scale;
    Int_t           ntuple_nuwgt;
    Double_t        ntuple_usr_wgts[MAXNUWEIGHT];   //[nuwgt]
    Char_t          ntuple_alphaspower;
    Char_t          ntuple_part[2];

    // Accessor methods
    Int_t           get_event_id() const { return ntuple_id; }
    Int_t           get_nparticle() const { return ntuple_nparticle; }
    const Float_t*  get_px() const { return ntuple_px; }
    const Float_t*  get_py() const { return ntuple_py; }
    const Float_t*  get_pz() const { return ntuple_pz; }
    const Float_t*  get_E() const { return ntuple_E; }
    Float_t         get_px(int i) const { return ntuple_px[i]; }
    Float_t         get_py(int i) const { return ntuple_py[i]; }
    Float_t         get_pz(int i) const { return ntuple_pz[i]; }
    Float_t         get_E(int i) const { return ntuple_E[i]; }
    Double_t        orig_alphas() const { return ntuple_alphas; }
    const Int_t*    get_kf() const { return ntuple_kf; }
    Int_t           get_kf(int i) const { return ntuple_kf[i]; }
    Double_t        orig_weight() const { return ntuple_weight; }
    Double_t        orig_weight2() const { return ntuple_weight2; }
    Double_t        orig_me_wgt() const { return ntuple_me_wgt; }
    Double_t        orig_me_wgt2() const { return ntuple_me_wgt2; }
    Double_t        get_x1() const { return ntuple_x1; }
    Double_t        get_x2() const { return ntuple_x2; }
    Double_t        get_x1p() const { return ntuple_x1p; }
    Double_t        get_x2p() const { return ntuple_x2p; }
    Int_t           get_id1() const { return ntuple_id1; }
    Int_t           get_id2() const { return ntuple_id2; }
    Double_t        orig_fac_scale() const { return ntuple_fac_scale; }
    Double_t        orig_ren_scale() const { return ntuple_ren_scale; }
    Int_t           get_nuwgt() const { return ntuple_nuwgt; }
    const Double_t* orig_usr_wgts() const { return ntuple_usr_wgts; }
    Double_t        orig_usr_wgts(int i) const { return ntuple_usr_wgts[i]; }
    Char_t          get_alphaspower() const { return ntuple_alphaspower; }
    const Char_t*   get_part() const { return ntuple_part; }
    Char_t          get_part(int i) const { return ntuple_part[i]; }

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
    TBranch        *b_alphaspower;   //!
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

    SelectorCommon(TTree* tree=0);
    ~SelectorCommon();

    // analysis
    Analysis* analysis;
    typedef std::vector<fastjet::PseudoJet> PseudoJetVector;

    fastjet::PseudoJet get_vec(int i) const;
    PseudoJetVector get_fjinput() const;
    std::vector<LSParticle> get_lsinput() const;
    static PseudoJetVector lsinput2fjinput(const std::vector<LSParticle>& in);


    enum {MODE_PLAIN=0, MODE_LOOPSIM};
    int analysis_mode;
    void process_single_event(bool do_reweight=true);

    // scale change
    typedef double (SelectorCommon::*RescalerType)(const double scale,
                                                   const PseudoJetVector& partons,
                                                   const PseudoJetVector& jets);
    double rescaler_multiplicative(const double scale,
                                   const PseudoJetVector& partons,
                                   const PseudoJetVector& jets);

    // pure jets scales
    double rescaler_ht(const double scale,
                       const PseudoJetVector& partons,
                       const PseudoJetVector& jets);
    double rescaler_hthat(const double scale,
                          const PseudoJetVector& partons,
                          const PseudoJetVector& jets);
    double rescaler_sumpt2(const double scale,
                           const PseudoJetVector& partons,
                           const PseudoJetVector& jets);
    double rescaler_sumpt2hat(const double scale,
                              const PseudoJetVector& partons,
                              const PseudoJetVector& jets);
    // two-nonqcd + jets scales
    double rescaler_maaht(const double scale,
                       const PseudoJetVector& input,
                       const PseudoJetVector& jets);
    double rescaler_maahthat(const double scale,
                          const PseudoJetVector& input,
                          const PseudoJetVector& jets);
    double rescaler_maa2sumpt2(const double scale,
                           const PseudoJetVector& input,
                           const PseudoJetVector& jets);
    double rescaler_maa2sumpt2hat(const double scale,
                              const PseudoJetVector& input,
                              const PseudoJetVector& jets);
    double rescaler_minlo(const double scale,
                          const PseudoJetVector& input,
                          const PseudoJetVector& jets);


    void setrescaler_none() { opt_rescaler = 0; }
    void setrescaler_multiplicative() { opt_rescaler = &SelectorCommon::rescaler_multiplicative; }
    void setrescaler_ht() { opt_rescaler = &SelectorCommon::rescaler_ht; }
    void setrescaler_hthat() { opt_rescaler = &SelectorCommon::rescaler_hthat; }
    void setrescaler_sumpt2() { opt_rescaler = &SelectorCommon::rescaler_sumpt2; }
    void setrescaler_sumpt2hat() { opt_rescaler = &SelectorCommon::rescaler_sumpt2hat; }
    void setrescaler_maaht() { opt_rescaler = &SelectorCommon::rescaler_maaht; }
    void setrescaler_maahthat() { opt_rescaler = &SelectorCommon::rescaler_maahthat; }
    void setrescaler_maa2sumpt2() { opt_rescaler = &SelectorCommon::rescaler_maa2sumpt2; }
    void setrescaler_maa2sumpt2hat() { opt_rescaler = &SelectorCommon::rescaler_maa2sumpt2hat; }
    void setrescaler_minlo() {
      clustering_def = fastjet::JetDefinition(new FlavourKTPlugin());
      clustering_def.delete_plugin_when_unused();
      lambda = LambdaQCD();
      std::cout << "Set LambdaQCD = " << lambda << std::endl;
      opt_rescaler = &SelectorCommon::rescaler_minlo;
    }

    // static functions
    static const double Nf;
    static const double CA;
    static const double CF;
    static const double b0;
    static const double b1;

    static int pdg2lha(int pdgnum);
    static double adim(Int_t flav);
    double Deltaf(double Q0sq, double Qsq, int flav);
    double Deltaf1(double Q0sq, double Qsq, int flav);
    double LambdaQCD(double muR=91.188, double aS=-1.) const;

    // member variables
    double opt_rescale_factor;
    int opt_rescale_n;
    RescalerType opt_rescaler;

    double fac_scalefactor;
    double ren_scalefactor;
    double alphafactor;
    // minlo-specific
    fastjet::JetDefinition clustering_def;
    double lambda;
    std::vector<double> minlo_scales;

    int opt_born_alphaspower;  // used only for APPLgrid init
    int event_order;
    int coll_weights_count;
    double coll_weights[8];
    double pdf_f1[4];
    double pdf_f2[4];
    double naked_weight;
    double event_weight;
    double get_ren_scale() const { return orig_ren_scale()*ren_scalefactor; }
    double get_fac_scale() const { return orig_fac_scale()*fac_scalefactor; }
    double get_event_weight() const { return event_weight; }
    int get_event_order() const { return event_order; }

    // LoopSim parameters
    double opt_loopsim_R;
    int opt_loopsim_nborn;

    // qfilter
    int opt_filter_inq;
    int opt_filter_nq;

    // pdf sets
    int opt_frompdf, opt_topdf;
    int lhaid1, lhaid2;
    double x1r, x2r;
    int get_lhaid1() const { return lhaid1; }
    int get_lhaid2() const { return lhaid2; }
    double get_x1r() const { return x1r; }
    double get_x2r() const { return x2r; }

    double pdfx1[13];
    double pdfx2[13];
    double pdfx1p[13];
    double pdfx2p[13];

    // alphas
    bool use_sherpa_alphas;
    void initAS(const int order, const double asMZ, const double mZ2,
                const std::vector<double>& qmasses);
    SHERPA::One_Running_AlphaS* sherpa_alphas;

    // beta0 workaround
    int opt_beta0fix, opt_cdr2fdhfix, opt_pi2o12fix;
    static double beta0pole2(int id1_, int id2_, int n_, const Int_t* kf_);
    static double pole2(int id1_, int id2_, int n_, const Int_t* kf_);
    static double cdr2fdh(int id1_, int id2_, int n_, const Int_t* kf_);

    // Q2, x1, x2 stats
    void stat_report();
    double stat_Q2_min, stat_Q2_max;
    double stat_x1_min, stat_x1_max;
    double stat_x2_min, stat_x2_max;

    // eventoscope
    int opt_stat_step;
    double xsval_cur, xserr_cur;
    std::vector<double> xsvals;
    std::vector<double> xserrs;

    // stats and warnings
    long print_event_step;
    long pdf_warning_limit;
    long pdf_warning_count;
    double pdf_warning_thresh;

  protected:
    double get_alphas(double mur);
    void prepare_event();
    void reweight(const PseudoJetVector& input,
                  const PseudoJetVector& jets);
    void reweight_pdfcheck();
    double reweight_applyfixes(double new_me_wgt, double log_r);
    void stat_update(double fac_scale, double ren_scale);
};

#if defined(__MAKECINT__)
#pragma link C++ class fastjet::JetDefinition;
#endif

#endif
