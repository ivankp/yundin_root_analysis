
#ifndef SELECTOR_COMMON_H
#define SELECTOR_COMMON_H

#include <Rtypes.h>

// Stuff
#include <fastjet/ClusterSequence.hh>

// Math
#include <inttypes.h>
#include <cstdlib>
#include <cmath>
using std::abs;
using std::pow;
using std::sqrt;

#include "SherpaAlphaS.h"

#ifndef DISABLE_LOOPSIM
#include <LoopSim.hh>
#endif // DISABLE_LOOPSIM

namespace RootAnalysis {
class NTupleEvent;
}

class SelectorReader;
class Analysis;

class SelectorCommon
{
  public :
    // Declaration of leaf types
    const Int_t*          input_id;
    const Int_t*          input_nparticle;
    const Float_t*        input_px;   //[nparticle]
    const Float_t*        input_py;   //[nparticle]
    const Float_t*        input_pz;   //[nparticle]
    const Float_t*        input_E;    //[nparticle]
    const Double_t*       input_alphas;
    const Int_t*          input_kf;   //[nparticle]
    const Double_t*       input_weight;
    const Double_t*       input_weight2;
    const Double_t*       input_me_wgt;
    const Double_t*       input_me_wgt2;
    const Double_t*       input_x1;
    const Double_t*       input_x2;
    const Double_t*       input_x1p;
    const Double_t*       input_x2p;
    const Int_t*          input_id1;
    const Int_t*          input_id2;
    const Double_t*       input_fac_scale;
    const Double_t*       input_ren_scale;
    const Int_t*          input_nuwgt;
    const Double_t*       input_usr_wgts;   //[nuwgt]
    const Char_t*         input_part;
    const Short_t*        input_alphaspower;

    // Accessor methods
    Int_t           get_event_id() const { return *input_id; }
    Int_t           get_nparticle() const { return *input_nparticle; }
    const Float_t*  get_px() const { return input_px; }
    const Float_t*  get_py() const { return input_py; }
    const Float_t*  get_pz() const { return input_pz; }
    const Float_t*  get_E() const { return input_E; }
    Float_t         get_px(int i) const { return input_px[i]; }
    Float_t         get_py(int i) const { return input_py[i]; }
    Float_t         get_pz(int i) const { return input_pz[i]; }
    Float_t         get_E(int i) const { return input_E[i]; }
    Double_t        orig_alphas() const { return *input_alphas; }
    const Int_t*    get_kf() const { return input_kf; }
    Int_t           get_kf(int i) const { return input_kf[i]; }
    Double_t        orig_weight() const { return *input_weight; }
    Double_t        orig_weight2() const { return *input_weight2; }
    Double_t        orig_me_wgt() const { return *input_me_wgt; }
    Double_t        orig_me_wgt2() const { return *input_me_wgt2; }
    Double_t        get_x1() const { return *input_x1; }
    Double_t        get_x2() const { return *input_x2; }
    Double_t        get_x1p() const { return *input_x1p; }
    Double_t        get_x2p() const { return *input_x2p; }
    Int_t           get_id1() const { return *input_id1; }
    Int_t           get_id2() const { return *input_id2; }
    Double_t        orig_fac_scale() const { return *input_fac_scale; }
    Double_t        orig_ren_scale() const { return *input_ren_scale; }
    Int_t           get_nuwgt() const { return *input_nuwgt; }
    const Double_t* orig_usr_wgts() const { return input_usr_wgts; }
    Double_t        orig_usr_wgts(int i) const { return input_usr_wgts[i]; }
    Int_t           get_alphaspower() const { return Int_t(*input_alphaspower); }
    const Char_t*   get_part() const { return input_part; }
    Char_t          get_part(int i) const { return input_part[i]; }

    void Init(const SelectorReader* reader);
    void Init(const RootAnalysis::NTupleEvent* event);
    bool Process();
    void SlaveBegin();
    void SlaveTerminate();

    SelectorCommon();
    ~SelectorCommon();

    // analysis
    Analysis* analysis;
    typedef std::vector<fastjet::PseudoJet> PseudoJetVector;

    fastjet::PseudoJet get_vec(int i) const;
    PseudoJetVector get_fjinput() const;
#ifndef DISABLE_LOOPSIM
    std::vector<LSParticle> get_lsinput() const;
    static PseudoJetVector lsinput2fjinput(const std::vector<LSParticle>& in);
#endif


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
    double rescaler_maa(const double scale,
                       const PseudoJetVector& input,
                       const PseudoJetVector& jets);
    double rescaler_maaht(const double scale,
                       const PseudoJetVector& input,
                       const PseudoJetVector& jets);
    double rescaler_maahthat(const double scale,
                          const PseudoJetVector& input,
                          const PseudoJetVector& jets);
    double rescaler_mwhthat(const double scale,
                       const PseudoJetVector& input,
                       const PseudoJetVector& jets);
    double rescaler_mwFhthat(const double scale,
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
    void setrescaler_mwhthat() { opt_rescaler = &SelectorCommon::rescaler_mwhthat; }
    void setrescaler_mwFhthat() { opt_rescaler = &SelectorCommon::rescaler_mwFhthat; }
    void setrescaler_minlo();

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

    // exact reweighting
    double event_trials;

  protected:
    double get_alphas(double mur);
    void prepare_event();
    void reweight(const PseudoJetVector& input,
                  const PseudoJetVector& jets);
    void reweight_pdfcheck();
    double reweight_applyfixes(double new_me_wgt, double log_r);
    void stat_update(double fac_scale, double ren_scale);

  public:
    intptr_t this_to_int() const { return reinterpret_cast<intptr_t>(this); }
};

#if defined(__MAKECINT__)
#pragma link C++ class fastjet::JetDefinition;
#endif

#endif
