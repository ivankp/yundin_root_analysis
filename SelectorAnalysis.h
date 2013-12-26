
#ifndef SELECTOR_ANALYSIS_H
#define SELECTOR_ANALYSIS_H

#include "SelectorCommon.h"

class SelectorCommon;

class Analysis
{
  public:
    typedef std::vector<fastjet::PseudoJet> PseudoJetVector;
    typedef std::vector<std::vector<Histogram*> > HistogramVector;
    static Analysis* create() { return new Analysis(); }

    Analysis();
    virtual ~Analysis();

    virtual bool check_cuts(SelectorCommon* event);
    virtual void analysis_bin(SelectorCommon* event);
    virtual void analysis_finalize();

    PseudoJetVector input;
    PseudoJetVector jets;

    void setJetNumber(const unsigned n);
    void setAntiKt(double R);

    fastjet::JetDefinition jet_definition;
    unsigned jet_number;
    double jet_ptmin;
    double jet_etamax;

    double event_count;
    double event_binned;

    LinearHistogram* jet_exclusive;
    LinearHistogram* jet_inclusive;
    std::vector<Histogram*> scale_wgt;
    std::vector<Histogram*> scale_nowgt;
    std::vector<std::vector<Histogram*> > jet_pt_n;
    std::vector<std::vector<Histogram*> > jet_eta_n;

    Grid* g_jet_inclusive;
    std::vector<Grid*> g_jet_pt_n;
    std::vector<Grid*> g_jet_eta_n;

    TString runname;

    virtual void reset();

  protected:
    std::set<TString> outputfiles;

    template <typename T>
    void clear_var(T*& var);

    void append_output_filename(const TString& name);

    void clear_histvec(std::vector<Histogram*>& histvec);
    void bin_histvec(const std::vector<Histogram*>& histvec,
                     int nextevt, double x, double w);
    void output_histvec(const std::vector<Histogram*>& histvec,
                        const TString& filename, std::ofstream& stream,
                        bool dryrun);

    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
    virtual void output_grids();
    virtual void clear();

    void fill_grid(Grid* grid, int nextevt, double x, double w, SelectorCommon* event);

    bool photonIsolation(const SelectorCommon* event, double photon_R,
                         double photon_n, double photon_eps) const;
};

class JetAnalysis : public Analysis
{
  public:
    static JetAnalysis* create() { return new JetAnalysis(); }

    JetAnalysis();

    virtual bool check_cuts(SelectorCommon* event);
    virtual void analysis_bin(SelectorCommon* event);

    double jet_pt1min;
    double jet_pt2min;
    double jet_pt3min;

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
};

class JetMAnalysis : public JetAnalysis
{
  public:
    typedef JetAnalysis BaseClass;

    static JetMAnalysis* create() { return new JetMAnalysis(); }

    JetMAnalysis();

    virtual bool check_cuts(SelectorCommon* event);
    virtual void analysis_bin(SelectorCommon* event);

    double ystar_min;
    double ystar_max;

    std::vector<Histogram*> jet_mass_jjj;

    Grid* g_jet_mass_jjj;

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
    virtual void output_grids();
    virtual void clear();
};

class PhotonJetAnalysis : public Analysis
{
  public:
    static PhotonJetAnalysis* create() { return new PhotonJetAnalysis(); }

    PhotonJetAnalysis();

    virtual bool check_cuts(SelectorCommon* event);
    virtual void analysis_bin(SelectorCommon* event);

    double jet_pt1min;  // extra cut on leading jet-pt
    double photon_R;
    double photon_n;
    double photon_eps;
    double photon_ptmin;
    double photon_etamax;
    double photon_jet_Rsep;

    std::vector<Histogram*> photon_pt;
    std::vector<Histogram*> photon_eta;
    std::vector<Histogram*> photon_jet_R11;
    std::vector<Histogram*> jet_jet_phi12;

    Grid* g_photon_pt;
    Grid* g_photon_eta;
    Grid* g_photon_jet_R11;
    Grid* g_jet_jet_phi12;

    virtual void reset();

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
    virtual void output_grids();
    virtual void clear();
};

class VJetAnalysis : public Analysis
{
  public:
    static VJetAnalysis* create() { return new VJetAnalysis(); }

    VJetAnalysis();

    virtual bool check_cuts(SelectorCommon* event);
    virtual void analysis_bin(SelectorCommon* event);

    double lepton_ptmin;
    double lepton_etamax;
    double lepton_etagap_min;
    double lepton_etagap_max;
    double etmiss_min;
    double vboson_mass_min;
    double vboson_mass_max;
    double lepton_jet_Rsep;
    double lepton_lepton_Rsep;
    double vboson_onshell_mass;

    std::vector<Histogram*> vboson_pt;
    std::vector<Histogram*> vboson_eta;

    Grid* g_vboson_pt;
    Grid* g_vboson_eta;

    virtual void reset();

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
    virtual void output_grids();
    virtual void clear();
};

class DiPhotonAnalysis : public Analysis
{
  public:
    static DiPhotonAnalysis* create() { return new DiPhotonAnalysis(); }

    DiPhotonAnalysis();

    virtual bool check_cuts(SelectorCommon* event);
    virtual void analysis_bin(SelectorCommon* event);

    double jet_pt1min;  // extra cut on leading jet-pt

    double photon_R;
    double photon_n;
    double photon_eps;
    double photon_pt1min;
    double photon_pt2min;
    double photon_etamax;

    double photon_photon_Rsep;
    double photon_jet_Rsep;

    std::vector<Histogram*> photon_mass;
    std::vector<Histogram*> photon_pt;
    std::vector<Histogram*> photon_eta;
    std::vector<Histogram*> photon_jet_R11;
    std::vector<Histogram*> jet_jet_phi12;
    std::vector<Histogram*> jet_jet_mass;
    std::vector<Histogram*> jet_jet_eta12;
    std::vector<Histogram*> diphoton_dijet_phi;
    std::vector<Histogram*> diphoton_dijet_ystar;

    Grid* g_photon_mass;
    Grid* g_photon_pt;
    Grid* g_photon_eta;
    Grid* g_photon_jet_R11;
    Grid* g_jet_jet_phi12;

    virtual void reset();

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
    virtual void output_grids();
    virtual void clear();
};

class DiPhotonAnalysisBH : public DiPhotonAnalysis
{
  public:
    static DiPhotonAnalysisBH* create() { return new DiPhotonAnalysisBH(); }

    virtual bool check_cuts(SelectorCommon* event);
};

#if defined(__MAKECINT__)
#pragma link C++ class Analysis;
#pragma link C++ class JetAnalysis;
#pragma link C++ class JetMAnalysis;
#pragma link C++ class VJetAnalysis;
#pragma link C++ class PhotonJetAnalysis;
#pragma link C++ class DiPhotonAnalysis;
#pragma link C++ class DiPhotonAnalysisBH;
#pragma link C++ class std::vector<Histogram*>;
#pragma link C++ class std::vector<std::vector<Histogram*> >;
#pragma link C++ class std::vector<Grid*>;
#endif

#endif // SELECTOR_ANALYSIS_H
