
#ifndef SELECTOR_ANALYSIS_H
#define SELECTOR_ANALYSIS_H

#include "SelectorCommon.h"
#include "SelectorHistograms.h"

class SelectorCommon;

class Analysis
{
  public:
    typedef std::vector<fastjet::PseudoJet> PseudoJetVector;
    typedef std::vector<std::vector<HistogramBase*> > HistogramVector;
    static Analysis* create() { return new Analysis(); }

    Analysis();
    virtual ~Analysis();

    void set_input(PseudoJetVector newinput);
    virtual bool check_cuts(const SelectorCommon* event);
    virtual void analysis_bin(const SelectorCommon* event);
    virtual void analysis_finalize(const SelectorCommon* event);

    PseudoJetVector input;
    PseudoJetVector jets;

    void setJetNumber(const unsigned n);
    void setAntiKt(double R);

    fastjet::JetDefinition jet_definition;
    unsigned jet_number;
    double jet_ptmin;
    double jet_etamax;
    double jet_ymax;

    double call_count;
    double event_count;
    double event_binned;

    LinearHistogram* jet_exclusive;
    LinearHistogram* jet_inclusive;
    std::vector<HistogramBase*> scale_wgt;
    std::vector<HistogramBase*> scale_nowgt;
    std::vector<HistogramBase*> jet_ht;
    std::vector<std::vector<HistogramBase*> > jet_pt_n;
    std::vector<std::vector<HistogramBase*> > jet_eta_n;
    std::vector<std::vector<HistogramBase*> > jet_y_n;

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

    void clear_histvec(std::vector<HistogramBase*>& histvec);
    void clear_histvec(std::vector<std::vector<HistogramBase*> >& histvecvec);
    void output_histvec(const std::vector<HistogramBase*>& histvec,
                        const TString& filename, std::ofstream& stream,
                        bool dryrun);

    void bin_histvec(const std::vector<HistogramBase*>& histvec,
                     int nextevt, double x, double w);
    void bin_histvec(const std::vector<HistogramBase*>& histvec,
                     int nextevt, double x, double y, double w);


    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
    virtual void output_grids();
    virtual void clear();
    virtual void finalize_stat(std::ostream& stream, const SelectorCommon* event);

    void fill_grid(Grid* grid, int nextevt, double x, double w, const SelectorCommon* event);

    bool photonIsolation(const SelectorCommon* event, double photon_R,
                         double photon_n, double photon_eps) const;
};

class JetAnalysis : public Analysis
{
  public:
    static JetAnalysis* create() { return new JetAnalysis(); }

    JetAnalysis();

    virtual bool check_cuts(const SelectorCommon* event);
    virtual void analysis_bin(const SelectorCommon* event);

    double jet_pt1min;
    double jet_pt2min;
    double jet_pt3min;

    double jet_ht2min;

    std::vector<HistogramBase*> jet_pt12ave;

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
    virtual void clear();
};

class Jet3Analysis : public JetAnalysis
{
  public:
    typedef JetAnalysis BaseClass;

    static Jet3Analysis* create() { return new Jet3Analysis(); }

    Jet3Analysis();

    virtual bool check_cuts(const SelectorCommon* event);
    virtual void analysis_bin(const SelectorCommon* event);

    double jet_eta1max, jet_eta2max, jet_eta2min;
    double jet_jet_DR23min, jet_jet_DR23max, jet_jet_M12min;

    std::vector<HistogramBase*> jet_jet_eta23;
    std::vector<HistogramBase*> jet_jet_phi23;
    std::vector<HistogramBase*> jet_jet_beta23;

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
    virtual void clear();
};

class FourJetMPIAnalysis : public JetAnalysis
{
  public:
    typedef JetAnalysis BaseClass;

    static FourJetMPIAnalysis* create() { return new FourJetMPIAnalysis(); }

    FourJetMPIAnalysis();
    int mpivars_d12_bins;
    double mpivars_d12_bin_low;
    double mpivars_d12_bin_high;

    virtual bool check_cuts(const SelectorCommon* event);
    virtual void analysis_bin(const SelectorCommon* event);

    std::vector<HistogramBase*> jets_d12_d34;

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
    virtual void clear();
};

class JetMAnalysis : public JetAnalysis
{
  public:
    typedef JetAnalysis BaseClass;

    static JetMAnalysis* create() { return new JetMAnalysis(); }

    JetMAnalysis();

    virtual bool check_cuts(const SelectorCommon* event);
    virtual void analysis_bin(const SelectorCommon* event);

    double ystar_min;
    double ystar_max;

    std::vector<HistogramBase*> jet_mass_jjj;

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

    virtual bool check_cuts(const SelectorCommon* event);
    virtual void analysis_bin(const SelectorCommon* event);

    double jet_pt1min;  // extra cut on leading jet-pt
    double photon_R;
    double photon_n;
    double photon_eps;
    double photon_ptmin;
    double photon_etamax;
    double photon_jet_Rsep;

    std::vector<HistogramBase*> photon_pt;
    std::vector<HistogramBase*> photon_eta;
    std::vector<HistogramBase*> photon_jet_R11;
    std::vector<HistogramBase*> jet_jet_phi12;

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

    virtual bool check_cuts(const SelectorCommon* event);
    virtual void analysis_bin(const SelectorCommon* event);

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

    fastjet::PseudoJet vboson;  // store vector boson momentum

    std::vector<HistogramBase*> vboson_pt;
    std::vector<HistogramBase*> vboson_eta;

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

    virtual bool check_cuts(const SelectorCommon* event);
    virtual void analysis_bin(const SelectorCommon* event);

    double jet_pt1min;  // extra cut on leading jet-pt

    double photon_R;
    double photon_n;
    double photon_eps;
    double photon_pt1min;
    double photon_pt2min;
    double photon_etamax;

    double photon_photon_Rsep;
    double photon_jet_Rsep;

    std::vector<HistogramBase*> photon_mass;
    std::vector<HistogramBase*> photon_pt;
    std::vector<HistogramBase*> photon_eta;
    std::vector<HistogramBase*> photon_jet_R11;
    std::vector<HistogramBase*> jet_jet_phi12;
    std::vector<HistogramBase*> jet_jet_mass;
    std::vector<HistogramBase*> jet_jet_eta12;
    std::vector<HistogramBase*> diphoton_dijet_phi;
    std::vector<HistogramBase*> diphoton_dijet_ystar;

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

    virtual bool check_cuts(const SelectorCommon* event);
};

#if defined(__MAKECINT__)
#pragma link C++ class Analysis;
#pragma link C++ class JetAnalysis;
#pragma link C++ class Jet3Analysis;
#pragma link C++ class FourJetMPIAnalysis;
#pragma link C++ class JetMAnalysis;
#pragma link C++ class VJetAnalysis;
#pragma link C++ class PhotonJetAnalysis;
#pragma link C++ class DiPhotonAnalysis;
#pragma link C++ class DiPhotonAnalysisBH;
#pragma link C++ class std::vector<HistogramBase*>;
#pragma link C++ class std::vector<std::vector<HistogramBase*> >;
#pragma link C++ class std::vector<Grid*>;
#endif

#endif // SELECTOR_ANALYSIS_H
