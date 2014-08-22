
#ifndef SELECTOR_ANALYSIS_H
#define SELECTOR_ANALYSIS_H

#include "SelectorCommon.h"
#include "SelectorHistograms.h"

class SelectorCommon;

class Analysis
{
  public:
    typedef std::vector<fastjet::PseudoJet> PseudoJetVector;
    typedef std::vector<std::vector<Histogram*> > HistogramVector;
    static Analysis* create() { return new Analysis(); }

    Analysis();
    virtual ~Analysis();

    void set_input(PseudoJetVector newinput);
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

    void bin_histvec(const std::vector<LinearHistogram2D*>& histvec,
                     int nextevt, double x, double y, double w);
    void output_histvec(const std::vector<LinearHistogram2D*>& histvec,
                        const TString& filename, std::ofstream& stream,
                        bool dryrun);


    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
    virtual void output_grids();
    virtual void clear();

    void fill_grid(Grid* grid, int nextevt, double x, double w, SelectorCommon* event);
};

class JetAnalysis : public Analysis
{
  public:
    static JetAnalysis* create() { return new JetAnalysis(); }

    JetAnalysis();

    virtual bool check_cuts(SelectorCommon* event);
    virtual void analysis_bin(SelectorCommon* event);

    double jet_pt1min;
    double jet_ht2min;

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
};

class Jet3Analysis : public JetAnalysis
{
  public:
    static Jet3Analysis* create() { return new Jet3Analysis(); }

    Jet3Analysis();

    virtual bool check_cuts(SelectorCommon* event);
    virtual void analysis_bin(SelectorCommon* event);

    double jet_eta1max, jet_eta2max, jet_eta2min;
    double jet_jet_DR23min, jet_jet_DR23max, jet_jet_M12min;

    std::vector<Histogram*> jet_jet_eta23;
    std::vector<Histogram*> jet_jet_phi23;
    std::vector<Histogram*> jet_jet_beta23;

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
};

class FourJetMPIAnalysis : public JetAnalysis
{
  public:
    static FourJetMPIAnalysis* create() { return new FourJetMPIAnalysis(); }

    FourJetMPIAnalysis();
    int mpivars_d12_bins;
    double mpivars_d12_bin_low;
    double mpivars_d12_bin_high;

    virtual bool check_cuts(SelectorCommon* event);
    virtual void analysis_bin(SelectorCommon* event);

    std::vector<LinearHistogram2D*> jets_d12_d34;

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream,
                                   bool dryrun);
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

class DiPhotonAnalysis : public Analysis
{
  public:
    static DiPhotonAnalysis* create() { return new DiPhotonAnalysis(); }

    DiPhotonAnalysis();

    virtual bool check_cuts(SelectorCommon* event);
    virtual void analysis_bin(SelectorCommon* event);

    double photon_R;
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

#if defined(__MAKECINT__)
#pragma link C++ class Analysis;
#pragma link C++ class JetAnalysis;
#pragma link C++ class Jet3Analysis;
#pragma link C++ class FourJetMPIAnalysis;
#pragma link C++ class PhotonJetAnalysis;
#pragma link C++ class DiPhotonAnalysis;
#pragma link C++ class std::vector<Histogram*>;
#pragma link C++ class std::vector<LinearHistogram2D*>;
#pragma link C++ class std::vector<std::vector<Histogram*> >;
#pragma link C++ class std::vector<Grid*>;
#endif

#endif // SELECTOR_ANALYSIS_H
