
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
    std::vector<std::vector<Histogram*> > jet_pt_n;
    std::vector<std::vector<Histogram*> > jet_eta_n;

    Grid* g_jet_inclusive;

    TString runname;

    virtual void reset();

    void addPtLinearHistograms(TString filename, int nbins,
                               std::vector<double>* ptlowlimits=0, std::vector<double>* pthighlimits=0);
    void addPtQuadraticHistograms(TString filename, int nbins, double f,
                                  std::vector<double>* ptlowlimits=0, std::vector<double>* pthighlimits=0);

    void addPtSmearedLinearHistograms(TString filename, int nbins, double s,
                                      std::vector<double>* ptlowlimits=0, std::vector<double>* pthighlimits=0);
    void addPtSmearedQuadraticHistograms(TString filename, int nbins, double f, double s,
                                         std::vector<double>* ptlowlimits=0, std::vector<double>* pthighlimits=0);

    void addEtaLinearHistograms(TString filename, int nbins,
                                std::vector<double>* etalowlimits=0, std::vector<double>* etahighlimits=0);
    void addEtaQuadraticHistograms(TString filename, int nbins, double f,
                                   std::vector<double>* etalowlimits=0, std::vector<double>* etahighlimits=0);

    void addEtaSmearedLinearHistograms(TString filename, int nbins, double s,
                                       std::vector<double>* etalowlimits=0, std::vector<double>* etahighlimits=0);
    void addEtaSmearedQuadraticHistograms(TString filename, int nbins, double f, double s,
                                          std::vector<double>* etalowlimits=0, std::vector<double>* etahighlimits=0);

  protected:
    std::set<TString> outputfiles;

    template <typename T>
    void clear_var(T*& var);

    virtual void output_histograms(const TString& filename, std::ofstream& stream);
    virtual void clear();

    void fill_grid(Grid* grid, int nextevt, double x, double w, SelectorCommon* event);

    template <typename T>
    void addPtHistograms(TString filename, int nbins,
                          double param1, double param2, double param3,
                          double low, double high,
                          std::vector<double>* lowlimits=0,
                          std::vector<double>* highlimits=0);
    template <typename T>
    void addEtaHistograms(TString filename, int nbins,
                          double param1, double param2, double param3,
                          double low, double high,
                          std::vector<double>* lowlimits=0,
                          std::vector<double>* highlimits=0);
};

class JetAnalysis : public Analysis
{
  public:
    static JetAnalysis* create() { return new JetAnalysis(); }

    JetAnalysis();

    virtual bool check_cuts(SelectorCommon* event);
    virtual void analysis_bin(SelectorCommon* event);

    double jet_pt1min;

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream);
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

    LinearHistogram* photon_mass;
    LinearHistogram* photon_jet_R11;
    LinearHistogram* jet_jet_phi12;

    virtual void reset();

  protected:
    virtual void output_histograms(const TString& filename, std::ofstream& stream);

    virtual void clear();

    template <class T>
    void addEtaHistograms(TString filename, int nbins,
                          double param1, double param2, double param3,
                          double low, double high,
                          std::vector<double>* lowlimits=0,
                          std::vector<double>* highlimits=0);
};

#if defined(__MAKECINT__)
#pragma link C++ class Analysis;
#pragma link C++ class JetAnalysis;
#pragma link C++ class DiPhotonAnalysis;
#endif

#endif // SELECTOR_ANALYSIS_H
