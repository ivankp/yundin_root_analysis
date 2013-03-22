//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 15 14:34:28 2013 by ROOT version 5.34/05
// from TTree T3/Reconst ntuple
// found on file: RSevent1.root
//////////////////////////////////////////////////////////

#ifndef T3selector_h
#define T3selector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Stuff

#include <TRegexp.h>
#include <fastjet/ClusterSequence.hh>
#include <fstream>

#include "SherpaAlphaS.h"

class Histogram
{
  public:
    Histogram(const TString& filename_, const TString& name_,
              int nbin_, double x1_, double x2_);
    virtual ~Histogram() {};

    virtual void bin(int nextevt, double x, double w) = 0;
    void print(std::ostream& stream, const TString& runname, double count, bool unweight=true);

    TString getFile() const;

  protected:
    void flush();
    void setedges();
    void fill(int evt, int n, double w);

    TString filename;
    TString name;

    double x1, x2, x12;
    const int nbin;
    int prevevt, lastidx;
    std::vector<int> curidx;
    std::vector<int> events;
    std::vector<double> curwgt;
    std::vector<double> wgt;
    std::vector<double> wgt2;
    std::vector<double> bwidth;
    std::vector<double> edge;
};

class LinearHistogram : public Histogram
{
  public:
    LinearHistogram(const TString& filename_, const TString& name_,
                    int nbin_, double x1_, double x2_);
    void bin(int nextevt, double x, double w);
};


class QuadraticHistogram : public Histogram
{
  public:
    QuadraticHistogram(const TString& filename_, const TString& name_,
                       int nbin_, double x1_, double x2_, double f);
    void bin(int nextevt, double x, double w);

  private:
    const double step;
    const double slope;
};

class T3analysis
{
  public:
    T3analysis();
    void setN(const unsigned njet_);
    void init(const TString& opt);

    void analysis_bin(const Int_t id, const Double_t weight,
                      const std::vector<fastjet::PseudoJet>& jets);
    void analysis_finalize();

    unsigned njet;
    fastjet::JetAlgorithm jetalgo;
    fastjet::JetDefinition jetdef;
    double jetR;
    double ptmin;
    double jet1ptmin;
    double etacut;

    double event_count;

    LinearHistogram* jet_exclusive;
    LinearHistogram* jet_inclusive;
    std::vector<std::vector<Histogram*> > jet_pt_n;
    std::vector<std::vector<Histogram*> > jet_eta_n;

    TString runname;

    void reset();

    void addPtLinearHistograms(TString filename, int nbins, std::vector<int> ptlimits);
    void addPtQuadraticHistograms(TString filename, int nbins, double f, std::vector<int> ptlimits);
    void addEtaLinearHistograms(TString filename, int nbins);
    void addEtaQuadraticHistograms(TString filename, int nbins, double f);

  private:
    double OptF(const TString subopt);
    int OptI(const TString subopt);

    std::set<TString> outputfiles;

    void clear();
};

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

#define MAXNPARTICLE 16
#define MAXNUWEIGHT 128

class T3selector : public TSelector
{
  public :
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

    T3selector(TTree * /*tree*/ =0)
    : fChain(0),
      use_sherpa_alphas(false), sherpa_alphas(0), beta0fix(0)
    { }

    virtual ~T3selector() { }
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

    ClassDef(T3selector, 0);

    // analysis stuff
    T3analysis analysis;

    int alphapower;
    double scalefactor;

    int FROMPDF, TOPDF;

    bool use_sherpa_alphas;
    void initAS(const int order, const double asMZ, const double mZ2,
                const std::vector<double>& qmasses);
    SHERPA::One_Running_AlphaS* sherpa_alphas;

    int beta0fix;

  protected:
    void reweight();
    double beta0pole2();
};

#if defined(__MAKECINT__)
#pragma link C++ class fastjet::JetDefinition;
#endif

#endif

#ifdef T3selector_cxx
void T3selector::Init(TTree *tree)
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

Bool_t T3selector::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef T3selector_cxx
