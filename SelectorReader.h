
#ifndef SELECTOR_READER_H
#define SELECTOR_READER_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Stuff

#define MAXNPARTICLE 100
#define MAXNUWEIGHT 32

class SelectorCommon;

class SelectorReader : public TSelector
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
    Float_t         ntuple_E[MAXNPARTICLE];    //[nparticle]
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
    Char_t          ntuple_part[2];
    Char_t          ntuple_alphaspower;

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
    TBranch        *b_part;   //!
    TBranch        *b_alphaspower;   //!

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

    ClassDef(SelectorReader, 0);

    // -----------------------------------------------------------------------
    // ROOT stuff END           ROOT stuff END          ROOT stuff END
    // -----------------------------------------------------------------------

    SelectorReader(TTree* tree=0);
    ~SelectorReader();

    void addSelector(SelectorCommon* selector);

  protected:
    std::vector<SelectorCommon*> selectors;
};

#endif
