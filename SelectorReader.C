
#include <TH2.h>
#include <TStyle.h>

#ifndef NDEBUG
  #include <iostream>
#endif

#include "SelectorCommon.h"
#include "SelectorReader.h"

// --------------------------------------------------------------------------- //
// Selector
// --------------------------------------------------------------------------- //

SelectorReader::SelectorReader(TTree* /*tree*/)
  : fChain(0) // ROOT
{

}


SelectorReader::~SelectorReader()
{

}

void SelectorReader::Init(TTree *tree)
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

  fChain->SetBranchAddress("id", &ntuple_id, &b_id);
  fChain->SetBranchAddress("nparticle", &ntuple_nparticle, &b_nparticle);
  fChain->SetBranchAddress("px", ntuple_px, &b_px);
  fChain->SetBranchAddress("py", ntuple_py, &b_py);
  fChain->SetBranchAddress("pz", ntuple_pz, &b_pz);
  fChain->SetBranchAddress("E", ntuple_E, &b_E);
  fChain->SetBranchAddress("alphas", &ntuple_alphas, &b_alphas);
  fChain->SetBranchAddress("kf", ntuple_kf, &b_kf);
  fChain->SetBranchAddress("weight", &ntuple_weight, &b_weight);
  fChain->SetBranchAddress("weight2", &ntuple_weight2, &b_weight2);
  fChain->SetBranchAddress("me_wgt", &ntuple_me_wgt, &b_me_wtg);
  fChain->SetBranchAddress("me_wgt2", &ntuple_me_wgt2, &b_me_wtg2);
  fChain->SetBranchAddress("x1", &ntuple_x1, &b_x1);
  fChain->SetBranchAddress("x2", &ntuple_x2, &b_x2);
  fChain->SetBranchAddress("x1p", &ntuple_x1p, &b_x1p);
  fChain->SetBranchAddress("x2p", &ntuple_x2p, &b_x2p);
  fChain->SetBranchAddress("id1", &ntuple_id1, &b_id1);
  fChain->SetBranchAddress("id2", &ntuple_id2, &b_id2);
  fChain->SetBranchAddress("fac_scale", &ntuple_fac_scale, &b_fac_scale);
  fChain->SetBranchAddress("ren_scale", &ntuple_ren_scale, &b_ren_scale);
  fChain->SetBranchAddress("nuwgt", &ntuple_nuwgt, &b_nuwgt);
  fChain->SetBranchAddress("usr_wgts", ntuple_usr_wgts, &b_usr_wgts);
  fChain->SetBranchAddress("alphasPower", &ntuple_alphaspower, &b_alphaspower);
  fChain->SetBranchAddress("part", ntuple_part, &b_part);
}

Bool_t SelectorReader::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  if (fChain && fChain->GetCurrentFile()) {
    std::cout << "File: " << fChain->GetCurrentFile()->GetName() << std::endl;
  }

  return kTRUE;
}

void SelectorReader::Begin(TTree* /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

}

void SelectorReader::SlaveBegin(TTree* /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  fastjet::ClusterSequence::set_fastjet_banner_stream(0); // silence fastjet
}


Bool_t SelectorReader::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either SelectorReader::GetEntry() or TBranch::GetEntry()
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

  bool rval;
  for (unsigned i = 0; i < selectors.size(); i++) {
    rval = selectors[i]->Process();
    if (not rval) {
      std::cerr << "Selector " << i << " failed to Process() event" << std::endl;
      Abort("Analysis failed");
    }
  }

  return true;
}

void SelectorReader::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  for (unsigned i = 0; i < selectors.size(); i++) {
    selectors[i]->SlaveTerminate();
  }
}

void SelectorReader::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}

void SelectorReader::addSelector(SelectorCommon* selector)
{
  selector->Init(this);
  selectors.push_back(selector);
}
