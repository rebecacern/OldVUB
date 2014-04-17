//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb  3 15:13:08 2012 by ROOT version 5.32/00
// from TChain myTree/myTree
//////////////////////////////////////////////////////////

#ifndef freyaLooper_h
#define freyaLooper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>


// Fixed size dimensions of array or collections stored in the TTree if any.

class freyaLooper {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        xlWeight;
   Double_t        lum;
   Int_t           npu;
   Int_t           nvertex;
   Double_t        metPt;
   Double_t        metPx;
   Double_t        metPy;
   std::vector<int>     *flavLepton;
   std::vector<double>  *ecalRelIsoLepton;
   std::vector<double>  *trkRelIsoLepton;
   std::vector<double>  *hcalRelIsoLepton;
   std::vector<double>  *totalRelIsoLepton;
   std::vector<double>  *d0Lepton;
   std::vector<double>  *etaLepton;
   std::vector<double>  *phiLepton;
   std::vector<double>  *eOverpLepton;
   std::vector<double>  *sigmaIetaIetaLepton;
   std::vector<double>  *deltaEtaSCTrkAtVtxLepton;
   std::vector<double>  *deltaPhiSCTrkAtVtxLepton;
   std::vector<double>  *convDistLepton;
   std::vector<double>  *convDcotLepton;
   std::vector<int>     *numPixHitsLepton;
   std::vector<double>  *fbremLepton;
   std::vector<double>  *ptLepton;
   std::vector<double>  *pxLepton;
   std::vector<double>  *pyLepton;
   std::vector<double>  *pzLepton;
   std::vector<double>  *eLepton;
   std::vector<double>  *qLepton;
   std::vector<double>  *ptJet;
   std::vector<double>  *pxJet;
   std::vector<double>  *pyJet;
   std::vector<double>  *pzJet;
   std::vector<double>  *eJet;
   std::vector<double>  *qJet;
   std::vector<double>  *btTCHPJet;
   std::vector<double>  *btTCHEJet;
   std::vector<double>  *btSSVHPJet;
   std::vector<double>  *btSSVHEJet;

   // List of branches
   TBranch        *b_xlWeight;   //!
   TBranch        *b_lum;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_nvertex;   //!
   TBranch        *b_metPt;   //!
   TBranch        *b_metPx;   //!
   TBranch        *b_metPy;   //!
   TBranch        *b_flavLepton;   //!
   TBranch        *b_ecalRelIsoLepton;   //!
   TBranch        *b_trkRelIsoLepton;   //!
   TBranch        *b_hcalRelIsoLepton;   //!
   TBranch        *b_totalRelIsoLepton;   //!
   TBranch        *b_d0Lepton;   //!
   TBranch        *b_etaLepton;   //!
   TBranch        *b_phiLepton;   //!
   TBranch        *b_eOverpLepton;   //!
   TBranch        *b_sigmaIetaIetaLepton;   //!
   TBranch        *b_deltaEtaSCTrkAtVtxLepton;   //!
   TBranch        *b_deltaPhiSCTrkAtVtxLepton;   //!
   TBranch        *b_convDistLepton;   //!
   TBranch        *b_convDcotLepton;   //!
   TBranch        *b_numPixHitsLepton;   //!
   TBranch        *b_fbremLepton;   //!
   TBranch        *b_ptLepton;   //!
   TBranch        *b_pxLepton;   //!
   TBranch        *b_pyLepton;   //!
   TBranch        *b_pzLepton;   //!
   TBranch        *b_eLepton;   //!
   TBranch        *b_qLepton;   //!
   TBranch        *b_ptJet;   //!
   TBranch        *b_pxJet;   //!
   TBranch        *b_pyJet;   //!
   TBranch        *b_pzJet;   //!
   TBranch        *b_eJet;   //!
   TBranch        *b_qJet;   //!
   TBranch        *b_btTCHPJet;   //!
   TBranch        *b_btTCHEJet;   //!
   TBranch        *b_btSSVHPJet;   //!
   TBranch        *b_btSSVHEJet;   //!

   freyaLooper(TTree *tree=0);
   virtual ~freyaLooper();
   virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   void             myLoop(int nsel, int mode, bool silent);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

freyaLooper::freyaLooper(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("outputs/EleStudies_2_data.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("outputs/EleStudies_2_data.root");
      }
      f->GetObject("myTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("myTree","myTree");
      chain->Add("outputs/EleStudies_0_data.root/myTree");
      chain->Add("outputs/EleStudies_1_data.root/myTree");
      chain->Add("outputs/EleStudies_2_data.root/myTree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

freyaLooper::~freyaLooper()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t freyaLooper::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t freyaLooper::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void freyaLooper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   flavLepton = 0;
   ecalRelIsoLepton = 0;
   trkRelIsoLepton = 0;
   hcalRelIsoLepton = 0;
   totalRelIsoLepton = 0;
   d0Lepton = 0;
   etaLepton = 0;
   phiLepton = 0;
   eOverpLepton = 0;
   sigmaIetaIetaLepton = 0;
   deltaEtaSCTrkAtVtxLepton = 0;
   deltaPhiSCTrkAtVtxLepton = 0;
   convDistLepton = 0;
   convDcotLepton = 0;
   numPixHitsLepton = 0;
   fbremLepton = 0;
   ptLepton = 0;
   pxLepton = 0;
   pyLepton = 0;
   pzLepton = 0;
   eLepton = 0;
   qLepton = 0;
   ptJet = 0;
   pxJet = 0;
   pyJet = 0;
   pzJet = 0;
   eJet = 0;
   qJet = 0;
   btTCHPJet = 0;
   btTCHEJet = 0;
   btSSVHPJet = 0;
   btSSVHEJet = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("xlWeight", &xlWeight, &b_xlWeight);
   fChain->SetBranchAddress("lum", &lum, &b_lum);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("nvertex", &nvertex, &b_nvertex);
   fChain->SetBranchAddress("metPt", &metPt, &b_metPt);
   fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
   fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
   fChain->SetBranchAddress("flavLepton", &flavLepton, &b_flavLepton);
   fChain->SetBranchAddress("ecalRelIsoLepton", &ecalRelIsoLepton, &b_ecalRelIsoLepton);
   fChain->SetBranchAddress("trkRelIsoLepton", &trkRelIsoLepton, &b_trkRelIsoLepton);
   fChain->SetBranchAddress("hcalRelIsoLepton", &hcalRelIsoLepton, &b_hcalRelIsoLepton);
   fChain->SetBranchAddress("totalRelIsoLepton", &totalRelIsoLepton, &b_totalRelIsoLepton);
   fChain->SetBranchAddress("d0Lepton", &d0Lepton, &b_d0Lepton);
   fChain->SetBranchAddress("etaLepton", &etaLepton, &b_etaLepton);
   fChain->SetBranchAddress("phiLepton", &phiLepton, &b_phiLepton);
   fChain->SetBranchAddress("eOverpLepton", &eOverpLepton, &b_eOverpLepton);
   fChain->SetBranchAddress("sigmaIetaIetaLepton", &sigmaIetaIetaLepton, &b_sigmaIetaIetaLepton);
   fChain->SetBranchAddress("deltaEtaSCTrkAtVtxLepton", &deltaEtaSCTrkAtVtxLepton, &b_deltaEtaSCTrkAtVtxLepton);
   fChain->SetBranchAddress("deltaPhiSCTrkAtVtxLepton", &deltaPhiSCTrkAtVtxLepton, &b_deltaPhiSCTrkAtVtxLepton);
   fChain->SetBranchAddress("convDistLepton", &convDistLepton, &b_convDistLepton);
   fChain->SetBranchAddress("convDcotLepton", &convDcotLepton, &b_convDcotLepton);
   fChain->SetBranchAddress("numPixHitsLepton", &numPixHitsLepton, &b_numPixHitsLepton);
   fChain->SetBranchAddress("fbremLepton", &fbremLepton, &b_fbremLepton);
   fChain->SetBranchAddress("ptLepton", &ptLepton, &b_ptLepton);
   fChain->SetBranchAddress("pxLepton", &pxLepton, &b_pxLepton);
   fChain->SetBranchAddress("pyLepton", &pyLepton, &b_pyLepton);
   fChain->SetBranchAddress("pzLepton", &pzLepton, &b_pzLepton);
   fChain->SetBranchAddress("eLepton", &eLepton, &b_eLepton);
   fChain->SetBranchAddress("qLepton", &qLepton, &b_qLepton);
   fChain->SetBranchAddress("ptJet", &ptJet, &b_ptJet);
   fChain->SetBranchAddress("pxJet", &pxJet, &b_pxJet);
   fChain->SetBranchAddress("pyJet", &pyJet, &b_pyJet);
   fChain->SetBranchAddress("pzJet", &pzJet, &b_pzJet);
   fChain->SetBranchAddress("eJet", &eJet, &b_eJet);
   fChain->SetBranchAddress("qJet", &qJet, &b_qJet);
   fChain->SetBranchAddress("btTCHPJet", &btTCHPJet, &b_btTCHPJet);
   fChain->SetBranchAddress("btTCHEJet", &btTCHEJet, &b_btTCHEJet);
   fChain->SetBranchAddress("btSSVHPJet", &btSSVHPJet, &b_btSSVHPJet);
   fChain->SetBranchAddress("btSSVHEJet", &btSSVHEJet, &b_btSSVHEJet);
   Notify();
}

Bool_t freyaLooper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void freyaLooper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t freyaLooper::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif
