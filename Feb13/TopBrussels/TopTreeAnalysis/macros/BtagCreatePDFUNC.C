/* Some Usage info

Setup under CMSSW_4_2_6

ln -s /afs/cern.ch/cms/slc5_amd64_gcc434/external/lhapdf/5.6.0-cms4/full/lib/libLHAPDF*.so .                                                                                                          
ln -s /afs/cern.ch/cms/slc5_amd64_gcc434/external/lhapdf/5.6.0-cms4/include/LHAPDF LHAPDF                                                                                                         
ln -s /afs/cern.ch/cms/slc5_amd64_gcc434/external/lhapdf/5.6.0-cms4/share/lhapdf/PDFsets PDFsets        

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib/


g++ -g -L ~/lib -I . -L . -l BtagAnalysis41 -l TopTreeAnaContent41 -l TopTreeAna41 -l LHAPDF -I `root-config --incdir` `root-config --libs` BtagCreatePDFUNC.C -o BtagCreatePDFUNC

*/

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// Root stuff
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TBranch.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TDirectory.h"

//pdf-systematic stuff
#include "LHAPDF/LHAPDF.h"

#include "../Content/interface/Dataset.h"
#include "../BtagEffAnalysis/interface/TRootNTuple.h"

using namespace std;

int main (int argc, char *argv[]) {

  // init PDF stuff
  LHAPDF::initPDFSet(1, "cteq66.LHgrid");
  cout << " Initialized PDF stuff" << endl;

  vector<string> datasets;
  datasets.push_back("BtagTrees/BtagTree_TTbarJets_SemiMuon");
  datasets.push_back("BtagTrees/BtagTree_TTbarJets_Other");
  datasets.push_back("BtagTrees/BtagTree_WJets");
  datasets.push_back("BtagTrees/BtagTree_ZJets");
  datasets.push_back("BtagTrees/BtagTree_QCD");
  datasets.push_back("BtagTrees/BtagTree_ST_tChannel_t");
  datasets.push_back("BtagTrees/BtagTree_ST_tChannel_tbar");
  datasets.push_back("BtagTrees/BtagTree_ST_tWChannel_t");
  datasets.push_back("BtagTrees/BtagTree_ST_tWChannel_tbar");
  
  /*ifstream file_to_read("data.txt");

  while(!file_to_read.eof()) {

    float v1;
    float v2;
    float v3;

    file_to_read >> v1 >> v2 >> v3;

    cout << v1 << " " << v2 << " " << v3 << endl;

    }*/

  for (int d=0; d<datasets.size();d++) {
    //for (int d=0; d<1;d++) {

    TFile* inFile = new TFile((datasets[d]+".root").c_str(),"READ");

    string outFile=datasets[d]+"_PDFWeights.txt";

    TRootNTuple *NTuple;
    TTree *tree_ = (TTree*)inFile->Get("tree");
    NTuple = new TRootNTuple();
    
    TBranch *branch = tree_->GetBranch("TheNTuple");
    branch->SetAddress(&NTuple);
    
    Int_t nEvent = tree_->GetEntries();
    cout << "+----> DataSet "<< datasets[d] << " nEvent " << nEvent << endl;

    cout << "PDF stuff" << endl;
    vector< vector<float> > pdfInfo;
    vector< float > mlb;
    for(unsigned int iEvt=0; iEvt<nEvent; iEvt++) {
      pdfInfo.push_back(vector<float>(45,0));
      mlb.push_back(0);
    }

    for(int i=0; i <=44; ++i) {
      cout << "PDF-stuff: " << i << endl;
      LHAPDF::usePDFMember(1,i);
      for(unsigned int iEvt=0; iEvt<nEvent; iEvt++) {
	tree_->GetEvent(iEvt);

	mlb[iEvt] = NTuple->mlj();
	
	float scale = NTuple->getq2();
	
	pdfInfo[iEvt][i] = LHAPDF::xfx(1, NTuple->getx1(), scale, NTuple->getid1()) * LHAPDF::xfx(1, NTuple->getx2(), scale, NTuple->getid2());
	//cout << NTuple->getx1() << " " << NTuple->getid1() << " " << NTuple->getx2() << " " << NTuple->getid2() << " " << scale << endl;
	////cout << "i: " << i << "  iEvt: " << iEvt << "  pdfWeight: " << pdfInfo[iEvt].weight[i] << endl;//pdfInfo[iEvt][0] << endl;
	
	//cout << LHAPDF::xfx(1, NTuple->getx1(), scale, NTuple->getid1()) << endl;
      }
    }

    cout << "Finished running dataset " << datasets[d] << endl;

    cout << "Writing output to file " << outFile << endl;
    
    fstream out(outFile.c_str(), ios::out | ios::trunc);

    for(unsigned int iEvt=0; iEvt<nEvent; iEvt++) {
      out << iEvt << " " << mlb[iEvt];
      for(int i=0; i <=44; ++i) 
	out << " " << pdfInfo[iEvt][i];//pdfInfo[iEvt][0];
      out << "\n";
    }
  }

  
}
