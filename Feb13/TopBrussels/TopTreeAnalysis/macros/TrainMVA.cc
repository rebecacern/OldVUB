#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib> 
#include <time.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/MultiCutPlot.h"
#include "../Tools/interface/CutImpactEvaluation.h"
#include "../Tools/interface/TemplateComparator.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/MVATrainer.h"
#include "../Tools/interface/MVAComputer.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/MCExpectation.h"
#include "../Content/interface/MCObsExpectation.h"
#include "../Content/interface/Dataset.h"
#include "../MCInformation/interface/MCWeighter.h"
#include "../Selection/interface/ElectronPlotter.h"
#include "../Selection/interface/MuonPlotter.h"
#include "../Selection/interface/JetPlotter.h"
#include "../Selection/interface/VertexPlotter.h"

#include "Style.C"


#define PI 3.14159265

using namespace std;
using namespace TopTree;

int main (int argc, char *argv[])
{

  // init clock

  clock_t start = clock();

  cout << "*******************************************" << endl;
  cout << " Begining of the program for MVA TRAINING ! " << endl;
  cout << "*******************************************" << endl;

  //SetStyle if needed
  setTDRStyle(); 
  //setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////
  ////////////////////////////////////////

  //xml file
  string xmlFileName ="../config/myJESconfig.xml";
  if (argc >= 2)
    xmlFileName=string(argv[1]);
  const char *xmlfile = xmlFileName.c_str();
  //Output files
  string rootFileName ("Output.root");

  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
 
  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<"Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;

  /////////////////////
  // Load MVATrainer //
  /////////////////////

  MVATrainer* trainer = new MVATrainer("TestTrainer","TMVA_output.root");

  trainer->bookInputVar("pt");
  trainer->bookInputVar("eta");

  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(datasets[i]);
  /////////////////////
  
  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex > vertex;
  vector < TRootMuon > muons;
  vector < TRootElectron > electrons;
  vector < TRootJet > jets;
  vector < TRootCaloJet > CaloJets;
  vector < TRootPFJet > PFJets;
  vector < TRootMET > mets;
  
  //////////////////////////
  // CREATE OUTPUT FILE
  //////////////////////////

  TFile *fout = new TFile (rootFileName.c_str (), "RECREATE");
  //Global variable
  TRootEvent* event = 0;

  //nof selected events
  float NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];

  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++)
  {

    if (verbose > 1)
      cout << "    * Dataset " << d << ": " << datasets[d].Name () << "/ title : " << datasets[d].Title () << endl;

    //open files and load
    cout<<"     - LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    //cout<<"LoadEvent"<<endl;
 

    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////

    int eventsSelected = 0;

    nEvents[d] = 0;
    if (verbose > 1)
      cout << "     - Loop over events " << endl;

    for (unsigned int ievt = 0; ievt < datasets[d].NofEvtsToRunOver(); ievt++) {

      nEvents[d]++;

      double perc = round(((double)ievt/(double)datasets[d].NofEvtsToRunOver())*100);

      int nBlocks = (int)perc/5;

      string blocks = "";
 
      for (int i=0; i<nBlocks; i++)
	blocks += "=";

      blocks += ">";

      for (int i=0; i<20-nBlocks; i++)
	blocks += " ";

      if(ievt%5000 == 0 || ievt == datasets[d].NofEvtsToRunOver()-2 ) std::cout<< "      -> Processing the "<<ievt<<"th event, [" << perc << "% " << blocks << "]. Time elapsed: " << (((double)clock() - start) / CLOCKS_PER_SEC)/60 << " minutes." << flush<<"\r";
     
      //load event
      if(anaEnv.JetType == 0 ) event = treeLoader.LoadEvent (ievt, &datasets[d], vertex, muons, electrons, jets, mets);
      if(anaEnv.JetType == 1 ) event = treeLoader.LoadEvent (ievt, &datasets[d], vertex, muons, electrons, CaloJets, mets);
      if(anaEnv.JetType == 2 ) event = treeLoader.LoadEvent (ievt, &datasets[d], vertex, muons, electrons, PFJets, mets);

      for (unsigned int i = 0; i<CaloJets.size(); i++) {
	  
	TRootCaloJet* jet = (TRootCaloJet*) &CaloJets[i];
	
	trainer->Fill("S","pt",jet->Pt());
	trainer->Fill("S","eta",jet->Eta());

	trainer->Fill("B","pt",jet->Pt()*0.5);
	trainer->Fill("B","eta",jet->Eta()*0.5);
      
      }

      ///////////////////////////////////////
    
       
     
    }				//loop on events
    cout<<endl;

    //important: free memory
    treeLoader.UnLoadDataset();

    std::cout<<"      -> This sample contains, " << datasets[d].NofEvtsToRunOver() << " events." << endl;

  }				//loop on datasets

  //Once everything is filled ...
  if (verbose > 0)
    cout << " - We ran over all the data ;-)" << endl;
  
  trainer->TrainMVA("Likelihood","Block","",0,0,"",0,0);

  ///////////////////
  // Writting
  //////////////////
  if (verbose > 1)
  	cout << "Writting the plots in the output root-file ..." << endl;

   //add configuration info
  fout->cd();
  configTree->Fill();
  configTree->Write();

  
  fout->Close();

  if (verbose > 1)
    cout << " - Done with writing the module outputs in the ouput file ..." << endl;
    cout << " - Closing the output file now..." << endl;
//  fout->Write();

  //fout->Close();

  //delete
  delete fout;
  delete event;
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;

  delete trainer;


  cout << "**********************************************************************" << endl;
  cout << "           End of the program !!" << endl;
  cout << " 		hasn't crashed yet ;-) " << endl;
  cout << "**********************************************************************" << endl;

  return 0;
}
