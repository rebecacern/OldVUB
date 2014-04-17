
///////////////////////////
///// TODO & COMMENTS /////
/////////////////////////// 

#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/JetTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../MCInformation/interface/MCWeighter.h"
#include "../Selection/interface/ElectronPlotter.h"
#include "../Selection/interface/MuonPlotter.h"
#include "../Selection/interface/JetPlotter.h"
#include "../Selection/interface/VertexPlotter.h"
#include "../Tools/interface/MVATrainer.h"
#include "../Tools/interface/MVAComputer.h"
#include "../MCInformation/interface/ResolutionFit.h"
#include "../JESMeasurement/interface/FullKinFit.h"
#include "../JESMeasurement/interface/JetCombiner.h"
#include "../Reconstruction/interface/Observables.h"
#include "../Reconstruction/interface/PlotObservables.h"
#include "../Reconstruction/interface/TTreeObservables.h"
#include "Style.C"

using namespace std;
using namespace TopTree;

Float_t ChiSquareMatching (TRootMuon* &, std::vector < TRootJet* > &, UInt_t &, UInt_t *, Float_t &, Float_t &, Float_t &, Float_t &, Float_t &, Float_t &, Bool_t &);

int main(int argc, char *argv[]) {

  clock_t start = clock();
  Bool_t DEBUG = (argc > 2 ? (atoi (argv[2]) == 1) : 0);
  cout << "************************************************************************" << endl;
  cout << " Running the Treemaker based on the Lepton+Jets RefSel for the SyncEx ! " << endl;
  cout << "************************************************************************" << endl;
  
  //SetStyle if needed
  //setTDRStyle(); 
  //setMyStyle(); //thick lines for the histogram plots...

  // Jet Chi2 combination
  UInt_t Permutation[4] = { 0, 0, 0, 0 };
  Float_t WMass = 90.0;		// Derived from the SanityChecker (Truth module)
  Float_t TopMass = 174.5;
  Float_t blMass = 91.1;
  //Float_t LepTopMassResol = 15.6;
  Float_t LepblMassResol = 36.0;
  Float_t HadTopMassResol = 21.1;
  Float_t HadWMassResol = 12.7;
  UInt_t NjetsForComb = 7;
  Float_t ChiSq = 0;


  /////////////////////
  // Configuration
  /////////////////////

    unsigned int lepton = 0; // 0 = muon+jets 1= electron +jets //note: I deleted the block(s) code for electron+jets...
    bool doPF2PAT = false;

  //xml file
  //string xmlfile ="../config/myRefSelconfig.xml";
  string xmlfile = "../config/myNPconfig.xml";

  //Output ROOT file
  string rootFileName ("Esel.root");

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
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
  
  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile.c_str());
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);

  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
  float Luminosity = oldLuminosity;
  for (unsigned int d = 0; d < datasets.size (); d++)
    {
      //cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
      if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
      string dataSetName = datasets[d]->Name();
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	{
	  Luminosity = datasets[d]->EquivalentLumi();
	}
    }
  if(Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  
  
  Observables obs;
  vector <string > lstVar;
  lstVar = obs.ListOfVariables ();

  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootJet* > init_jets_corrected;
  vector < TRootMET* > mets; //maybe mets_corrected has to be made in the future as well!!
  vector < TRootMET* > mets_corrected;
  vector < TRootJet* > CombinedJets;

  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  //Global variable
  TRootEvent* event = 0;
  
  ////////////////////////////////////
  /// Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////

  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;


   ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTable;

  if (lepton == 0) {
    cout << " - Preparing the cutflow table for mu+jets refSelV4" << endl;
    CutsSelecTable.push_back(string("initial"));
    CutsSelecTable.push_back(string("trigged"));
    CutsSelecTable.push_back(string("Good PV"));
    CutsSelecTable.push_back(string("1 selected muon"));
    CutsSelecTable.push_back(string("Veto 2nd muon"));
    CutsSelecTable.push_back(string("Veto electron"));
    char LabelNJets[100];
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
    CutsSelecTable.push_back(string(LabelNJets));
  } else if (lepton == 1){
    cout << " - Preparing the cutflow table for e+jets refSelV4" << endl;
    CutsSelecTable.push_back(string("initial"));
    CutsSelecTable.push_back(string("trigged"));
    CutsSelecTable.push_back(string("Good PV"));
    CutsSelecTable.push_back(string("1 isolated electron"));
    CutsSelecTable.push_back(string("Loose muon veto"));
    CutsSelecTable.push_back(string("Z veto"));
    CutsSelecTable.push_back(string("Conversion rejection"));
    CutsSelecTable.push_back(string("Partnertrack veto"));
    char LabelNJets[100];
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
    CutsSelecTable.push_back(string(LabelNJets));

  } else {

    cout << "!!!!Lepton should equal 0 or 1, exiting" << endl;

    exit(1);

  }

  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ... size: " << CutsSelecTable.size() <<endl;
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity); //was commented...?
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;


  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
		
		//open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);

    selecTable.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
    
    string dataSetName = datasets[d]->Name();
    string dataSetTitle = datasets[d]->Title (); //now placed here


    //piece from TTbarJes.cc
    /////////////////////////////////////
    /// Initialize JEC factors
    /////////////////////////////////////
   	
    vector<JetCorrectorParameters> vCorrParam;
    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
    {
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L2Relative.txt");
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L3Absolute.txt");
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L2L3Residual.txt");
      vCorrParam.push_back(*L2JetCorPar);
      vCorrParam.push_back(*L3JetCorPar);
      vCorrParam.push_back(*ResJetCorPar);
    }
    else // MC!
    {
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L2Relative.txt");
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L3Absolute.txt");
      vCorrParam.push_back(*L2JetCorPar);
      vCorrParam.push_back(*L3JetCorPar);
    }
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/START38_V14_AK5PF_Uncertainty.txt");    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);


    vector<string> BranchesTtree;
    BranchesTtree = obs.ListOfVariables ();
    BranchesTtree.push_back("scaleFactor"); //not an observable...
    //cout<<"BranchesTtree.size() = "<<BranchesTtree.size()<<endl;
    //for(unsigned int i=0;i<BranchesTtree.size();i++) cout<<" "<<BranchesTtree[i]<<endl;
    TTreeObservables TtreeObs(BranchesTtree,dataSetTitle);
    
    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    //nEvents[d] = 0;
    int itrigger = -1, previousRun = -1;
    string previousFilename("");
    float minAllJetsPtMET = 9999;
    int iFile = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
//    for (unsigned int ievt = 0; ievt < 50000; ievt++)
    {

      selecTable.Fill(d,0,1.);

      
      //nEvents[d]++;

      //cout << "anaEnv.JetType -> " << anaEnv.JetType << " --- " <<  endl;
      //cout << "anaEnv.METType -> " << anaEnv.METType << " --- " <<  endl;


      if(ievt%500 == 0)
        std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
        
      //load event

      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      //if(!(event->runId()==1 && event->eventId()==1233785 && event->lumiBlockId()==3)) continue; 
      
      //cout<<"For jets of the event before all corrections (i.e. straight from the toptree)"<<endl;
      //for(unsigned int i=0; i<init_jets.size(); i++){
      //   cout<<"init_jets["<<i<<"]->Pt() = "<<init_jets[i]->Pt()<<endl;
      //} 
          
      //long piece from TTbarJes.cc
      vector<TRootGenJet*> genjets;
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {
        genjets = treeLoader.LoadGenJet(ievt);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      // scale factor for the event
      float scaleFactor = 1.;
      
      // Load the GenEvent and calculate the branching ratio correction
      //if(dataSetName.find("TTbarJets") == 0)
      if(dataSetName == "TTJets")
      {
        TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt);
        if( genEvt->isSemiLeptonic() )
	   scaleFactor *= (0.108*9.)*(0.676*1.5);
	else if( genEvt->isFullHadronic() )
	   scaleFactor *= (0.676*1.5)*(0.676*1.5);
	else if( genEvt->isFullLeptonic() )
  	   scaleFactor *= (0.108*9.)*(0.108*9.);
      }
      
      // Clone the init_jets vector, otherwise the corrections will be removed
      for(unsigned int i=0; i<init_jets_corrected.size(); i++)
        if(init_jets_corrected[i]) delete init_jets_corrected[i];
      init_jets_corrected.clear();

      for(unsigned int i=0; i<init_jets.size(); i++)
        init_jets_corrected.push_back( (TRootJet*) init_jets[i]->Clone() );
	
	
      // Clone the mets vector, otherwise the corrections will be removed //added this w.r.t. changes of Stijn
      for(unsigned int i=0; i<mets_corrected.size(); i++)
        if(mets_corrected[i]) delete mets_corrected[i];
      mets_corrected.clear();

      for(unsigned int i=0; i<mets.size(); i++)
        mets_corrected.push_back( (TRootMET*) mets[i]->Clone() );

      /*//testing
      cout<<"new event"<<endl;
      cout<<"  mets[0]->Px() = "<<mets[0]->Px()<<endl;
      cout<<"  mets[0]->Py() = "<<mets[0]->Py()<<endl;
      cout<<"  mets[0]->Pz() = "<<mets[0]->Pz()<<endl;
      cout<<"  mets[0]->E() = "<<mets[0]->E()<<endl;
      cout<<"  mets[0]->Et() = "<<mets[0]->Et()<<endl;
      cout<<" now setting Px = 1 and Py = 3"<<endl;
      mets[0]->SetPx(1.);
      mets[0]->SetPy(3.);
      cout<<" new values"<<endl;
      cout<<"  mets[0]->Px() = "<<mets[0]->Px()<<endl;
      cout<<"  mets[0]->Py() = "<<mets[0]->Py()<<endl;
      cout<<"  mets[0]->Pz() = "<<mets[0]->Pz()<<endl;
      cout<<"  mets[0]->E() = "<<mets[0]->E()<<endl;
      cout<<"  mets[0]->Et() = "<<mets[0]->Et()<<endl;
      cout<<" now setting E = sqrt(Px^2 + Py^2)"<<endl;
      mets[0]->SetE(sqrt(pow(mets[0]->Px(),2) + pow(mets[0]->Py(),2)));
      cout<<"  mets[0]->E() = "<<mets[0]->E()<<endl;
      cout<<"  mets[0]->Et() = "<<mets[0]->Et()<<endl;*/
      
      // Correct the JetId fractions for PFJets
      if (init_jets_corrected.size() > 0)
      {
        if(init_jets_corrected[0]->jetType() == 2 || doPF2PAT)
        {
          for(unsigned int j=0; j<init_jets_corrected.size(); j++)
          {
            TRootPFJet* pfJet = jetTools->convertToPFJets(init_jets_corrected[j]);
            pfJet->setChargedHadronEnergyFraction(pfJet->chargedHadronEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
            pfJet->setNeutralHadronEnergyFraction(pfJet->neutralHadronEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
            pfJet->setChargedEmEnergyFraction(pfJet->chargedEmEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
            pfJet->setChargedMuEnergyFraction(pfJet->chargedMuEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
            pfJet->setNeutralEmEnergyFraction(pfJet->neutralEmEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
          }
        }
      }

      /////////////
      // TRIGGER //
      /////////////
      
      bool trigged=false;
      
      if (lepton == 0) {
	//if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
	
/*	int currentRun = event->runId();
	//if(previousRun != currentRun) {
	//  previousRun = currentRun;
	//cout << currentRun << endl;
	itrigger = treeLoader.iTrigger ("HLT_Mu9", currentRun);
	//}
*/	
	string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
	if(previousFilename != currentFilename){
	     previousFilename = currentFilename;
	     iFile++;
	     cout<<"   iFile = "<<iFile<<endl;
	}
	  
	int currentRun = event->runId();
	  //cout<<", currentRun "<<currentRun<<endl;	  
      	if(previousRun != currentRun){	  
	     cout<<" currentRun "<<currentRun<<endl;
             previousRun = currentRun;
	     if(currentRun<147196){
	       cout<<"HLT_Mu9 should be used"<<endl;
	       itrigger = treeLoader.iTrigger (string ("HLT_Mu9"), currentRun, iFile);
	       cout<<" itrigger "<<itrigger<<endl;
	       }
	     else if(currentRun>=147196){
               cout<<"HLT_Mu15_v1 should be used"<<endl;
	       itrigger = treeLoader.iTrigger (string ("HLT_Mu15_v1"), currentRun, iFile); // HLT_Mu9 is prescaled from a certain run in data...
	       cout<<" itrigger "<<itrigger<<endl;
	     }
	}
	trigged = treeLoader.EventTrigged (itrigger);
	
	//} else 
	//trigged=true;
      }


      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
      {
        // Apply the scraping veto
        bool isBeamBG = true;
        if(event->nTracks() > 10)
        {
          if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
            isBeamBG = false;
        }
        if(isBeamBG) continue;
        
        // Apply the JSON
        // bool passedJSON = treeLoader.EventPassedJSON(datasets[d], event->runId(), event->lumiBlockId());
        // if(!passedJSON) continue;
      }

      // Apply Jet Corrections on-the-fly ((z) energy scale)
      jetTools->correctJets(init_jets_corrected);      
      // Apply Jet Corrections on the MET (z), is this correct?? probably... btw, init_jets (wrongly corrected) and init_jets_corrected do not change very much... in the 5th digit sometimes... suspicious??
      //cout<<"  With 'old' JES correction, but the original TRootMET: mets[0]->Et() = "<<mets[0]->Et()<<endl;
      //cout<<"  With 'old' JES correction: mets_corrected[0]->Et() = "<<mets_corrected[0]->Et()<<endl;
      for(unsigned int i=0; i<init_jets_corrected.size(); i++){
         //cout<<"   init_jets["<<i<<"]->Pt() = "<<init_jets[i]->Pt()<<", init_jets_corrected["<<i<<"]->Pt() = "<<init_jets_corrected[i]->Pt()<<endl;
         mets_corrected[0]->SetPx(mets_corrected[0]->Px() + init_jets[i]->Px() - init_jets_corrected[i]->Px());
         mets_corrected[0]->SetPy(mets_corrected[0]->Py() + init_jets[i]->Py() - init_jets_corrected[i]->Py());
	 mets_corrected[0]->SetE(sqrt(pow(mets_corrected[0]->Px(),2) + pow(mets_corrected[0]->Py(),2)));
      }      
      //cout<<"  With 'new' JES correction: mets_corrected[0]->Et() = "<<mets_corrected[0]->Et()<<endl;
      
      
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {
        // Correct for the difference in muon efficiency (HLT and Id) between Data and MC
        scaleFactor *= 0.965;
        
        // Correct jets and MET for JES uncertainy systematics
        //jetTools->correctJetJESUnc(init_jets_corrected, mets_corrected[0],"minus"); //other direction: "plus"
        //cout<<endl;
	// Match RecoJets with GenJets
        vector< pair<size_t, size_t> > indexVector; //first index = RecoJet, second index = GenJet
        vector<bool> mLock(genjets.size(),false);   // when locked, genJet is already matched to a recoJet
        for(size_t i=0; i<init_jets_corrected.size(); i++)
        {
          pair<size_t, size_t> tmpIndex;
          float minDR = 9999.;
          for(size_t j=0; j<genjets.size(); j++)
          {
            if( ! mLock[j] )
            {
              if( init_jets_corrected[i]->DeltaR(*genjets[j]) < 0.4 && init_jets_corrected[i]->DeltaR(*genjets[j]) < minDR )
              {
                minDR = init_jets_corrected[i]->DeltaR(*genjets[j]);
                tmpIndex = pair<size_t, size_t>(i,j);
              }
            }
          }
          if(minDR < 999.)
          {
            mLock[tmpIndex.second] = true;
            indexVector.push_back(tmpIndex);
          }
        }

        // Apply correction for jet energy resolution on-the-fly, only for recoJets matched with a genJet
        for(size_t i=0; i<indexVector.size(); i++)
        {
          if( genjets[indexVector[i].second]->Pt() < 15 ) continue;
	  //cout<<" mets_corrected[0]->Px() iter "<<i<<" (BEFORE new correction) = "<<mets_corrected[0]->Px()<<endl;
          //wait!! this is tricky... init_jets is not 'uncorrected', but rather 'wronly corrected'... twiki is VERY strange, what do they mean with 'uncorrected' jets? I use 'corrected' jets...
	  mets_corrected[0]->SetPx(mets_corrected[0]->Px() + init_jets_corrected[indexVector[i].first]->Px());
	  mets_corrected[0]->SetPy(mets_corrected[0]->Py() + init_jets_corrected[indexVector[i].first]->Py());
	  float corrFactor = 0.1; // factor is either 0.1 for bias correction, 0.0 for JER_minus and 0.2 for JER_plus
          float deltapt = ( init_jets_corrected[indexVector[i].first]->Pt() - genjets[indexVector[i].second]->Pt() ) * corrFactor;
          float ptscale = max(0.0, ( init_jets_corrected[indexVector[i].first]->Pt() + deltapt) / init_jets_corrected[indexVector[i].first]->Pt() );
          if(ptscale > 0.0){
            init_jets_corrected[indexVector[i].first]->SetPxPyPzE(init_jets_corrected[indexVector[i].first]->Px()*ptscale, init_jets_corrected[indexVector[i].first]->Py()*ptscale,
	       init_jets_corrected[indexVector[i].first]->Pz()*ptscale, init_jets_corrected[indexVector[i].first]->E()*ptscale);
        //    init_jets[indexVector[i].first]->SetPxPyPzE(init_jets[indexVector[i].first]->Px()*ptscale, init_jets[indexVector[i].first]->Py()*ptscale,
	//       init_jets[indexVector[i].first]->Pz()*ptscale, init_jets[indexVector[i].first]->E()*ptscale);
	  }
	  mets_corrected[0]->SetPx(mets_corrected[0]->Px() - init_jets_corrected[indexVector[i].first]->Px());
	  mets_corrected[0]->SetPy(mets_corrected[0]->Py() - init_jets_corrected[indexVector[i].first]->Py());
	  mets_corrected[0]->SetE(sqrt(pow(mets_corrected[0]->Px(),2) + pow(mets_corrected[0]->Py(),2)));
	  //cout<<" mets_corrected[0]->Px() iter "<<i<<" (AFTER new correction) = "<<mets_corrected[0]->Px()<<endl;

	}
        
        //Scale jets with a certain factor
        //jetTools->scaleJets(init_jets_corrected, 1.);
      }

      //cout<<"For jets of the event after all corrections, just before selection"<<endl;
      //for(unsigned int i=0; i<6; i++){ //init_jets_corrected.size()
      //   TRootPFJet* pfJet = jetTools->convertToPFJets(init_jets_corrected[i]);
      //   cout<<" pfJet "<<i<<":"<<endl;
	 //cout<<"   Pt() = "<<pfJet->Pt()<<endl;
	 //cout<<"   Eta() = "<<pfJet->Eta()<<endl;
	 //cout<<"   nOfDaughters() = "<<pfJet->nOfDaughters()<<endl;
	 //cout<<"   chargedEmEnergyFraction() = "<<pfJet->chargedEmEnergyFraction()<<endl;
	 //cout<<"   neutralHadronEnergyFraction() = "<<pfJet->neutralHadronEnergyFraction()<<endl;
	 //cout<<"   neutralEmEnergyFraction() = "<<pfJet->neutralEmEnergyFraction()<<endl;
	 //cout<<"   chargedHadronEnergyFraction() = "<<pfJet->chargedHadronEnergyFraction()<<endl;
	 //cout<<"   chargedMultiplicity()() = "<<pfJet->chargedMultiplicity()<<endl;
      //}
      /////////////////////////////
      //   Selection
      /////////////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets_corrected); //mets also corrected...
      
      bool isGoodPV = isGoodPV = selection.isPVSelected(vertex, 4,24,2.);  

      //*******************************//
      //***** Muon+Jets RefSel V4 *****//
      //*******************************//

      if (lepton == 0) {
    
	// FILL THE SELECTION TABLE //
	
	if(trigged){
	  selecTable.Fill(d,1,1.);
	  
	  if(isGoodPV){
	    selecTable.Fill(d,2,1.);
	    
	    vector<TRootJet*> selectedJets;
	    vector<TRootMuon*> selectedMuons;
	   
	    bool doPF2PAT = false; // not supported/synced atm...

	    if (init_jets_corrected.size() > 0) {
	      if (init_jets_corrected[0]->jetType() == 1 || doPF2PAT) { // calojets
		//cout << "Selecting for caloJets" << endl;
		selectedJets = selection.GetSelectedJets(true);
		selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	      }	else {
		//cout << "Selecting for PF/JPT jets" << endl;
		vector<TRootMuon*> overlapMuons = selection.GetSelectedMuons(vertex[0]);
		//selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1); // refSelV4 values
		selectedJets = selection.GetSelectedJets(overlapMuons,true);
		selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);	
	      }
	    }

	    vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
	    vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(false);

	    if(selectedMuons.size()==1){
	      selecTable.Fill(d,3,1.);
	      
	      if(vetoMuons.size()==1){
		selecTable.Fill(d,4,1.);
		
		if(vetoElectrons.size()==0) {
		  selecTable.Fill(d,5,1.);			    
		  
		  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-3) {
		    selecTable.Fill(d,6,1.);
		    
		    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2) {
		      selecTable.Fill(d,7,1.);
		      
		      if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1)
			{
			  selecTable.Fill(d,8,1.);
			  
			  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets)
			    {
			      selecTable.Fill(d,9,1.);

			      // EVENT IS SELECTED!!!!
//////			      cout<<event->runId()<<":"<<event->eventId()<<":"<<event->lumiBlockId()<<endl;
 			      //////////////////////
	                      //LOOP OVER VARIABLES
	  		      for (unsigned int v = 0; v < lstVar.size (); v++) {	     			
			//	if (!anaEnv.runOnObs ((int) v))
	        	//	  continue;
  
	     			float variable = -999;
	     			string variable_string;	        
	        		ChiSq = ChiSquareMatching (selectedMuons[0], selectedJets, NjetsForComb, Permutation, WMass, TopMass, blMass, HadWMassResol, HadTopMassResol, LepblMassResol, DEBUG);
	        		CombinedJets.clear ();
	   
 	        		if (selectedJets.size () >= 4)
	          		  for (Int_t i = 0; i < 4; i++) {
		     		    CombinedJets.push_back (selectedJets[Permutation[i]]);
	          		  }
	        		else
	          		  CombinedJets = selectedJets;
	        	        string a ("isInSR");
		
	        		//Observables obsEvt (selectedMuons[0], CombinedJets, selectedJets, mets[0], a, ChiSq);
		
				//now other arguments needed... this may be a temporary workaround in stead of adapting the whole observables class...
				vector<TRootJet> CombinedJets_nopointers;
				vector<TRootJet> selectedJets_nopointers;	
				for(unsigned int mmm=0;mmm<CombinedJets.size();mmm++){
				  CombinedJets_nopointers.push_back(*(CombinedJets[mmm]));
				}
				for(unsigned int mmm=0;mmm<CombinedJets.size();mmm++){
				  selectedJets_nopointers.push_back(*(selectedJets[mmm]));
				}	
				Observables obsEvt (*(selectedMuons[0]), CombinedJets_nopointers, selectedJets_nopointers, *(mets_corrected[0]), a, ChiSq);
				//the following is to set the scaleFactor of the event, to be able to store this in a branch in the end... the scaleFactor is not an observable of course, so in the future, a lot of care has to be taken w.r.t. e.g. the ObservablesRanker
				//cout<<"scaleFactor event = "<<scaleFactor<<endl;
				obsEvt.setEventscaleFactor(scaleFactor); //for data, scaleFactor should be 1
	        		
				//float thisAllJetsPtMET = obsEvt.Variable("AllJetsPtMET");
				//if(thisAllJetsPtMET<minAllJetsPtMET){				
				//  minAllJetsPtMET = thisAllJetsPtMET;
				//  cout<<"New minimum AllJetsPtMET = "<<minAllJetsPtMET<<endl;
				//}
				/*for (unsigned int kkk=0;kkk<obsEvt.Variables().size();kkk++){	  
      		  		  if ( lstVar[v]== obsEvt.Variables()[kkk].first) { 	     
		     			//cout<<"FOUND THE CORRECT ONEEE!!!!!!!!!    "<<lstVar[v]<< "------------->>>>>>>>>>>>>>>>  "<< obsEvt.Variables ()[kkk].second<<"  " <<obsEvt.Variables ()[kkk].first<<"   "<<kkk<<endl;
		     		    Index_=kkk;
		  		  }
	  			}*/
	        		//Index_=v;
				/*
				if (anaEnv.MCRound==0)
		  		  BinContentObsAll_.push_back( pair< string,float > (lstVar[v],obsEvt.Variables ()[Index_].second));
		  
	        		BinContentObsPerSample_.push_back( pair< string,float > (lstVar[v],obsEvt.Variables ()[Index_].second));
			////        if(lstVar[v]=="MET"){}
		             */
				bool fill;fill=0;
	        		if (v==lstVar.size()-1) 
		  		  fill=1;
                		TtreeObs.Fill (obsEvt,fill);
                		/*Plots.Fill(obsEvt,fillweight,fill);
	        		variable = obsEvt.Variables ()[Index_].second;
                		variable_string = lstVar[v];
	        		if (anaEnv.nPseudoExp!=0) {
	    	  		  container.Add(variable);
	  			}	

	        		ValuesInSR.push_back (variable);	    
	        		ValuesInSR_WithString.push_back (pair < string, float > (variable_string,variable));*/
	         	    /*
  	        		string dstitle=datasets[d]->Title();
				if(anaEnv.isMC == 0){
		   		  hRelIsoData->Fill(RelIso);
				}
				if(anaEnv.isMC == 1){
		   		  myMultiSamplePlot.Fill(RelIso,datasets[d],true,Luminosity);
				}	
    			    */
	  		      }			//loop over variables
 
 
			    }
			}
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      
    } // event loop
    
    if( anaEnv.nPseudoExp==0 || dataSetTitle=="Data")
      TtreeObs.Write (true); ////////////////write all tree OBS for each dataset
    
  } // dataset loop

  //Selection tables
  selecTable.TableCalculator(false, true, true, true);
  string selectiontable = "SelectionTable_SyncEx";
  selectiontable = selectiontable +".tex"; 	
  selecTable.Write(selectiontable.c_str());


  //TFile* fout = new TFile("out.root","RECREATE");

  fout->cd();

  TDirectory* th1dir = fout->mkdir("1D_histograms");

  th1dir->cd();

  fout->cd();

  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++) {
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    //temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    //temp->SetBinContent(N+1,0);
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( ("plots/"+it->first+".png").c_str() );
  }
  
  fout->Close();

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

}

Float_t
ChiSquareMatching (TRootMuon* & muon, std::vector < TRootJet* > &jets, UInt_t & n, UInt_t * Permutation, Float_t & WMass_, Float_t & TopMass_, Float_t & blMass_, Float_t & HadWMassResol_, Float_t & HadTopMassResol_, Float_t & LepblMassResol_,
		   Bool_t & DEBUG)
{

  if (jets.size () < 4)
    return -1;
  UInt_t NbOfJets = 0;
  (jets.size () < n ? NbOfJets = jets.size () : NbOfJets = n);
  UInt_t *numbers = new UInt_t[NbOfJets];
  for (UInt_t i = 0; i < NbOfJets; i++)
    numbers[i] = i;
  UInt_t *comb = new UInt_t[4];
  for (UInt_t i = 0; i < 4; i++)
    comb[i] = i;
  UInt_t NbOfPermutations = 0;
  Float_t ChiSquare = 0, ChiSquare_tmp = 0;
  Float_t HadWChiSquare = 0, HadTopChiSquare = 0, LepblChiSquare = 0;
  Float_t HadWCandMass_tmp = -1, HadTopCandMass_tmp = -1, LepBlCandMass_tmp = -1;
  Float_t HadWCandMass = -1, HadTopCandMass = -1, LepBlCandMass = -1;
  //DEBUG = true;
  ChiSquare = 9999999;
  do {
    if (DEBUG) {
      std::cout << "-- Enter loop for jet combination " << std::endl;
      std::cout << "-- Considered Jet combination : " << comb[0] << "/" << comb[1] << "/" << comb[2] << "/" << comb[3] << std::endl;
    }
    NbOfPermutations = 0;
    do {
      NbOfPermutations++;
      if (DEBUG)
	std::cout << "--- Current permutations : " << comb[0] << "/" << comb[1] << "/" << comb[2] << "/" << comb[3] << std::endl;
      HadWCandMass_tmp = (*jets[comb[0]] + *jets[comb[1]]).M ();
      HadTopCandMass_tmp = (*jets[comb[0]] + *jets[comb[1]] + *jets[comb[2]]).M ();
      LepBlCandMass_tmp = (*jets[comb[3]] + *muon).M ();

      HadWChiSquare = pow ((HadWCandMass_tmp - WMass_), 2) / pow (HadWMassResol_, 2);
      HadTopChiSquare = pow ((HadTopCandMass_tmp - TopMass_), 2) / pow (HadTopMassResol_, 2);
      LepblChiSquare = pow ((LepBlCandMass_tmp - blMass_), 2) / pow (LepblMassResol_, 2);

      ChiSquare_tmp = HadWChiSquare + HadTopChiSquare + LepblChiSquare;
      if (ChiSquare_tmp < ChiSquare) {
	//for(unsigned int j=0;j<4;j++) Permutation_tmp[i][j] = comb[j];
	//ChiSquare[i] = ChiSquare_tmp[i];
	for (UInt_t i = 0; i < 4; i++)
	  Permutation[i] = comb[i];
	ChiSquare = ChiSquare_tmp;
	HadWCandMass = HadWCandMass_tmp;
	HadTopCandMass = HadTopCandMass_tmp;
	LepBlCandMass = LepBlCandMass_tmp;
      }
    }
    while (next_permutation (comb, comb + 4));
    if (DEBUG) {
      std::cout << "--- ChiSquare matching : Nb of permutations found = " << NbOfPermutations << std::endl;
      std::cout << "--- ChiSquare matching : Min ChiSquare found = " << ChiSquare << std::endl;
      std::cout << "--- ChiSquare matching : associated to the permutation : " << Permutation[0] << "/" << Permutation[1] << "/" << Permutation[2] << "/" << Permutation[3] << std::endl;
    }
  }
  while (stdcomb::next_combination (numbers, numbers + NbOfJets, comb, comb + 4));
  if (DEBUG) {
    std::cout << "ChiSquare = " << ChiSquare << std::endl;
    std::cout << "Associated had W mass   = " << HadWCandMass << std::endl;
    std::cout << "Associated had top mass = " << HadTopCandMass << std::endl;
    std::cout << "Associated lep b+l mass = " << LepBlCandMass << std::endl;
  }
  return ChiSquare;
  delete comb;
  delete numbers;
}
 
