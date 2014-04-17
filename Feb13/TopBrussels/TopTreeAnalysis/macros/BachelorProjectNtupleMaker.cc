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
#include "../Tools/interface/JetTools.h"
#include "../MCInformation/interface/ResolutionFit.h"
#include "../MCInformation/interface/JetPartonMatching.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../BtagEffAnalysis/interface/TRootNTuple.h"
#include "Style.C"

using namespace std;
using namespace TopTree;

int main (int argc, char *argv[])
{
	
	bool doPF2PAT = false;
	
	
	clock_t start = clock();
	
	cout << "***********************************************************************" << endl;
	cout << " Beginning of the program for creating the Bachelor Projects NTuples ! " << endl;
	cout << "***********************************************************************" << endl;
	
	//SetStyle if needed
	//setTDRStyle(); 
	setMyStyle();
	
	/////////////////////
	// Configuration
	/////////////////////
	
	//xml file
	string xmlFileName ="../config/myBTAGconfig.xml";
	//string xmlFileName ="../config/myRefSelconfig.xml";
	if (argc >= 2)
		xmlFileName=string(argv[1]);
	const char *xmlfile = xmlFileName.c_str();
	
	cout << "used config file: " << xmlfile << endl;
	
	//Output ROOT file
	string rootFileName ("BtaggingOutput.root");
	
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
	float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
	
 	cout << "analysis environment luminosity for rescaling "<< oldLuminosity << endl;
	
	/////////////////////
	// Load Datasets
	/////////////////////
	
	TTreeLoader treeLoader;
	cout << " - Load datasets ..." << endl;
	vector < Dataset* > datasets;
	treeLoader.LoadDatasets (datasets, xmlfile);
	for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
	
	float Luminosity = oldLuminosity;
	bool isSemiMu = false;
	for (unsigned int d = 0; d < datasets.size (); d++)
    {
		//cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
		if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
		string dataSetName = datasets[d]->Name();
		if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
		{
			Luminosity = datasets[d]->EquivalentLumi();
			break;
		}
		if(dataSetName.find("TTbarJets_SemiMu") == 0) isSemiMu = true;
    }
	if(Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
	
	//vector of objects
	cout << " - Variable declaration ..." << endl;
	vector < TRootVertex* > vertex;
	vector < TRootMuon* > init_muons;
	vector < TRootElectron* > init_electrons;
	vector < TRootJet* > init_jets;
	vector < TRootMET* > mets;
	
	TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
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
	for (unsigned int d = 0; d < datasets.size (); d++) {
		
		int nSelected=0;
		if (verbose > 1)
			cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
		if (verbose > 1)
			std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
		
		//open files and load
		cout<<"LoadEvent"<<endl;
		treeLoader.LoadDataset (datasets[d], anaEnv);
		cout<<"LoadEvent"<<endl;
		
		string previousFilename = "";
		int iFile = -1;
		string dataSetName = datasets[d]->Name();
		
		/////////////////////////////////////
		/// Initialize JEC factors
		/////////////////////////////////////
		
		vector<JetCorrectorParameters> vCorrParam;
		if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
		{
			
			JetCorrectorParameters *L2JetCorPar;
			JetCorrectorParameters *L3JetCorPar;
			JetCorrectorParameters *ResJetCorPar;
			
			if(anaEnv.JetType == 1) {
				cout << "Loading Jet Correction factors for CaloJets" << endl;
				L2JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5Calo_L2Relative.txt");
				L3JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5Calo_L3Absolute.txt");
				ResJetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5Calo_L2L3Residual.txt");
			}
			
			if(anaEnv.JetType == 2) {
				cout << "Loading Jet Correction factors for PFJets" << endl;
				L2JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L2Relative.txt");
				L3JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L3Absolute.txt");
				ResJetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L2L3Residual.txt");
			}
			
			vCorrParam.push_back(*L2JetCorPar);
			vCorrParam.push_back(*L3JetCorPar);
			vCorrParam.push_back(*ResJetCorPar);
		}
		else // MC!
		{      JetCorrectorParameters *L2JetCorPar;
			JetCorrectorParameters *L3JetCorPar;
			
			if(anaEnv.JetType == 1) {
				cout << "Loading Jet Correction factors for CaloJets" << endl;
				L2JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5Calo_L2Relative.txt");
				L3JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5Calo_L3Absolute.txt");
			}
			
			if(anaEnv.JetType == 2) {
				cout << "Loading Jet Correction factors for PFJets" << endl;
				L2JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L2Relative.txt");
				L3JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L3Absolute.txt");
			}
			
			vCorrParam.push_back(*L2JetCorPar);
			vCorrParam.push_back(*L3JetCorPar);
		}
		
		JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/START38_V14_AK5PF_Uncertainty.txt");
		
		JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
		
		
		////////////////////////////////////////////////////////////
		// CREATE OUTPUT FILE AND TTREE FOR STORAGE OF THE NTUPLE //
		////////////////////////////////////////////////////////////
		
		string dir = "BachelorProjectNTuples/";
		
		mkdir(dir.c_str(),0777);
		
		if (argc > 2)
			dir += string(argv[2])+"/";
		
		mkdir(dir.c_str(),0777);
		
		
		string Title = dir+"NTuple_"+dataSetName+".root";
		
		cout << "INFO: creating  "+Title << endl;
		
		TFile* TreeFile = new TFile(Title.c_str(),"RECREATE");
		
		TTree* Tree= new TTree("tree","Tree containing the Bachelor Students NTuple");
		
		// Create the branches
		
		//-- Special Branches --//
		
		Dataset dataSetProcessed = *datasets[d];
		
		Tree->Branch("Dataset",&dataSetProcessed);    
		
		Int_t eventID, runID, lumiBlockID;
		
		Tree->Branch("eventID",&eventID);    
		Tree->Branch("runID",&runID);    
		Tree->Branch("lumiBlockID",&lumiBlockID);    
		
		Float_t event_scaleFactor = 0;
			
		Tree->Branch("eventWeight",&event_scaleFactor);    
		
		//-- MUON --//
		
		Int_t Muon_Charge = 0;
		Float_t Muon_RelIso = 0;
		TLorentzVector Muon;
		
		Tree->Branch("Muon_Charge",&Muon_Charge);
		Tree->Branch("Muon_RelIso",&Muon_RelIso);
		Tree->Branch("Muon_p4",&Muon);
		
		//-- MET --//
		
		TLorentzVector MET;
		
		Tree->Branch("MET_p4",&MET);

		//-- JETS --//
		
		TClonesArray* jets = new TClonesArray("TLorentzVector",50); 
		
		Tree->Branch("Jets_p4","TClonesArray",&jets);
		
		TArrayF* jets_pdgId = new TArrayF(50);
		TArrayF* jets_mpdgId = new TArrayF(50);
		//TArrayF* jets_pdgId_DeltaR = new TArrayF(50);
		TArrayF* jets_bTag = new TArrayF(50);
		
		Tree->Branch("Jets_pdgID","TArrayF",&jets_pdgId);
		Tree->Branch("Jets_Mother_pdgID","TArrayF",&jets_mpdgId);
		//Tree->Branch("Jets_pdgID_matchDeltaR","TArrayF",&jets_pdgId_DeltaR);

		Tree->Branch("Jets_bTagDiscriminator","TArrayF",&jets_bTag);
		
		////////////////////////////////////
		//	Loop on events
		////////////////////////////////////
		
		nEvents[d] = 0;
		int itrigger = -1, previousRun = -1;
		if (verbose > 1)
			cout << "	Loop over events " << endl;
		
		for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
		//for (unsigned int ievt = 0; ievt < 10000; ievt++)
		{
			nEvents[d]++;
			if(ievt%1000 == 0)
				std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << " +> # selected: " << nSelected << flush<<"\r";
			// std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
			
			//load event
			event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
			
			init_jets = init_jets;
			
			vector<TRootGenJet*> genjets;
			if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
			{
				genjets = treeLoader.LoadGenJet(ievt,false);
				sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
			}
			
			// scale factor for the event
			float scaleFactor = 1.;
			
			// Load the GenEvent and calculate the branching ratio correction
			if(dataSetName.find("TTbarJets") == 0)
			{
				//cout << "LOADING GenEvent" << endl;
				TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
				if( genEvt->isSemiLeptonic() )
					scaleFactor *= (0.108*9.)*(0.676*1.5);
				else if( genEvt->isFullHadronic() )
					scaleFactor *= (0.676*1.5)*(0.676*1.5);
				else if( genEvt->isFullLeptonic() )
					scaleFactor *= (0.108*9.)*(0.108*9.);
			}
			
			// Clone the init_jets vector, otherwise the corrections will be removed
			/*for(unsigned int i=0; i<init_jets.size(); i++)
			 if(init_jets[i]) delete init_jets[i];
			 init_jets.clear();
			 
			 for(unsigned int i=0; i<init_jets.size(); i++)
			 init_jets.push_back( (TRootJet*) init_jets[i]->Clone() );
			 */
			
			// Correct the JetId fractions for PFJets
			if (init_jets.size() > 0)
			{
				if(init_jets[0]->jetType() == 2 || doPF2PAT)
				{
					for(unsigned int j=0; j<init_jets.size(); j++)
					{
						TRootPFJet* pfJet = jetTools->convertToPFJets(init_jets[j]);
						pfJet->setChargedHadronEnergyFraction(pfJet->chargedHadronEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
						pfJet->setNeutralHadronEnergyFraction(pfJet->neutralHadronEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
						pfJet->setChargedEmEnergyFraction(pfJet->chargedEmEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
						pfJet->setChargedMuEnergyFraction(pfJet->chargedMuEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
						pfJet->setNeutralEmEnergyFraction(pfJet->neutralEmEnergyFraction()*pfJet->getJetCorrFactor("L1L2L3"));
					}
				}
			}
			
			// check which file in the dataset it is to have the HLTInfo right
			
			string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
			if(previousFilename != currentFilename){
				previousFilename = currentFilename;
				iFile++;
				cout<<"File changed!!! => iFile = "<<iFile<<endl;
			}
			
			int currentRun = event->runId();
			if(previousRun != currentRun)
			{
				previousRun = currentRun;
				
				if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
					
					//cout << "data" << endl;
					if (event->runId() < 147196)
						itrigger = treeLoader.iTrigger (string ("HLT_Mu9"), currentRun, iFile);
					else if (event->runId() >= 147196)
						itrigger = treeLoader.iTrigger (string ("HLT_Mu15_v1"), currentRun, iFile);
					
				} else {
					
					itrigger = treeLoader.iTrigger (string ("HLT_Mu9"), currentRun);
				}
			}
						
			// Apply Jet Corrections on-the-fly
			jetTools->correctJets(init_jets);
			
			float after = init_jets[0]->Pt();
			
			if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
			{
				//cout << "Do I pass by?" << endl;
				// Correct for the difference in muon efficiency (HLT and Id) between Data and MC
				scaleFactor *= 0.965;
				
				// Correct jets for JES uncertainy systematics
				//        jetTools->correctJetJESUnc(init_jets, "minus");
				
				// Match RecoJets with GenJets
				vector< pair<size_t, size_t> > indexVector; //first index = RecoJet, second index = GenJet
				vector<bool> mLock(genjets.size(),false);   // when locked, genJet is already matched to a recoJet
				for(size_t i=0; i<init_jets.size(); i++)
				{
					pair<size_t, size_t> tmpIndex;
					float minDR = 9999.;
					for(size_t j=0; j<genjets.size(); j++)
					{
						//cout << "Do I smell GenJets?" << endl;
						if( ! mLock[j] )
						{
							if( init_jets[i]->DeltaR(*genjets[j]) < 0.4 && init_jets[i]->DeltaR(*genjets[j]) < minDR )
							{
								minDR = init_jets[i]->DeltaR(*genjets[j]);
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
				
				//cout << indexVector.size() << endl;
				
				// Apply correction for jet energy resolution on-the-fly, only for recoJets matched with a genJet
				for(size_t i=0; i<indexVector.size(); i++)
				{
					//cout << "Smearing the shit out of the crap" << endl;
					if( genjets[indexVector[i].second]->Pt() < 15 ) continue;
					float corrFactor = 0.1; // factor is either 0.1 for bias correction, 0.0 for JER_minus and 0.2 for JER_plus
					float deltapt = ( init_jets[indexVector[i].first]->Pt() - genjets[indexVector[i].second]->Pt() ) * corrFactor;
					float ptscale = max(0.0, ( init_jets[indexVector[i].first]->Pt() + deltapt) / init_jets[indexVector[i].first]->Pt() );
					if(ptscale > 0.0)
						init_jets[indexVector[i].first]->SetPxPyPzE(init_jets[indexVector[i].first]->Px()*ptscale, init_jets[indexVector[i].first]->Py()*ptscale,
																	init_jets[indexVector[i].first]->Pz()*ptscale, init_jets[indexVector[i].first]->E()*ptscale);
				}
				
				//Scale jets with a certain factor
				//jetTools->scaleJets(init_jets, 1.);
			}
			
			//cout << scaleFactor << endl;
			
			/////////////////////////////
			//   Selection
			/////////////////////////////
			
			//Declare selection instance    
			Selection selection(init_jets, init_muons, init_electrons, mets);
			
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
			
			bool trigged = treeLoader.EventTrigged (itrigger);
			bool isGoodPV = false;
			//if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // temporary fix since there is a different cut on Data and MC!
			//  isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, 24, anaEnv.PVertexRhoCut);
			//else
			
			isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);
			
			vector<TRootJet*> selectedJets;
			vector<TRootMuon*> selectedMuons;
			vector<TRootMuon*> selectedLooseMuons = selection.GetSelectedLooseMuons();
			
			if (init_jets.size() > 0) {
				if (init_jets[0]->jetType() == 1 || doPF2PAT) { // calojets
					//cout << "Selecting for caloJets" << endl;
					selectedJets = selection.GetSelectedJets(true);
					selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
				} else {
					//cout << "Selecting for PF/JPT jets" << endl;
					vector<TRootMuon*> overlapMuons = selection.GetSelectedMuons(vertex[0]);
					//selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1); // refSelV4 values
					selectedJets = selection.GetSelectedJets(overlapMuons,true);
					selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);	
				}
			}
			
			vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(false);
			
			vector<TRootJet*> selectedLooseJets = selection.GetSelectedJets(20., anaEnv.JetsEtaCutSR, true);
			vector<TRootMCParticle*> mcParticles;
			
			if(dataSetName.find("TTbarJets_SemiMu") == 0)
			{
				treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
				sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
			}
			
			if( ! (trigged && isGoodPV) ) continue;
			// Don't look further at events with a bad primary vertex or not passing the trigger
			
			if(selectedMuons.size() != 1) continue;
			
			if(selectedJets.size() < 4) continue;
			
			bool eventselected = false;
			if (trigged && isGoodPV)
				if (selectedMuons.size() == 1)
					if (selectedLooseMuons.size() == 1)
						if (vetoElectrons.size() == 0)
							if (selectedJets.size() >= 4)
								eventselected = true;
			
			if (!eventselected) continue;
			
			nSelected++;
			
			//////////////////////////////////////////////////////////////////////////
			// the event is selected so now we will perform the jet-parton matching
			//////////////////////////////////////////////////////////////////////////
			
			//sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)     
			
			
			//vector<TRootMCParticle*> mcParticles;

			if(dataSetName != "Data" && dataSetName != "data" && dataSetName != "DATA") // NOT Data!
			{
				treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  // false to not reload jet coll
				sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
			}
			
			// FILL THE TUPLE
			
			eventID = event->eventId();
			runID = event->runId();
			lumiBlockID = event->lumiBlockId();
			
			event_scaleFactor = scaleFactor;
			
			Muon_Charge = selectedMuons[0]->charge();
			Muon_RelIso = selectedMuons[0]->relativeIso03();
			Muon = *selectedMuons[0];
			
			if (mets.size() > 0)
				MET = *mets[0];
			
			//if (selectedJets.size() > 8) 
			//cout << "Number of selected jets " << selectedJets.size() << endl;
			
			// MATCHING (overkill version)
			
			pair<unsigned int,unsigned int> leptonicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
			pair<unsigned int,unsigned int> hadronicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
			pair<unsigned int,unsigned int> hadronicWJet1_ = pair<unsigned int,unsigned int>(9999,9999);
			pair<unsigned int,unsigned int> hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
			
			vector<TRootMCParticle*> mcParticlesMatching_;
			
			vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
			
			if(dataSetName.find("TTbarJets_SemiMu") == 0)
			{
				vector<TLorentzVector> mcParticlesTLV, selectedLooseJetsTLV;
				
				bool muPlusFromTop = false, muMinusFromTop = false;
				int nTTbarQuarks = 0;
				for(unsigned int i=0; i<mcParticles.size(); i++)
				{
					if( mcParticles[i]->status() != 3) continue;
					
					if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 )
					{
						if(muMinusFromTop) cerr<<"muMinusFromTop was already true"<<endl;
						muMinusFromTop = true;
					}
					if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 )
					{
						if(muPlusFromTop) cerr<<"muPlusFromTop was already true"<<endl;
						muPlusFromTop = true;
					}
					
					if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 )
					{
						mcParticlesTLV.push_back(*mcParticles[i]);
						mcParticlesMatching_.push_back(mcParticles[i]);	
						
					}
				}
				
				if(muPlusFromTop && muMinusFromTop)
					cerr<<"muPlusFromTop and muMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;
				
				// take all the selectedLooseJets_ to study the radiation stuff, selectedLooseJets_ are already ordened in decreasing Pt()
				for(unsigned int i=0; i<selectedJets.size(); i++)
					selectedLooseJetsTLV.push_back(*selectedJets[i]);
				
				JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedLooseJetsTLV, 2, true, true, 0.3);
				
				if(matching.getNumberOfAvailableCombinations() != 1)
					cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
								
				for(unsigned int i=0; i<mcParticlesTLV.size(); i++)
				{
					int matchedJetNumber = matching.getMatchForParton(i, 0);
					if(matchedJetNumber != -1)
						JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
				}
				
				for(unsigned int i=0; i<JetPartonPair.size(); i++)
				{
					unsigned int j = JetPartonPair[i].second;
					
					//      cout<<"JetPartonPair["<<i<<"]: "<<endl;
					//      cout<<"Jet Pt: "<<selectedLooseJets_[JetPartonPair[i].first]->Pt()<<"  Eta: "<<selectedLooseJets_[JetPartonPair[i].first]->Eta()<<"  Phi: "<<selectedLooseJets_[JetPartonPair[i].first]->Phi()<<endl;
					//      cout<<"MCParticle Pt: "<<mcParticlesMatching_[JetPartonPair[i].second]->Pt()<<"  Eta: "<<mcParticlesMatching_[JetPartonPair[i].second]->Eta()<<"  Phi: "<<mcParticlesMatching_[JetPartonPair[i].second]->Phi()<<endl;
					//      cout<<"\ttype: "<<mcParticlesMatching_[JetPartonPair[i].second]->type()<<"  motherType: "<<mcParticlesMatching_[JetPartonPair[i].second]->motherType()<<"  grannyType: "<<mcParticlesMatching_[JetPartonPair[i].second]->grannyType()<<endl;
					//      cout<<"DR(MCParticle, jet) = "<<selectedLooseJets_[JetPartonPair[i].first]->DeltaR(*mcParticlesMatching_[JetPartonPair[i].second])<<endl;
					
					if( fabs(mcParticlesMatching_[j]->type()) < 6 )
					{
						if( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -6 )
						   || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == 6 ) )
						{
							if(hadronicWJet1_.first == 9999) 
								hadronicWJet1_ = JetPartonPair[i];
							else if(hadronicWJet2_.first == 9999) 
								hadronicWJet2_ = JetPartonPair[i];
							else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
						}
					}
					if( fabs(mcParticlesMatching_[j]->type()) == 5 )
					{
						if( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -6 )
						   || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 6 ) )
							hadronicBJet_ = JetPartonPair[i];
						else if( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == 6 )
								|| ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == -6 ) )
							leptonicBJet_ = JetPartonPair[i];
					}
					
				}
				
			} //if dataset Semi mu ttbar
			
			
			for (unsigned int p=0;p<selectedJets.size();p++) {
			
				new ((*jets)[p]) TLorentzVector (*selectedJets[p]);	
				
				Float_t tagger = selectedJets[p]->btag_trackCountingHighEffBJetTags();
				jets_bTag->AddAt(tagger,p);
				
				for(unsigned int i=0; i<JetPartonPair.size(); i++)
				{
				
					unsigned int j = JetPartonPair[i].first;
					
					if (p == j) {
						
						//cout << "Jet " << j << " matched with parton " << JetPartonPair[i].second << " of type " << mcParticlesMatching_[JetPartonPair[i].second]->type() << endl;
						
						jets_pdgId->AddAt(mcParticlesMatching_[JetPartonPair[i].second]->type(),p);
						jets_mpdgId->AddAt(mcParticlesMatching_[JetPartonPair[i].second]->motherType(),p);
						
					}
					
				}
				
			}
			
			Tree->Fill();
			
			jets->Clear();
			jets_pdgId->Reset();
			jets_mpdgId->Reset();
			jets_bTag->Reset();
			
			
			
			//if (ievt > 10)
			//ievt = 10000000;
		}
		
		Tree->Write();
		
		TreeFile->Close();
		
		cout << "!!! Selected " << nSelected << " events " << endl;
	} 			//loop on datasets
	
	//Once everything is filled ...
	if (verbose > 0)
		cout << " We ran over all the data ;-)" << endl;
	
	
	cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
	
	cout << "********************************************" << endl;
	cout << "           End of the program !!            " << endl;
	cout << "           hasn't crashed yet ;-)           " << endl;
	cout << "********************************************" << endl;
	
	return 0;
}
