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
#include "Style.C"

using namespace std;
using namespace TopTree;

int main() {
  
  clock_t start = clock();
  
  
  cout << "*************************************************" << endl;
  cout << " Welcome to the single top tW analysis "  << endl;
  cout << "*************************************************" << endl;
  
  setMyStyle();
  string xmlfile ="../config/tWconfig.xml";
  double lumi = 670;
  
  bool reweightPU = true;
  
  //modes: 0 emu, 1mumu, 2ee 
  int  mode = 0;  
  if (mode != 0 &&  mode !=1 && mode !=2) mode = 0;
  if      (mode == 0) 	cout << " Electron-Muon Mixed channel " << endl;
  else if (mode == 1) 	cout << " Di-Muon channel " << endl;
  else if (mode == 2) 	cout << " Di-Electron channel " << endl;
  cout << "*************************************************" << endl;
  
  int sample = 0; 
  char name[100];
  double xlweight;
  if (sample == 666){		sprintf(name, "data");			xlweight = 1;}
  else if (sample == 0){	sprintf(name, "tt");			xlweight = lumi*157.5/1164194;}
  else if (sample == 1){	sprintf(name, "tw");			xlweight = lumi*10.6/489412;}
  else if (sample == 2){	sprintf(name, "t");			xlweight = lumi*20.93/484030;}
  else if (sample == 3){	sprintf(name, "s");			xlweight = lumi*1.53/494944;}
  else if (sample == 4){	sprintf(name, "wjets"); 		xlweight = lumi*31314/15014623;}
  else if (sample == 5){	sprintf(name, "zjets"); 		xlweight = lumi*3048/2227482;} 
  else if (sample == 6){	sprintf(name, "vqq");			xlweight = lumi*36/738209;}
  else if (sample == 7){	sprintf(name, "ww");			xlweight = lumi*42.9/2061548;}
  else if (sample == 77){	sprintf(name, "ww2l");			xlweight = lumi*4.52/109935;}
  else if (sample == 8){	sprintf(name, "wz");			xlweight = lumi*18.3/2108184;}
  else if (sample == 9){	sprintf(name, "zz");			xlweight = lumi*7.67/2108322;}
  else if (sample == 10){	sprintf(name, "qcd_mu");		xlweight = lumi*84679.3/27113941;} 
  else if (sample == 11){	sprintf(name, "qcd_bc_2030");		xlweight = lumi*132160/1763439;}//
  else if (sample == 12){	sprintf(name, "qcd_bc_3080");  		xlweight = lumi*136804/735502;}//
  else if (sample == 13){	sprintf(name, "qcd_bc_80170"); 	 	xlweight = lumi*9360/899446;}//
  else if (sample == 14){	sprintf(name, "qcd_em_2030"); 	 	xlweight = lumi*2454400/30518227;}//
  else if (sample == 15){	sprintf(name, "qcd_em_3080"); 	 	xlweight = lumi*3866200/71748368;}//
  else if (sample == 16){	sprintf(name, "qcd_em_80170"); 	 	xlweight = lumi*139500/7635061;}//
  else {			sprintf(name, "test"); 	 		xlweight = 1;}
  
  
  char rootFileName[100];
  sprintf(rootFileName,"outputs/out_%d_%s.root", mode, name);
  
  cout << " " << lumi << " pb: " << name << " - " << xlweight << endl;
  
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
  
  AnalysisEnvironment anaEnv;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  
  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile.c_str());
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootMET* > mets;
  
  TFile *fout = new TFile (rootFileName, "RECREATE");
  
  TRootEvent* event = 0;
  
  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;

  
  TH1F* cutflow = new TH1F("cutflow", " ", 31,  -0.5, 30.5 );
  TH1F* cutflow_raw = new TH1F("cutflow_raw", " ", 31,  -0.5, 30.5 );
  
  cutflow->Sumw2();
  cutflow_raw->Sumw2();
  
  double xlWeight; 
  double lum;
  
  double metPt;
  double metPx;
  double metPy;
  
  std::vector<double> *ptLepton;
  std::vector<double> *pxLepton;
  std::vector<double> *pyLepton;
  std::vector<double> *pzLepton;
  std::vector<double> *eLepton;
  std::vector<double> *qLepton;
  
  std::vector<double> *ptJet;
  std::vector<double> *pxJet;
  std::vector<double> *pyJet;
  std::vector<double> *pzJet;
  std::vector<double> *eJet;
  std::vector<double> *qJet;
  std::vector<double> *btTCHPJet;
  std::vector<double> *btTCHEJet;
  std::vector<double> *btSSVHPJet;
  std::vector<double> *btSSVHEJet;
  
  TTree* myTree = new TTree("myTree", "   ");
  
  myTree->Branch("xlWeight", &xlWeight, "xlWeight/D");
  myTree->Branch("lum", &lum, "lum/D");
 
  myTree->Branch("metPt", &metPt, "metPt/D");
  myTree->Branch("metPx", &metPx, "metPx/D");
  myTree->Branch("metPy", &metPy, "metPy/D");
  
  myTree->Branch("ptLepton","std::vector<double>",&ptLepton);
  myTree->Branch("pxLepton","std::vector<double>",&pxLepton);
  myTree->Branch("pyLepton","std::vector<double>",&pyLepton);
  myTree->Branch("pzLepton","std::vector<double>",&pzLepton);
  myTree->Branch("eLepton","std::vector<double>",&eLepton);
  myTree->Branch("qLepton","std::vector<double>",&qLepton);
  
  myTree->Branch("ptJet","std::vector<double>",&ptJet);
  myTree->Branch("pxJet","std::vector<double>",&pxJet);
  myTree->Branch("pyJet","std::vector<double>",&pyJet);
  myTree->Branch("pzJet","std::vector<double>",&pzJet);
  myTree->Branch("eJet","std::vector<double>",&eJet);
  myTree->Branch("qJet","std::vector<double>",&qJet);
  myTree->Branch("btTCHPJet","std::vector<double>",&btTCHPJet);
  myTree->Branch("btTCHEJet","std::vector<double>",&btTCHEJet);
  myTree->Branch("btSSVHPJet","std::vector<double>",&btSSVHPJet);
  myTree->Branch("btSSVHEJet","std::vector<double>",&btSSVHEJet);
  

  MCWeighter* weighter = new MCWeighter();
  string var = "pileup";
  if (reweightPU){
    
    string dir = "ReweightHistos/pileup";
    TFile* fData = new TFile((dir+"/data.root").c_str(),"READ");
    if (fData->IsOpen())
      histo1D["PUReweigting_HistoData"]=(TH1F*) fData->Get(var.c_str())->Clone();
    if (!fData->IsOpen() || !histo1D["PUReweigting_HistoData"]) {
      cerr << "PUReweighting:: Please make shure you have the DATA histogram under " << dir << "/Data.root that contains " << var << endl;
      exit(1);
    } else {
      cout << "PUReweighting:: Loaded " << var << " from " << dir << "/data.root" << endl;
      histo1D["PUReweigting_HistoData"]->SetName("PUReweigting_HistoData");
    }
    
    for (unsigned int d = 0; d < datasets.size (); d++) {
      
      float dataLuminosity = lumi;
      
      float dataSetLuminosity = datasets[d]->EquivalentLumi();
      float dataSetXSection = datasets[d]->Xsection();
      string dataSetName = datasets[d]->Name();
      
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") continue;
      string puFileName = dir+"/"+dataSetName+".root";
      TFile* fMC = new TFile(puFileName.c_str(),"READ");

      if (fMC->IsOpen())
	histo1D["PUReweigting_HistoMC_"+dataSetName]=(TH1F*) fMC->Get(var.c_str())->Clone();
      if (!fMC->IsOpen() || !histo1D["PUReweigting_HistoMC_"+dataSetName]) {
	cerr << "PUReweighting:: Error loading Reweighting file " << puFileName << " or the file does not contain " << var << endl;
	exit(1);
      } else {
	cout << "PUReweighting:: Loaded " << var << " from " << puFileName << endl;
	histo1D["PUReweigting_HistoMC_"+dataSetName]->SetName(("PUReweigting_HistoMC_"+dataSetName).c_str());
	cout << "Lumi: " << dataLuminosity << " XS: " << dataSetXSection << " nExp: " << dataSetXSection*dataLuminosity << endl;
	histo1D["PUReweigting_HistoMC_"+dataSetName]->Scale((dataSetXSection*dataLuminosity)/histo1D["PUReweigting_HistoMC_"+dataSetName]->Integral());
	if (histo1D.find("PUReweigting_HistoMC") == histo1D.end()) {
	histo1D["PUReweigting_HistoMC"] = (TH1F*)histo1D["PUReweigting_HistoMC_"+dataSetName]->Clone();	
	histo1D["PUReweigting_HistoMC"]->SetName("PUReweigting_HistoMC");
	} else {
	  histo1D["PUReweigting_HistoMC"]->Add(histo1D["PUReweigting_HistoMC_"+dataSetName]);
	}
	
      }
      
      float scaleMC=histo1D["PUReweigting_HistoMC"]->Integral();
      float scaleData=histo1D["PUReweigting_HistoData"]->Integral();
      histo1D["PUReweigting_HistoMC"]->Scale(1./scaleMC);
      histo1D["PUReweigting_HistoData"]->Scale(1./scaleData);
      
      weighter->LoadVariable(var, true, histo1D["PUReweigting_HistoMC"], histo1D["PUReweigting_HistoData"], NULL, NULL);
      weighter->Write(new TFile("MCWeighterMonitor.root","RECREATE"));
      
      histo1D["PUReweigting_HistoMC"]->Scale(scaleMC);
      histo1D["PUReweigting_HistoData"]->Scale(scaleData);
    }
  }
  
  
  
  for (unsigned int d = 0; d < datasets.size (); d++)
    {
      treeLoader.LoadDataset (datasets[d], anaEnv);
      string dataSetName = datasets[d]->Name();
      
      cout << datasets[d]->NofEvtsToRunOver() << " total events" << endl;
      for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
	{
	  
	  if(ievt%500 == 0) std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
	  event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
	  
	  double weight = xlweight;
	  
	  if (reweightPU){
	    if(dataSetName != "Data" && dataSetName != "data" && dataSetName != "DATA")
	      weight *=weighter->GetWeight(var,event->nPu(0));
	  }
	  
	  /////////////
	  // TRIGGER //
	  /////////////
	  
	  int currentRun = event->runId();
	  int itrigger = -5;
	  int isecondtrigger = -5;
	  
	  char triggername[100];
	  char triggerbase[100];
	  char secondtriggerbase[100];
	  
	  if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
	    if (mode == 0) sprintf(triggerbase,"HLT_Mu8_Ele17_CaloIdL");
	    else if (mode == 1) sprintf(triggerbase,"HLT_DoubleMu7");
	    else sprintf(triggerbase,"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL");
	    
	    if (mode == 0) sprintf(secondtriggerbase,"HLT_Mu17_Ele8_CaloIdL");
	    
	    for (int i = 1; i < 16; i++){
	      sprintf(triggername,"%s_v%d", triggerbase, i);
	      itrigger = treeLoader.iTrigger (string (triggername), currentRun);
	      if (itrigger != 9999) i = 100;
	      else if (i == 15) cout << "NO VALID TRIGGER FOUND FOR THIS RUN " << event->runId() << endl;
	    }
	    if (mode == 0){
	      for (int i = 1; i < 16; i++){
	        sprintf(triggername,"%s_v%d", secondtriggerbase, i);
	        isecondtrigger = treeLoader.iTrigger (string (triggername), currentRun);
	        if (isecondtrigger != 9999) i = 100;
	        else if (i == 15) cout << "NO VALID TRIGGER FOUND FOR THIS RUN " << event->runId() << endl;
	      }
	    }
          } else {
	    if (mode == 0) sprintf(triggerbase,"HLT_Mu5_Ele17");
	    else if (mode == 1) sprintf(triggerbase,"HLT_DoubleMu5");
	    else sprintf(triggerbase,"HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R");
	    
	    if (mode == 0) sprintf(secondtriggerbase,"HLT_Mu11_Ele8");
	    
	    for (int i = 1; i < 16; i++){
	      sprintf(triggername,"%s_v%d", triggerbase, i);
	      itrigger = treeLoader.iTrigger (string (triggername), currentRun);
	      if (itrigger != 9999) i = 100;
	      else if (i == 15) cout << "NO VALID TRIGGER FOUND FOR THIS RUN " << event->runId() << endl;
	    }
	    if (mode == 0){
	      for (int i = 1; i < 16; i++){
	        sprintf(triggername,"%s_v%d", secondtriggerbase, i);
	        isecondtrigger = treeLoader.iTrigger (string (triggername), currentRun);
	        if (isecondtrigger != 9999) i = 100;
	        else if (i == 15) cout << "NO VALID TRIGGER FOUND FOR THIS RUN " << event->runId() << endl;
	      }
	    }
	    
	  }
	  
	  
	  bool trigged = false;
	  if (mode == 0) trigged = itrigger + isecondtrigger;
	  else trigged = itrigger;
	  
	  
	  Selection selection(init_jets, init_muons, init_electrons, mets);
	  
	  bool isGoodPV = isGoodPV = selection.isPVSelected(vertex, 4,24,2.);  
	  
	  cutflow->Fill(1, weight);
	  cutflow_raw->Fill(1);
	  if(trigged){
	    cutflow->Fill(2, weight);
	    cutflow_raw->Fill(2);
	    if(isGoodPV){
	      cutflow->Fill(3, weight);
	      cutflow_raw->Fill(3);
	      
	      vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(vertex[0]);
	      vector<TRootMuon*> selectedMuons = selection.GetSelectedDiMuons();
	      vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseElectrons(true);
	      vector<TRootMuon*> looseMuons = selection.GetSelectedLooseMuons();
	      vector<TRootJet*> selectedJets = selection.GetSelectedJets(15,2.4,selectedElectrons,true);
	      
	      bool leptonSelection = false;
	      if 	(mode == 0 && selectedElectrons.size()== 1 && selectedMuons.size()== 1) leptonSelection = true;
	      else if 	(mode == 1 && selectedElectrons.size()== 0 && selectedMuons.size()== 2) leptonSelection = true;
	      else if 	(mode == 2 && selectedElectrons.size()== 2 && selectedMuons.size()== 0) leptonSelection = true;
	      
	      if (leptonSelection) {
		
		bool charge = false;
		double q0, q1;
		TLorentzVector lepton0, lepton1;
		
		if (mode == 0){
		  TRootElectron* electron = (TRootElectron*) selectedElectrons[0];
		  TRootMuon* muon = (TRootMuon*) selectedMuons[0];
		  if (electron->charge()*muon->charge() < 0) charge = true;
		  if (electron->Pt() > muon->Pt()){
		    lepton0.SetPxPyPzE(electron->Px(), electron->Py(), electron->Pz(), electron->Energy());
		    lepton1.SetPxPyPzE(muon->Px(), muon->Py(), muon->Pz(), muon->Energy());
		    q0 = electron->charge();
		    q1 = muon->charge();
		  } else {
		    lepton0.SetPxPyPzE(muon->Px(), muon->Py(), muon->Pz(), muon->Energy());
		    lepton1.SetPxPyPzE(electron->Px(), electron->Py(), electron->Pz(), electron->Energy());
		    q0 = muon->charge();
		    q1 = electron->charge();
		  }
		} else if (mode == 1){
		  TRootMuon* muon0 = (TRootMuon*) selectedMuons[0];
		  TRootMuon* muon1 = (TRootMuon*) selectedMuons[1];
		  if (muon0->charge()*muon1->charge() < 0) charge = true;
		  lepton0.SetPxPyPzE(muon0->Px(), muon0->Py(), muon0->Pz(), muon0->Energy());
		  lepton1.SetPxPyPzE(muon1->Px(), muon1->Py(), muon1->Pz(), muon1->Energy());
		  q0 = muon0->charge();
		  q1 = muon1->charge();
		} else {
		  TRootElectron* electron0 = (TRootElectron*) selectedElectrons[0];
		  TRootElectron* electron1 = (TRootElectron*) selectedElectrons[1];
		  if (electron0->charge()*electron1->charge() < 0) charge = true;
		  lepton0.SetPxPyPzE(electron0->Px(), electron0->Py(), electron0->Pz(), electron0->Energy());
		  lepton1.SetPxPyPzE(electron1->Px(), electron1->Py(), electron1->Pz(), electron1->Energy());
		  q0 = electron0->charge();
		  q1 = electron1->charge();
		}
		
		if (charge){
		  cutflow->Fill(4, weight);
		  cutflow_raw->Fill(4);
		  
		  bool leptonVeto = false;
		  if 	  (mode == 0 && looseMuons.size()== 1 && looseElectrons.size() == 1) leptonVeto = true;
		  else if (mode == 1 && looseMuons.size()== 2 && looseElectrons.size() == 0) leptonVeto = true;
		  else if (mode == 2 && looseMuons.size()== 0 && looseElectrons.size() == 2) leptonVeto = true;
		  
		  if (leptonVeto) {
		    cutflow->Fill(5, weight);
		    cutflow_raw->Fill(5);
		    
		    TLorentzVector pair = lepton0 + lepton1;   
		    if ((mode !=0 && pair.M() > 20) || mode == 0){
		      lum = lumi;
		      
		      xlWeight = weight;
		      
		      metPt = mets[0]->Pt();
		      metPx = mets[0]->Px();
		      metPy = mets[0]->Py();
		      
		      ptLepton = new std::vector<double>; 
		      pxLepton = new std::vector<double>; 
		      pyLepton = new std::vector<double>; 
		      pzLepton = new std::vector<double>; 
		      eLepton = new std::vector<double>; 
		      qLepton = new std::vector<double>;
		      
		      ptJet = new std::vector<double>; 
		      pxJet = new std::vector<double>; 
		      pyJet = new std::vector<double>; 
		      pzJet = new std::vector<double>; 
		      eJet = new std::vector<double>; 
		      qJet = new std::vector<double>; 
		      btTCHPJet = new std::vector<double>; 
		      btTCHEJet = new std::vector<double>; 
		      btSSVHPJet = new std::vector<double>;
		      btSSVHEJet = new std::vector<double>; 
		    
		      ptLepton->push_back(lepton0.Pt());
		      ptLepton->push_back(lepton1.Pt());
		      
		      pxLepton->push_back(lepton0.Px());
		      pxLepton->push_back(lepton1.Px());
		      
		      pyLepton->push_back(lepton0.Py());
		      pyLepton->push_back(lepton1.Py());
		      
		      pzLepton->push_back(lepton0.Pz());
		      pzLepton->push_back(lepton1.Pz());
		    
		      eLepton->push_back(lepton0.Energy());
		      eLepton->push_back(lepton1.Energy());
		      
		      qLepton->push_back(q0);
		      qLepton->push_back(q1);
		      
		      for (int i =0; i < selectedJets.size(); i ++){
			TRootJet* tempJet = (TRootJet*) selectedJets[i];
			ptJet->push_back(tempJet->Pt());
			pxJet->push_back(tempJet->Px());
			pyJet->push_back(tempJet->Py());
			pzJet->push_back(tempJet->Pz());
			eJet->push_back(tempJet->Energy());
			qJet->push_back(tempJet->charge());
			btTCHPJet->push_back(tempJet->btag_trackCountingHighPurBJetTags() );
			btTCHEJet->push_back(tempJet->btag_trackCountingHighEffBJetTags() );
			btSSVHPJet->push_back(tempJet->btag_simpleSecondaryVertexHighPurBJetTags() );
			btSSVHEJet->push_back(tempJet->btag_simpleSecondaryVertexHighEffBJetTags() );  
		      }
		    
		      myTree->Fill();
		      
		      delete ptLepton;
		      delete pxLepton;
		      delete pyLepton;
		      delete pzLepton;
		      delete eLepton;
		      delete qLepton;
		      
		      delete ptJet;
		      delete pxJet;
		      delete pyJet;
		      delete pzJet;
		      delete eJet;
		      delete qJet;
		      delete btTCHPJet;
		      delete btTCHEJet;
		      delete btSSVHPJet;
		      delete btSSVHEJet;
		      
		      if (pair.M() > 101 || pair.M() < 81 || mode == 0){
			cutflow->Fill(6, weight);
			cutflow_raw->Fill(6);
			if (mets[0]->Pt() > 30 || mode == 0){
			  cutflow->Fill(7, weight);
			  cutflow_raw->Fill(7);
			  
			  int nJetsBT = 0;
			  int nJets = 0;
			  bool bTagged = false;
			  int iJet = -5;
			  int iSF;
			  for (int i =0; i < selectedJets.size(); i ++){
			    TRootJet* tempJet = (TRootJet*) selectedJets[i];
			    TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
			    if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
			      iSF = rand() % 101;
			      if (iSF < 96 || sample == 666) nJetsBT++;
			    }
			    if (tempJet->Pt() > 30 && TMath::Min(fabs(lepton0.DeltaR(tJet)), fabs(lepton1.DeltaR(tJet))) > 0.3) {
			      nJets++;
			      iJet = i;
			      if (tempJet->btag_simpleSecondaryVertexHighEffBJetTags() > 1.74){
				iSF = rand() % 101;
				if (iSF < 96 || sample == 666) bTagged = true;
			      }
			    }
			  }
			  if(nJets == 1){
			    TRootJet* jet = (TRootJet*) selectedJets[iJet];
			    cutflow->Fill(8, weight);
			    cutflow_raw->Fill(8);
			    if (bTagged && nJetsBT == 1){
			      cutflow->Fill(9,weight);
			      cutflow_raw->Fill(9);
			      
			      double ptSysPx = lepton0.Px() + lepton1.Px() + jet->Px() + mets[0]->Px();
			      double ptSysPy = lepton0.Py() + lepton1.Py() + jet->Py() + mets[0]->Py();
			      double ptSys = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
			      double Ht = lepton0.Pt() + lepton1.Pt() + jet->Pt() + mets[0]->Pt(); 
			      
			      if (ptSys < 60){
				cutflow->Fill(10, weight);
				cutflow_raw->Fill(10);
				if (Ht > 160 || mode != 0){
				  cutflow->Fill(11, weight);
				  cutflow_raw->Fill(11);
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
	    }
	  
	  }
	} // event loop
    } // dataset loop
  
  
  cout << "*************************************************" << endl;
  cout << "Results Normalized: " <<  endl;
  cout << "*************************************************" << endl;
  cout << "All:       " <<  cutflow->GetBinContent(2) << " +/- "  << cutflow->GetBinError(2) << endl;
  cout << "HLT:       " <<  cutflow->GetBinContent(3) << " +/- "  << cutflow->GetBinError(3) << endl;
  cout << "PV:        " <<  cutflow->GetBinContent(4) << " +/- "  << cutflow->GetBinError(4) << endl;
  cout << "Lep. Sel:  " <<  cutflow->GetBinContent(5) << " +/- "  << cutflow->GetBinError(5) << endl;
  cout << "Lep. Veto: " <<  cutflow->GetBinContent(6) << " +/- "  << cutflow->GetBinError(6) << endl;
  cout << "mll:       " <<  cutflow->GetBinContent(7) << " +/- "  << cutflow->GetBinError(7) << endl;
  cout << "MET:       " <<  cutflow->GetBinContent(8) << " +/- "  << cutflow->GetBinError(8) << endl;
  cout << "1 jet:     " <<  cutflow->GetBinContent(9) << " +/- "  << cutflow->GetBinError(9) << endl;
  cout << "1 jet BT:  " <<  cutflow->GetBinContent(10) << " +/- " << cutflow->GetBinError(10) << endl;
  cout << "Pt tW:     " <<  cutflow->GetBinContent(11) << " +/- " << cutflow->GetBinError(11) << endl;
  cout << "Ht:        " <<  cutflow->GetBinContent(12) << " +/- " << cutflow->GetBinError(12) << endl;
  cout << "*************************************************" << endl; 

  fout->Write();
  fout->Close();
  
  cout << "It took you " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
}

