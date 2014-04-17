#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>

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
#include "../MCInformation/interface/JetPartonMatching.h"
#include "../Tools/interface/MVATrainer.h"
#include "../Tools/interface/MVAComputer.h"
#include "../MCInformation/interface/ResolutionFit.h"
#include "../JESMeasurement/interface/FullKinFit.h"
#include "Style.C"

using namespace std;
using namespace TopTree;

int main (int argc, char *argv[])
{

  clock_t start = clock();

  cout << "*****************************************************" << endl;
  cout << " Beginning of the program for TTbar JES Application! " << endl;
  cout << "*****************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

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
  string rootFileName ("TTbarJES_Application.root");

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
  float Luminosity = anaEnv.Luminosity;	// in 1/pb

  ////////////////////////////////////////////////////
  // open the ttbarJES output to get the corrections
  ////////////////////////////////////////////////////

  std::map<double,double> BCorrections,LightCorrections;

  TFile* JES = new TFile("TTbarJES.root","READ");

  TGraphErrors* CorrB = (TGraphErrors*) JES->Get("EstimatedDEb_VS_X");

  if (!CorrB)
    cout << "Error: the B energy scale corrections could not be loaded" << endl;
  else {

    cout << "JES: loading B energy scale corrections" << endl;
    for (unsigned int i=0; i<CorrB->GetN(); i++) {
  
      double x = 0;
      double y = 0;
      
      CorrB->GetPoint(i,x,y);
      
      cout << "B correction: Contour X=" << x << " DEl=" << y << endl;
      
      BCorrections[x]=y;
      
    }
  }
  
  TGraphErrors* CorrL = (TGraphErrors*) JES->Get("EstimatedDEl_VS_X");

  if (!CorrL)
    cout << "Error: the light energy scale corrections could not be loaded" << endl;
  else {

    cout << "JES: loading light energy scale corrections" << endl;
    for (unsigned int i=0; i<CorrL->GetN(); i++) {
  
      double x = 0;
      double y = 0;
      
      CorrL->GetPoint(i,x,y);
      
      cout << "Light correction: Contour X=" << x << " DEl=" << y << endl;
      
      LightCorrections[x]=y;
      
    }
  }

  JES->Close();
  delete JES;

  delete CorrB;
  delete CorrL;

  ////////////
  // HISTOS 
  ///////////

  std::map<string,TH1F*> histo1D;
  
  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(datasets[i]);
  /////////////////////

  ////////////////////////////////////////
  
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
      cout << "   Dataset " << d << ": " << datasets[d].Name () << "/ title : " << datasets[d].Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d].NofEvtsToRunOver() << " events." << endl;
    
    //open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"LoadEvent"<<endl;
    
    //vector of objects
    cout << " - Variable declaration ..." << endl;
    vector < TRootVertex > vertex;
    vector < TRootMuon > init_muons;
    vector < TRootElectron > init_electrons;
    vector < TRootJet > init_jets;
    vector < TRootCaloJet > init_Calojets;
    vector < TRootPFJet > init_PFjets;
    vector < TRootMET > mets;
    
    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    nEvents[d] = 0;
    int itrigger = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    for (unsigned int ievt = 0; ievt < datasets[d].NofEvtsToRunOver(); ievt++) {
      //for (unsigned int ievt = 0; ievt < 10000; ievt++) {
	nEvents[d]++;
	if(ievt%500 == 0)
	  std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
	
	float lumi = 50; // scale plots to this lumi (in pb-1)
      
	//load event
	if(anaEnv.JetType == 0 ) event = treeLoader.LoadEvent (ievt, &datasets[d], vertex, init_muons, init_electrons, init_jets, mets);
	if(anaEnv.JetType == 1 ) event = treeLoader.LoadEvent (ievt, &datasets[d], vertex, init_muons, init_electrons, init_Calojets, mets);
	if(anaEnv.JetType == 2 ) event = treeLoader.LoadEvent (ievt, &datasets[d], vertex, init_muons, init_electrons, init_PFjets, mets);
	
	int currentRun = event->runId();
	if(previousRun != currentRun)
	  {
	    previousRun = currentRun;
	    itrigger = treeLoader.iTrigger (&datasets[d], string ("HLT_Mu9"), currentRun);
	  }
	

	/////////////////////////////
	//   Selection
	/////////////////////////////
	//Declare selection instance    
	Selection selection(anaEnv.JetType, init_jets, init_Calojets, init_PFjets, init_muons, init_electrons, mets);
	selection.SetConfiguration(&anaEnv);

	//      if(datasets[d].Name() == "Data" || datasets[d].Name() == "data" || datasets[d].Name() == "DATA")
//      {
//        bool passedJSON = treeLoader.EventPassedJSON(datasets[d], event->runId(), event->lumiBlockId());
//        if(!passedJSON) continue;
//      }

	bool trigged = treeLoader.EventTrigged (itrigger);
	bool isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);

	vector<TRootJet> selectedJets = selection.GetSelectedJets(anaEnv.JetsPtCutSR,anaEnv.JetsEtaCutSR, anaEnv.applyJetID);
	vector<TRootMuon> selectedMuons = selection.GetSelectedMuons(anaEnv.MuonPtCutSR,anaEnv.MuonEtaCutSR,anaEnv.MuonRelIsoCutSR, selectedJets);
	vector<TRootMuon> vetoMuons = selection.GetSelectedLooseMuons(anaEnv.MuonPtCutVetoSR,anaEnv.MuonEtaCutVetoSR,anaEnv.MuonRelIsoCutVetoSR);
	vector<TRootElectron> vetoElectrons = selection.GetSelectedLooseElectrons(anaEnv.ElectronPtCut,anaEnv.ElectronEtaCut,anaEnv.ElectronRelIsoCut);
	
	vector<TRootJet> selectedLooseJets = selection.GetSelectedJets(20., anaEnv.JetsEtaCutSR, true);
	
	vector<TRootMCParticle> mcParticles;
	
	if(datasets[d].Name().find("TTbarJets_SemiMu") == 0)
	  {
	    treeLoader.LoadMCEvent(ievt, &datasets[d], 0, 0, mcParticles);
	    
	    sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
	  }
	
	bool eventSelected = false;
	bool all4JetsMatched = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-mu partons
	bool all4PartonsMatched = false; // True if the 4 ttbar semi-mu partons are matched to 4 jets (not necessarily the 4 highest pt jets)
	
	if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-3){
	  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2){
	    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1){
	      if(selectedJets.size()>=(unsigned int)anaEnv.NofJets){
		eventSelected = true;
		if(verbose>2) cout << event->runId() << ":" << event->eventId() << ":" << event->lumiBlockId()  << ":" << setprecision(8) << selectedMuons[0].Pt() << endl;
	      }
	    }
	  }
	}

	if( ! (trigged && isGoodPV) ) continue;
	// Don't look further at events with a bad primary vertex or not passing the trigger
	
	sort(selectedLooseJets.begin(),selectedLooseJets.end(),HighestPt()); // HighestPt() is included from the Selection class
     
	pair<unsigned int, unsigned int> leptonicBJet, hadronicBJet, hadronicWJet1, hadronicWJet2;
      leptonicBJet = hadronicBJet = hadronicWJet1 = hadronicWJet2 = std::pair<unsigned int,unsigned int>(9999,9999);
        
      
      // Jet-MCParticle matching, only for TTbar semi-mu
      string dataSetName = datasets[d].Name();
      if(dataSetName.find("TTbarJets_SemiMu") == 0)
      {
      
        vector<TLorentzVector> mcParticlesTLV, selectedLooseJetsTLV;
        vector<TRootMCParticle> mcParticlesMatching; // MCParticles used for the matching
        
        bool muPlusFromTop = false, muMinusFromTop = false;
        int nTTbarQuarks = 0;
        for(unsigned int i=0; i<mcParticles.size(); i++)
        {
          if( mcParticles[i].status() != 3) continue;
          
          if( mcParticles[i].type() == 13 && mcParticles[i].motherType() == -24 && mcParticles[i].grannyType() == -6 )
          {
            if(muMinusFromTop) cerr<<"muMinusFromTop was already true"<<endl;
            muMinusFromTop = true;
          }
          if( mcParticles[i].type() == -13 && mcParticles[i].motherType() == 24 && mcParticles[i].grannyType() == 6 )
          {
            if(muPlusFromTop) cerr<<"muPlusFromTop was already true"<<endl;
            muPlusFromTop = true;
          }
          
          if( abs(mcParticles[i].type()) < 6 || abs(mcParticles[i].type()) == 21 )
          {
            mcParticlesTLV.push_back(mcParticles[i]);
            mcParticlesMatching.push_back(mcParticles[i]);
            
            if( fabs(mcParticles[i].motherType()) == 6 || fabs(mcParticles[i].grannyType()) == 6 )
            {
              nTTbarQuarks++;
            }
	  }
        }
        if(muPlusFromTop && muMinusFromTop)
          cerr<<"muPlusFromTop and muMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;

        // take all the selectedLooseJets to study the radiation stuff, selectedLooseJets are already ordened in decreasing Pt()
        if(selectedLooseJets.size() < 4) continue;
        
        for(unsigned int i=0; i<selectedLooseJets.size(); i++)
          selectedLooseJetsTLV.push_back(selectedLooseJets[i]);
        
        JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedLooseJetsTLV, 2, true, true, 0.3);
        
        if(matching.getNumberOfAvailableCombinations() != 1)
          cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;

        vector< pair<unsigned int, unsigned int> > JetPartonPair, ISRJetPartonPair; // First one is jet number, second one is mcParticle number
        
        for(unsigned int i=0; i<mcParticlesTLV.size(); i++)
        {
          int matchedJetNumber = matching.getMatchForParton(i, 0);
          if(matchedJetNumber != -1)
            JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
        }
       
        for(unsigned int i=0; i<JetPartonPair.size(); i++)
        {
	        unsigned int j = JetPartonPair[i].second;

//          cout<<"JetPartonPair["<<i<<"]: "<<endl;
//          cout<<"Jet Pt: "<<selectedLooseJets[JetPartonPair[i].first].Pt()<<"  Eta: "<<selectedLooseJets[JetPartonPair[i].first].Eta()<<"  Phi: "<<selectedLooseJets[JetPartonPair[i].first].Phi()<<endl;
//          cout<<"MCParticle Pt: "<<mcParticlesMatching[JetPartonPair[i].second].Pt()<<"  Eta: "<<mcParticlesMatching[JetPartonPair[i].second].Eta()<<"  Phi: "<<mcParticlesMatching[JetPartonPair[i].second].Phi()<<endl;
//          cout<<"\ttype: "<<mcParticlesMatching[JetPartonPair[i].second].type()<<"  motherType: "<<mcParticlesMatching[JetPartonPair[i].second].motherType()<<"  grannyType: "<<mcParticlesMatching[JetPartonPair[i].second].grannyType()<<endl;
//          cout<<"DR(MCParticle, jet) = "<<selectedLooseJets[JetPartonPair[i].first].DeltaR(mcParticlesMatching[JetPartonPair[i].second])<<endl;
          
          if( fabs(mcParticlesMatching[j].type()) < 6 )
          {
            if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == -24 && mcParticlesMatching[j].grannyType() == -6 )
              || ( muMinusFromTop && mcParticlesMatching[j].motherType() == 24 && mcParticlesMatching[j].grannyType() == 6 ) )
            {
              if(hadronicWJet1.first == 9999) 
		            hadronicWJet1 = JetPartonPair[i];
              else if(hadronicWJet2.first == 9999) 
		            hadronicWJet2 = JetPartonPair[i];
              else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
            }
          }
          if( fabs(mcParticlesMatching[j].type()) == 5 )
          {
            if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == -6 )
              || ( muMinusFromTop && mcParticlesMatching[j].motherType() == 6 ) )
              hadronicBJet = JetPartonPair[i];
            else if( ( muPlusFromTop && mcParticlesMatching[j].motherType() == 6 )
              || ( muMinusFromTop && mcParticlesMatching[j].motherType() == -6 ) )
              leptonicBJet = JetPartonPair[i];
          }
          
          // look for ISR stuff
          if( fabs(mcParticlesMatching[j].type()) != 6 && fabs(mcParticlesMatching[j].motherType()) != 24 && fabs(mcParticlesMatching[j].motherType()) != 6 &&
            fabs(mcParticlesMatching[j].grannyType()) != 24 && fabs(mcParticlesMatching[j].grannyType()) != 6 )
          {
            ISRJetPartonPair.push_back(JetPartonPair[i]);
          }
        }

        if(hadronicWJet1.first != 9999 && hadronicWJet2.first != 9999 && hadronicBJet.first != 9999 && leptonicBJet.first != 9999)
        {
          all4PartonsMatched = true;
          if(hadronicWJet1.first < 4 && hadronicWJet2.first < 4 && hadronicBJet.first < 4 && leptonicBJet.first < 4)
            all4JetsMatched = true;
        }

	//cout << "ievt " << ievt << " all matched " << all4JetsMatched << endl;

	// DRAW THE W and TOP masses.

	if (eventSelected && all4JetsMatched) {
	  
	  // uncorrected

	  int i = hadronicWJet1.first;
	  int j = hadronicWJet2.first;
	  int k = hadronicBJet.first;

	  TLorentzVector WBoson = selectedJets[i]+selectedJets[j];
	  TLorentzVector HadBQuark = selectedJets[k];
	  
	  TLorentzVector TopQuark = WBoson+HadBQuark;
	  
	  if (histo1D.find("WMass") == histo1D.end()) {
	    
	    histo1D["WMass"] = new TH1F("WMass_Uncorrected","W mass distribution without corrections;m_{W},Number of events",200,0,200);
	    
	  }
	  
	  //cout << WBoson.M() << endl;
	  histo1D["WMass"]->Fill(WBoson.M());
	  
	  if (histo1D.find("TopMass") == histo1D.end()) {
	    
	    histo1D["TopMass"] = new TH1F("TopMass_Uncorrected","Top mass distribution without corrections",400,0,400);
	    
	  }
	  
	  histo1D["TopMass"]->Fill(TopQuark.M());

	  // apply corrections

	  for (std::map<double,double>::const_iterator it=LightCorrections.begin(); it != LightCorrections.end() ; ++it) {

	    stringstream s; s << it->first;
	    
	    string titleW = "WMass_LightCorr_X"+s.str();

	    if (histo1D.find(titleW) == histo1D.end()) {
	    
	      histo1D[titleW] = new TH1F(titleW.c_str(),"W mass distribution with corrections",200,0,200);
	    
	    }

	    /*TLorentzVector Jet1 = selectedJets[i];
	    TLorentzVector CorrJet1 = selectedJets[i]+(selectedJets[i]*it->second);
	    TLorentzVector Jet2 = selectedJets[j];
	    TLorentzVector CorrJet2 = selectedJets[j]+(selectedJets[j]*it->second);

	    if (!strcmp(s.str().c_str(),"0.3")) {
	    cout << "---------------------------------" << it->first << "---------------------------" << endl;
	    cout << "Correction: " << it->second << endl; 
	    cout << "Jet1 (E,px,py,pz): " << Jet1.E() << " " << Jet1.Px() << " " << Jet1.Py() << " " << Jet1.Pz() << endl;
	    cout << "CorrJet1 (E,px,py,pz): " << CorrJet1.E() << " " << CorrJet1.Px() << " " << CorrJet1.Py() << " " << CorrJet1.Pz() << endl;
	    cout << "Diff1 (E,px,py,pz): " << (-Jet1.E()+CorrJet1.E())/CorrJet1.E() << " " << (-Jet1.Px()+CorrJet1.Px())/CorrJet1.Px() << " " << (-Jet1.Py()+CorrJet1.Py())/CorrJet1.Py() << " " << (-Jet1.Pz()+CorrJet1.Pz())/CorrJet1.Pz() << endl; 
	    cout << "Jet2 (E,px,py,pz): " << Jet2.E() << " " << Jet2.Px() << " " << Jet2.Py() << " " << Jet2.Pz() << endl;
	    cout << "Diff2 (E,px,py,pz): " << (-Jet2.E()+CorrJet2.E())/CorrJet2.E() << " " << (-Jet2.Px()+CorrJet2.Px())/CorrJet2.Px() << " " << (-Jet2.Py()+CorrJet2.Py())/CorrJet2.Py() << " " << (-Jet2.Pz()+CorrJet2.Pz())/CorrJet2.Pz() << endl;
	    cout << "CorrJet2 (E,px,py,pz): " << CorrJet2.E() << " " << CorrJet2.Px() << " " << CorrJet2.Py() << " " << CorrJet2.Pz() << endl;
	    cout << "WMass before " << (Jet1+Jet2).M() << endl; 
	    cout << "WMass after " << (CorrJet1+CorrJet2).M() << endl;
	    cout << "Diff Mass" << ((Jet1+Jet2).M()-(CorrJet1+CorrJet2).M())/(Jet1+Jet2).M()<< endl;
	    cout << "Abs Diff Mass" << ((Jet1+Jet2).M()-(CorrJet1+CorrJet2).M())<< endl;
	    }*/
	    TLorentzVector CorrectedWBoson = (selectedJets[i]+(selectedJets[i]*it->second))+(selectedJets[j]+(selectedJets[j]*it->second));

	    histo1D[titleW]->Fill(CorrectedWBoson.M());

	  }

	  for (std::map<double,double>::const_iterator it=BCorrections.begin(); it != BCorrections.end() ; ++it) {

	    stringstream s; s << it->first;
	    
	    string title = "TopMass_Light_B_Corr_X"+s.str();

	    if (histo1D.find(title) == histo1D.end()) {
	    
	      histo1D[title] = new TH1F(title.c_str(),"Top mass distribution with corrections",400,0,400);
	    
	    }

	    TLorentzVector CorrectedWBoson = (selectedJets[i]+(selectedJets[i]*LightCorrections[it->first]))+(selectedJets[j]+(selectedJets[j]*LightCorrections[it->first]));
	    TLorentzVector CorrectedTopQuark = CorrectedWBoson+(selectedJets[k]+(selectedJets[k]*it->second));

	    histo1D[title]->Fill(CorrectedTopQuark.M());

	  }
	}
      } 
    }
  }

  // WRITING PLOTS

  string pathPNG = "PlotsJES/";

  fout->cd();

  // mW

  vector<TH1F*> listOfWHistos; listOfWHistos.push_back(histo1D["WMass"]);

  TF1* fitUncorr = new TF1("Uncorrected W mass fit","gaus");
  
  double left = histo1D["WMass"]->GetMean()-histo1D["WMass"]->GetRMS();
  double right = histo1D["WMass"]->GetMean()+histo1D["WMass"]->GetRMS();

  fitUncorr->SetRange(left,right);
    
  histo1D["WMass"]->Fit(fitUncorr,"RQNO");
  
  fitUncorr->SetLineColor(1);

  // mtop

  vector<TH1F*> listOfTopHistos; listOfTopHistos.push_back(histo1D["TopMass"]);

  TF1* mTopfitUncorr = new TF1("Uncorrected Top mass fit","gaus");
  
  double left2 = histo1D["TopMass"]->GetMean()-histo1D["TopMass"]->GetRMS();
  double right2 = histo1D["TopMass"]->GetMean()+histo1D["TopMass"]->GetRMS();

  mTopfitUncorr->SetRange(left2,right2);
    
  histo1D["TopMass"]->Fit(mTopfitUncorr,"RQNO");
  
  mTopfitUncorr->SetLineColor(1);

  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++) {
      TH1F *temp = it->second;
      int N = temp->GetNbinsX();
      //temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
      //temp->SetBinContent(N+1,0);
      temp->Write();
      TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
      tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );

      if (strstr(it->first.c_str(),"WMass_")) {
	listOfWHistos.push_back(it->second);
	//TCanvas* tempCanvas = TCanvasCreator(listOfWHistos,"lp","WMass Correction",1);
	//tempCanvas->SaveAs( (pathPNG+"CombinedWMass"+it->first+".png").c_str() );
	//tempCanvas->Write();

	TLegend* leg = new TLegend(0.75,0.8,0.99,1.);
	leg->SetFillColor(0);
	leg->SetShadowColor(0);
	
	TCanvas *c1 = new TCanvas("WMass","WMass");
	c1->cd();
	float ymax = 1.1*YMax(listOfWHistos);

	for(unsigned int i=0;i<listOfWHistos.size();i++){

	  if (i > 0) {
	    
	    TF1* fitCorr = new TF1("Corrected W mass fit","gaus");
	    
	    double left = listOfWHistos[i]->GetMean()-listOfWHistos[i]->GetRMS();
	    double right = listOfWHistos[i]->GetMean()+listOfWHistos[i]->GetRMS();
	    
	    fitCorr->SetRange(left,right);
	    
	    listOfWHistos[i]->Fit(fitCorr,"RQNO");
	    
	    fitCorr->SetLineColor(i+1);
	    
	    fitCorr->Draw("same");

	    stringstream s,t; s<<fitCorr->GetParameter(1); t << fitCorr->GetParError(1);
	    string mass = "m_{W}^{Fit-Corrected} = "+s.str()+" #pm "+t.str()+" GeV/c^{2}";
	    
	    double y = 0.8*ymax;
	    TLatex* text1 = new TLatex(0,y,mass.c_str());
	    text1->SetTextColor(i+1);
	    text1->Draw();

	    //delete text1;
	  }
	  
	  fitUncorr->Draw("same");

	  stringstream s,t; s<<fitUncorr->GetParameter(1); t << fitUncorr->GetParError(1);
	  string massUC = "m_{W}^{Fit-Uncorrected} = "+s.str()+" #pm "+t.str()+" GeV/c^{2}";

	  double y = 0.9*ymax;
	  TLatex* text1 = new TLatex(0,y,massUC.c_str());
	  text1->Draw();

	  //delete text1;

	  listOfWHistos[i]->SetLineColor(i+1);
	  //listOfWHistos[i]->SetTitle("WMass");

	  listOfWHistos[i]->GetXaxis()->SetTitle("m_{W} [GeV/c^{2}]");
	  listOfWHistos[i]->GetYaxis()->SetTitle("Number of events");
	  listOfWHistos[i]->SetTitle("m_{W}");

	  if(i==0){
	    listOfWHistos[i]->GetYaxis()->SetRangeUser(0.,ymax);
	    listOfWHistos[i]->Draw();
	  }
	  else listOfWHistos[i]->Draw("same");

	  if (i==0)
	    leg->AddEntry(listOfWHistos[i],"WMass NO extra corrections" , "lp");
	  else
	  leg->AddEntry(listOfWHistos[i],it->first.c_str() , "lp");

	}
	leg->Draw();

	listOfWHistos.pop_back();
	//	delete tempCanvas;

	c1->Write();
	c1->SaveAs( (pathPNG+"CombinedWMass"+it->first+".png").c_str() );
	delete c1;
	delete leg;

      }

      // mtop

      if (strstr(it->first.c_str(),"TopMass_")) {
	listOfTopHistos.push_back(it->second);
	//TCanvas* tempCanvas = TCanvasCreator(listOfWHistos,"lp","WMass Correction",1);
	//tempCanvas->SaveAs( (pathPNG+"CombinedWMass"+it->first+".png").c_str() );
	//tempCanvas->Write();

	TLegend* leg = new TLegend(0.75,0.8,0.99,1.);
	leg->SetFillColor(0);
	leg->SetShadowColor(0);
	
	TCanvas *c1 = new TCanvas("TopMass","TopMass");
	c1->cd();
	float ymax = 1.1*YMax(listOfWHistos);

	for(unsigned int i=0;i<listOfTopHistos.size();i++){

	  if (i > 0) {
	    
	    TF1* fitCorr = new TF1("Corrected Top mass fit","gaus");
	    
	    double left = listOfTopHistos[i]->GetMean()-listOfTopHistos[i]->GetRMS();
	    double right = listOfTopHistos[i]->GetMean()+listOfTopHistos[i]->GetRMS();
	    
	    fitCorr->SetRange(left,right);
	    
	    listOfTopHistos[i]->Fit(fitCorr,"RQNO");
	    
	    fitCorr->SetLineColor(i+1);

	    fitCorr->Draw("same");

	    stringstream s,t; s<<fitCorr->GetParameter(1); t << fitCorr->GetParError(1);
	    string mass = "m_{Top}^{Fit-Corrected} = "+s.str()+" #pm "+t.str()+" GeV/c^{2}";
	    
	    double y = 0.8*ymax;
	    TLatex* text1 = new TLatex(0,y,mass.c_str());

	    text1->SetTextColor(i+1);

	    text1->Draw();

	    //delete text1;
	  }
	  
	  mTopfitUncorr->Draw("same");

	  stringstream s,t; s<<mTopfitUncorr->GetParameter(1); t << mTopfitUncorr->GetParError(1);
	  string massUC = "m_{Top}^{Fit-Uncorrected} = "+s.str()+" #pm "+t.str()+" GeV/c^{2}";

	  double y = 0.9*ymax;
	  TLatex* text1 = new TLatex(0,y,massUC.c_str());
	  text1->Draw();

	  //delete text1;

	  listOfTopHistos[i]->SetLineColor(i+1);
	  //listOfTopHistos[i]->SetTitle("WMass");

	  listOfTopHistos[i]->GetXaxis()->SetTitle("m_{t} [GeV/c^{2}]");
	  listOfTopHistos[i]->GetYaxis()->SetTitle("Number of events");
	  listOfTopHistos[i]->SetTitle("m_{t}");

	  if(i==0){
	    listOfTopHistos[i]->GetYaxis()->SetRangeUser(0.,ymax);
	    listOfTopHistos[i]->Draw();
	  }
	  else listOfTopHistos[i]->Draw("same");

	  if (i==0)
	    leg->AddEntry(listOfTopHistos[i],"WMass NO extra corrections" , "lp");
	  else
	  leg->AddEntry(listOfTopHistos[i],it->first.c_str() , "lp");

	}
	leg->Draw();

	listOfTopHistos.pop_back();
	//	delete tempCanvas;

	c1->Write();
	c1->SaveAs( (pathPNG+"CombinedTopMass"+it->first+".png").c_str() );
	delete c1;
	delete leg;

      }
  }
  
  //TCanvas* tempCanvas = TCanvasCreator(listOfWHistos,"lp","WMass Correction",0);
  //tempCanvas->SaveAs( (pathPNG+"WMassJESApplied.png").c_str() );
  //tempCanvas->Write();

  fout->Close();

  delete fout;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;

  return 0;

}
