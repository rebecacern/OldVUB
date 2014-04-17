// simple file to execute the b-tag efficiency measurement function in myNTupleAnalyzer.cc

#include <TString.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>
#include "../interface/myNTupleAnalyzer.h"
using namespace std;

int main () {
	
	int nBins = 50;

	float min = -10;
	float max = 30;
	
	bool fixSize = false;
	
	if (fixSize) {
		
		float diff = max-min;

		float binSize = diff/nBins;
		
		cout << "Fix Binning: double rangesbTag[" << nBins+1 << "] = {";
		
		for (unsigned int i=0; i<nBins;i++) {
			
			cout << min+(i*binSize) << ",";
			
		}cout << max << "};" << endl;
	
		return 0;
		
	}

	std::vector<std::string> samples;

	samples.push_back("TTrees/BtagTree_TTbarJets_SemiMuon.root");	
	samples.push_back("TTrees/BtagTree_TTbarJets_Other.root");
	samples.push_back("TTrees/BtagTree_ST_sChannel.root");
	samples.push_back("TTrees/BtagTree_ST_tChannel.root");
	samples.push_back("TTrees/BtagTree_ST_tWChannel.root");
	samples.push_back("TTrees/BtagTree_WJets.root");
	samples.push_back("TTrees/BtagTree_ZJets.root");
	//samples.push_back("TTrees/BtagTree_QCD_Mu15.root");

	//samples.push_back("TTrees/BtagTree_Data.root");

	std::map<int,std::vector<float> > bTag;
	
	for (std::vector<std::string>::const_iterator it=samples.begin(); it != samples.end(); ++it) {
	
		TFile* f = new TFile((*it).c_str(),"READ");
		
		TRootNTuple *NTuple;
		TTree *tree_ = (TTree*)f->Get("tree");
		NTuple = new TRootNTuple();
		
		tree_->Branch("TheNTuple","TRootNTuple",&NTuple);	     	    
		
		TBranch *branch = tree_->GetBranch("TheNTuple");
		branch->SetAddress(&NTuple);
		
		Int_t nEvent = tree_->GetEntries();
		
		cout<<"+> Start event loop on " << nEvent << " events" << endl;
		for(int ev=0; ev < nEvent; ev++){   
			
			if((int) ev/100 == (double) ev/100)  {
					cout << "+>>> analyzing event " << ev << " of " << nEvent <<  flush << "\r" ;//<< endl;
					
			}    
			tree_->GetEvent(ev); 

			if (NTuple->btag_trackCountingHighEffBJetTags() > min && NTuple->btag_trackCountingHighEffBJetTags() < max)
			
			  bTag[0].push_back(NTuple->btag_trackCountingHighEffBJetTags());
			
			if (NTuple->btag_trackCountingHighPurBJetTags() > min && NTuple->btag_trackCountingHighPurBJetTags() < max)
			
				bTag[1].push_back(NTuple->btag_trackCountingHighPurBJetTags());
			/*bTag[2].push_back(NTuple->btag_combinedSecondaryVertexBJetTags());
			bTag[3].push_back(NTuple->btag_combinedSecondaryVertexMVABJetTags());
			bTag[4].push_back(NTuple->btag_impactParameterMVABJetTags());
			bTag[5].push_back(NTuple->btag_jetBProbabilityBJetTags());
			bTag[6].push_back(NTuple->btag_jetProbabilityBJetTags());
			bTag[7].push_back(NTuple->btag_simpleSecondaryVertexBJetTags());
			bTag[8].push_back(NTuple->btag_softElectronBJetTags());
			bTag[9].push_back(NTuple->btag_softMuonBJetTags());
			bTag[10].push_back(NTuple->btag_softMuonNoIPBJetTags());*/
			
			
		}
		delete NTuple;
		delete tree_;
		f->Close();
		
	}
	cout << endl;
	// calculate binning for btaggers
	
	for (std::map<int,std::vector<float> >::iterator it=bTag.begin(); it != bTag.end(); ++it) {
	
	  std::vector<float> bTagVals, binLowEdges;
		
		sort(it->second.begin(),it->second.end());
		
		int len = it->second.size();
		
		int nPerBin = len/nBins;
		
		//cout << nPerBin << endl;
		
		//cout << it->second[0] << endl;
		
		for (int i=0; i<it->second.size(); i++) 
			if (it->second[i] > -100)
				bTagVals.push_back(it->second[i]);
		
		//cout << "NEW " << bTagVals[0] << endl;

		for (int i=0; i<bTagVals.size()-1; i++) {

		  if (bTagVals[i] > bTagVals[i+1])
		    cout << "ERROR " << i << " " << i+1 << endl;
		  
		}cout << endl;
						
		for (int i=0; i<nBins; i++) {
		
			float min = 0;
			
			if (i*nPerBin < len) {
				//cout << i*nPerBin << endl;

				min = bTagVals[i*nPerBin];
			}
						
			if (min > -100)
			  binLowEdges.push_back(min);
				
			//cout << min << endl;
		}

		string lowBounds = "{";

		for (int i=0; i<binLowEdges.size()-1; i++) {

		  stringstream s; s << binLowEdges[i];

		  if (binLowEdges[i] > binLowEdges[i+1])
		    cout << "BinEdges ERROR " << i << " " << i+1 << endl;
		  else
		    lowBounds += s.str()+",";

		}
		
		cout << binLowEdges.size() << endl;
		stringstream t; t << binLowEdges[binLowEdges.size()-1];
		stringstream m; m << max;
			
		lowBounds += t.str()+","+m.str()+".}";
		
		cout << "bTagger #" << it->first << " Binning: double rangesbTag[" << nBins+1 << "] = " << lowBounds << ";" << endl;
		

	}

return 0;

}
