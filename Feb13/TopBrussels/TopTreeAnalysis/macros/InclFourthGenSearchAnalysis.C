//compile using:
//g++   -I `root-config --incdir` `root-config --libs` macroname.C -o macroname
//

// general
#include "TH1.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TF1.h"
#include "TKey.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFormula.h"
#include "TAxis.h"
#include "TRandom3.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "tdrstyle.C"


TFile *file, *JESPlus_file, *JESMinus_file, *JERPlus_file, *JERMinus_file;

const int nPE = 10000;
const int nmasses = 6;
const int nboxes = 6;
const int ndatasetnamesSM = 7;
const int nCKM_A_values = 10;	

bool doSystematics = true; // if this is false, the next booleans will be ignored
bool doLumiUnc = false; 
float lumisyst = 0.06; // integrated luminosity uncertainty is taken to be 6%
bool doJESUnc = false;
bool doJERUnc = false;

bool doXSTTJetsUnc = false;
float XSTTJets = 0.10;
bool doXSWJetsUnc = true;
float XSWJets = 0.30;


std::string IntToStr( int n )
{
        std::ostringstream result;
        result << n;
        return result.str();
}


void SortVectorAscending(std::vector<float> &pvalues){
	std::sort(pvalues.begin(),pvalues.end());
}

//////////////////////////////////////////
// Function to calculate the percentile //
//////////////////////////////////////////
// calculate the percentile (related to the specified second argument) of e.g. p-values;
// the 50th percentile is the median, and in principle, the 95th percentile is p95
float GetProbaPercentile(float pvalues[nmasses][nPE], double percent,  int mass, int nofPE){
	std::vector<float> vectorpvalues;
	for(unsigned int i=0;i<nofPE;i++){
		vectorpvalues.push_back(pvalues[mass][i]);
	}
	SortVectorAscending(vectorpvalues);
	double pvaluesArray[nPE];  // needs to be an array of doubles (Double_t) to make use of TMath::Quantiles();
	for(unsigned int i=0;i<nofPE;i++){
		pvaluesArray[i] = vectorpvalues[i];
	}
	vectorpvalues.clear();
	double probability[1], quantile[1];
	probability[0] = percent/100;
	TMath::Quantiles(nofPE,1,pvaluesArray,quantile,probability);  //see ROOT documentation for non-default calculation of quantiles
	return (float) quantile[0];
}

//////////////////////////////////////////////////////////////////////////
// Function to calculate the p value corresponding to the given Q value //
//////////////////////////////////////////////////////////////////////////
float Proba_;
float CalculateProbability(float Qvalues_MC[nmasses][nPE], float Qvalue_data, int mass, int nofPE){
	int verbosity = 0;

	int integral_larger = 0;
	int integral_total = 0;
	for(unsigned int i=0; i<nofPE; i++){
		if(!(Qvalues_MC[mass][i] != Qvalues_MC[mass][i]) && !(Qvalue_data != Qvalue_data)){  // to exclude nan
			if(Qvalues_MC[mass][i]>Qvalue_data) 
				integral_larger++;
			integral_total++;
		}
		else std::cout << "Warning: nan Q value" << std::endl;
	}
	if(verbosity>0) 
		std::cout << "Qvalue_data: " << Qvalue_data << ", integral_larger: " << integral_larger << ", integral_total: " << integral_total << std::endl;
	if(integral_total!=0){
		Proba_ = float(integral_larger)/float(integral_total);
	}
	else 
		std::cout << "Warning: integral Monte Carlo Q distribution = 0" << std::endl;
	return Proba_;	
}






int main()
{
  using namespace std;

  setTDRStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gSystem->mkdir("SomePlots");


	string postfix = "";
	if(doSystematics && doLumiUnc) postfix = postfix+"_LumiUnc";
	if(doSystematics && doJESUnc) postfix = postfix+"_JESUnc";
	if(doSystematics && doJERUnc) postfix = postfix+"_JERUnc";
	if(doSystematics && doXSTTJetsUnc) postfix = postfix+"_TTJetsUnc";
	if(doSystematics && doXSWJetsUnc) postfix = postfix+"_WJetsUnc";

	string outfilename = "InclFourthGenSearch_Qvalues"+postfix+".root";
  	TFile *fout = new TFile (outfilename.c_str(), "RECREATE");

	cout << "DO SYSTEMATICS " << doSystematics << endl;	
	if(doSystematics){
		cout << "	do lumi syst " << doLumiUnc << " with uncertainty " << lumisyst*100 << "%" << endl;	
		cout << "	do JES syst " << doJESUnc << endl;	
		cout << "	do JER syst " << doJERUnc << endl;	
		cout << "	do XS TTJets syst " << doXSTTJetsUnc << " with uncertainty " << XSTTJets*100 << "%" << endl;	
		cout << "	do XS WJets syst " << doXSWJetsUnc << " with uncertainty " << XSWJets*100 << "%" << endl;	
  	}
	
	TString fileName = "InclFourthGenSearch_500InvPb.root";
	file = TFile::Open(fileName);
	
	TString JESPlus_fileName = "InclFourthGenSearch_JESPlus.root";
	if (doSystematics && doJESUnc) JESPlus_file = TFile::Open(JESPlus_fileName);
	TString JESMinus_fileName = "InclFourthGenSearch_JESMinus.root";
	if (doSystematics && doJESUnc) JESMinus_file = TFile::Open(JESMinus_fileName);
	
	TString JERPlus_fileName = "InclFourthGenSearch_JERPlus.root";
	if (doSystematics && doJERUnc) JERPlus_file = TFile::Open(JERPlus_fileName);
	TString JERMinus_fileName = "InclFourthGenSearch_JERMinus.root";
	if (doSystematics && doJERUnc) JERMinus_file = TFile::Open(JERMinus_fileName);


	if (!file) abort();
	if (doSystematics && doJESUnc && (!JESPlus_file || !JESMinus_file)) abort();
	if (doSystematics && doJERUnc && (!JERPlus_file || !JERMinus_file)) abort();
	
	//make sure the 3 single top processes are the 3 last in the array!!!
	string datasetnamesSM[ndatasetnamesSM] = {"TTJets_SemiMuon","TTJets_Other","WJets","ZJets","ST_tChannel","ST_tWChannel","ST_sChannel"};
	string tprime[nmasses] = {"Tprime350","Tprime400","Tprime450","Tprime500","Tprime550","Tprime600"};
	string stprime[nmasses] = {"STprime350","STprime400","STprime450","STprime500","STprime550","STprime600"};
	string bprime[nmasses] = {"Bprime350","Bprime400","Bprime450","Bprime500","Bprime550","Bprime600"};
	string sbprime[nmasses] = {"SBprime350","SBprime400","SBprime450","SBprime500","SBprime550","SBprime600"};
	string data = "Data";

	string boxes[nboxes] = {"1D_histograms/HT2_1B_0W","1D_histograms/HT2_1B_1W","1D_histograms/HT2_2B_0W","1D_histograms/HT2_2B_1W","1D_histograms/HT2_2B_2W","1D_histograms/HT2_2B_3W"};


	//HT cuts in different boxes
	//int bincuts[nboxes] = {4,4,5,4,1,1};
	int bincuts[nboxes] = {0,0,0,0,0,0};
	if(bincuts[0] != 0) postfix = postfix+"_HTCUTS";
	

	///////////////////////////////////////
	// get data histograms for all boxes //
	///////////////////////////////////////
	TH1F * Data_boxes[nboxes];

	for(int box_i = 0; box_i < nboxes; box_i++){
		string histoname_data = boxes[box_i].c_str();
		histoname_data = histoname_data + data.c_str();
		Data_boxes[box_i] = (TH1F*) file->Get(histoname_data.c_str()); 
	}
	cout << "retrieved histograms for the data for all boxes" << endl;

	/////////////////////////////////////////////////////////////
	// get signal histograms for all boxes and all mass points //
	/////////////////////////////////////////////////////////////
	TH1F * Tprime_boxes_masses[nboxes][nmasses], * STprime_boxes_masses[nboxes][nmasses], * Bprime_boxes_masses[nboxes][nmasses], * SBprime_boxes_masses[nboxes][nmasses];
	TH1F * JESPlus_Tprime_boxes_masses[nboxes][nmasses], * JESPlus_STprime_boxes_masses[nboxes][nmasses], * JESPlus_Bprime_boxes_masses[nboxes][nmasses], * JESPlus_SBprime_boxes_masses[nboxes][nmasses];
	TH1F * JESMinus_Tprime_boxes_masses[nboxes][nmasses], * JESMinus_STprime_boxes_masses[nboxes][nmasses], * JESMinus_Bprime_boxes_masses[nboxes][nmasses], * JESMinus_SBprime_boxes_masses[nboxes][nmasses];
	TH1F * JERPlus_Tprime_boxes_masses[nboxes][nmasses], * JERPlus_STprime_boxes_masses[nboxes][nmasses], * JERPlus_Bprime_boxes_masses[nboxes][nmasses], * JERPlus_SBprime_boxes_masses[nboxes][nmasses];
	TH1F * JERMinus_Tprime_boxes_masses[nboxes][nmasses], * JERMinus_STprime_boxes_masses[nboxes][nmasses], * JERMinus_Bprime_boxes_masses[nboxes][nmasses], * JERMinus_SBprime_boxes_masses[nboxes][nmasses];

	for(int box_i = 0; box_i < nboxes; box_i++){

		for(int mass_j = 0; mass_j < nmasses; mass_j++){
			string histoname_tprime = boxes[box_i].c_str();
			histoname_tprime = histoname_tprime + tprime[mass_j].c_str();
			string histoname_stprime = boxes[box_i].c_str();
			histoname_stprime = histoname_stprime + stprime[mass_j].c_str();
			string histoname_bprime = boxes[box_i].c_str();
			histoname_bprime = histoname_bprime + bprime[mass_j].c_str();
			string histoname_sbprime = boxes[box_i].c_str();
			histoname_sbprime = histoname_sbprime + sbprime[mass_j].c_str();
			Tprime_boxes_masses[box_i][mass_j] = (TH1F*) file->Get(histoname_tprime.c_str()); 
			STprime_boxes_masses[box_i][mass_j] = (TH1F*) file->Get(histoname_stprime.c_str());
			Bprime_boxes_masses[box_i][mass_j] = (TH1F*) file->Get(histoname_bprime.c_str());
			SBprime_boxes_masses[box_i][mass_j] = (TH1F*) file->Get(histoname_sbprime.c_str());
			if(doSystematics && doJESUnc){
				JESPlus_Tprime_boxes_masses[box_i][mass_j] = (TH1F*) JESPlus_file->Get(histoname_tprime.c_str()); 
				JESPlus_STprime_boxes_masses[box_i][mass_j] = (TH1F*) JESPlus_file->Get(histoname_stprime.c_str());
				JESPlus_Bprime_boxes_masses[box_i][mass_j] = (TH1F*) JESPlus_file->Get(histoname_bprime.c_str());
				JESPlus_SBprime_boxes_masses[box_i][mass_j] = (TH1F*) JESPlus_file->Get(histoname_sbprime.c_str());
				JESMinus_Tprime_boxes_masses[box_i][mass_j] = (TH1F*) JESMinus_file->Get(histoname_tprime.c_str()); 
				JESMinus_STprime_boxes_masses[box_i][mass_j] = (TH1F*) JESMinus_file->Get(histoname_stprime.c_str());
				JESMinus_Bprime_boxes_masses[box_i][mass_j] = (TH1F*) JESMinus_file->Get(histoname_bprime.c_str());
				JESMinus_SBprime_boxes_masses[box_i][mass_j] = (TH1F*) JESMinus_file->Get(histoname_sbprime.c_str());
			}
			if(doSystematics && doJERUnc){
				JERPlus_Tprime_boxes_masses[box_i][mass_j] = (TH1F*) JERPlus_file->Get(histoname_tprime.c_str()); 
				JERPlus_STprime_boxes_masses[box_i][mass_j] = (TH1F*) JERPlus_file->Get(histoname_stprime.c_str());
				JERPlus_Bprime_boxes_masses[box_i][mass_j] = (TH1F*) JERPlus_file->Get(histoname_bprime.c_str());
				JERPlus_SBprime_boxes_masses[box_i][mass_j] = (TH1F*) JERPlus_file->Get(histoname_sbprime.c_str());
				JERMinus_Tprime_boxes_masses[box_i][mass_j] = (TH1F*) JERMinus_file->Get(histoname_tprime.c_str()); 
				JERMinus_STprime_boxes_masses[box_i][mass_j] = (TH1F*) JERMinus_file->Get(histoname_stprime.c_str());
				JERMinus_Bprime_boxes_masses[box_i][mass_j] = (TH1F*) JERMinus_file->Get(histoname_bprime.c_str());
				JERMinus_SBprime_boxes_masses[box_i][mass_j] = (TH1F*) JERMinus_file->Get(histoname_sbprime.c_str());
			}
		}
		
	}
	cout << "retrieved histograms for the signal processes for all boxes and masses" << endl;

	///////////////////////////////////////////////////
	// get histograms for SM processes for all boxes //
	///////////////////////////////////////////////////
	TH1F * SM_boxes_processes[nboxes][ndatasetnamesSM];
	TH1F * JESPlus_SM_boxes_processes[nboxes][ndatasetnamesSM];
	TH1F * JESMinus_SM_boxes_processes[nboxes][ndatasetnamesSM];
	TH1F * JERPlus_SM_boxes_processes[nboxes][ndatasetnamesSM];
	TH1F * JERMinus_SM_boxes_processes[nboxes][ndatasetnamesSM];
	
	for(int box_i = 0; box_i < nboxes; box_i++){
		
		for(int name_j= 0; name_j < ndatasetnamesSM; name_j++){
			string histoname = boxes[box_i].c_str();
			histoname = histoname + datasetnamesSM[name_j].c_str();
			SM_boxes_processes[box_i][name_j] = (TH1F*) file->Get(histoname.c_str()); 
			if(doSystematics && doJESUnc){
				JESPlus_SM_boxes_processes[box_i][name_j] = (TH1F*) JESPlus_file->Get(histoname.c_str()); 
				JESMinus_SM_boxes_processes[box_i][name_j] = (TH1F*) JESMinus_file->Get(histoname.c_str());
			} 
			if(doSystematics && doJERUnc){
				JERPlus_SM_boxes_processes[box_i][name_j] = (TH1F*) JERPlus_file->Get(histoname.c_str()); 
				JERMinus_SM_boxes_processes[box_i][name_j] = (TH1F*) JERMinus_file->Get(histoname.c_str());
			}
		}
	}
	cout << "retrieved histograms for SM processes for all boxes" << endl;


	/////////////////////////////////////////////////////////////////////////////////////////
	// combine the SM processes that do not depend on CKM_A  and the 3 single top channels //
	/////////////////////////////////////////////////////////////////////////////////////////
	TH1F * SM_boxes[nboxes];
	TH1F * SingleTop_boxes[nboxes];
	TH1F * JESPlus_SM_boxes[nboxes];
	TH1F * JESPlus_SingleTop_boxes[nboxes];
	TH1F * JESMinus_SM_boxes[nboxes];
	TH1F * JESMinus_SingleTop_boxes[nboxes];
	TH1F * JERPlus_SM_boxes[nboxes];
	TH1F * JERPlus_SingleTop_boxes[nboxes];
	TH1F * JERMinus_SM_boxes[nboxes];
	TH1F * JERMinus_SingleTop_boxes[nboxes];
	TH1F * XSTTJetsPlus_SM_boxes[nboxes]; //combine the SM processes with a raised ttjets XS
	TH1F * XSTTJetsMinus_SM_boxes[nboxes];//combine the SM processes with a lowered ttjets XS
	TH1F * XSWJetsPlus_SM_boxes[nboxes]; //combine the SM processes with a raised Wjets XS
	TH1F * XSWJetsMinus_SM_boxes[nboxes];//combine the SM processes with a lowered Wjets XS
	
	for(int box_i = 0; box_i < nboxes; box_i++){
		///////// combine SM processes that are not single top
		SM_boxes[box_i] = (TH1F*) SM_boxes_processes[box_i][0]->Clone();
		if(doSystematics && doJESUnc){
			JESPlus_SM_boxes[box_i] = (TH1F*) JESPlus_SM_boxes_processes[box_i][0]->Clone();
			JESMinus_SM_boxes[box_i] = (TH1F*) JESMinus_SM_boxes_processes[box_i][0]->Clone();
		}
		if(doSystematics && doJERUnc){
			JERPlus_SM_boxes[box_i] = (TH1F*) JERPlus_SM_boxes_processes[box_i][0]->Clone();
			JERMinus_SM_boxes[box_i] = (TH1F*) JERMinus_SM_boxes_processes[box_i][0]->Clone();
		}
		if(doSystematics && doXSTTJetsUnc){
			XSTTJetsPlus_SM_boxes[box_i] = (TH1F*) SM_boxes_processes[box_i][0]->Clone();
			XSTTJetsMinus_SM_boxes[box_i] = (TH1F*) SM_boxes_processes[box_i][0]->Clone();
		}
		if(doSystematics && doXSWJetsUnc){
			XSWJetsPlus_SM_boxes[box_i] = (TH1F*) SM_boxes_processes[box_i][0]->Clone();
			XSWJetsMinus_SM_boxes[box_i] = (TH1F*) SM_boxes_processes[box_i][0]->Clone();
		}
		
		
		//{"TTJets_SemiMuon","TTJets_Other","WJets","ZJets","ST_tChannel","ST_tWChannel","ST_sChannel"};
		for(int process_j = 1; process_j <ndatasetnamesSM -3; process_j++){ //7 is size of array with SMprocesses, the 3 last ones are single top processes
			SM_boxes[box_i] -> Add(SM_boxes_processes[box_i][process_j]);
			if(doSystematics && doJESUnc){
				JESPlus_SM_boxes[box_i] -> Add(JESPlus_SM_boxes_processes[box_i][process_j]);
				JESMinus_SM_boxes[box_i] -> Add(JESMinus_SM_boxes_processes[box_i][process_j]);
			}
			if(doSystematics && doJERUnc){
				JERPlus_SM_boxes[box_i] -> Add(JERPlus_SM_boxes_processes[box_i][process_j]);
				JERMinus_SM_boxes[box_i] -> Add(JERMinus_SM_boxes_processes[box_i][process_j]);
			}
			if(doSystematics && doXSTTJetsUnc){
				if(process_j == 1) {//ttother histo
					XSTTJetsPlus_SM_boxes[box_i]->Add(SM_boxes_processes[box_i][process_j]); //add ttother to ttsemimu
					XSTTJetsPlus_SM_boxes[box_i]->Scale(1.+XSTTJets);//scale the combination of ttsemimu and ttother
					XSTTJetsMinus_SM_boxes[box_i]->Add(SM_boxes_processes[box_i][process_j]); //add ttother to ttsemimu
					XSTTJetsMinus_SM_boxes[box_i]->Scale(1.-XSTTJets);//scale the combination of ttsemimu and ttother
				}
				XSTTJetsPlus_SM_boxes[box_i]->Add(SM_boxes_processes[box_i][process_j]); 
				XSTTJetsMinus_SM_boxes[box_i]->Add(SM_boxes_processes[box_i][process_j]); 
			}
			if(doSystematics && doXSWJetsUnc){
				if(process_j == 2) { //WJets histo
					TH1F * WJetsHistoPlus = (TH1F*) SM_boxes_processes[box_i][process_j]->Clone();
					WJetsHistoPlus->Scale(1.+XSWJets);
					XSWJetsPlus_SM_boxes[box_i]->Add(WJetsHistoPlus); //add WJets
					TH1F * WJetsHistoMinus = (TH1F*) SM_boxes_processes[box_i][process_j]->Clone();
					WJetsHistoMinus->Scale(1.-XSWJets);
					XSWJetsMinus_SM_boxes[box_i]->Add(WJetsHistoMinus); //add WJets
				}
				XSWJetsPlus_SM_boxes[box_i]->Add(SM_boxes_processes[box_i][process_j]); 
				XSWJetsMinus_SM_boxes[box_i]->Add(SM_boxes_processes[box_i][process_j]); 
			}
			
		}
		
		
		
	cout << "will combine single top processes now" << endl;		
		
		
		///////// combine single top processes
		SingleTop_boxes[box_i] = (TH1F*) SM_boxes_processes[box_i][ndatasetnamesSM-1]->Clone();
		if(doSystematics && doJESUnc){
			JESPlus_SingleTop_boxes[box_i] = (TH1F*) JESPlus_SM_boxes_processes[box_i][ndatasetnamesSM-1]->Clone();
			JESMinus_SingleTop_boxes[box_i] = (TH1F*) JESMinus_SM_boxes_processes[box_i][ndatasetnamesSM-1]->Clone();
		}
		if(doSystematics && doJERUnc){
			JERPlus_SingleTop_boxes[box_i] = (TH1F*) JERPlus_SM_boxes_processes[box_i][ndatasetnamesSM-1]->Clone();
			JERMinus_SingleTop_boxes[box_i] = (TH1F*) JERMinus_SM_boxes_processes[box_i][ndatasetnamesSM-1]->Clone();
		}
			
		//not needed for XSTTJetsUnc or XSWJetsUnc to add these again, just use the standard SingleTop boxes
		for(int process_j = ndatasetnamesSM-1-1; process_j > ndatasetnamesSM-1-3; process_j--){// start counting from 6, because array filled from 0->6
			SingleTop_boxes[box_i] -> Add(SM_boxes_processes[box_i][process_j]);
			if(doSystematics && doJESUnc){
				JESPlus_SingleTop_boxes[box_i] -> Add(JESPlus_SM_boxes_processes[box_i][process_j]);
				JESMinus_SingleTop_boxes[box_i] -> Add(JESMinus_SM_boxes_processes[box_i][process_j]);
			}
			if(doSystematics && doJERUnc){
				JERPlus_SingleTop_boxes[box_i] -> Add(JERPlus_SM_boxes_processes[box_i][process_j]);
				JERMinus_SingleTop_boxes[box_i] -> Add(JERMinus_SM_boxes_processes[box_i][process_j]);
			}
		}
	}
	cout << "non-EWK SM processes are combined into a single histogram" << endl;





TH1F * HIST_Q_value_NEWPHYS[nCKM_A_values+2][nmasses], * HIST_Q_value_BACKGROUND[nCKM_A_values+2][nmasses];
string masses[nmasses] = {"350","400","450","500","550","600"};
float CLs_median[nCKM_A_values+2][nmasses];
float CLs_16th[nCKM_A_values+2][nmasses];
float CLs_84th[nCKM_A_values+2][nmasses];
float CLs_data[nCKM_A_values+2][nmasses];

for(int ckm_a = 0; ckm_a<=nCKM_A_values+1; ckm_a++){
	for(int mass_j = 0; mass_j < nmasses; mass_j++){
		string histname_qvalues_NP = "HIST_Q_value_NEWPHYS_";	
		string histname_qvalues_BKG = "HIST_Q_value_BACKGROUND_";	
		histname_qvalues_NP = histname_qvalues_NP+masses[mass_j].c_str();
		histname_qvalues_BKG = histname_qvalues_BKG+masses[mass_j].c_str();
		histname_qvalues_NP = histname_qvalues_NP+"_"+IntToStr(ckm_a).c_str();
		histname_qvalues_BKG = histname_qvalues_BKG+"_"+IntToStr(ckm_a).c_str();
		HIST_Q_value_NEWPHYS[ckm_a][mass_j] = new TH1F(histname_qvalues_NP.c_str(),histname_qvalues_NP.c_str(),1000,-50,50);
		HIST_Q_value_BACKGROUND[ckm_a][mass_j] = new TH1F(histname_qvalues_BKG.c_str(),histname_qvalues_BKG.c_str(),1000,-50,50);
	}
	/////////////////////////////////////////////////////
	// change the cross sections of the EWK processes //	
	////////////////////////////////////////////////////		
	float CKM_A, CKM_B;
	if(ckm_a==nCKM_A_values+1){// special case, because non-unitary matrix: maximal cross section for single top and single t' / b'
		CKM_A = 1;
		CKM_B = 1;
	}else{
		CKM_A = 1 - ((float) ckm_a / (float) nCKM_A_values);
		CKM_B = 1 - CKM_A;		
	}
	//	ckm_a == 0 -> single top XS maximal, single t' / b' XS = 0
	//	ckm_a == CKM_A_values -> single top XS = 0, single t' / b' XS maximal
	//	ckm_a == CKM_A_values+1 -> single top XS maximal, single t' / b' XS maximal
	
	TH1F * SingleTop_boxes_clone[nboxes];
	TH1F * STprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * SBprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JESPlus_SingleTop_boxes_clone[nboxes];
	TH1F * JESPlus_STprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JESPlus_SBprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JESMinus_SingleTop_boxes_clone[nboxes];
	TH1F * JESMinus_STprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JESMinus_SBprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JERPlus_SingleTop_boxes_clone[nboxes];
	TH1F * JERPlus_STprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JERPlus_SBprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JERMinus_SingleTop_boxes_clone[nboxes];
	TH1F * JERMinus_STprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JERMinus_SBprime_boxes_masses_clone[nboxes][nmasses];
	// for single top, multiply the cross section (scale the histogram) by CKM_A
	for(int box_i = 0; box_i < nboxes; box_i++){
		SingleTop_boxes_clone[box_i] = (TH1F*) SingleTop_boxes[box_i]->Clone();
		SingleTop_boxes_clone[box_i] -> Scale(CKM_A);
		if(doSystematics && doJESUnc){
			JESPlus_SingleTop_boxes_clone[box_i] = (TH1F*) JESPlus_SingleTop_boxes[box_i]->Clone();
			JESPlus_SingleTop_boxes_clone[box_i] -> Scale(CKM_A);
			JESMinus_SingleTop_boxes_clone[box_i] = (TH1F*) JESMinus_SingleTop_boxes[box_i]->Clone();
			JESMinus_SingleTop_boxes_clone[box_i] -> Scale(CKM_A);
		}
		if(doSystematics && doJERUnc){
			JERPlus_SingleTop_boxes_clone[box_i] = (TH1F*) JERPlus_SingleTop_boxes[box_i]->Clone();
			JERPlus_SingleTop_boxes_clone[box_i] -> Scale(CKM_A);
			JERMinus_SingleTop_boxes_clone[box_i] = (TH1F*) JERMinus_SingleTop_boxes[box_i]->Clone();
			JERMinus_SingleTop_boxes_clone[box_i] -> Scale(CKM_A);
		}
	}
	
	// for stprime and sbprime, multiply the cross section (scale the histogram) by (1-CKM_A)
	for(int box_i = 0; box_i < nboxes; box_i++){
		for(int mass_j = 0; mass_j < nmasses; mass_j++){
			STprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) STprime_boxes_masses[box_i][mass_j]->Clone();
			SBprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) SBprime_boxes_masses[box_i][mass_j]->Clone();
			STprime_boxes_masses_clone[box_i][mass_j] -> Scale(CKM_B);
			SBprime_boxes_masses_clone[box_i][mass_j] -> Scale(CKM_B);
			if(doSystematics && doJESUnc){
				JESPlus_STprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JESPlus_STprime_boxes_masses[box_i][mass_j]->Clone();
				JESPlus_SBprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JESPlus_SBprime_boxes_masses[box_i][mass_j]->Clone();
				JESPlus_STprime_boxes_masses_clone[box_i][mass_j] -> Scale(CKM_B);
				JESPlus_SBprime_boxes_masses_clone[box_i][mass_j] -> Scale(CKM_B);
				JESMinus_STprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JESMinus_STprime_boxes_masses[box_i][mass_j]->Clone();
				JESMinus_SBprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JESMinus_SBprime_boxes_masses[box_i][mass_j]->Clone();
				JESMinus_STprime_boxes_masses_clone[box_i][mass_j] -> Scale(CKM_B);
				JESMinus_SBprime_boxes_masses_clone[box_i][mass_j] -> Scale(CKM_B);
			}
			if(doSystematics && doJERUnc){
				JERPlus_STprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JERPlus_STprime_boxes_masses[box_i][mass_j]->Clone();
				JERPlus_SBprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JERPlus_SBprime_boxes_masses[box_i][mass_j]->Clone();
				JERPlus_STprime_boxes_masses_clone[box_i][mass_j] -> Scale(CKM_B);
				JERPlus_SBprime_boxes_masses_clone[box_i][mass_j] -> Scale(CKM_B);
				JERMinus_STprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JERMinus_STprime_boxes_masses[box_i][mass_j]->Clone();
				JERMinus_SBprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JERMinus_SBprime_boxes_masses[box_i][mass_j]->Clone();
				JERMinus_STprime_boxes_masses_clone[box_i][mass_j] -> Scale(CKM_B);
				JERMinus_SBprime_boxes_masses_clone[box_i][mass_j] -> Scale(CKM_B);
			}
		}
	}
	cout << "changed the cross sections of the EWK processes according to CKM_A = " << CKM_A << " and 1-CKM_A =  " << CKM_B << endl;

	///////////////////////////////////////////////////////////////
	// combine the signal processes and the background processes //
	///////////////////////////////////////////////////////////////
	TH1F * SIGNAL_boxes_masses[nboxes][nmasses];
	TH1F * BACKGROUND_boxes[nboxes];
	TH1F * Tprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * Bprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JESPlus_SIGNAL_boxes_masses[nboxes][nmasses];	
	TH1F * JESPlus_BACKGROUND_boxes[nboxes];	
	TH1F * JESPlus_Tprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JESPlus_Bprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JESMinus_SIGNAL_boxes_masses[nboxes][nmasses];	
	TH1F * JESMinus_BACKGROUND_boxes[nboxes];	
	TH1F * JESMinus_Tprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JESMinus_Bprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JERPlus_SIGNAL_boxes_masses[nboxes][nmasses];	
	TH1F * JERPlus_BACKGROUND_boxes[nboxes];	
	TH1F * JERPlus_Tprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JERPlus_Bprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JERMinus_SIGNAL_boxes_masses[nboxes][nmasses];	
	TH1F * JERMinus_BACKGROUND_boxes[nboxes];	
	TH1F * JERMinus_Tprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * JERMinus_Bprime_boxes_masses_clone[nboxes][nmasses];
	TH1F * XSTTJetsPlus_BACKGROUND_boxes[nboxes];
	TH1F * XSTTJetsMinus_BACKGROUND_boxes[nboxes];
	TH1F * XSWJetsPlus_BACKGROUND_boxes[nboxes];
	TH1F * XSWJetsMinus_BACKGROUND_boxes[nboxes];
	
	for(int box_i = 0; box_i < nboxes; box_i++){
		BACKGROUND_boxes[box_i] = (TH1F*) SM_boxes[box_i]->Clone();
		BACKGROUND_boxes[box_i] -> Add(SingleTop_boxes_clone[box_i]);
		if(doSystematics && doJESUnc){
			JESPlus_BACKGROUND_boxes[box_i] = (TH1F*) JESPlus_SM_boxes[box_i]->Clone();
			JESPlus_BACKGROUND_boxes[box_i] -> Add(JESPlus_SingleTop_boxes_clone[box_i]);
			JESMinus_BACKGROUND_boxes[box_i] = (TH1F*) JESMinus_SM_boxes[box_i]->Clone();
			JESMinus_BACKGROUND_boxes[box_i] -> Add(JESMinus_SingleTop_boxes_clone[box_i]);
		}	
		if(doSystematics && doJERUnc){
			JERPlus_BACKGROUND_boxes[box_i] = (TH1F*) JERPlus_SM_boxes[box_i]->Clone();
			JERPlus_BACKGROUND_boxes[box_i] -> Add(JERPlus_SingleTop_boxes_clone[box_i]);
			JERMinus_BACKGROUND_boxes[box_i] = (TH1F*) JERMinus_SM_boxes[box_i]->Clone();
			JERMinus_BACKGROUND_boxes[box_i] -> Add(JERMinus_SingleTop_boxes_clone[box_i]);
		}	
		if(doSystematics && doXSTTJetsUnc){
			XSTTJetsPlus_BACKGROUND_boxes[box_i] = (TH1F*) XSTTJetsPlus_SM_boxes[box_i]->Clone();
			XSTTJetsPlus_BACKGROUND_boxes[box_i] -> Add(SingleTop_boxes_clone[box_i]);
			XSTTJetsMinus_BACKGROUND_boxes[box_i] = (TH1F*) XSTTJetsMinus_SM_boxes[box_i]->Clone();		
			XSTTJetsMinus_BACKGROUND_boxes[box_i] -> Add(SingleTop_boxes_clone[box_i]);		
		}
		if(doSystematics && doXSWJetsUnc){
			XSWJetsPlus_BACKGROUND_boxes[box_i] = (TH1F*) XSWJetsPlus_SM_boxes[box_i]->Clone();
			XSWJetsPlus_BACKGROUND_boxes[box_i] -> Add(SingleTop_boxes_clone[box_i]);
			XSWJetsMinus_BACKGROUND_boxes[box_i] = (TH1F*) XSWJetsMinus_SM_boxes[box_i]->Clone();		
			XSWJetsMinus_BACKGROUND_boxes[box_i] -> Add(SingleTop_boxes_clone[box_i]);		
		}
		for(int mass_j = 0; mass_j < nmasses; mass_j++){
			Tprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) Tprime_boxes_masses[box_i][mass_j]->Clone();
			Bprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) Bprime_boxes_masses[box_i][mass_j]->Clone();
			SIGNAL_boxes_masses[box_i][mass_j] = (TH1F*) Tprime_boxes_masses_clone[box_i][mass_j]->Clone();
			SIGNAL_boxes_masses[box_i][mass_j] -> Add(Bprime_boxes_masses_clone[box_i][mass_j]);
			SIGNAL_boxes_masses[box_i][mass_j] -> Add(STprime_boxes_masses_clone[box_i][mass_j]);
			SIGNAL_boxes_masses[box_i][mass_j] -> Add(SBprime_boxes_masses_clone[box_i][mass_j]);
			if(doSystematics && doJESUnc){
				JESPlus_Tprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JESPlus_Tprime_boxes_masses[box_i][mass_j]->Clone();
				JESPlus_Bprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JESPlus_Bprime_boxes_masses[box_i][mass_j]->Clone();
				JESPlus_SIGNAL_boxes_masses[box_i][mass_j] = (TH1F*) JESPlus_Tprime_boxes_masses_clone[box_i][mass_j]->Clone();
				JESPlus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JESPlus_Bprime_boxes_masses_clone[box_i][mass_j]);
				JESPlus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JESPlus_STprime_boxes_masses_clone[box_i][mass_j]);
				JESPlus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JESPlus_SBprime_boxes_masses_clone[box_i][mass_j]);
				JESMinus_Tprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JESMinus_Tprime_boxes_masses[box_i][mass_j]->Clone();
				JESMinus_Bprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JESMinus_Bprime_boxes_masses[box_i][mass_j]->Clone();
				JESMinus_SIGNAL_boxes_masses[box_i][mass_j] = (TH1F*) JESMinus_Tprime_boxes_masses_clone[box_i][mass_j]->Clone();
				JESMinus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JESMinus_Bprime_boxes_masses_clone[box_i][mass_j]);
				JESMinus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JESMinus_STprime_boxes_masses_clone[box_i][mass_j]);
				JESMinus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JESMinus_SBprime_boxes_masses_clone[box_i][mass_j]);
			}
			if(doSystematics && doJERUnc){
				JERPlus_Tprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JERPlus_Tprime_boxes_masses[box_i][mass_j]->Clone();
				JERPlus_Bprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JERPlus_Bprime_boxes_masses[box_i][mass_j]->Clone();
				JERPlus_SIGNAL_boxes_masses[box_i][mass_j] = (TH1F*) JERPlus_Tprime_boxes_masses_clone[box_i][mass_j]->Clone();
				JERPlus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JERPlus_Bprime_boxes_masses_clone[box_i][mass_j]);
				JERPlus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JERPlus_STprime_boxes_masses_clone[box_i][mass_j]);
				JERPlus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JERPlus_SBprime_boxes_masses_clone[box_i][mass_j]);
				JERMinus_Tprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JERMinus_Tprime_boxes_masses[box_i][mass_j]->Clone();
				JERMinus_Bprime_boxes_masses_clone[box_i][mass_j] = (TH1F*) JERMinus_Bprime_boxes_masses[box_i][mass_j]->Clone();
				JERMinus_SIGNAL_boxes_masses[box_i][mass_j] = (TH1F*) JERMinus_Tprime_boxes_masses_clone[box_i][mass_j]->Clone();
				JERMinus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JERMinus_Bprime_boxes_masses_clone[box_i][mass_j]);
				JERMinus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JERMinus_STprime_boxes_masses_clone[box_i][mass_j]);
				JERMinus_SIGNAL_boxes_masses[box_i][mass_j] -> Add(JERMinus_SBprime_boxes_masses_clone[box_i][mass_j]);
			}
		}
	}
	cout << "different signal processes and background processes are combined into 2 separate histograms" << endl;

	////////////////////////////////////////////////////
	// Q values for data considering different masses //
	////////////////////////////////////////////////////
	float Q_value_masses[nmasses];

	for(int mass_j = 0; mass_j < nmasses; mass_j++){
		Q_value_masses[mass_j]= 0;
		for(int box_i = 0; box_i < nboxes; box_i++){
//			for(int bin_k = 0; bin_k < Data_boxes[box_i] -> GetNbinsX()+1; bin_k++){	
			for(int bin_k = bincuts[box_i]; bin_k < Data_boxes[box_i] -> GetNbinsX()+1; bin_k++){	
				Q_value_masses[mass_j]= Q_value_masses[mass_j] 
					- log(
					pow((SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k) + BACKGROUND_boxes[box_i]->GetBinContent(bin_k))
					/BACKGROUND_boxes[box_i]->GetBinContent(bin_k),Data_boxes[box_i]->GetBinContent(bin_k) )
					)
					+ SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k);
			}
		}
		cout << "Q_value for data with Vtb = sqrt("<< CKM_A << ") and Vt'b = Vb't = sqrt(" << CKM_B << ") and for mt' = mb' = " << 350+(50*mass_j) << " GeV  is : " << Q_value_masses[mass_j] << endl;	
	}
	
	///////////////////////////////////////////////////////////////////////
	// Q values for SM and SM+NP hypotheses considering different masses //
	///////////////////////////////////////////////////////////////////////
	float array_Qvalue_NEWPHYS[nmasses][nPE];
	float array_Qvalue_BACKGROUND[nmasses][nPE];

	TRandom3 * lumismearing = new TRandom3(0); 
	TRandom3 * randomJESUnc = new TRandom3(0); 
	TRandom3 * randomJERUnc = new TRandom3(0); 
	TRandom3 * randomXSTTJetsUnc = new TRandom3(0); 
	TRandom3 * randomXSWJetsUnc = new TRandom3(0); 
	TRandom3 * rand_BKG = new TRandom3(0);
	TRandom3 * rand_NP = new TRandom3(0);
	
	float Q_value_NP, Q_value_BKG;
	float pseudodata_BKG, pseudodata_NP;
	for(int mass_j = 0; mass_j < nmasses; mass_j++){
		
		cout << endl;
		cout << "For mt' = mb' = "<< 350+(50*mass_j) << " GeV" <<endl;

		for(int PE=0; PE< nPE; PE++){
			double lumisyst = lumismearing->Gaus(0,lumisyst); // Normal distribution with mean = 0 and sigma = luminosity uncertainty
			double JESUnc = randomJESUnc->Gaus(0,1);
			double JERUnc = randomJERUnc->Gaus(0,1);
			double XSTTUnc = randomXSTTJetsUnc->Gaus(0,1);
			double XSWUnc = randomXSWJetsUnc->Gaus(0,1);
			
			//cout << "generated random numbers if necessary " << endl;
			
			if(PE%1000 == 0)
			  std::cout<<"			performing the "<<PE<<"th pseudo-experiment" <<flush<<"\r";
			Q_value_NP = 0;
			Q_value_BKG = 0;
			
			for(int box_i = 0; box_i < nboxes; box_i++){
		
//				for(int bin_k = 0; bin_k < BACKGROUND_boxes[box_i] -> GetNbinsX()+1; bin_k++){	
				for(int bin_k = bincuts[box_i]; bin_k < BACKGROUND_boxes[box_i] -> GetNbinsX()+1; bin_k++){	
		
					pseudodata_BKG = 0;
					pseudodata_NP = 0;

	pseudodata_BKG = BACKGROUND_boxes[box_i]->GetBinContent(bin_k);
	pseudodata_NP = SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k) + BACKGROUND_boxes[box_i]->GetBinContent(bin_k);

if(doSystematics){
	if(doLumiUnc){
		if(lumismearing>0){//scale all bincontents up
			pseudodata_BKG = pseudodata_BKG + BACKGROUND_boxes[box_i]->GetBinContent(bin_k)*fabs(lumisyst);
			pseudodata_NP = pseudodata_NP + (SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k)+BACKGROUND_boxes[box_i]->GetBinContent(bin_k))*fabs(lumisyst);							}else{//scale all bincontents down
			pseudodata_BKG = pseudodata_BKG - BACKGROUND_boxes[box_i]->GetBinContent(bin_k)*fabs(lumisyst);
			pseudodata_NP = pseudodata_NP - (SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k)+BACKGROUND_boxes[box_i]->GetBinContent(bin_k))*fabs(lumisyst);
		}
	}
	
	if(doJESUnc){
		float JESDeltaBinContent_BKG = 0;
		float JESDeltaBinContent_NP = 0;
		if(randomJESUnc>0){//use JESPlus histograms
			JESDeltaBinContent_BKG = JESPlus_BACKGROUND_boxes[box_i]->GetBinContent(bin_k)-BACKGROUND_boxes[box_i]->GetBinContent(bin_k);
			JESDeltaBinContent_NP = JESDeltaBinContent_BKG+JESPlus_SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k)-SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k);
		}else{//use JESMinus histograms
			JESDeltaBinContent_BKG = JESMinus_BACKGROUND_boxes[box_i]->GetBinContent(bin_k)-BACKGROUND_boxes[box_i]->GetBinContent(bin_k);
			JESDeltaBinContent_NP = JESDeltaBinContent_BKG+JESMinus_SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k)-SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k);	
		}

		pseudodata_BKG = pseudodata_BKG + (JESDeltaBinContent_BKG*fabs(JESUnc));
		pseudodata_NP = pseudodata_NP + (JESDeltaBinContent_NP*fabs(JESUnc));

		if(pseudodata_BKG < 0) pseudodata_BKG = 0;
		if(pseudodata_NP < 0) pseudodata_NP = 0;		
	}
	
	if(doJERUnc){
		float JERDeltaBinContent_BKG = 0;
		float JERDeltaBinContent_NP = 0;
		if(randomJERUnc>0){//use JERPlus histograms
			JERDeltaBinContent_BKG = JERPlus_BACKGROUND_boxes[box_i]->GetBinContent(bin_k)-BACKGROUND_boxes[box_i]->GetBinContent(bin_k);
			JERDeltaBinContent_NP = JERDeltaBinContent_BKG+JERPlus_SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k)-SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k);
		}else{//use JERMinus histograms
			JERDeltaBinContent_BKG = JERMinus_BACKGROUND_boxes[box_i]->GetBinContent(bin_k)-BACKGROUND_boxes[box_i]->GetBinContent(bin_k);
			JERDeltaBinContent_NP = JERDeltaBinContent_BKG+JERMinus_SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k)-SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k);	
		}

		pseudodata_BKG = pseudodata_BKG + (JERDeltaBinContent_BKG*fabs(JERUnc));
		pseudodata_NP = pseudodata_NP + (JERDeltaBinContent_NP*fabs(JERUnc));

		if(pseudodata_BKG < 0) pseudodata_BKG = 0;
		if(pseudodata_NP < 0) pseudodata_NP = 0;		
	}
	
	if(doXSTTJetsUnc){
		float XSTTJetsDeltaBinContent_BKG = 0;
		float XSTTJetsDeltaBinContent_NP = 0;
		if(randomXSTTJetsUnc>0){//use XSTTJetsPlus histograms
			XSTTJetsDeltaBinContent_BKG = XSTTJetsPlus_BACKGROUND_boxes[box_i]->GetBinContent(bin_k)-BACKGROUND_boxes[box_i]->GetBinContent(bin_k);
			XSTTJetsDeltaBinContent_NP = XSTTJetsDeltaBinContent_BKG;//only the background changes, signal does not contribute!
		}else{//use XSTTJetsMinus histograms
			XSTTJetsDeltaBinContent_BKG = XSTTJetsMinus_BACKGROUND_boxes[box_i]->GetBinContent(bin_k)-BACKGROUND_boxes[box_i]->GetBinContent(bin_k);
			XSTTJetsDeltaBinContent_NP = XSTTJetsDeltaBinContent_BKG;	
		}

		pseudodata_BKG = pseudodata_BKG + (XSTTJetsDeltaBinContent_BKG*fabs(XSTTUnc));
		pseudodata_NP = pseudodata_NP + (XSTTJetsDeltaBinContent_NP*fabs(XSTTUnc));//only the background changes, signal does not contribute!

		if(pseudodata_BKG < 0) pseudodata_BKG = 0;
		if(pseudodata_NP < 0) pseudodata_NP = 0;		
	
	}

	if(doXSWJetsUnc){
		float XSWJetsDeltaBinContent_BKG = 0;
		float XSWJetsDeltaBinContent_NP = 0;
		if(randomXSWJetsUnc>0){//use XSWJetsPlus histograms
			XSWJetsDeltaBinContent_BKG = XSWJetsPlus_BACKGROUND_boxes[box_i]->GetBinContent(bin_k)-BACKGROUND_boxes[box_i]->GetBinContent(bin_k);
			XSWJetsDeltaBinContent_NP = XSWJetsDeltaBinContent_BKG;//only the background changes, signal does not contribute!
		}else{//use XSWJetsMinus histograms
			XSWJetsDeltaBinContent_BKG = XSWJetsMinus_BACKGROUND_boxes[box_i]->GetBinContent(bin_k)-BACKGROUND_boxes[box_i]->GetBinContent(bin_k);
			XSWJetsDeltaBinContent_NP = XSWJetsDeltaBinContent_BKG;	
		}

		pseudodata_BKG = pseudodata_BKG + (XSWJetsDeltaBinContent_BKG*fabs(XSWUnc));
		pseudodata_NP = pseudodata_NP + (XSWJetsDeltaBinContent_NP*fabs(XSWUnc));//only the background changes, signal does not contribute!

		if(pseudodata_BKG < 0) pseudodata_BKG = 0;
		if(pseudodata_NP < 0) pseudodata_NP = 0;		
	
	}

	
}
	//cout << "pseudodata_BKG " << pseudodata_BKG << endl;
	//cout << "pseudodata_NP " << pseudodata_NP << endl;

					Q_value_BKG = Q_value_BKG
						- log(
						pow((SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k) + BACKGROUND_boxes[box_i]->GetBinContent(bin_k)) 
						/BACKGROUND_boxes[box_i]->GetBinContent(bin_k), rand_BKG->Poisson(pseudodata_BKG) )
						)
						+ SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k);
					
					Q_value_NP = Q_value_NP
						- log(
						pow((SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k) + BACKGROUND_boxes[box_i]->GetBinContent(bin_k)) 
						/ BACKGROUND_boxes[box_i]->GetBinContent(bin_k), rand_NP->Poisson(pseudodata_NP ) )
						)
						+ SIGNAL_boxes_masses[box_i][mass_j]->GetBinContent(bin_k);


					
				} // bin loop
		
			} // box loop		
			
			//cout << "for pseudo-experiment " << PE << " Q_value in case of SM-only hypothesis is " << Q_value_BKG << " and in case of SM+NP " << Q_value_NP <<endl;
			array_Qvalue_NEWPHYS[mass_j][PE] = Q_value_NP;
			array_Qvalue_BACKGROUND[mass_j][PE] = Q_value_BKG;
			HIST_Q_value_NEWPHYS[ckm_a][mass_j]->Fill(Q_value_NP);
			HIST_Q_value_BACKGROUND[ckm_a][mass_j]->Fill(Q_value_BKG);

//			delete lumismearing;
//			delete randomJESUnc;
		}// PE loop

		fout->cd();
		HIST_Q_value_NEWPHYS[ckm_a][mass_j]->Write();
		HIST_Q_value_BACKGROUND[ckm_a][mass_j]->Write();

	}// mass loop
	delete rand_BKG;
	delete rand_NP;
	delete lumismearing;
	delete randomJESUnc;
	delete randomJERUnc;


//////////////////////////
// Calculate CLs limits //
//////////////////////////
	cout << endl;
	cout << "calculating the CLs limits" << endl;
	for(int mass_j = 0; mass_j < nmasses; mass_j++){
		cout << "expected CLs for "<< 350+(mass_j*50) << " GeV: "<< endl;
		//median
		float Q_SMonlyMedian =0;
		Q_SMonlyMedian = GetProbaPercentile(array_Qvalue_BACKGROUND,50, mass_j, nPE);
		cout << "Q_SMonlyMedian: " << Q_SMonlyMedian << endl; 
		float CLbpluss_median = CalculateProbability(array_Qvalue_NEWPHYS, Q_SMonlyMedian, mass_j, nPE);
		cout << "CLbpluss_median: " << CLbpluss_median << endl; 
		float CLb_median = CalculateProbability(array_Qvalue_BACKGROUND, Q_SMonlyMedian, mass_j, nPE);
		cout << "CLb_median: " << CLb_median << endl; 
		CLs_median[ckm_a][mass_j] = CLbpluss_median / CLb_median;
		cout << "CLs_median: " << CLs_median[ckm_a][mass_j] << endl; 

	  	//16th percentile
	  	float Q_SMonly16thPerc = 0;
	  	Q_SMonly16thPerc = GetProbaPercentile(array_Qvalue_BACKGROUND,84, mass_j, nPE);
		  float CLbpluss_16th = CalculateProbability(array_Qvalue_NEWPHYS, Q_SMonly16thPerc, mass_j, nPE);
		  float CLb_16th = CalculateProbability(array_Qvalue_BACKGROUND, Q_SMonly16thPerc, mass_j, nPE);
	  	CLs_16th[ckm_a][mass_j] = CLbpluss_16th / CLb_16th;

	  	//84th percentile
	  	float Q_SMonly84thPerc = 0;
	  	Q_SMonly84thPerc = GetProbaPercentile(array_Qvalue_BACKGROUND,16, mass_j, nPE);
	  	float CLbpluss_84th = CalculateProbability( array_Qvalue_NEWPHYS, Q_SMonly84thPerc, mass_j, nPE);
	  	float CLb_84th = CalculateProbability(array_Qvalue_BACKGROUND, Q_SMonly84thPerc, mass_j, nPE);
	  	CLs_84th[ckm_a][mass_j] = CLbpluss_84th / CLb_84th;

		cout << "	expected CLs " << CLs_median[ckm_a][mass_j] << endl;
		cout << "		CLs 16th percentile " << CLs_16th[ckm_a][mass_j] << " and 84th percentile " << CLs_84th[ckm_a][mass_j] << endl;
	
		float CLbpluss_data = CalculateProbability(array_Qvalue_NEWPHYS,Q_value_masses[mass_j], mass_j, nPE);
		float CLb_data = CalculateProbability(array_Qvalue_BACKGROUND,Q_value_masses[mass_j], mass_j, nPE);
		CLs_data[ckm_a][mass_j] = CLbpluss_data / CLb_data;

		cout << "	observed CLs "<< CLs_data[ckm_a][mass_j] << endl;

	}// mass loop	

	
}//CKM_A loop	

fout->cd();
float massvalues[nmasses] = {350,400,450,500,550,600};

float mass_exp[nCKM_A_values+2][2], mass_obs[nCKM_A_values+2][2],mass_exp_16th[nCKM_A_values+2][2], mass_exp_84th[nCKM_A_values+2][2];
float CLs_exp[nCKM_A_values+2][2], CLs_obs[nCKM_A_values+2][2],CLs_exp_16th[nCKM_A_values+2][2],CLs_exp_84th[nCKM_A_values+2][2];
float masslimit_exp[nCKM_A_values+2], masslimit_obs[nCKM_A_values+2],masslimit_exp_16th[nCKM_A_values+2],masslimit_exp_84th[nCKM_A_values+2];

TGraph *CLs_median_versus_mass[nCKM_A_values+2];
TGraph *CLs_16th_versus_mass[nCKM_A_values+2];
TGraph *CLs_84th_versus_mass[nCKM_A_values+2];
TGraph *CLs_data_versus_mass[nCKM_A_values+2];

for(int ckm_a = 0; ckm_a<=nCKM_A_values+1; ckm_a++){
	//cout << "Vtb = sqrt("<< CKM_A << ") and Vt'b = Vb't = sqrt(" << CKM_B << "): " << endl;	
	mass_exp[ckm_a][0] = 0; mass_exp[ckm_a][1] = 0; CLs_exp[ckm_a][0] = 0; CLs_exp[ckm_a][1] = 0;
	mass_obs[ckm_a][0] = 0; mass_obs[ckm_a][1] = 0; CLs_obs[ckm_a][0] = 0; CLs_obs[ckm_a][1] = 0;
	mass_exp_16th[ckm_a][0] = 0; mass_exp_16th[ckm_a][1] = 0; CLs_exp_16th[ckm_a][0] = 0; CLs_exp_16th[ckm_a][1] = 0;
	mass_exp_84th[ckm_a][0] = 0; mass_exp_84th[ckm_a][1] = 0; CLs_exp_84th[ckm_a][0] = 0; CLs_exp_84th[ckm_a][1] = 0;
	masslimit_exp[ckm_a]= 0; masslimit_exp_84th[ckm_a]= 0; masslimit_exp_16th[ckm_a]= 0; masslimit_obs[ckm_a]= 0; 
	for(int i = 0 ; i < nmasses-1; i++){
		if(CLs_median[ckm_a][i]<0.05 && CLs_median[ckm_a][i+1]>0.05){
			cout << "	expected masslimit is between " << massvalues[i] << " and " << massvalues[i+1] << " GeV" << endl;
			mass_exp[ckm_a][0] = massvalues[i];
			mass_exp[ckm_a][1] = massvalues[i+1];
			CLs_exp[ckm_a][0] = CLs_median[ckm_a][i];
			CLs_exp[ckm_a][1] = CLs_median[ckm_a][i+1];
		}
		if(CLs_data[ckm_a][i]<0.05 && CLs_data[ckm_a][i+1]>0.05){
			cout << "	observed masslimit is between " << massvalues[i] << " and " << massvalues[i+1] << " GeV" << endl;
			mass_obs[ckm_a][0] = massvalues[i];
			mass_obs[ckm_a][1] = massvalues[i+1];
			CLs_obs[ckm_a][0] = CLs_data[ckm_a][i];
			CLs_obs[ckm_a][1] = CLs_data[ckm_a][i+1];
		}
		if(CLs_16th[ckm_a][i]<0.05 && CLs_16th[ckm_a][i+1]>0.05){
			cout << "	16th masslimit is between " << massvalues[i] << " and " << massvalues[i+1] << " GeV" << endl;
			mass_exp_16th[ckm_a][0] = massvalues[i];
			mass_exp_16th[ckm_a][1] = massvalues[i+1];
			CLs_exp_16th[ckm_a][0] = CLs_16th[ckm_a][i];
			CLs_exp_16th[ckm_a][1] = CLs_16th[ckm_a][i+1];
		}
		if(CLs_84th[ckm_a][i]<0.05 && CLs_84th[ckm_a][i+1]>0.05){
			cout << "	84th masslimit is between " << massvalues[i] << " and " << massvalues[i+1] << " GeV" << endl;
			mass_exp_84th[ckm_a][0] = massvalues[i];
			mass_exp_84th[ckm_a][1] = massvalues[i+1];
			CLs_exp_84th[ckm_a][0] = CLs_84th[ckm_a][i];
			CLs_exp_84th[ckm_a][1] = CLs_84th[ckm_a][i+1];
		
		}
		
	}
	if(mass_exp[ckm_a][0] == 0 || mass_exp[ckm_a][1] == 0){
		cout << "	expected masslimit is out of the scanned mass range!" << endl;
	}
	if(mass_exp_16th[ckm_a][0] == 0 || mass_exp_16th[ckm_a][1] == 0){
		cout << "	expected 16th masslimit is out of the scanned mass range!" << endl;
	}
	if(mass_exp_84th[ckm_a][0] == 0 || mass_exp_84th[ckm_a][1] == 0){
		cout << "	expected 84th masslimit is out of the scanned mass range!" << endl;
	}
	if(mass_obs[ckm_a][0] == 0 || mass_obs[ckm_a][1] == 0){
		cout << "	observed masslimit is out of the scanned mass range!" << endl;
	}

  	double CLs95CL = 0.05;

	TGraph mass_versus_CLs_exp(2,CLs_exp[ckm_a],mass_exp[ckm_a]);
	TF1 f1("f1","pol1");
  	if(mass_exp[ckm_a][0] != 0 && mass_exp[ckm_a][1] != 0) mass_versus_CLs_exp.Fit("f1","RQ0");    
	if(mass_exp[ckm_a][0] != 0 && mass_exp[ckm_a][1] != 0) masslimit_exp[ckm_a] = f1.Eval(CLs95CL);

	TGraph mass_versus_CLs_obs(2,CLs_obs[ckm_a],mass_obs[ckm_a]);
	TF1 f2("f2","pol1");
  	if(mass_obs[ckm_a][0] != 0 && mass_obs[ckm_a][1] != 0) mass_versus_CLs_obs.Fit("f2","RQ0");    
	if(mass_obs[ckm_a][0] != 0 && mass_obs[ckm_a][1] != 0) masslimit_obs[ckm_a] = f2.Eval(CLs95CL);

	TGraph mass_versus_CLs_exp_16th(2,CLs_exp_16th[ckm_a],mass_exp_16th[ckm_a]);
	TF1 f3("f3","pol1");
  	if(mass_exp_16th[ckm_a][0] != 0 && mass_exp_16th[ckm_a][1] != 0) mass_versus_CLs_exp_16th.Fit("f3","RQ0");    
	if(mass_exp_16th[ckm_a][0] != 0 && mass_exp_16th[ckm_a][1] != 0) masslimit_exp_16th[ckm_a] = f3.Eval(CLs95CL);

	TGraph mass_versus_CLs_exp_84th(2,CLs_exp_84th[ckm_a],mass_exp_84th[ckm_a]);
	TF1 f4("f4","pol1");
  	if(mass_exp_84th[ckm_a][0] != 0 && mass_exp_84th[ckm_a][1] != 0) mass_versus_CLs_exp_84th.Fit("f4","RQ0");    
	if(mass_exp_84th[ckm_a][0] != 0 && mass_exp_84th[ckm_a][1] != 0) masslimit_exp_84th[ckm_a] = f4.Eval(CLs95CL);

	cout << "masses below " << masslimit_exp[ckm_a] << " are expected to be excluded" << endl;
	cout << "masses below " << masslimit_obs[ckm_a] << " are observed to be excluded" << endl;
	cout << "uncertainty on expected masslimit (lower) " << masslimit_exp_16th[ckm_a] << endl;
	cout << "uncertainty on expected masslimit (upper) " << masslimit_exp_84th[ckm_a] << endl;

	cout << "putting plots in nice format" << endl;
		
	CLs_median_versus_mass[ckm_a] = new TGraph(nmasses,massvalues,CLs_median[ckm_a]);
	CLs_16th_versus_mass[ckm_a] = new TGraph(nmasses,massvalues,CLs_16th[ckm_a]);
	CLs_84th_versus_mass[ckm_a] = new TGraph(nmasses,massvalues,CLs_84th[ckm_a]);
	CLs_data_versus_mass[ckm_a] = new TGraph(nmasses,massvalues,CLs_data[ckm_a]);

	fout->cd();
	TString strCLs_median_versus_mass = "CLs_median_versus_mass_";
	strCLs_median_versus_mass += IntToStr(ckm_a).c_str();
	CLs_median_versus_mass[ckm_a]->SetTitle(strCLs_median_versus_mass);
	TString strCLs_16th_versus_mass = "CLs_16th_versus_mass_";
	strCLs_16th_versus_mass += IntToStr(ckm_a).c_str();
	CLs_16th_versus_mass[ckm_a]->SetTitle(strCLs_16th_versus_mass);
	TString strCLs_84th_versus_mass = "CLs_84th_versus_mass_";
	strCLs_84th_versus_mass += IntToStr(ckm_a).c_str();
	CLs_84th_versus_mass[ckm_a]->SetTitle(strCLs_84th_versus_mass);
	TString strCLs_data_versus_mass = "CLs_data_versus_mass_";
	strCLs_data_versus_mass += IntToStr(ckm_a).c_str();
	CLs_data_versus_mass[ckm_a]->SetTitle(strCLs_data_versus_mass);
	
	CLs_median_versus_mass[ckm_a]->Write();
	CLs_16th_versus_mass[ckm_a]->Write();
	CLs_84th_versus_mass[ckm_a]->Write();
	CLs_data_versus_mass[ckm_a]->Write();
	

 	TString strplotsA = "SomePlots/CLsLimit_versus_mass_10k_6boxes_";
	strplotsA += IntToStr(ckm_a).c_str();
	strplotsA += postfix;
	strplotsA += ".pdf";
 	TCanvas * PlotsA = new TCanvas("PlotsA","");
	PlotsA->SetGridx();
 	PlotsA->SetGridy();
	PlotsA->SetTickx();
	PlotsA->SetTicky();
	CLs_data_versus_mass[ckm_a] -> SetTitle("");
  	CLs_data_versus_mass[ckm_a] -> GetXaxis()-> SetTitle("m_{t'} = m_{b'} (GeV/c^{2})");
  	CLs_data_versus_mass[ckm_a] -> GetXaxis()->SetTitleOffset(1.2);
  	CLs_data_versus_mass[ckm_a] -> GetYaxis()-> SetTitle("CLs");
  	CLs_data_versus_mass[ckm_a]-> GetYaxis()->SetTitleOffset(1.45);
	CLs_median_versus_mass[ckm_a]->SetLineColor(4);
	CLs_median_versus_mass[ckm_a]->SetLineWidth(2);
	CLs_16th_versus_mass[ckm_a]->SetLineColor(4);
	CLs_16th_versus_mass[ckm_a]->SetLineStyle(2);
	CLs_16th_versus_mass[ckm_a]->SetLineWidth(2);
	CLs_84th_versus_mass[ckm_a]->SetLineColor(4);
	CLs_84th_versus_mass[ckm_a]->SetLineStyle(2);
	CLs_84th_versus_mass[ckm_a]->SetLineWidth(2);
	CLs_data_versus_mass[ckm_a]->SetLineColor(1);
	CLs_data_versus_mass[ckm_a]->SetLineWidth(2);
	CLs_data_versus_mass[ckm_a]->Draw("ALp");
	CLs_median_versus_mass[ckm_a]->Draw("Lp");
	CLs_16th_versus_mass[ckm_a]->Draw("Lp");
	CLs_84th_versus_mass[ckm_a]->Draw("Lp");
	TLegend *leg=new TLegend(0.2,0.8,0.45,0.95);
	leg->SetFillColor(0);
	leg->AddEntry(CLs_median_versus_mass[ckm_a],"expected","l");
	leg->AddEntry(CLs_data_versus_mass[ckm_a],"observed","l");
	//leg->Draw("Same");
	if(ckm_a==nCKM_A_values+1) PlotsA->Print(strplotsA);
	
	delete PlotsA;

  	TLine l1(330,CLs95CL,620,CLs95CL);
  	l1.SetLineWidth(2);
  	l1.SetLineColor(2);

	TString strplotsB = "SomePlots/CLsLimit_expectedlog_versus_mass_10k_6boxes_";
	strplotsB += IntToStr(ckm_a).c_str();
	strplotsB += postfix;
	strplotsB += ".pdf";
	TCanvas * plotsB = new TCanvas("plotsB","");
	plotsB->SetGridx();
  	plotsB->SetGridy();
	plotsB->SetTickx();
	plotsB->SetTicky();
	plotsB->SetLogy();
	CLs_median_versus_mass[ckm_a] -> SetTitle("");
  	CLs_median_versus_mass[ckm_a] -> GetXaxis()-> SetTitle("m_{t'} = m_{b'} (GeV/c^{2})");
  	CLs_median_versus_mass[ckm_a] -> GetXaxis()->SetTitleOffset(1.2);
  	CLs_median_versus_mass[ckm_a] -> GetYaxis()-> SetTitle("CLs");
  	CLs_median_versus_mass[ckm_a] -> GetYaxis()->SetTitleOffset(1.45);
	CLs_median_versus_mass[ckm_a] -> SetMaximum(1.0);
	CLs_median_versus_mass[ckm_a] -> SetMinimum(0.00001);
	CLs_median_versus_mass[ckm_a]->SetLineColor(4);
	CLs_median_versus_mass[ckm_a]->SetLineWidth(2);
	CLs_16th_versus_mass[ckm_a]->SetLineColor(4);
	CLs_16th_versus_mass[ckm_a]->SetLineStyle(2);
	CLs_16th_versus_mass[ckm_a]->SetLineWidth(2);
	CLs_84th_versus_mass[ckm_a]->SetLineColor(4);
	CLs_84th_versus_mass[ckm_a]->SetLineStyle(2);
	CLs_84th_versus_mass[ckm_a]->SetLineWidth(2);
	CLs_median_versus_mass[ckm_a]->Draw("ALp");
	CLs_16th_versus_mass[ckm_a]->Draw("Lp");
	CLs_84th_versus_mass[ckm_a]->Draw("Lp");
	l1.Draw("C");
	if(ckm_a==nCKM_A_values+1) plotsB->Print(strplotsB);
	
	delete plotsB;
}//CKM_A loop	



cout << "making the CKM_A versus mass plot" << endl;
// needed to store the masslimits in a different array with one place less, because the last value is the hypothetical case with CKM_A = 1-CKM_A = 1;
float exp_masslimits[nCKM_A_values+1], obs_masslimits[nCKM_A_values+1],exp16th_masslimits[nCKM_A_values+1],exp84th_masslimits[nCKM_A_values+1];
float CKMA[nCKM_A_values+1];
int n = nCKM_A_values+1;
for(int ckm_a = 0; ckm_a<nCKM_A_values+1; ckm_a++){
	CKMA[ckm_a] = 1 - ((float) ckm_a / (float) nCKM_A_values);
	exp_masslimits[ckm_a] = masslimit_exp[ckm_a]; 
	obs_masslimits[ckm_a] = masslimit_obs[ckm_a];
	exp16th_masslimits[ckm_a] = masslimit_exp_16th[ckm_a]; 
	exp84th_masslimits[ckm_a] = masslimit_exp_84th[ckm_a]; 
	cout << "for CKMA " << CKMA[ckm_a] << " exp_masslimit =  " << exp_masslimits[ckm_a] << endl;
	cout << "for CKMA " << CKMA[ckm_a] << " obs_masslimit =  " << obs_masslimits[ckm_a] << endl;
	cout << "for CKMA " << CKMA[ckm_a] << " exp16th_masslimit =  " << exp16th_masslimits[ckm_a] << endl;
	cout << "for CKMA " << CKMA[ckm_a] << " exp84th_masslimit =  " << exp84th_masslimits[ckm_a] << endl;
}

TGraph * CKM_A_versus_mass_expected = new TGraph(n,CKMA,exp_masslimits);
TGraph * CKM_A_versus_mass_observed = new TGraph(n,CKMA,obs_masslimits);
TGraph * CKM_A_versus_mass_expected_16th = new TGraph(n,CKMA,exp16th_masslimits);
TGraph * CKM_A_versus_mass_expected_84th = new TGraph(n,CKMA,exp84th_masslimits);

fout->cd();
CKM_A_versus_mass_expected -> Write();
CKM_A_versus_mass_observed -> Write();
fout->Close();

TCanvas * Plots = new TCanvas("Plots","");
Plots->SetGridx();
Plots->SetGridy();
Plots->SetTickx();
Plots->SetTicky();
CKM_A_versus_mass_expected -> SetMinimum(340);
CKM_A_versus_mass_expected -> SetMaximum(610);
CKM_A_versus_mass_expected -> GetXaxis()-> SetRangeUser(0.,1.);
CKM_A_versus_mass_expected -> SetTitle("");
CKM_A_versus_mass_expected -> GetXaxis()-> SetTitle("A");
CKM_A_versus_mass_expected -> GetXaxis()-> SetTitleOffset(1.2);
CKM_A_versus_mass_expected -> GetYaxis()-> SetTitle("m_{t'} = m_{b'} (GeV/c^{2})");
CKM_A_versus_mass_expected -> GetYaxis()-> SetTitleOffset(1.45);
CKM_A_versus_mass_expected->SetLineColor(4);
CKM_A_versus_mass_expected->SetLineWidth(2);
CKM_A_versus_mass_observed->SetLineColor(1);
CKM_A_versus_mass_observed->SetLineWidth(2);
CKM_A_versus_mass_expected_16th->SetLineStyle(2);
CKM_A_versus_mass_expected_16th->SetLineColor(4);
CKM_A_versus_mass_expected_16th->SetLineWidth(2);
CKM_A_versus_mass_expected_84th->SetLineStyle(2);
CKM_A_versus_mass_expected_84th->SetLineColor(4);
CKM_A_versus_mass_expected_84th->SetLineWidth(2);
CKM_A_versus_mass_expected->Draw("AL");
CKM_A_versus_mass_observed->Draw("L");
CKM_A_versus_mass_expected_16th->Draw("L");
CKM_A_versus_mass_expected_84th->Draw("L");
TLegend *leg=new TLegend(0.7,0.78,0.95,0.9);
leg->SetFillColor(0);
leg->AddEntry(CKM_A_versus_mass_expected,"expected","l");
leg->AddEntry(CKM_A_versus_mass_observed,"observed","l");
leg->Draw("Same");
TString strplots = "SomePlots/CKM_A_versus_mass";
strplots += postfix;
strplots += ".pdf";
Plots->Print(strplots);
delete Plots;

cout << "DONE" << endl;

  return 0;

}
