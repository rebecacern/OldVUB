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


std::string IntToStr( int n )
{
        std::ostringstream result;
        result << n;
        return result.str();
}


//////////////////////////////////////////
// Function to calculate the percentile //
//////////////////////////////////////////
// calculate the percentile (related to the specified second argument) of e.g. p-values;
// the 50th percentile is the median, and in principle, the 95th percentile is p95
float GetProbaPercentile(std::vector<float> vectorpvalues, double percent, int nofPE){
	std::sort(vectorpvalues.begin(),vectorpvalues.end());
	double pvaluesArray[nofPE];  // needs to be an array of doubles (Double_t) to make use of TMath::Quantiles();
	for(unsigned int i=0;i<nofPE;i++){
		pvaluesArray[i] = vectorpvalues[i];
	}
	vectorpvalues.clear();
	double probability[1], quantile[1];
	probability[0] = percent/100;
	TMath::Quantiles(nofPE,1,pvaluesArray,quantile,probability);  //see ROOT documentation for non-default calculation of quantiles
	return (float) quantile[0];
}

int main()
{
  using namespace std;

	const int nCKM_A_values = 7; //10 is default
	const int massvalues = 6;
	int n = nCKM_A_values+1;
	const int ckm_a_start = 0;

	float masses[massvalues] = {350,400,450,500,550,600};

	TString directories[nCKM_A_values+1] = {"CKMA0/","CKMA1/","CKMA2/","CKMA3/","CKMA4/","CKMA5/","CKMA6/","CKMA7/"};
	TString filenames[massvalues] = {"mass350/incl4gen_limit_plr_mujets_expected_test.ascii","mass400/incl4gen_limit_plr_mujets_expected_test.ascii","mass450/incl4gen_limit_plr_mujets_expected_test.ascii","mass500/incl4gen_limit_plr_mujets_expected_test.ascii","mass550/incl4gen_limit_plr_mujets_expected_test.ascii","mass600/incl4gen_limit_plr_mujets_expected_test.ascii"};
	TString filenames_observed[massvalues] = {"mass350/incl4gen_limit_plr_mujets_observed_test.ascii","mass400/incl4gen_limit_plr_mujets_observed_test.ascii","mass450/incl4gen_limit_plr_mujets_observed_test.ascii","mass500/incl4gen_limit_plr_mujets_observed_test.ascii","mass550/incl4gen_limit_plr_mujets_observed_test.ascii","mass600/incl4gen_limit_plr_mujets_observed_test.ascii"};
	
	float expected_limit[nCKM_A_values+1-ckm_a_start][massvalues];
	float expected_limit_plussigma[nCKM_A_values+1-ckm_a_start][massvalues];
	float expected_limit_minussigma[nCKM_A_values+1-ckm_a_start][massvalues];
	float observed_limit[nCKM_A_values+1-ckm_a_start][massvalues];

	for(int ckm_a = ckm_a_start; ckm_a < nCKM_A_values+1; ckm_a++){
			cout << "ckm_a "<< ckm_a << endl;
		for(int mass = 0; mass < massvalues; mass++){
			cout << "	massvalue "<< 350 + (mass*50) << endl;
			
			int atboundary = 0;
			
			TString filename_tmp = directories[ckm_a]+filenames[mass];
			ifstream file_to_read(filename_tmp);

  		char line[256];
  		vector<TString> upperlimits_string;
 			vector<float> upperlimits;
  		while(!file_to_read.eof()){
  			file_to_read.getline(line,256);
  			upperlimits_string.push_back(line);
  		}
			
  		for(int i = 0; i < upperlimits_string.size()-1; i++){
				//cout << "upperlimit " << i << ": " << upperlimits_string[i] << endl;
				
				if( atof(upperlimits_string[i])== 0 ||  atof(upperlimits_string[i])== 5) atboundary++;
				if( atof(upperlimits_string[i])!= 0 && atof(upperlimits_string[i])!= 5)
					upperlimits.push_back(atof(upperlimits_string[i]));
				//cout << "upperlimit in float format " << i << ": " << upperlimits[i] << endl;
		
				//if(upperlimits[i] == 0 || upperlimits[i] == 5) atboundary++;
			}
			int nPE = upperlimits.size();
			cout << "		fraction of pseudo-experiments at boundary " << (float) atboundary/nPE << endl;
			
			expected_limit[ckm_a][mass] = GetProbaPercentile(upperlimits,50,nPE);
			cout << "		expected limit " << expected_limit[ckm_a][mass] << endl;
			expected_limit_plussigma[ckm_a][mass] = GetProbaPercentile(upperlimits,84,nPE);
			cout << "		expected limit plussigma " << expected_limit_plussigma[ckm_a][mass] << endl;
			expected_limit_minussigma[ckm_a][mass] = GetProbaPercentile(upperlimits,16,nPE);
			cout << "		expected limit minussigma " << expected_limit_minussigma[ckm_a][mass] << endl;

			upperlimits.clear();
			upperlimits_string.clear();

			TString filename_obs_tmp = directories[ckm_a]+filenames_observed[mass];
			ifstream file_to_read_obs(filename_obs_tmp);

  		char number[256];
			file_to_read_obs.getline(number,256);
			observed_limit[ckm_a][mass] = atof(number);

		}
	
	}

  setTDRStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gSystem->mkdir("Plots");
	
	TGraph *limit_median_versus_mass[nCKM_A_values+1-ckm_a_start];
	TGraph *limit_plussigma_versus_mass[nCKM_A_values+1-ckm_a_start];
	TGraph *limit_minussigma_versus_mass[nCKM_A_values+1-ckm_a_start];
	TGraph *limit_data_versus_mass[nCKM_A_values+1-ckm_a_start];


	float mass_exp[nCKM_A_values+1-ckm_a_start][2],mass_obs[nCKM_A_values+1-ckm_a_start][2],mass_exp_plussigma[nCKM_A_values+1-ckm_a_start][2], mass_exp_minussigma[nCKM_A_values+1-ckm_a_start][2];
	float limit_exp[nCKM_A_values+2][2],limit_obs[nCKM_A_values+1-ckm_a_start][2],limit_exp_plussigma[nCKM_A_values+1-ckm_a_start][2],limit_exp_minussigma[nCKM_A_values+1-ckm_a_start][2];
	float masslimit_exp[nCKM_A_values+1-ckm_a_start],masslimit_obs[nCKM_A_values+1-ckm_a_start],masslimit_exp_plussigma[nCKM_A_values+1-ckm_a_start],masslimit_exp_minussigma[nCKM_A_values+1-ckm_a_start];
	
	float CKMA[nCKM_A_values+1-ckm_a_start];
	
for(int ckm_a = ckm_a_start; ckm_a<nCKM_A_values+1; ckm_a++){
	CKMA[ckm_a] = 1 - ((float) ckm_a / (float) 10);
	cout << "CKMA[ckm_a]" << CKMA[ckm_a] << endl;

	limit_median_versus_mass[ckm_a] = new TGraph(massvalues,masses,expected_limit[ckm_a]);
	limit_plussigma_versus_mass[ckm_a] = new TGraph(massvalues,masses,expected_limit_plussigma[ckm_a]);
	limit_minussigma_versus_mass[ckm_a] = new TGraph(massvalues,masses,expected_limit_minussigma[ckm_a]);
	limit_data_versus_mass[ckm_a] = new TGraph(massvalues,masses,observed_limit[ckm_a]);

	TString strlimit_median_versus_mass = "limit_median_versus_mass_";
	strlimit_median_versus_mass += IntToStr(ckm_a).c_str();
	limit_median_versus_mass[ckm_a]->SetTitle(strlimit_median_versus_mass);
	TString strlimit_plussigma_versus_mass = "limit_plussigma_versus_mass_";
	strlimit_plussigma_versus_mass += IntToStr(ckm_a).c_str();
	limit_plussigma_versus_mass[ckm_a]->SetTitle(strlimit_plussigma_versus_mass);
	TString strlimit_minussigma_versus_mass = "limit_minussigma_versus_mass_";
	strlimit_minussigma_versus_mass += IntToStr(ckm_a).c_str();
	limit_minussigma_versus_mass[ckm_a]->SetTitle(strlimit_minussigma_versus_mass);
	TString strlimit_data_versus_mass = "limit_data_versus_mass_";
	strlimit_data_versus_mass += IntToStr(ckm_a).c_str();
	limit_data_versus_mass[ckm_a]->SetTitle(strlimit_data_versus_mass);
	

	mass_exp[ckm_a][0] = 0; mass_exp[ckm_a][1] = 0; limit_exp[ckm_a][0] = 0; limit_exp[ckm_a][1] = 0;
	mass_obs[ckm_a][0] = 0; mass_obs[ckm_a][1] = 0; limit_obs[ckm_a][0] = 0; limit_obs[ckm_a][1] = 0;
	mass_exp_plussigma[ckm_a][0] = 0; mass_exp_plussigma[ckm_a][1] = 0; limit_exp_plussigma[ckm_a][0] = 0; limit_exp_plussigma[ckm_a][1] = 0;
	mass_exp_minussigma[ckm_a][0] = 0; mass_exp_minussigma[ckm_a][1] = 0; limit_exp_minussigma[ckm_a][0] = 0; limit_exp_minussigma[ckm_a][1] = 0;
	masslimit_exp[ckm_a]= 0; masslimit_exp_minussigma[ckm_a]= 0; masslimit_exp_plussigma[ckm_a]= 0; masslimit_obs[ckm_a]= 0; 

	for(int i = 0 ; i < massvalues-1; i++){
		if(expected_limit[ckm_a][i]<1 && expected_limit[ckm_a][i+1]>1){
			cout << "	expected masslimit is between " << masses[i] << " and " << masses[i+1] << " GeV" << endl;
			mass_exp[ckm_a][0] = masses[i];
			mass_exp[ckm_a][1] = masses[i+1];
			limit_exp[ckm_a][0] = expected_limit[ckm_a][i];
			limit_exp[ckm_a][1] = expected_limit[ckm_a][i+1];
		}
		if(observed_limit[ckm_a][i]<1 && observed_limit[ckm_a][i+1]>1){
			cout << "	observed masslimit is between " << masses[i] << " and " << masses[i+1] << " GeV" << endl;
			mass_obs[ckm_a][0] = masses[i];
			mass_obs[ckm_a][1] = masses[i+1];
			limit_obs[ckm_a][0] = observed_limit[ckm_a][i];
			limit_obs[ckm_a][1] = observed_limit[ckm_a][i+1];
		}
		if(expected_limit_plussigma[ckm_a][i]<1 && expected_limit_plussigma[ckm_a][i+1]>1){
			cout << "	plussigma masslimit is between " << masses[i] << " and " << masses[i+1] << " GeV" << endl;
			mass_exp_plussigma[ckm_a][0] = masses[i];
			mass_exp_plussigma[ckm_a][1] = masses[i+1];
			limit_exp_plussigma[ckm_a][0] = expected_limit_plussigma[ckm_a][i];
			limit_exp_plussigma[ckm_a][1] = expected_limit_plussigma[ckm_a][i+1];
		}
		if(expected_limit_minussigma[ckm_a][i]<1 && expected_limit_minussigma[ckm_a][i+1]>1){
			cout << "	minussigma masslimit is between " << masses[i] << " and " << masses[i+1] << " GeV" << endl;
			mass_exp_minussigma[ckm_a][0] = masses[i];
			mass_exp_minussigma[ckm_a][1] = masses[i+1];
			limit_exp_minussigma[ckm_a][0] = expected_limit_minussigma[ckm_a][i];
			limit_exp_minussigma[ckm_a][1] = expected_limit_minussigma[ckm_a][i+1];
		
		}
		
	}
	
	if(mass_exp[ckm_a][0] == 0 && mass_exp[ckm_a][0] == 0 ){
		cout << "	expected masslimit < 350 GeV!" << endl;
		if(expected_limit[ckm_a][0]>1) masslimit_exp[ckm_a] = 345;
		if(expected_limit[ckm_a][5]<1) masslimit_exp[ckm_a] = 605;			
	}
	
	if(mass_exp_plussigma[ckm_a][0] == 0 && mass_exp_plussigma[ckm_a][1] == 0){
		cout << "	expected plussigma masslimit  < 350 GeV!" << endl;
		if(expected_limit_plussigma[ckm_a][0]>1) masslimit_exp_plussigma[ckm_a] = 345;	
		if(expected_limit_plussigma[ckm_a][1]<1) masslimit_exp_plussigma[ckm_a] = 605;	
	}
	if(mass_exp_minussigma[ckm_a][0] == 0 && mass_exp_minussigma[ckm_a][1] == 0){
		cout << "	expected minussigma masslimit  < 350 GeV!" << endl;
		if(expected_limit_minussigma[ckm_a][0]>1) masslimit_exp_minussigma[ckm_a] = 345;	
		if(expected_limit_minussigma[ckm_a][1]<1) masslimit_exp_minussigma[ckm_a] = 605;	
	}
	if(mass_obs[ckm_a][0] == 0 && mass_obs[ckm_a][1] == 0){
		cout << "	observed masslimit  < 350 GeV!" << endl;
		if(observed_limit[ckm_a][0]>1) masslimit_obs[ckm_a] = 345;
		if(observed_limit[ckm_a][1]<1) masslimit_obs[ckm_a] = 605;
	}
	
	float limitpoint = 1.;
	TGraph mass_versus_limit_exp(2,limit_exp[ckm_a],mass_exp[ckm_a]);
	TF1 f1("f1","pol1");
  if(mass_exp[ckm_a][0] != 0 && mass_exp[ckm_a][1] != 0) mass_versus_limit_exp.Fit("f1","Q");    
	if(mass_exp[ckm_a][0] != 0 && mass_exp[ckm_a][1] != 0) masslimit_exp[ckm_a] = f1.Eval(limitpoint);
	
	TGraph mass_versus_limit_obs(2,limit_obs[ckm_a],mass_obs[ckm_a]);
	TF1 f2("f2","pol1");
  if(mass_obs[ckm_a][0] != 0 && mass_obs[ckm_a][1] != 0) mass_versus_limit_obs.Fit("f2","Q");    
	if(mass_obs[ckm_a][0] != 0 && mass_obs[ckm_a][1] != 0) masslimit_obs[ckm_a] = f2.Eval(limitpoint);

	TGraph mass_versus_limit_exp_plussigma(2,limit_exp_plussigma[ckm_a],mass_exp_plussigma[ckm_a]);
	TF1 f3("f3","pol1");
  if(mass_exp_plussigma[ckm_a][0] != 0 && mass_exp_plussigma[ckm_a][1] != 0) mass_versus_limit_exp_plussigma.Fit("f3","Q");    
	if(mass_exp_plussigma[ckm_a][0] != 0 && mass_exp_plussigma[ckm_a][1] != 0) masslimit_exp_plussigma[ckm_a] = f3.Eval(limitpoint);

	TGraph mass_versus_limit_exp_minussigma(2,limit_exp_minussigma[ckm_a],mass_exp_minussigma[ckm_a]);
	TF1 f4("f4","pol1");
  if(mass_exp_minussigma[ckm_a][0] != 0 && mass_exp_minussigma[ckm_a][1] != 0) mass_versus_limit_exp_minussigma.Fit("f4","Q");    
	if(mass_exp_minussigma[ckm_a][0] != 0 && mass_exp_minussigma[ckm_a][1] != 0) masslimit_exp_minussigma[ckm_a] = f4.Eval(limitpoint);

	cout << "masses below " << masslimit_exp[ckm_a] << " are expected to be excluded" << endl;
	cout << "masses below " << masslimit_obs[ckm_a] << " are observed to be excluded" << endl;
	cout << "plussigma expected masslimit " << masslimit_exp_plussigma[ckm_a] << endl;
	cout << "minussigma on expected masslimit " << masslimit_exp_minussigma[ckm_a] << endl;
}

	TGraph * CKM_A_versus_mass_expected = new TGraph(n,CKMA,masslimit_exp);
	TGraph * CKM_A_versus_mass_observed = new TGraph(n,CKMA,masslimit_obs);

	TGraph * CKM_A_versus_mass_expected_1sigma = new TGraph(2*n);
	for(int i = 0; i<nCKM_A_values+1; i++){
		CKM_A_versus_mass_expected_1sigma->SetPoint(i,CKMA[i],masslimit_exp_plussigma[i]);
		CKM_A_versus_mass_expected_1sigma->SetPoint(n+i,CKMA[n-i-1],masslimit_exp_minussigma[n-i-1]);
	}
	
	TCanvas * Plots = new TCanvas("Plots","");
	Plots->SetGridx();
	Plots->SetGridy();
	Plots->SetTickx();
	Plots->SetTicky();
	CKM_A_versus_mass_expected_1sigma -> SetMinimum(350);
	CKM_A_versus_mass_expected_1sigma -> SetMaximum(600);
	CKM_A_versus_mass_expected_1sigma -> GetXaxis()-> SetRangeUser(0.3,1.);
	CKM_A_versus_mass_expected_1sigma -> SetTitle("");	
	CKM_A_versus_mass_expected_1sigma -> GetXaxis()-> SetTitle("A");
	CKM_A_versus_mass_expected_1sigma -> GetXaxis()-> SetTitleOffset(1.2);
	CKM_A_versus_mass_expected_1sigma -> GetYaxis()-> SetTitle("m_{t'} = m_{b'} (GeV/c^{2})");
	CKM_A_versus_mass_expected_1sigma -> GetYaxis()-> SetTitleOffset(1.45);
	CKM_A_versus_mass_expected->SetLineColor(4); //green
	CKM_A_versus_mass_expected->SetLineWidth(2);
	CKM_A_versus_mass_observed->SetLineColor(1); //black
	CKM_A_versus_mass_observed->SetLineWidth(2);
	CKM_A_versus_mass_expected_1sigma->SetFillColor(8); //green
	CKM_A_versus_mass_expected_1sigma->Draw("FA");
	CKM_A_versus_mass_observed->Draw("L");
	CKM_A_versus_mass_expected->Draw("L"); 
	TLegend *leg=new TLegend(0.2,0.2,0.5,0.4);
	leg->SetFillColor(0);
	leg->AddEntry(CKM_A_versus_mass_expected,"expected","l");
	leg->AddEntry(CKM_A_versus_mass_observed,"observed","l");
	leg->AddEntry(CKM_A_versus_mass_expected_1sigma,"68% C.L.","f");
	leg->Draw("Same");
	TString strplots = "Plots/CKM_A_versus_mass";
	strplots += ".pdf";
	Plots->Print(strplots);
	delete Plots;


}
