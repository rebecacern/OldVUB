//Rebeca Gonzalez Suarez
//rebeca@cern.ch

#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "inputs.h"

using namespace std;

void datacardmaker(){

  bool verbose = false;
  
  //choose regions
  const int SRbin = 28;
  const int CR1bin = 30;
  const int CR2bin = 31;
  
  //choosing precision
  int pre = 4;
  
  ofstream datacard("singletop_tW_5fb.txt"); 
 
  char myRootFile[300];
  sprintf(myRootFile,"results/Histos_cutbased_full.root");
  
  TFile *_file0 = TFile::Open(myRootFile);
  
  const int nProcess = 4;
  TString processName[nProcess+1] =  { "data", "twdr", "tt", "zjets", "others"};
  
  TH1F*  hdata [3];
  TH1F*  hnom [3][3];
  
  for (int process = 0; process < nProcess; process++){
    for (int i = 0; i < 3; i++){
      
      int mode = 0;
      if (i < 2) mode = i+1;
      char modeName[300];
      if (i == 0) sprintf(modeName,"emu");
      else if (i == 1) sprintf(modeName,"mumu");
      else sprintf(modeName,"ee");
      char title[300];
      sprintf(title,"R_%s_",modeName);
      if (verbose) cout << "Reading " << title + processName[process] << ", ordered at " << mode << endl;
      if (process == 0) hdata[mode] = (TH1F*) (_file0->Get(title + processName[process]))->Clone();
      else {
	hnom[mode][process-1] = (TH1F*) (_file0->Get(title + processName[process]))->Clone();
	if (process == 3){
	  TH1F* h1 = (TH1F*) (_file0->Get(title + processName[process+1]))->Clone();
	  if (verbose) cout << "Adding " << title + processName[process+1] << ", ordered at " << mode << endl;
	  hnom[mode][process-1]->Add(h1);
	}
      }    
    }
  } 
  if (verbose) cout << endl;
  
  datacard << setprecision(pre) << "# this is the version with *exclusive* jet / tag bins " << endl;
  datacard << setprecision(pre) << "# based on 4.9/fb " << endl;
  datacard << setprecision(pre) << "imax 9 # number of bins " << endl;
  datacard << setprecision(pre) << "jmax 2 # number of processes - 1 " << endl;
  datacard << setprecision(pre) << "kmax * # number of uncertainties " << endl;
  datacard << setprecision(pre) << "------------ " << endl;
  datacard << setprecision(pre) << "bin \t\tee1j1t\t emu1j1t\t mumu1j1t\t ee2j1t\t emu2j1t\t mumu2j1t\t ee2j2t\t emu2j2t\t mumu2j2t " << endl;
  datacard << setprecision(pre) << "observation" ;
  
  for (int i = 0; i < 3; i++){
    datacard << setprecision(pre) << "\t " << hdata[i]->GetBinContent(SRbin) ;
  }
  for (int i = 0; i < 3; i++){
    datacard << setprecision(pre) << "\t " << hdata[i]->GetBinContent(CR1bin) ;
  }
  for (int i = 0; i < 3; i++){
    datacard << setprecision(pre) << "\t " << hdata[i]->GetBinContent(CR2bin) ;
  }
  
  datacard << setprecision(pre) << endl;
  datacard << setprecision(pre) << "------------ " << endl;
  datacard << setprecision(pre) << "bin        \tee1j1t\t ee1j1t\t ee1j1t\t emu1j1t\t emu1j1t\t emu1j1t\t mumu1j1t\t mumu1j1t\t mumu1j1t\t ee2j1t\t ee2j1t \t ee2j1t\t emu2j1t\t emu2j1t\t emu2j1t\t mumu2j1t\t mumu2j1t\t mumu2j1t\t ee2j2t\t ee2j2t\t ee2j2t\t emu2j2t\t emu2j2t\t emu2j2t\t mumu2j2t\t mumu2j2t\t mumu2j2t" << endl;
  datacard << setprecision(pre) << "process    \tst\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt    \t other\t st\t tt\t other" << endl;
  datacard << setprecision(pre) << "process    \t0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2" << endl;
  datacard << setprecision(pre) << "rate       \t" ;
  
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << hnom[i][j]->GetBinContent(SRbin) << "\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << hnom[i][j]->GetBinContent(CR1bin) << "\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << hnom[i][j]->GetBinContent(CR2bin) << "\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "---------------------- " << endl;
  datacard << setprecision(pre) <<  "lumi      lnN\t1.022\t 1.022\t -\t1.022\t 1.022\t 1.022\t 1.022\t 1.022\t -\t 1.022\t 1.022\t 1.022\t 1.022\t 1.022\t 1.022\t1.022\t 1.022\t 1.022	1.022\t1.022\t 1.022    1.022\t 1.022\t 1.022    1.022\t 1.022\t 1.022 " << endl;
  datacard << setprecision(pre) <<  "hlte      lnN\t1.015\t 1.015\t -\t1.011\t 1.011\t 1.011\t -\t  -\t  -\t 1.015\t 1.015\t 1.015\t 1.011\t 1.011\t 1.011\t-\t -\t -\t1.015\t1.015\t 1.015    1.011\t 1.011\t 1.011	-\t	-\t	- " << endl;
  datacard << setprecision(pre) <<  "hltmu     lnN\t-\t  -\t  -\t1.011\t 1.011\t 1.011\t 1.015\t 1.015\t -\t -\t	-\t	-\t1.011\t 1.011\t 1.011	1.015\t 1.015\t 1.015\t -\t	-\t  -\t 1.011\t 1.011\t 1.011    1.015\t 1.015\t 1.015 " << endl;
  datacard << setprecision(pre) <<  "ele       lnN\t1.02\t1.02\t-\t1.02\t1.02\t 1.02\t -\t  -\t -\t1.02\t 1.02\t1.02\t1.02\t1.02\t 1.02\t -\t -\t -\t1.02\t 1.02\t 1.02\t 1.02\t 1.02\t1.02\t -\t -\t - " << endl;
  datacard << setprecision(pre) <<  "mu        lnN\t-\t  -\t  -\t1.01\t1.01\t 1.01\t 1.01\t1.01\t-\t  -\t  -\t -\t1.01\t1.01\t 1.01\t 1.01\t1.01\t1.01	-\t -\t  -\t 1.01\t 1.01\t1.01	 1.01\t1.01\t1.01 " << endl;
  
 
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;  
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
    
    for (int j = 0; j < 3; j++){
      if (j == 0){ 
        char title[300];
        sprintf(title,"R_%s_tw_sup",modeName);
	if (verbose) cout << "Reading " << modeName << ", " << title << ", ordered at " << mode << endl;
        hup[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
        sprintf(title,"R_%s_tw_sdo",modeName);
	if (verbose) cout << "Reading " << modeName << ", " << title << ", ordered at " << mode << endl;
	hdown[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
      } else if (j == 1) {
        char title[300];
        sprintf(title,"R_%s_tt_scaleup",modeName);
	if (verbose) cout << "Reading " << modeName << ", " << title << ", ordered at " << mode << endl;
        hup[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
        sprintf(title,"R_%s_tt_scaledown",modeName);
	if (verbose) cout << "Reading " << modeName << ", " << title << ", ordered at " << mode << endl;
	hdown[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
      }
    }
  } 
  if (verbose) cout << endl;
  
  datacard << setprecision(pre) << "ttscale   lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
        double avn = (hnom[0][j]->GetBinContent(SRbin)+hnom[1][j]->GetBinContent(SRbin)+hnom[2][j]->GetBinContent(SRbin));
        double avup = (hup[0][j]->GetBinContent(SRbin)+hup[1][j]->GetBinContent(SRbin)+hup[2][j]->GetBinContent(SRbin));
	double avdown = (hdown[0][j]->GetBinContent(SRbin)+hdown[1][j]->GetBinContent(SRbin)+hdown[2][j]->GetBinContent(SRbin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
        double avn = (hnom[0][j]->GetBinContent(CR1bin)+hnom[1][j]->GetBinContent(CR1bin)+hnom[2][j]->GetBinContent(CR1bin));
        double avup = (hup[0][j]->GetBinContent(CR1bin)+hup[1][j]->GetBinContent(CR1bin)+hup[2][j]->GetBinContent(CR1bin));
	double avdown = (hdown[0][j]->GetBinContent(CR1bin)+hdown[1][j]->GetBinContent(CR1bin)+hdown[2][j]->GetBinContent(CR1bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
        double avn = (hnom[0][j]->GetBinContent(CR2bin)+hnom[1][j]->GetBinContent(CR2bin)+hnom[2][j]->GetBinContent(CR2bin));
        double avup = (hup[0][j]->GetBinContent(CR2bin)+hup[1][j]->GetBinContent(CR2bin)+hup[2][j]->GetBinContent(CR2bin));
	double avdown = (hdown[0][j]->GetBinContent(CR2bin)+hdown[1][j]->GetBinContent(CR2bin)+hdown[2][j]->GetBinContent(CR2bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  
  datacard << setprecision(pre) << "twscale   lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double avn = (hnom[0][j]->GetBinContent(SRbin)+hnom[1][j]->GetBinContent(SRbin)+hnom[2][j]->GetBinContent(SRbin));
        double avup = (hup[0][j]->GetBinContent(SRbin)+hup[1][j]->GetBinContent(SRbin)+hup[2][j]->GetBinContent(SRbin));
	double avdown = (hdown[0][j]->GetBinContent(SRbin)+hdown[1][j]->GetBinContent(SRbin)+hdown[2][j]->GetBinContent(SRbin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  double maxv = TMath::Max(fabs(upv), fabs(downv)); 
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
        double avn = (hnom[0][j]->GetBinContent(CR1bin)+hnom[1][j]->GetBinContent(CR1bin)+hnom[2][j]->GetBinContent(CR1bin));
        double avup = (hup[0][j]->GetBinContent(CR1bin)+hup[1][j]->GetBinContent(CR1bin)+hup[2][j]->GetBinContent(CR1bin));
	double avdown = (hdown[0][j]->GetBinContent(CR1bin)+hdown[1][j]->GetBinContent(CR1bin)+hdown[2][j]->GetBinContent(CR1bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
        double avn = (hnom[0][j]->GetBinContent(CR2bin)+hnom[1][j]->GetBinContent(CR2bin)+hnom[2][j]->GetBinContent(CR2bin));
        double avup = (hup[0][j]->GetBinContent(CR2bin)+hup[1][j]->GetBinContent(CR2bin)+hup[2][j]->GetBinContent(CR2bin));
	double avdown = (hdown[0][j]->GetBinContent(CR2bin)+hdown[1][j]->GetBinContent(CR2bin)+hdown[2][j]->GetBinContent(CR2bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;  
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else  sprintf(modeName,"ee");
    
    for (int j = 0; j < 3; j++){
      if (j == 1) {
        char title[300];
        sprintf(title,"R_%s_tt_matchingup",modeName);
	if (verbose) cout << "Reading " << modeName << ", " << title << ", ordered at " << mode << endl;
        hup[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
        sprintf(title,"R_%s_tt_matchingdown",modeName);
	if (verbose) cout << "Reading " << modeName << ", " << title << ", ordered at " << mode << endl;
	hdown[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
      }
    }
  } 
  if (verbose) cout << endl;
 
  
  datacard << setprecision(pre) << "ttmatch   lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
        double avn = (hnom[0][j]->GetBinContent(SRbin)+hnom[1][j]->GetBinContent(SRbin)+hnom[2][j]->GetBinContent(SRbin));
        double avup = (hup[0][j]->GetBinContent(SRbin)+hup[1][j]->GetBinContent(SRbin)+hup[2][j]->GetBinContent(SRbin));
	double avdown = (hdown[0][j]->GetBinContent(SRbin)+hdown[1][j]->GetBinContent(SRbin)+hdown[2][j]->GetBinContent(SRbin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
        double avn = (hnom[0][j]->GetBinContent(CR1bin)+hnom[1][j]->GetBinContent(CR1bin)+hnom[2][j]->GetBinContent(CR1bin));
        double avup = (hup[0][j]->GetBinContent(CR1bin)+hup[1][j]->GetBinContent(CR1bin)+hup[2][j]->GetBinContent(CR1bin));
	double avdown = (hdown[0][j]->GetBinContent(CR1bin)+hdown[1][j]->GetBinContent(CR1bin)+hdown[2][j]->GetBinContent(CR1bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
        double avn = (hnom[0][j]->GetBinContent(CR2bin)+hnom[1][j]->GetBinContent(CR2bin)+hnom[2][j]->GetBinContent(CR2bin));
        double avup = (hup[0][j]->GetBinContent(CR2bin)+hup[1][j]->GetBinContent(CR2bin)+hup[2][j]->GetBinContent(CR2bin));
	double avdown = (hdown[0][j]->GetBinContent(CR2bin)+hdown[1][j]->GetBinContent(CR2bin)+hdown[2][j]->GetBinContent(CR2bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  
  TH1F*  h [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;  
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
     
    for (int j = 0; j < 3; j++){
      if (j == 0){ 
	char title[300];
	sprintf(title,"R_%s_twds",modeName);
	if (verbose) cout << "Reading " << modeName << ", " << title << ", ordered at " << mode << endl;
	h[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
      }
    }
  } 
  if (verbose) cout << endl;
   
   
  datacard << setprecision(pre) << "twdrds    lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = hnom[0][j]->GetBinContent(SRbin) + hnom[1][j]->GetBinContent(SRbin) + hnom[2][j]->GetBinContent(SRbin);
	if (average !=0) datacard << setprecision(pre) << 1 + ((h[0][j]->GetBinContent(SRbin) + h[1][j]->GetBinContent(SRbin) + h[2][j]->GetBinContent(SRbin)) - average)/average << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      } 
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = hnom[0][j]->GetBinContent(CR1bin) + hnom[1][j]->GetBinContent(CR1bin) + hnom[2][j]->GetBinContent(CR1bin);
	if (average !=0) datacard << setprecision(pre) << 1 + ((h[0][j]->GetBinContent(CR1bin) + h[1][j]->GetBinContent(CR1bin) + h[2][j]->GetBinContent(CR1bin)) - average)/average << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = hnom[0][j]->GetBinContent(CR2bin) + hnom[1][j]->GetBinContent(CR2bin) + hnom[2][j]->GetBinContent(CR2bin);
	if (average !=0) datacard << setprecision(pre) << 1 + ((h[0][j]->GetBinContent(CR2bin) + h[1][j]->GetBinContent(CR2bin) + h[2][j]->GetBinContent(CR2bin)) - average)/average << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  TString SystName = "PU";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
     
    for (int j = 0; j < 3; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      if (verbose) cout << "Reading " << modeName << ", " << title + processName[j+1]+ "__" + SystName + "__plus"  << ", ordered at " << mode << endl;
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      if (j == 2){
        if (verbose) cout << "Addding " << modeName << ", " << title + processName[j+2]+ "__" + SystName + "__plus"  << ", ordered at " << mode << endl;
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__plus"))->Clone();
	hup[mode][j]->Add(h1);
      }
      if (verbose) cout << "Reading " << modeName << ", " << title + processName[j+1]+ "__" + SystName + "__minus"  << ", ordered at " << mode << endl;
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
      if (j == 2){
        if (verbose) cout << "Addding " << modeName << ", " << title + processName[j+2]+ "__" + SystName + "__minus" << ", ordered at " << mode << endl;
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__minus"))->Clone();
	hdown[mode][j]->Add(h1);
      } 
    }
  } 
  if (verbose) cout << endl;
   
  datacard << setprecision(pre) << "pu        lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(SRbin) !=0){
	  double upv = (hup[i][j]->GetBinContent(SRbin) - hnom[i][j]->GetBinContent(SRbin))/hnom[i][j]->GetBinContent(SRbin);
	  double downv =  (hdown[i][j]->GetBinContent(SRbin) - hnom[i][j]->GetBinContent(SRbin))/hnom[i][j]->GetBinContent(SRbin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0){
	  double upv = (hup[i][j]->GetBinContent(CR1bin) - hnom[i][j]->GetBinContent(CR1bin))/hnom[i][j]->GetBinContent(CR1bin);
	  double downv =  (hdown[i][j]->GetBinContent(CR1bin) - hnom[i][j]->GetBinContent(CR1bin))/hnom[i][j]->GetBinContent(CR1bin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0){
	  double upv = (hup[i][j]->GetBinContent(CR2bin) - hnom[i][j]->GetBinContent(CR2bin))/hnom[i][j]->GetBinContent(CR2bin);
	  double downv =  (hdown[i][j]->GetBinContent(CR2bin) - hnom[i][j]->GetBinContent(CR2bin))/hnom[i][j]->GetBinContent(CR2bin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
   
  TString SystName = "JES";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
     
    for (int j = 0; j < 3; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      if (verbose) cout << "Reading " << modeName << ", " << title + processName[j+1]+ "__" + SystName + "__plus"  << ", ordered at " << mode << endl;
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      if (j == 2){
        if (verbose) cout << "Adding " << modeName << ", " << title + processName[j+2]+ "__" + SystName + "__plus"  << ", ordered at " << mode << endl;
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__plus"))->Clone();
	hup[mode][j]->Add(h1);
      }
      if (verbose) cout << "Reading " << modeName << ", " << title + processName[j+1]+ "__" + SystName + "__minus"  << ", ordered at " << mode << endl;
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
      if (j == 2){
        if (verbose) cout << "Adding " << modeName << ", " << title + processName[j+2]+ "__" + SystName + "__minus"  << ", ordered at " << mode << endl;
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__minus"))->Clone();
	hdown[mode][j]->Add(h1);
      }  
    }
  } 
  if (verbose) cout << endl;
   
  datacard << setprecision(pre) << "jes       lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(SRbin) !=0){
	  double upv = (hup[i][j]->GetBinContent(SRbin) - hnom[i][j]->GetBinContent(SRbin))/hnom[i][j]->GetBinContent(SRbin);
	  double downv =  (hdown[i][j]->GetBinContent(SRbin) - hnom[i][j]->GetBinContent(SRbin))/hnom[i][j]->GetBinContent(SRbin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0){
	  double upv = (hup[i][j]->GetBinContent(CR1bin) - hnom[i][j]->GetBinContent(CR1bin))/hnom[i][j]->GetBinContent(CR1bin);
          double downv =  (hdown[i][j]->GetBinContent(CR1bin) - hnom[i][j]->GetBinContent(CR1bin))/hnom[i][j]->GetBinContent(CR1bin);
          double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0){
	  double upv = (hup[i][j]->GetBinContent(CR2bin) - hnom[i][j]->GetBinContent(CR2bin))/hnom[i][j]->GetBinContent(CR2bin);
	  double downv =  (hdown[i][j]->GetBinContent(CR2bin) - hnom[i][j]->GetBinContent(CR2bin))/hnom[i][j]->GetBinContent(CR2bin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
   
  TString SystName = "btag";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
     
    for (int j = 0; j < 3; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      if (verbose) cout << "Reading " << modeName << ", " <<title + processName[j+1]+ "__" + SystName + "__plus"  << ", ordered at " << mode << endl;
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      if (j == 2){
        if (verbose) cout << "Adding " << modeName << ", " << title + processName[j+2]+ "__" + SystName + "__plus"  << ", ordered at " << mode << endl;
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__plus"))->Clone();
	hup[mode][j]->Add(h1);
      }
      if (verbose) cout << "Reading " << modeName << ", " << title + processName[j+1]+ "__" + SystName + "__minus"  << ", ordered at " << mode << endl;
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
      if (j == 2){
        if (verbose) cout << "Adding " << modeName << ", " << title + processName[j+2]+ "__" + SystName + "__minus" << ", ordered at " << mode << endl;
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__minus"))->Clone();
	hdown[mode][j]->Add(h1);
      } 
    }
  } 
  if (verbose) cout << endl;
   
  datacard << setprecision(pre) << "btag      lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(SRbin) !=0){
	  double upv = (hup[i][j]->GetBinContent(SRbin) - hnom[i][j]->GetBinContent(SRbin))/hnom[i][j]->GetBinContent(SRbin);
	  double downv =  (hdown[i][j]->GetBinContent(SRbin) - hnom[i][j]->GetBinContent(SRbin))/hnom[i][j]->GetBinContent(SRbin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0){
	  double upv = (hup[i][j]->GetBinContent(CR1bin) - hnom[i][j]->GetBinContent(CR1bin))/hnom[i][j]->GetBinContent(CR1bin);
          double downv =  (hdown[i][j]->GetBinContent(CR1bin) - hnom[i][j]->GetBinContent(CR1bin))/hnom[i][j]->GetBinContent(CR1bin);
          double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0){
	  double upv = (hup[i][j]->GetBinContent(CR2bin) - hnom[i][j]->GetBinContent(CR2bin))/hnom[i][j]->GetBinContent(CR2bin);
	  double downv =  (hdown[i][j]->GetBinContent(CR2bin) - hnom[i][j]->GetBinContent(CR2bin))/hnom[i][j]->GetBinContent(CR2bin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  } datacard << setprecision(pre) << endl;  
  
  TString SystName = "JER";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
     
    for (int j = 0; j < 3; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      if (verbose) cout << "Reading " << modeName << ", " << title + processName[j+1]+ "__" + SystName + "__plus" << ", ordered at " << mode << endl;
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      if (j == 2){
        if (verbose) cout << "Adding " << modeName << ", " << title + processName[j+2]+ "__" + SystName + "__plus" << ", ordered at " << mode << endl;
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__plus"))->Clone();
	hup[mode][j]->Add(h1);
      }
      if (verbose) cout << "Reading " << modeName << ", " << title + processName[j+1]+ "__" + SystName + "__minus" << ", ordered at " << mode << endl;
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
      if (j == 2){
        if (verbose) cout << "Adding " << modeName << ", " << title + processName[j+2]+ "__" + SystName + "__minus" << ", ordered at " << mode << endl;
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__minus"))->Clone();
	hdown[mode][j]->Add(h1);
      } 
    }
  } 
  if (verbose) cout << endl;
  
  datacard << setprecision(pre) << "jer       lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(SRbin) !=0){
	  double upv = (hup[i][j]->GetBinContent(SRbin) - hnom[i][j]->GetBinContent(SRbin))/hnom[i][j]->GetBinContent(SRbin);
	  double downv =  (hdown[i][j]->GetBinContent(SRbin) - hnom[i][j]->GetBinContent(SRbin))/hnom[i][j]->GetBinContent(SRbin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0){
	  double upv = (hup[i][j]->GetBinContent(CR1bin) - hnom[i][j]->GetBinContent(CR1bin))/hnom[i][j]->GetBinContent(CR1bin);
	  double downv =  (hdown[i][j]->GetBinContent(CR1bin) - hnom[i][j]->GetBinContent(CR1bin))/hnom[i][j]->GetBinContent(CR1bin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0){
	  double upv = (hup[i][j]->GetBinContent(CR2bin) - hnom[i][j]->GetBinContent(CR2bin))/hnom[i][j]->GetBinContent(CR2bin);
	  double downv =  (hdown[i][j]->GetBinContent(CR2bin) - hnom[i][j]->GetBinContent(CR2bin))/hnom[i][j]->GetBinContent(CR2bin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
   
  TString SystName = "MET";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
    
    for (int j = 0; j < 3; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      if (verbose) cout << "Reading " << modeName << ", " << title + processName[j+1]+ "__" + SystName + "__plus" << ", ordered at " << mode << endl;
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      if (j == 2){
        if (verbose) cout << "Adding " << modeName << ", " << title + processName[j+2]+ "__" + SystName + "__plus" << ", ordered at " << mode << endl;
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__plus"))->Clone();
	hup[mode][j]->Add(h1);
      }
      if (verbose) cout << "Reading " << modeName << ", " << title + processName[j+1]+ "__" + SystName + "__minus" << ", ordered at " << mode << endl;
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
      if (j == 2){
        if (verbose) cout << "Adding " << modeName << ", " << title + processName[j+2]+ "__" + SystName + "__minus" << ", ordered at " << mode << endl;
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__minus"))->Clone();
	hdown[mode][j]->Add(h1);
      }
    }
  } 
  if (verbose) cout << endl;
   
  datacard << setprecision(pre) << "met       lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(SRbin) !=0){
	  double upv = (hup[i][j]->GetBinContent(SRbin) - hnom[i][j]->GetBinContent(SRbin))/hnom[i][j]->GetBinContent(SRbin);
	  double downv =  (hdown[i][j]->GetBinContent(SRbin) - hnom[i][j]->GetBinContent(SRbin))/hnom[i][j]->GetBinContent(SRbin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0){
	  double upv = (hup[i][j]->GetBinContent(CR1bin) - hnom[i][j]->GetBinContent(CR1bin))/hnom[i][j]->GetBinContent(CR1bin);
	  double downv =  (hdown[i][j]->GetBinContent(CR1bin) - hnom[i][j]->GetBinContent(CR1bin))/hnom[i][j]->GetBinContent(CR1bin);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
      	if (hnom[i][j]->GetBinContent(CR2bin) !=0){
	  double upv = (hup[i][j]->GetBinContent(CR2bin) - hnom[i][j]->GetBinContent(CR2bin))/hnom[i][j]->GetBinContent(CR2bin);
          double downv =  (hdown[i][j]->GetBinContent(CR2bin) - hnom[i][j]->GetBinContent(CR2bin))/hnom[i][j]->GetBinContent(CR2bin);
          double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(pre) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(pre) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(pre) << 1+ upv<< "\t ";
	  else datacard << setprecision(pre) << "-\t ";
	} else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  datacard << setprecision(pre) << "pdf       lnN\t1.020/0.978     1.024/0.975     -     1.018/0.980      1.024/0.975     -     1.018/0.98     1.024/0.975     -    1.025/0.974    1.024/0.975    -     1.021/0.977     1.023/0.976    -     1.019/0.978     1.023/0.976    -    1.019/0.979     1.022/0.976     -     1.021/0.978     1.024/0.975    -     1.02/0.978      1.023/0.976     -" << endl;
   
  datacard << setprecision(pre) << "dynorm    lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==0 && j == 2)datacard << setprecision(pre) << "1.34\t ";
      if(i ==1 && j == 2)datacard << setprecision(pre) << "1.28\t ";
      if(i ==2 && j == 2)datacard << setprecision(pre) << "1.27\t ";
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==0 && j == 2)datacard << setprecision(pre) << "1.32\t ";
      if(i ==1 && j == 2)datacard << setprecision(pre) << "1.30\t ";
      if(i ==2 && j == 2)datacard << setprecision(pre) << "1.28\t ";
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==0 && j == 2)datacard << setprecision(pre) << "1.27\t ";
      if(i ==1 && j == 2)datacard << setprecision(pre) << "1.28\t ";
      if(i ==2 && j == 2)datacard << setprecision(pre) << "1.33\t ";
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  datacard << "ttxs      lnN\t -\t 1.06\t -\t -\t1.06\t -\t -\t 1.06\t -\t -\t 1.06\t -\t -\t 1.06\t -\t -\t 1.06\t - \t -\t 1.06 \t - \t -\t 1.06 \t - \t- \t1.06 \t -" << endl;
  datacard << "dyxs      lnN\t -\t -\t 1.05\t  -\t -\t1.05\t -\t -\t 1.05\t -\t -\t 1.05\t -\t -\t 1.05\t -\t -\t 1.05\t - \t -\t 1.05 \t - \t -\t 1.05 \t - \t- \t1.05 " << endl;
   
  datacard << setprecision(pre) << "# mc statistics for signal:" << endl;
  datacard << setprecision(pre) << "mcstatst1 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==0 && j == 0){
	if (hnom[i][j]->GetBinContent(SRbin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(SRbin)/hnom[i][j]->GetBinContent(SRbin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatst2 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==0 && j == 0){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR1bin)/hnom[i][j]->GetBinContent(CR1bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatst3 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==0 && j == 0){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR2bin)/hnom[i][j]->GetBinContent(CR2bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatst4 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==1 && j == 0){
	if (hnom[i][j]->GetBinContent(SRbin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(SRbin)/hnom[i][j]->GetBinContent(SRbin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){  
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  datacard << setprecision(pre) << "mcstatst5 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==1 && j == 0){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR1bin)/hnom[i][j]->GetBinContent(CR1bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){  
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  datacard << setprecision(pre) << "mcstatst6 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){  
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==1 && j == 0){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR2bin)/hnom[i][j]->GetBinContent(CR2bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
   
  datacard << setprecision(pre) << "mcstatst7 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==2 && j == 0){
	if (hnom[i][j]->GetBinContent(SRbin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(SRbin)/hnom[i][j]->GetBinContent(SRbin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
   
  datacard << setprecision(pre) << "mcstatst8 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==2 && j == 0){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR1bin)/hnom[i][j]->GetBinContent(CR1bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
   
  datacard << setprecision(pre) << "mcstatst9 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==2 && j == 0){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR2bin)/hnom[i][j]->GetBinContent(CR2bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  datacard << setprecision(pre) << "# mc statitics for background" << endl;
  datacard << setprecision(pre) << "mcstatot1 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 1){
	if (hnom[i][j]->GetBinContent(SRbin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(SRbin)/hnom[i][j]->GetBinContent(SRbin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  datacard << setprecision(pre) << "mcstatot2 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 1){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR1bin)/hnom[i][j]->GetBinContent(CR1bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot3 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 1){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR2bin)/hnom[i][j]->GetBinContent(CR2bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot4 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 1){
	if (hnom[i][j]->GetBinContent(SRbin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(SRbin)/hnom[i][j]->GetBinContent(SRbin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot5 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 1){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR1bin)/hnom[i][j]->GetBinContent(CR1bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot6 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 1){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR2bin)/hnom[i][j]->GetBinContent(CR2bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot7 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 1){
	if (hnom[i][j]->GetBinContent(SRbin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(SRbin)/hnom[i][j]->GetBinContent(SRbin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot8 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 1){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR1bin)/hnom[i][j]->GetBinContent(CR1bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot9 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 1){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR2bin)/hnom[i][j]->GetBinContent(CR2bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot10 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 2){
	if (hnom[i][j]->GetBinContent(SRbin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(SRbin)/hnom[i][j]->GetBinContent(SRbin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot11 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 2){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR1bin)/hnom[i][j]->GetBinContent(CR1bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot12 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 2){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR2bin)/hnom[i][j]->GetBinContent(CR2bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot13 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 2){
	if (hnom[i][j]->GetBinContent(SRbin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(SRbin)/hnom[i][j]->GetBinContent(SRbin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot14 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 2){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR1bin)/hnom[i][j]->GetBinContent(CR1bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  datacard << setprecision(pre) << "mcstatot15 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 2){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR2bin)/hnom[i][j]->GetBinContent(CR2bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  datacard << setprecision(pre) << "mcstatot16 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 2){
	if (hnom[i][j]->GetBinContent(SRbin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(SRbin)/hnom[i][j]->GetBinContent(SRbin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  datacard << setprecision(pre) << "mcstatot17 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 2){
	if (hnom[i][j]->GetBinContent(CR1bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR1bin)/hnom[i][j]->GetBinContent(CR1bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
   
  datacard << setprecision(pre) << "mcstatot18 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(pre) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 2){
	if (hnom[i][j]->GetBinContent(CR2bin) !=0) datacard << setprecision(pre) << 1 + (hnom[i][j]->GetBinError(CR2bin)/hnom[i][j]->GetBinContent(CR2bin)) << "\t ";
	else datacard << setprecision(pre) << "-\t ";
      }
      else datacard << setprecision(pre) << "-\t ";
    }
  }
  datacard << setprecision(pre) << endl;
  
  
  
}



