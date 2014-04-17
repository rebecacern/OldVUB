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

void syst(){

  bool verbose = false;
  
  //choose regions
  const int SRbin = 28;
  const int CR1bin = 30;
  const int CR2bin = 31;
  
  

  char myRootFile[300];
  sprintf(myRootFile,"results/Histos_cutbased_full.root");
  
  TFile *_file0 = TFile::Open(myRootFile);
  
  const int nProcess = 5;
  TString processName[nProcess] =  { "data", "twdr", "tt", "zjets", "others"};
  
  TH1F*  hdata [3];
  TH1F*  hnom [3][4];
  
  for (int process = 0; process < nProcess; process++){
    for (int i = 0; i < 3; i++){
      
      int mode = i;
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
	
      }    
    }
  } 
  if (verbose) cout << endl;
  
  cout << "Breakdown of the systematics " << endl;
  cout << "* Lumi " << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      2.2%" << endl;
  }
  cout << "* HLT (electron):" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      1.1%" << endl;
    else if (i == 1) cout << "mode:" << i << "      -" << endl;
    else if (i == 2) cout << "mode:" << i << "      1.5%" << endl;
  }
  cout << "* HLT (muon):" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      1.1%" << endl;
    else if (i == 1) cout << "mode:" << i << "      1.15%" << endl;
    else if (i == 2) cout << "mode:" << i << "      -" << endl;
  }
  cout << "* Electron ID:" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      2%" << endl;
    else if (i == 1) cout << "mode:" << i << "      -" << endl;
    else if (i == 2) cout << "mode:" << i << "      2%" << endl;
  }
  cout << "* Muon ID:" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      1%" << endl;
    else if (i == 1) cout << "mode:" << i << "      1%" << endl;
    else if (i == 2) cout << "mode:" << i << "      -" << endl;
  }
  
  cout.precision(2);
  
  TH1F*  hup [3][4];
  TH1F*  hdown [3][4];
  for (int i = 0; i < 3; i++){
  
    int mode = i;
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
  
  cout << "* Factorization/Normalization Scale " << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <3; j++){ 
      if(j < 2){
        double avn = (hnom[0][j]->GetBinContent(SRbin)+hnom[1][j]->GetBinContent(SRbin)+hnom[2][j]->GetBinContent(SRbin));
        double avup = (hup[0][j]->GetBinContent(SRbin)+hup[1][j]->GetBinContent(SRbin)+hup[2][j]->GetBinContent(SRbin));
	double avdown = (hdown[0][j]->GetBinContent(SRbin)+hdown[1][j]->GetBinContent(SRbin)+hdown[2][j]->GetBinContent(SRbin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << "\t[2j1t]" ;
    for (int j = 0; j <3; j++){ 
      if(j < 2){
        double avn = (hnom[0][j]->GetBinContent(CR1bin)+hnom[1][j]->GetBinContent(CR1bin)+hnom[2][j]->GetBinContent(CR1bin));
        double avup = (hup[0][j]->GetBinContent(CR1bin)+hup[1][j]->GetBinContent(CR1bin)+hup[2][j]->GetBinContent(CR1bin));
	double avdown = (hdown[0][j]->GetBinContent(CR1bin)+hdown[1][j]->GetBinContent(CR1bin)+hdown[2][j]->GetBinContent(CR1bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    } cout << "\t[2j2t]" ;
    for (int j = 0; j <3; j++){ 
      if(j < 2){
        double avn = (hnom[0][j]->GetBinContent(CR2bin)+hnom[1][j]->GetBinContent(CR2bin)+hnom[2][j]->GetBinContent(CR2bin));
        double avup = (hup[0][j]->GetBinContent(CR2bin)+hup[1][j]->GetBinContent(CR2bin)+hup[2][j]->GetBinContent(CR2bin));
	double avdown = (hdown[0][j]->GetBinContent(CR2bin)+hdown[1][j]->GetBinContent(CR2bin)+hdown[2][j]->GetBinContent(CR2bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << endl;
  }
  cout << endl;
  
  for (int i = 0; i < 3; i++){
    
    int mode = i;
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
 
  
  cout << "* ME/PS matching thresholds" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <3; j++){ 
      if(j == 1){
        double avn = (hnom[0][j]->GetBinContent(SRbin)+hnom[1][j]->GetBinContent(SRbin)+hnom[2][j]->GetBinContent(SRbin));
        double avup = (hup[0][j]->GetBinContent(SRbin)+hup[1][j]->GetBinContent(SRbin)+hup[2][j]->GetBinContent(SRbin));
	double avdown = (hdown[0][j]->GetBinContent(SRbin)+hdown[1][j]->GetBinContent(SRbin)+hdown[2][j]->GetBinContent(SRbin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << "\t[2j1t]" ;
    for (int j = 0; j <3; j++){ 
      if(j == 1){
        double avn = (hnom[0][j]->GetBinContent(CR1bin)+hnom[1][j]->GetBinContent(CR1bin)+hnom[2][j]->GetBinContent(CR1bin));
        double avup = (hup[0][j]->GetBinContent(CR1bin)+hup[1][j]->GetBinContent(CR1bin)+hup[2][j]->GetBinContent(CR1bin));
	double avdown = (hdown[0][j]->GetBinContent(CR1bin)+hdown[1][j]->GetBinContent(CR1bin)+hdown[2][j]->GetBinContent(CR1bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    } cout << "\t[2j2t]" ;
    for (int j = 0; j <3; j++){ 
      if(j == 1){
        double avn = (hnom[0][j]->GetBinContent(CR2bin)+hnom[1][j]->GetBinContent(CR2bin)+hnom[2][j]->GetBinContent(CR2bin));
        double avup = (hup[0][j]->GetBinContent(CR2bin)+hup[1][j]->GetBinContent(CR2bin)+hup[2][j]->GetBinContent(CR2bin));
	double avdown = (hdown[0][j]->GetBinContent(CR2bin)+hdown[1][j]->GetBinContent(CR2bin)+hdown[2][j]->GetBinContent(CR2bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << endl;
  }
  cout << endl;
  
  TH1F*  h [3][3];
  for (int i = 0; i < 3; i++){
    int mode = i;
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
   
   
  cout << "* DRDS scheme" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <3; j++){ 
      if(j == 0){
        double avn = (hnom[0][j]->GetBinContent(SRbin)+hnom[1][j]->GetBinContent(SRbin)+hnom[2][j]->GetBinContent(SRbin));
        double avup = (h[0][j]->GetBinContent(SRbin)+h[1][j]->GetBinContent(SRbin)+h[2][j]->GetBinContent(SRbin));
	double upv = (avup - avn)/avn;
	if (avn !=0){
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%\t ";
	} 
      }
    }
    cout << "\t[2j1t]" ;
    for (int j = 0; j <3; j++){ 
      if(j == 0){
        double avn = (hnom[0][j]->GetBinContent(CR1bin)+hnom[1][j]->GetBinContent(CR1bin)+hnom[2][j]->GetBinContent(CR1bin));
        double avup = (h[0][j]->GetBinContent(CR1bin)+h[1][j]->GetBinContent(CR1bin)+h[2][j]->GetBinContent(CR1bin));
	double upv = (avup - avn)/avn;
	if (avn !=0){
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%\t ";
	
	} 
      }
    } cout << "\t[2j2t]" ;
    for (int j = 0; j <3; j++){ 
      if(j == 0){
        double avn = (hnom[0][j]->GetBinContent(CR2bin)+hnom[1][j]->GetBinContent(CR2bin)+hnom[2][j]->GetBinContent(CR2bin));
        double avup = (h[0][j]->GetBinContent(CR2bin)+h[1][j]->GetBinContent(CR2bin)+h[2][j]->GetBinContent(CR2bin));
	double upv = (avup - avn)/avn;
	if (avn !=0){
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%\t ";
	
	} 
      }
    } 
    cout << endl;
  }
  cout << endl;
   
  TString SystName = "PU";
  
  for (int i = 0; i < 3; i++){
    int mode = i;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
   
    for (int j = 0; j < 4; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
    }
     
  } 
  if (verbose) cout << endl;
   
  cout << "* PU" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
        double avn = (hnom[i][j]->GetBinContent(SRbin));
        double avup = (hup[i][j]->GetBinContent(SRbin));
	double avdown = (hdown[i][j]->GetBinContent(SRbin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
	double avn = (hnom[i][j]->GetBinContent(CR1bin));
        double avup = (hup[i][j]->GetBinContent(CR1bin));
	double avdown = (hdown[i][j]->GetBinContent(CR1bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    } cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
	double avn = (hnom[i][j]->GetBinContent(CR2bin));
	double avup = (hup[i][j]->GetBinContent(CR2bin));
	double avdown = (hdown[i][j]->GetBinContent(CR2bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << endl;
  }
  cout << endl;

   
   
  cout << "* tt cross-section:" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      6%" << endl;
    else if (i == 1) cout << "mode:" << i << "      6%" << endl;
    else if (i == 2) cout << "mode:" << i << "      6" << endl;
  }
  
   
  TString SystName = "JES";
  
  for (int i = 0; i < 3; i++){
    int mode = i;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
   
    for (int j = 0; j < 4; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
    }
     
  } 
  if (verbose) cout << endl;
   
   
  cout << "* JES" << endl; 
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
        double avn = (hnom[i][j]->GetBinContent(SRbin));
        double avup = (hup[i][j]->GetBinContent(SRbin));
	double avdown = (hdown[i][j]->GetBinContent(SRbin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
	double avn = (hnom[i][j]->GetBinContent(CR1bin));
        double avup = (hup[i][j]->GetBinContent(CR1bin));
	double avdown = (hdown[i][j]->GetBinContent(CR1bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    } cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
	double avn = (hnom[i][j]->GetBinContent(CR2bin));
	double avup = (hup[i][j]->GetBinContent(CR2bin));
	double avdown = (hdown[i][j]->GetBinContent(CR2bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << endl;
  }
  cout << endl;

   
   
  TString SystName = "btag";
  
  for (int i = 0; i < 3; i++){
    int mode = i;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
   
    for (int j = 0; j < 4; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
    }
     
  } 
  if (verbose) cout << endl;
   
  cout << "* B-tagging" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
        double avn = (hnom[i][j]->GetBinContent(SRbin));
        double avup = (hup[i][j]->GetBinContent(SRbin));
	double avdown = (hdown[i][j]->GetBinContent(SRbin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
	double avn = (hnom[i][j]->GetBinContent(CR1bin));
        double avup = (hup[i][j]->GetBinContent(CR1bin));
	double avdown = (hdown[i][j]->GetBinContent(CR1bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    } cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
	double avn = (hnom[i][j]->GetBinContent(CR2bin));
	double avup = (hup[i][j]->GetBinContent(CR2bin));
	double avdown = (hdown[i][j]->GetBinContent(CR2bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << endl;
  }
  cout << endl;

  
  TString SystName = "JER";
  
  for (int i = 0; i < 3; i++){
    int mode = i;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
   
    for (int j = 0; j < 4; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
    }
     
  } 
  if (verbose) cout << endl;
   
  cout << "* JER" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
        double avn = (hnom[i][j]->GetBinContent(SRbin));
        double avup = (hup[i][j]->GetBinContent(SRbin));
	double avdown = (hdown[i][j]->GetBinContent(SRbin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
	double avn = (hnom[i][j]->GetBinContent(CR1bin));
        double avup = (hup[i][j]->GetBinContent(CR1bin));
	double avdown = (hdown[i][j]->GetBinContent(CR1bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    } cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
	double avn = (hnom[i][j]->GetBinContent(CR2bin));
	double avup = (hup[i][j]->GetBinContent(CR2bin));
	double avdown = (hdown[i][j]->GetBinContent(CR2bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << endl;
  }
  cout << endl;

   
  TString SystName = "MET";
  
  for (int i = 0; i < 3; i++){
    int mode = i;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
   
    for (int j = 0; j < 4; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
    }
     
  } 
  if (verbose) cout << endl;
   
  cout << "* Missing ET (unclustered)" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
        double avn = (hnom[i][j]->GetBinContent(SRbin));
        double avup = (hup[i][j]->GetBinContent(SRbin));
	double avdown = (hdown[i][j]->GetBinContent(SRbin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
	double avn = (hnom[i][j]->GetBinContent(CR1bin));
        double avup = (hup[i][j]->GetBinContent(CR1bin));
	double avdown = (hdown[i][j]->GetBinContent(CR1bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    } cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
	double avn = (hnom[i][j]->GetBinContent(CR2bin));
	double avup = (hup[i][j]->GetBinContent(CR2bin));
	double avdown = (hdown[i][j]->GetBinContent(CR2bin));
	if (avn !=0){
	  double upv = (avup - avn)/avn;
	  double downv =  (avdown - avn)/avn;
	  cout << processName[j+1] << ":" ;
	  cout <<  upv*100 << "%/" <<  downv*100<< "%\t ";
	} 
      }
    }
    cout << endl;
  }
  cout << endl;

  
  cout << "* pdf (signal only, for all of them check https://fblekman.web.cern.ch/fblekman/documents/pdf_systs_cut-and-count.txt):" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      -2.0/1.8%" << endl;
    else if (i == 1) cout << "mode:" << i << "      -2.0/1.8%" << endl;
    else if (i == 2) cout << "mode:" << i << "      -2.2/2.0%" << endl;
  }
  
  
  cout << "*DY re-weighted/non re-weighted" << endl;
  for (int i = 0; i < 3; i++){
   
    cout << "mode:" << i << " ";
    if(i ==2)cout <<  "34%\t ";
    if(i ==0)cout <<  "28%\t ";
    if(i ==1)cout <<  "27%\t ";
  
    cout << "\t[2j1t]" ;

    if(i ==2)cout <<"32%\t";
    if(i ==0)cout <<"30%\t";
    if(i ==1)cout <<"28%\t";
   
    cout << "\t[2j2t]" ;
   
    if(i ==2)cout <<  "27%\t ";
    if(i ==0)cout <<  "28%\t ";
    if(i ==1)cout <<  "33%\t ";
    cout << endl;
  }
  cout << endl;
  
  
  cout << "* MC statistics:" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
        double avn = (hnom[i][j]->GetBinContent(SRbin));
	if (avn !=0){
	  cout << processName[j+1] << ":" ;
	  cout <<  hnom[i][j]->GetBinError(SRbin)*100/avn << "%\t ";
	} 
      }
    }
    cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
        double avn = (hnom[i][j]->GetBinContent(CR1bin));
	if (avn !=0){
	  cout << processName[j+1] << ":" ;
	  cout <<  hnom[i][j]->GetBinError(CR1bin)*100/avn << "%\t ";
	} 
      }
    } cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if(j < 4){
        double avn = (hnom[i][j]->GetBinContent(CR2bin));
	if (avn !=0){
	  cout << processName[j+1] << ":" ;
	  cout <<  hnom[i][j]->GetBinError(CR2bin)*100/avn << "%\t ";
	} 
      }
    } 
	
    cout << endl;
  }
  cout << endl;

  cout << endl;
}
