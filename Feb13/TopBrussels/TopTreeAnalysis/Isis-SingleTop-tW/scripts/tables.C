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
void tables(int mode = 0){

  double lumi = luminosity;
  
  if (mode == 0 )        lumi = 5085.246;
  else if ( mode == 1)   lumi = 4919.924;
  else if ( mode == 2)   lumi = 4895.249;
  
  
  char myTexFile[300];
  sprintf(myTexFile,"tables/table_%d_%dpb.tex", mode, lumi);
  ofstream salida(myTexFile); 
  
  char myRootFile[300];
  sprintf(myRootFile,"results/an_%dpb_%d.root", lumi, mode);
  TFile *_file0 = TFile::Open(myRootFile);
  
  const int nProcess = 10;
  TString processName[nProcess] = { "twdr", "tt", "di", "wjets", "zjets", "st", "qcd_mu", "others", "data", "mc"};
  TString processLabel[nProcess] = { "\\textbf{$tW$}", "\\textbf{$t \\bar{t}$}", "\\textbf{$WW$, $WZ$, $ZZ$}", "\\textbf{$W+jets$}", "\\textbf{$Z+jets$}",
			       "\\textbf{$single$ $top$}","\\textbf{$QCD$}", "\\textbf{$Other$}", "\\textbf{$MuEG$}", "\\textbf{$MC$}"};
  
  if (mode == 1)        processLabel[8] = "\\textbf{$DoubleMu$}";
  else if (mode == 2)   processLabel[8] = "\\textbf{$DoubleElectron$}";
  
  TString cutLabel[8] = { "blank0", "blank1", "Lepton Sel.", "Inv. Mass", "$E_{T}^{miss}$", "1 Jet", "b-tagging", "$H_{T}$"};
  
  TH1F*  h [nProcess];
  for(int i=0; i<10; i++){
    if (i == 9) {
     h[i] =  (TH1F*)h[0]->Clone();
     h[i]->Add(h[1]);
    // h[i]->Add(h[2]);
     //h[i]->Add(h[3]);
     h[i]->Add(h[4]);
     //h[i]->Add(h[5]);
    // h[i]->Add(h[6]);
      h[i]->Add(h[7]);
    }
    else h[i] = (TH1F*) _file0->Get("cuts_"+processName[i]);
    // Lepton ID and HLT SF
    // if (mode == 0 && i != 8) h[i]->Scale(0.97713);
    // if (mode == 1 && i != 8) h[i]->Scale(0.910067);
    // if (mode == 2 && i != 8) h[i]->Scale(0.945736);
  }
  
  double vectorValue[nProcess][17][4];
  for (int i = 0; i < 16; i++){
    for (int j = 0; j < 10; j++){
      vectorValue[j][i][0] = h[j]->GetBinContent(i);
      vectorValue[j][i][1] = precision(h[j]->GetBinError(i));
      vectorValue[j][i][2] = h[j]->GetBinError(i);
    }  
  }

  salida << "\\documentclass{cmspaper}" << endl;
  salida << "\\begin{document}" << endl;
  salida << endl;
  salida << endl;
  
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|c|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 0; i < 7; i++){
    salida << " & " << processLabel[i] ;
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
  
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 7; j++){
      if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] != 0){
        salida << " & $\\leq$ " << setprecision(vectorValue[j][i-1][1]) << 2*vectorValue[j][i-1][2];
      } else if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] == 0){
        salida << " & $-$ " ;
      } else {
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
      }
    }
    salida <<  " \\\\  " << endl; 
  }
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;

  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 0; i < 9; i++){
    if (i == 0 || i == 1 || i ==4 || i == 7) salida << " & " << processLabel[i] ;
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
  
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 8; j++){
      if (j == 2) j = 4;
      if (j == 5) j = 7;
      if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] != 0){
        salida << " & $\\leq$ " << setprecision(vectorValue[j][i-1][1]) << 2*vectorValue[j][i-1][2];
      } else if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] == 0){
        salida << " & $-$ " ;
      } else {
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][2];
	
      } 
    }
    salida <<  " \\\\  " << endl;
  }
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 8; i < 10; i++){
    salida << " & " << processLabel[i] ;
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
  
  for (int i=2; i < 8; i++){
    salida << cutLabel[i] << "	";
    for (int j = 8; j < 10; j++){
      if (j == 8){
	if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] != 0){
	  salida << " & $\\leq$ " << setprecision(vectorValue[j][i-1][1]) << 2*vectorValue[j][i-1][2];
	} else if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] == 0){
	  salida << " & $-$ " ;
	} else {
	  salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][0] ; 
	} 
      } else {
	if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] != 0){
	  salida << " & $\\leq$ " << setprecision(vectorValue[j][i-1][1]) << 2*vectorValue[j][i-1][2];
	} else if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] == 0){
	  salida << " & $-$ " ;
	} else {
	  salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][0] ; 
	  salida << " $\\pm $"  << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][2];
	} 
      }
      
    }
    salida <<  " \\\\  " << endl;
  }
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
  
  salida << "\\end{document}" << endl;
  
}


double normalization(double nevents, double xsec, double lumi){
  
  if (nevents !=0) return lumi*xsec/nevents;
  else return 1;
  
}

double precision(double error){
  
  int precisionValue;
  double factErr = 0; 
  int iN = 0;
  if (error == 0 || error >= 1) precisionValue = 0;
  else if (error < 1) {
    iN = 0;
    factErr = 0; 
    while (factErr < 1){
      factErr = error*(10**iN);
      iN++;  
    }
    precisionValue = iN-1;
  }
  
  if (factErr > 9.5) precisionValue-=1;
  
  return precisionValue;
  
}

double errorEfi (double efi, double errTotal, double totalEvents, double error, double number){
  
  if (totalEvents !=0 && number != 0) return efi*((errTotal/totalEvents)+(error/number));
  else return 0;
  
}

double efficiency (double finalevents, double totalevents){
  
  if (totalevents !=0) return finalevents*100/totalevents;
  else return 0;
  
}
