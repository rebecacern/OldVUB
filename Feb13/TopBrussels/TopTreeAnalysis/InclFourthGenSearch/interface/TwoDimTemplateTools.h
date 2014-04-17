#ifndef TwoDimTemplateTools_h
#define TwoDimTemplateTools_h
#include "TopTreeAnalysis/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysis/Content/interface/Dataset.h"
#include "TopTreeAnalysis/Tools/interface/MultiSamplePlot.h"

#include <iostream>
#include <iomanip>
#include <string.h>
#include <sstream>
#include <map>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

using namespace std;

class TwoDimTemplateTools
{
 public:
   TwoDimTemplateTools();
   TwoDimTemplateTools(string identifier, string xvariable, int nbinsxvariable, string yvariable, int nbinsyvariable);
   TwoDimTemplateTools(const TwoDimTemplateTools & t);
   ~TwoDimTemplateTools();
   
   void SetDatasets(vector < Dataset* > datasets);
   void LoadTwoDimBinning(const string binningFileName); /**Loading of the binning for 2D distribution*/
   void Fill_for2DBinning(double xvariable_value,float yvariable_value,float fillweight); /** for the step when the 2D binning is to be created*/
   void Fill(double xvariable_value,float yvariable_value,float fillweight,unsigned int datasetindex); /** for the step after the 2D binning is already created*/
   void Write_for2DBinning(const string binningFileName); /** for the step when the 2D binning is to be created*/
   void Write(TFile* fout,TDirectory* th1dir); /** for the step after the 2D binning is already created*/
   void Convert2Dto1D(string postfix);
 
 private:
   string xvariable_;
   string yvariable_;
   int nbinsxvariable_;
   int nbinsyvariable_;
   string xvariable_XaxisLabel_;
   string xvariable_YaxisLabel_;
   int nfixedbins_;
   
   double *xbins;
   
   string binningFileName_;
   string identifier_;
   
   map< int, MultiSamplePlot* > MSPlot_xvariableBins; /** the integer indicates the HT bin, and the plot that corresponds to it is the Mtop distribution in that HT bin*/
   map< int, vector < pair < TH1F*, Dataset* > > > Histos_xvariableBins; /**the integer indicates the HT bin, the histos are Mtop distributions per dataset*/
   map<string,vector< float > > VariableValuesMap; /** the string indicates the variable, the vector is a list of the values of this variables for all events*/
   vector< float > eventweightvector; /** a vector of the event weights, needed to get flat binning*/
   map< string, int > nbinsMap; /** the string indicates the variable, the integer the number of bins in the 2D distribution*/
   
   map<string,pair<float,float > > historanges;   
   
   vector < Dataset* > datasets_;
   
   vector< pair < TH1F*, Dataset* > > templates1D; /** all templates*/
   vector< pair < TH1F*, Dataset* > > templates1D_forMS; /** templates for the multisampleplot that is going to be written*/
   
};

#endif
