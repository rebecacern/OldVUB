// simple file to execute the b-tag efficiency measurement function in myNTupleAnalyzer.cc

#include <TString.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>
#include "../interface/myNTupleAnalyzer.h"
using namespace std;

//--------- Configuration
//name of the outrootfiles and outrootpath should be defined here

TString *inName = new TString("./TTrees/");

//fixed for all samples
//double chisqCut[24]={1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,15.,20.,30.,50.,9999999999.};
//double chisqCut[28]={1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,15.,20.,30.,50.,9999999999.,100.,500.,1000.,5000.,};
//double chisqCut[24]={1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.,20.,25.,30.,40.,50.,60.,70.,85.,100.,500.,1000.,5000.,9999999999.};
//double chisqCut[24]={2.,3.,4.,5.,6.,7.,8.,9.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,70.,80.,90.,100.,9999999999.};
double chisqCut[24]={30.,40.,50.,60.,70.,80.,125.,100.,110.,120.,130,140,150,160,170,180,190,200,220,240,260,300,350,9999999999.};

//double chisqCut[24]={30.,40.,50.,60.,70.,80.,90.,100.,150.,200.,250,300,350,400,450,500,600,700,800,900,1000,2000,3000,5000};

int verbosity=0; //default

int nSystematic=0; //default (nominal samples) 

// GOOD
/*int leftlimit = 50;
 int centerleftlimit = 140;
 int centerrightlimit = 140;
 int rightlimit = 240;*/

 /*int leftlimit = 70;
  int centerleftlimit = 170;
  int centerrightlimit = 170;
  int rightlimit = 320;*/ 

int leftlimit = 70;
int centerleftlimit = 170;
int centerrightlimit = 170;
int rightlimit = 300;


/*int leftlimit = 50;
 int centerleftlimit = 140;
 int centerrightlimit = 160;
 int rightlimit = 260;
 */

//good for e+jets
//int leftlimit = 80;
//int centerleftlimit = 160;
//int centerrightlimit = 160;
//int rightlimit = 220;

//////////////////////////////////////////////////
/////
///// Settings for running the btag efficiency with fixed bins (remark, the boolean doVarBins=false in myNtupleAnalyzer.cc)
/////
//////////////////////////////////////////////////

//TString *outName = new TString("TTfastall_ST_W_Z_1fb_novarbinbtag_novarbinscrw_mlb90_300_chisqchange_"); //for this one not yet in btag eff. in pteta-bins
/*TString *outName = new TString("TTfastall_ST_W_Z_1fb_novarbinbtag_mlb90_300_chisqchange_"); //for this one not yet in btag eff. in pteta-bins
int runSamples[6]={16,22,23,24,26,28};         // to run on fastsim ttbar, single top and W/Zjets madgraph samples
////int runSamples[6]={0,22,23,24,26,28};      // to run on fullsim ttbar, single top and W/Zjets madgraph samples
int nRunSamples=6;
int doJESsample=6;*/

//to make the root files per sample
//TString *outName = new TString("Sample0_1fb_novarbinbtag_mlb90_300_chisqchange_");
//TString *outName = new TString("Sample1_1fb_novarbinbtag_mlb90_300_chisqchange_");
//TString *outName = new TString("Sample2_1fb_novarbinbtag_mlb90_300_chisqchange_");
//TString *outName = new TString("Sample3_1fb_novarbinbtag_mlb90_300_chisqchange_");
//TString *outName = new TString("Sample4_1fb_novarbinbtag_mlb90_300_chisqchange_");
//TString *outName = new TString("Sample5_1fb_novarbinbtag_mlb90_300_chisqchange_");
/*TString *outName = new TString("Sample6_1fb_novarbinbtag_mlb90_300_chisqchange_");
int runSamples[7]={16,22,23,24,26,28,61};
int nRunSamples=7;
int doJESsample=6;*/

/////
//////////////////////////////////////////////////

//////////////////////////////////////////////////
/////
///// Settings for running the btag efficiency with variable bins (remark, the boolean doVarBins=true in myNtupleAnalyzer.cc)
/////
//////////////////////////////////////////////////

TString *outName = new TString("DATA_"); //for this one not yet in btag eff. in pteta-bin

int runSamples[15]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};         // to run on fastsim ttbar, single top and W/Zjets madgraph samples/
int nRunSamples=10;

//int runSamples[1]={7};      // to run on fullsim ttbar, single top and W/Zjets madgraph samples
//int nRunSamples=1;

int doJESsample=6;

//bool doFfromMC = true; 
bool doFfromMC = false; 
// If this boolean is true, I should make sure I have 20 bins in pt and eta so that I can obtain the fit function for the MC F values.
// This boolean needs to be false when I have 5 bins, and I obtain the 'data' value of \epsilon_{b}


/////
//////////////////////////////////////////////////

//////////////////////////////////////////////////
/////
///// Settings for running the btag efficiency with fixed bins for correlation plots
/////
//////////////////////////////////////////////////

/*TString *outName = new TString("TTfastall_ST_W_Z_1fb_novarbinbtag_mlbnolimits_5_10_chisqchange_"); //for this one not yet in btag eff. in pteta-bins
int runSamples[6]={16,22,23,24,26,28};         // to run on fastsim ttbar, single top and W/Zjets madgraph samples
int nRunSamples=6;
int doJESsample=6;
int leftlimit = 0;
int centerleftlimit = 500;
int centerrightlimit = 500;
int rightlimit = 500;
*/
/////
//////////////////////////////////////////////////


//////////////////////////////////////////////////
/////
///// Settings for the making the pseudo-experiments
/////
//////////////////////////////////////////////////

/*//TString *outName = new TString("TTfastall_ST_W_Z_1fb_mlb90_300_PseudoExpRR2_alter_chisqchange_"); //for this one not yet in btag eff. in pteta-bins
TString *outName = new TString("TTfastall_ST_W_Z_1fb_mlb90_300_PseudoExpRR2_correct_chisqchange_"); //for this one not yet in btag eff. in pteta-bins
//TString *outName = new TString("TTfastall_ST_W_Z_1fb_mlb90_300_PseudoExpRR_Fnorw_chisqchange_"); //for this one not yet in btag eff. in pteta-bins
int runSamples[6]={16,22,23,24,26,28};
int nRunSamples=6;
int doJESsample=6;
*/

float desiredIntLum=30000.0; 

//float desiredIntLum=4568.68; 

int nPseudoExp=1;//= is there to run on all data
bool doPseudoExp=false; //to do the pseudo-exps
int nBtag=-1;

/*double desiredIntLum=1000;
int nPseudoExp=3000;
bool doPseudoExp=true; //to do the pseudo-exps
*/
/////
//////////////////////////////////////////////////


//////////////////////////////////////////////////
/////
///// Settings for the making the JES
/////
//////////////////////////////////////////////////

bool doJESchange=false;
//bool doJESchange=true;
double JESfactor[13]={0.8,0.85,0.9,0.925,0.95,0.975,1,1.025,1.05,1.075,1.1,1.15,1.2};
////int doJESsample=-1;

//////////////////////////////////////////////////


//////////////////////////////////////////////////
/////
///// Settings for the making of Fbias vs effBias.
/////
//////////////////////////////////////////////////

//TString *outName = new TString("x_TTfastall_Fbias_pteta_ST_W_Z_1fb_5pt5etabins_chisqchange_23_Fbias_");
//int runSamples[6]={16,22,23,24,26,28};// to run only on fastsim ttbar, single top and W/Zjets madgraph samples
//int nRunSamples=6;// to run only on fullsim ttbar, single top and W/Zjets madgraph samples
//int doJESsample=6;

//variables for introducing an arbitrary bias on F
int inFbias=0;
int nBiases=11;
double Biases[11]={0.85,0.88,0.91,0.94,0.97,1,1.03,1.06,1.09,1.12,1.15};

//variables for introducing an arbitrary deficit on the amount of non b jets in the sample
bool doBackgroundFraction=false;
int inBackgroundFraction=0;
int nbackgroundFraction=7;
double backgroundFraction[7]={0.4,0.5,0.6,0.7,0.8,0.9,1};

/////
//////////////////////////////////////////////////





//////////////////////////////////////////////////
/////
///// Settings for the making the mtop
/////
//////////////////////////////////////////////////

/*//TString *outName = new TString("x_TTfastall_totsample_mtop168_SWrwchange_pteta_ST_W_Z_1fb_5pt5etabins_chisqchange_");
//int runSamples[1]={17};
TString *outName = new TString("x_TTfastall_totsample_Wmore_SWrwchange_pteta_ST_W_Z_1fb_5pt5etabins_chisqchange_");
int runSamples[6]={16,22,23,24,26,28};
int nRunSamples=6;
int doJESsample=6;*/

/////
//////////////////////////////////////////////////



//int leftlimit = 0;
//int centerleftlimit = 160;
//int centerrightlimit = 160;
//int rightlimit = 500;

/*double leftlimitperc=0.220497; //obtained with 90/160/210 with {0,22,23,24,26,28}
double centerlimitperc=0.592576;
double rightlimitperc=0.754151;*/

double leftlimitperc=0.225444; //obtained with 90/160/210 with {08}
double centerlimitperc=0.612447; // not used
double rightlimitperc=0.77487; // not used


/*int leftlimit = 80;
int centerleftlimit = 155;
int centerrightlimit = 155;
int rightlimit = 190;
*/



bool doSCreweigh=true; //to use the reweighted parameters for the non b ratio between left and right coming from the mlb distribution in control sample
bool doTwoLights=true; //to use the two light jets to fill mlb

bool useFit=true; // to use the fitted ratio between left and right pt to reweigh the btag distro (false uses non-fitted ratio)
//bool useFit=false; // to use the fitted ratio between left and right pt to reweigh the btag distro (false uses non-fitted ratio)
bool do2D=false; //to tod do the reweighing in 2D for the left-right pt/eta reweighing (useFit should be switched to false (don't know what happens if it's true))

//bool do2Dcontrol=false; //to do the reweighing in 2D for the signal-control pt/eta reweighing 
bool do2Dcontrol=true; //to do the reweighing in 2D for the signal-control pt/eta reweighing 

//bool doPtEtaBin=false; //to do the measurement in bins of pt/eta
bool doPtEtaBin=false; //to do the measurement in bins of pt/eta
double chisqCutApplied=6;

//-------- execution of the module



double chisqCutNumber=-1;


int inBin;
int outBin;

int inRunSample=-1;



int main(int argc, char* argv[]){

  time_t curr=time(0);
  cout << "current start time is: " << ctime(&curr) <<endl;

  //if(argc!=10) cout << "ERROR: you provided " << argc-1 << " argument(s) while 10 argument is needed" << endl;;

  verbosity=atoi(argv[1]);

  chisqCutNumber=chisqCutApplied;
  chisqCutNumber=atoi(argv[2]);
  //doJESsample=atoi(argv[3]);
    int fitMode=atoi(argv[3]);
  inFbias=atoi(argv[4]);
  inBackgroundFraction=atoi(argv[5]);
  inRunSample=atoi(argv[6]);
  inBin=atoi(argv[7]);
  outBin=atoi(argv[8]);
	 
 	if (argc == 10)
		nSystematic = (int)atoi(argv[9]);
	
	else if (argc > 10) {
		
		doPseudoExp = (bool)atoi(argv[9]);
		
		if (doPseudoExp) {
			desiredIntLum = atof(argv[10]);
			nPseudoExp = atoi(argv[11]);
		}
        
        if (doPseudoExp && argc > 12)
            nBtag = atoi(argv[12]);
	}
	
  cout << "chiSquare Cut Number: " << chisqCutNumber << endl;
  cout << "chiSquare Cut : " << chisqCut[(int) chisqCutNumber] << endl;

  bool doFbiasTest=false;
  if(inFbias>-1) doFbiasTest=true;
  cout << "doFbiasTest : " << doFbiasTest << endl;
  cout << "doFbiasTest in Number : " << inFbias << endl;
  if(doFbiasTest) cout << "doFbiasTest Value : " << Biases[inFbias] << endl;

  if(inBackgroundFraction>-1) doBackgroundFraction=true;
  cout << "doBackgroundFraction : " << doBackgroundFraction << endl;
  cout << "doBackgroundFraction : " << inBackgroundFraction << endl;
  if(doBackgroundFraction) cout << "doBackgroundFraction : " << backgroundFraction[inBackgroundFraction] << endl;

  cout << "desiredIntLum : " << desiredIntLum << endl;

  cout << "outName: " << *outName << endl;
  
  std::string chisqCutName;
  std::string FbiasName;
  std::stringstream out;
  
  if(!doFbiasTest){
    out << chisqCutNumber;
    chisqCutName = out.str();
    cout << "chisqCutName " <<chisqCutName << endl;
    outName->Append(chisqCutName);
    outName->Append("/");
  } else {
    out << inFbias;
    FbiasName = out.str();
    cout << "FbiasName " << FbiasName<< endl;
    outName->Append(FbiasName);
    outName->Append("/");
  }

  
  cout << "inBin " << inBin << endl;
  cout << "outBin " << outBin << endl;
  
  std::string inBinName;
  std::stringstream out_inBin;
  out_inBin << inBin;
  inBinName = out_inBin.str();

  std::string outBinName;
  std::stringstream out_outBin;
  out_outBin << outBin;
  outBinName = out_outBin.str();


  //outName->Append(inBinName);
  ////outName->Append(outBinName);
  //outName->Append("/");
 
  cout << "outName: " << *outName << endl;

  int tmp_nRunSamples=0;
  
  if(inRunSample>-1){
    runSamples[0]=runSamples[inRunSample];
    tmp_nRunSamples=1;
  } else {
    tmp_nRunSamples=nRunSamples;
  }

  for(int i=0; i<tmp_nRunSamples; i++){
    cout << "#RunSamples: " << tmp_nRunSamples << " -- sample to run: " << runSamples[i] << endl;
  }

	cout << "desired lumi!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << desiredIntLum << endl;
  myNTupleAnalyzer myAnalyzer(inName,outName,runSamples,tmp_nRunSamples,chisqCut[(int) chisqCutNumber],inBin,outBin,desiredIntLum,nPseudoExp,doPseudoExp);
    
    chisqCut[(int) chisqCutNumber]=175;
    myNTupleAnalyzer myAnalyzer2(inName,outName,runSamples,tmp_nRunSamples,chisqCut[(int) chisqCutNumber],inBin,outBin,desiredIntLum,nPseudoExp,doPseudoExp);
    //myNTupleAnalyzer myAnalyzer2(inName,outName,runSamples,tmp_nRunSamples,750.0,inBin,outBin,desiredIntLum,nPseudoExp,doPseudoExp);

  if(doFbiasTest){ myAnalyzer.setFMCBias(Biases[inFbias]); }
  if(doBackgroundFraction) myAnalyzer.setBackgroundFraction(backgroundFraction[inBackgroundFraction]);
  
  //double FMCBiases[20];
  //double effMCBiases[3][20];

  int nWP=18;
  double *eff_ = new double[18];
  double *effMCVal_ = new double[18];
  double *effRRMCVal_ = new double[18];
  double *effVal_ = new double[18];
  double *effRRVal_ = new double[18];

  //double *effMCValBias_ = new double[18];
  //double *effRRMCValBias_ = new double[18];

  double FMC=0;  
  double FMCBias_=0;  

  // fout = new TFile(outRootName,"RECREATE");
  //cout << "}-- Created output root file: " << outRootName << endl; 

  double tempFactor=-1;
  if(doFbiasTest) {tempFactor=1;} 
  if(!doJESchange){
    tempFactor=1;
  } else {
    tempFactor=JESfactor[doJESsample];
  }

  //cout << "JES-factor: " << tempFactor << endl;
  
 
  bool doNewF=false;
  if(doFfromMC) doNewF=false;
  //myAnalyzer.run(verbosity, leftlimit, centerleftlimit, centerrightlimit, rightlimit, doSCreweigh, doTwoLights,useFit, do2D, do2Dcontrol,doPtEtaBin,doJESchange,tempFactor,doNewF);

    //int decay=1;
    //fitMode=0; // choose fit variable -> 0: m_lj 1: M3 2: 2D fit of (m_lj,M3)
 
    /*if (decay==1) {
        
        cout << "Warning: E+jets selected, changing left and right defs" << endl;
        
        leftlimit = 80;
        centerleftlimit = 160;
        centerrightlimit = 160;
        rightlimit = 220;
        
    }*/
    
  //int *percentiles = new int(3);
    //mu+jets
    myAnalyzer.run(verbosity, leftlimit, centerleftlimit, centerrightlimit, rightlimit, doSCreweigh, doTwoLights,useFit, do2D, do2Dcontrol,doPtEtaBin,doJESchange,tempFactor,doNewF,leftlimitperc,centerlimitperc,rightlimitperc,1.,doFfromMC,nSystematic,0,fitMode,nBtag);
    
    //cout << "Warning: E+jets selected, changing left and right defs" << endl;

    //leftlimit = 70;
    //centerleftlimit = 150;
    //centerrightlimit = 150;
    //rightlimit = 250;
    //e+jets
    
    //myAnalyzer2.run(verbosity, leftlimit, centerleftlimit, centerrightlimit, rightlimit, doSCreweigh, doTwoLights,useFit, do2D, do2Dcontrol,doPtEtaBin,doJESchange,tempFactor,doNewF,leftlimitperc,centerlimitperc,rightlimitperc,1.,doFfromMC,nSystematic,1,fitMode,nBtag);
  
  //  myAnalyzer.getPercentiles(percentiles);
  //cout << percentiles[0] << " " << percentiles[1] << " " << percentiles [2]<< endl;

  //myAnalyzer.run(verbosity, percentiles[0], percentiles[1], percentiles[1], percentiles[2], doSCreweigh, doTwoLights,useFit, do2D, do2Dcontrol,doPtEtaBin,doJESchange,tempFactor,doNewF,-1,-1,-1,2.,doFfromMC);
  
  
  time_t currEnd=time(0);
  cout << "current end time is: " << ctime(&currEnd) <<endl;
  
  time_t currDiff=currEnd-curr;
  cout << "bTag measurement calculation time is: " << currDiff << " seconds = " << (double) currDiff/60 << " minutes" <<endl;
  
}




//possible extentions
/* int nBinsVarLR[2]={50,20};
 double lowRangVarLR[2]={0,0};
 double upRangeVarLR[2]={500,2.4};
 int nBinsVarSC[2]={50,10}; //for 1D, then the 10 is not taken into account
// int nBinsVarSC[2]={20,10};
 double lowRangVarSC[2]={0,0};
 double upRangeVarSC[2]={500,2.4};

//bool doVarBins=false;
bool doVarBins=true;
bool doShift=false; //warning: when switching this to true make sure the correct shift is applied!
//bool doShift=true;

bool dobValueBdisc=true;

bool doCoutBinning = false;


//int lumi
//int nPseudoExp
bool doPseudoExp=false; //to do the pseudo-exps
*/
