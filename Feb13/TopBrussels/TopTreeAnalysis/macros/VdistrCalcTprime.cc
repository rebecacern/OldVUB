#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "TF1.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TObject.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <fstream>
#include <string>
#include <fstream>
#include "TLatex.h"
#include "TString.h"
#include "TInterpreter.h"
#include <fstream>
#include "TH1.h"
#include "TGraphSmooth.h"
#include "TSystem.h"
//user code
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../StatProcedure/interface/WeightProbaCalculator.h"

//#include "Style.C"

using namespace std;

void DrawSmooth(Int_t pad, const char *title, const char *xt, const char *yt);

int main(int argc, char*argv[]){
    
    cout << "	" << "========================================================" << endl;
    cout << "	" << "TOWARDS EXCLUSION PLOT FOR TPRIME (calculating 1-beta, ...)" << endl;
      
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetOptFit(1111);
   
    TCanvas *c = new TCanvas("c","Exclusion limit",0,0,700,700);
    TCanvas *cc = new TCanvas("cc","Exclusion limit2",0,0,700,700);
    
    TGraph *grout;
    TGraph *grintest = new TGraph(4);
   
    //xml file
    const char *xmlfile = "../config/myTprimeconfig.xml";
    string thisROUNDdir = "./ROUND13/"; //directory for (most of??) the files of this round...
    AnalysisEnvironment anaEnv;
    cout << "Loading environment ..." << endl;
    AnalysisEnvironmentLoader anaLoad (anaEnv, xmlfile);
    vector<float> FractionHWEvts = anaEnv.FractionHWEvts;
    string SignificanceCombination("SUM_LRatioMinded"); // to combine squared bin significances to get a combined weight for the event
    //string SignificanceCombination("SUM");
    bool applyondata = true;
    bool likelihood = true; //if true: inverses how the p-values are calculated... (test statistic S+B less than for B on average). For V-values this should be false.
    
    //float xsectrange[] = {0.1,6}; //15 //25 //e.g. {0.5,4} for 200/pb (ROUND23)
    float xsectrange[] = {0.1,10}; //{0.1,15}
    float xsectstep = 0.3; //e.g. 0.3 for/200pb (ROUND23)
    //float xsect = xsectrange[0], xsectstep = 1;
    //float massestprimeHardCoded[6] = {350,400,450,500,550,600};
    //float xsection_theory[] = {2.94,1.30,0.617,0.310,0.162,0.088}; //NLO cross sections bprime pair production, from paper EXO-(...) but should be also found in another paper (because they just use it, they didn't obtain it)
    float massestprimeHardCoded[4] = {350,400,450,500};
//    float massestprimeHardCoded[2] = {350,400};
    float xsection_theory[] = {2.94,1.30,0.617,0.310};
//    float xsection_theory[] = {2.94,1.30};
    vector<float> massestprime;
    vector<float> crosssectionstprime; //in pb
    vector<float> NofEventsTprimeTopTree; //= {22166,20701,17817,20601}; //still correctly handled???
    unsigned int numberofmasses = 4; //4
    float xaxisarray[numberofmasses]; //masses tprime
    //float NofEventsTopTreetprime[4] = {22166,20701,17817,20601};
    //float NofEventsTTreetprime[4] = {1844,1818,1682,1922};
//    float RefSelEff_tprime[4] = {0.08319,0.08782,0.09440,0.093296}; //RefSelV4 of top+lepton group (mu-channel)
    float RefSelEff_tprime[4] = {1,1,1,1}; //dummy
    
    /////////////////////
    // Load Datasets
    /////////////////////
    TTreeLoader treeLoader;
    cout << " - Load datasets ..." << endl;
    vector < Dataset* > datasets;
    vector < Dataset* > datasets_allconfigMC;
    treeLoader.LoadDatasets (datasets, xmlfile);
    
    int massiter = 0;
    for (unsigned int d = 0; d < datasets.size (); d++) {
      string dataSetTitle = datasets[d]->Title();
      string dataSetName = datasets[d]->Name();
      if(dataSetName.find("prime")<=dataSetName.size() || dataSetTitle.find("prime")<=dataSetTitle.size()){
	massestprime.push_back(datasets[d]->Mass());
	xaxisarray[massiter] = datasets[d]->Mass();
	massiter++;
	//crosssectionstprime.push_back(datasets[d]->Xsection());
	NofEventsTprimeTopTree.push_back(datasets[d]->eventTree()->GetEntriesFast());
	//cout<<"     nofevents t' samples? "<<datasets[d]->eventTree()->GetEntriesFast()<<endl;
      }
    }
    
    float xsect = xsectrange[0];
    crosssectionstprime.clear();
    while(xsect>=xsectrange[0] && xsect<=xsectrange[1]){//for example same as in GOFanalysis //25 for lumi = 36.1, I have 10 for 250 lumi
       crosssectionstprime.push_back(xsect);
       xsect = xsect + xsectstep;
    }    
    
    float alfa = 0.05; //determines the confidence level for exclusion...
    float cl = (1 - alfa)*100; //confidence level (in %)

    unsigned int nofFractionsHWEvts = FractionHWEvts.size();
    float FractionHWEvtsArray[nofFractionsHWEvts]; //just in ARRAY format
    TH2F *hpoweroftest_tprime[nofFractionsHWEvts]; //= new TH2F ("1-#beta", "#sigma (pb)", m0Nbins, m0min, m0max, m12Nbins, m12min, m12max);
    TH2F *hpvaluedata_tprime[nofFractionsHWEvts];
    TH2F *hNofSelEvts_tprime;
    TH2F *hSignSelEvts_tprime;
    
    //attempt to use CLs method (reference: A.L. Read, Modified Frequentist analysis of Search Results (The CLs method))    
    TH2F *hCLs_ExpMedian_tprime[nofFractionsHWEvts]; 
    TH2F *hCLs_Exp16thPerc_tprime[nofFractionsHWEvts];    
    TH2F *hCLs_Exp84thPerc_tprime[nofFractionsHWEvts]; 
    TH2F *hCLs_Data_tprime[nofFractionsHWEvts];
    
    TFile *_output;
    TFile* new_output ;

    //for the plots...
    float massmin = massestprime[0] - (massestprime[1] - massestprime[0])/2; //300 - (350-300)/2 = 275    
    float massmax = massestprime[massestprime.size()-1] + (massestprime[1] - massestprime[0])/2; //450 + (350-300)/2 = 475
    //int massNbins = 4;
    float xsectionmin = xsectrange[0];
    float xsectionmax = xsectrange[1]; //26
    int xsectionNbins = int((xsectrange[1] - xsectrange[0])/xsectstep); //10; // 25 //dangerous, when not ints... to be checked in future!!
    cout<<"massmin = "<<massmin<<", massmax = "<<massmax<<", numberofmasses = "<<numberofmasses<<", xsectionmin = "<<xsectionmin<<", xsectionmax = "<<xsectionmax<<", xsectionNbins = "<<xsectionNbins<<endl;

    float xsection_power50[nofFractionsHWEvts][numberofmasses], xsection_power16[nofFractionsHWEvts][numberofmasses], xsection_power84[nofFractionsHWEvts][numberofmasses]; //interpolated xsections to give the 1-beta = 50%, 16%, 84% respectively. One array (due to the fractions x) for each mass.
    float xsection_alfadata[nofFractionsHWEvts][numberofmasses];
    float xsection_counting2sigma[numberofmasses];
    float xsection_Median_ClsAlpha[nofFractionsHWEvts][numberofmasses],xsection_16thPerc_ClsAlpha[nofFractionsHWEvts][numberofmasses],xsection_84thPerc_ClsAlpha[nofFractionsHWEvts][numberofmasses],xsection_Data_ClsAlpha[nofFractionsHWEvts][numberofmasses];    

////    TGraph *grin[nofFractionsHWEvts];
    TGraphAsymmErrors *grin[nofFractionsHWEvts];
    TGraph *grin_up[nofFractionsHWEvts];
    TGraph *grin_down[nofFractionsHWEvts];
    TGraph *grin_data[nofFractionsHWEvts];
    
    TGraphAsymmErrors *grin_ExpExclusionCls[nofFractionsHWEvts];
    TGraph *grin_dataCls[nofFractionsHWEvts];
    
    TGraph *gr_plus = new TGraph(10);
    TGraph *gr_ = new TGraph(10);

    TGraph *gr2 = new TGraph(10);
    TGraph *gr2_plus = new TGraph(10);
    TGraph *gr2_ = new TGraph(10);

    ostringstream Lumistrstream;
    Lumistrstream << anaEnv.Luminosity;
    for(unsigned int k=0;k<nofFractionsHWEvts;k++){
    	FractionHWEvtsArray[k] = FractionHWEvts[k];
	ostringstream Fractionstrstream;
     	Fractionstrstream << FractionHWEvts[k];
	hpoweroftest_tprime[k] = new TH2F ("1-#beta for x="+TString(Fractionstrstream.str())+" (Lumi = "+TString(Lumistrstream.str())+"/pb)", "1-#beta vs tprime space (Lumi = "+TString(Lumistrstream.str())+"/pb)", numberofmasses, massmin, massmax, xsectionNbins, xsectionmin, xsectionmax);
        hpoweroftest_tprime[k]->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
	hpoweroftest_tprime[k]->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");//"#sigma t' pair production (pb)"
	
	hpvaluedata_tprime[k] = new TH2F ("p-value of data for x="+TString(Fractionstrstream.str())+" (Lumi = "+TString(Lumistrstream.str())+"/pb)", "p-value of data vs tprime space (Lumi = "+TString(Lumistrstream.str())+"/pb)", numberofmasses, massmin, massmax, xsectionNbins, xsectionmin, xsectionmax);
        hpvaluedata_tprime[k]->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
	hpvaluedata_tprime[k]->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");
    
        hCLs_ExpMedian_tprime[k] = new TH2F ("CLs Median for x="+TString(Fractionstrstream.str())+" (Lumi = "+TString(Lumistrstream.str())+"/pb)", "CLs Median vs tprime space (Lumi = "+TString(Lumistrstream.str())+"/pb)", numberofmasses, massmin, massmax, xsectionNbins, xsectionmin, xsectionmax);
        hCLs_ExpMedian_tprime[k]->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
	hCLs_ExpMedian_tprime[k]->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");

        hCLs_Exp16thPerc_tprime[k] = new TH2F ("CLs 16th Percentile for x="+TString(Fractionstrstream.str())+" (Lumi = "+TString(Lumistrstream.str())+"/pb)", "CLs 16th Percentile vs tprime space (Lumi = "+TString(Lumistrstream.str())+"/pb)", numberofmasses, massmin, massmax, xsectionNbins, xsectionmin, xsectionmax);
        hCLs_Exp16thPerc_tprime[k]->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
	hCLs_Exp16thPerc_tprime[k]->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");

        hCLs_Exp84thPerc_tprime[k] = new TH2F ("CLs 84th Percentile for x="+TString(Fractionstrstream.str())+" (Lumi = "+TString(Lumistrstream.str())+"/pb)", "CLs 84th Percentile vs tprime space (Lumi = "+TString(Lumistrstream.str())+"/pb)", numberofmasses, massmin, massmax, xsectionNbins, xsectionmin, xsectionmax);
        hCLs_Exp84thPerc_tprime[k]->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
	hCLs_Exp84thPerc_tprime[k]->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");
	
	hCLs_Data_tprime[k] = new TH2F ("CLs of data for x="+TString(Fractionstrstream.str())+" (Lumi = "+TString(Lumistrstream.str())+"/pb)", "CLs of data vs tprime space (Lumi = "+TString(Lumistrstream.str())+"/pb)", numberofmasses, massmin, massmax, xsectionNbins, xsectionmin, xsectionmax);
        hCLs_Data_tprime[k]->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
	hCLs_Data_tprime[k]->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");
    }
    hNofSelEvts_tprime = new TH2F ("Number of selected t' pair events (Lumi = "+TString(Lumistrstream.str())+"/pb)", "Number of selected t' pair events (Lumi = "+TString(Lumistrstream.str())+"/pb)", numberofmasses, massmin, massmax, xsectionNbins, xsectionmin, xsectionmax);
    hNofSelEvts_tprime->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
    hNofSelEvts_tprime->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");
    hSignSelEvts_tprime = new TH2F ("Selected t' pair events significance s/sqrt{b} (Lumi = "+TString(Lumistrstream.str())+"/pb)", "Selected t' pair events significance s/sqrt{b} (Lumi = "+TString(Lumistrstream.str())+"/pb)", numberofmasses, massmin, massmax, xsectionNbins, xsectionmin, xsectionmax);
    hSignSelEvts_tprime->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
    hSignSelEvts_tprime->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");

    string rootFileNameNullHypo;
    string rootFileNameAltHypo; //SM pseudos in case of exclusion
    string rootFileNameData;
    string SaveFile;


    cout<<"...Starting scan..."<<endl;
    for (unsigned int miter=0;miter<massestprime.size();miter++){
     
     vector<vector<float> > myPowerOfTestValues, myPValuesData; //for a certain mass, one element of this vector is a vector of the power of tests or p-values for the different fractions x, and you have a number of these vectors since you have a number of cross sections...
     vector<vector<float> > myCLsValues_ExpMedian, myCLsValues_Exp16thPerc, myCLsValues_Exp84thPerc, myCLsValues_Data;
     myPowerOfTestValues.clear();
     myPValuesData.clear();
     myCLsValues_ExpMedian.clear();
     myCLsValues_Exp16thPerc.clear();
     myCLsValues_Exp84thPerc.clear();
     float PowerOfTestArray[nofFractionsHWEvts], PValueDataArray[nofFractionsHWEvts];
     float CLsArray_ExpMedian[nofFractionsHWEvts], CLsArray_Exp16thPerc[nofFractionsHWEvts],CLsArray_Exp84thPerc[nofFractionsHWEvts],CLsArray_Data[nofFractionsHWEvts];
     vector<float> PowerOfTestVector, PValueDataVector;
     vector<float> CLsVector_ExpMedian, CLsVector_Exp16thPerc, CLsVector_Exp84thPerc, CLsVector_Data;     
     vector<float> mySignSelEvtsValues;
     
     for (unsigned int ixsect=0;ixsect<crosssectionstprime.size();ixsect++){ //CHANGE (to use crosssections vector)
	cout<<"  -->  At point ("<<massestprime[miter]<<" GeV, "<<crosssectionstprime[ixsect]<<" pb)"<<endl;
	
	ostringstream masstprime_sstream, crosssections_sstream, lumi_sstream; //nPseudoExp_sstream;
        masstprime_sstream << massestprime[miter];
        crosssections_sstream << crosssectionstprime[ixsect];
        lumi_sstream << anaEnv.Luminosity;
	
	float NofSelEvts_tprime = crosssectionstprime[ixsect] * anaEnv.Luminosity * RefSelEff_tprime[miter];
	hNofSelEvts_tprime->Fill(massestprime[miter],crosssectionstprime[ixsect],NofSelEvts_tprime);
	float SignSelEvts_tprime = NofSelEvts_tprime / sqrt(428*anaEnv.Luminosity/36.1); // the number of background events is an estimation based on what what obtained for 36.1/pb, namely 428 Standard Model MC events
	hSignSelEvts_tprime->Fill(massestprime[miter],crosssectionstprime[ixsect],SignSelEvts_tprime);
	mySignSelEvtsValues.push_back(SignSelEvts_tprime);
	
	PowerOfTestVector.clear();
	PValueDataVector.clear();
	CLsVector_ExpMedian.clear();
	CLsVector_Exp16thPerc.clear();
	CLsVector_Exp84thPerc.clear();
	CLsVector_Data.clear();

        //from now on: filenames should follow convention of output files from GOFanalysis.cc
	rootFileNameNullHypo = thisROUNDdir + "SMTprime_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_" + SignificanceCombination + "_Vtree.root";
	if(SignificanceCombination == "SUM_LRatioMinded") rootFileNameAltHypo = thisROUNDdir + "TpAnalysis_SMonly_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Lumi" + lumi_sstream.str() + "_" + SignificanceCombination + "_Vtree.root";
	else rootFileNameAltHypo = thisROUNDdir + "TpAnalysis_SMonly_Lumi" + lumi_sstream.str() + "_" + SignificanceCombination + "_Vtree.root";	
	
	if(!applyondata) SaveFile = thisROUNDdir + "Results_exclusionTprime_Lumi" + lumi_sstream.str() + "_Nosystematics.root";	
	if(applyondata){
	   SaveFile = thisROUNDdir + "Results_exclusionTprime_Lumi" + lumi_sstream.str() + "_WithData_JESsystematics.root";
	   rootFileNameData = thisROUNDdir + "Data_Lumi" +lumi_sstream.str() +"_" + masstprime_sstream.str() + "GeV_" + crosssections_sstream.str() + "pb_Vtree.root";	
	}


  	WeightProbaCalculator myWeightProbaCalculator;	

       	 //final calculations related to V_alfa and 1-beta (power of the test)
         //vector<float> PowerOfTests; // 1 - beta
        for(int k=0;k<nofFractionsHWEvts;k++){
	  /*if(0.1 != FractionHWEvts[k]){	    
	    //cout<<"Skipping calculation for fraction "<<FractionHWEvts[k]<<endl;
	    continue;	    
	  }*/
	  	  
	  ostringstream Fractionstrstream, alfastrstream;
     	  Fractionstrstream << FractionHWEvts[k];  //if you make a choice here, you can take e.g. k = 3
  	  vector<float> NullHypo_Vvect;
  	  vector<float> AltHypo_Vvect;
	  NullHypo_Vvect.clear();
	  AltHypo_Vvect.clear();
	  float CLs, CLb, CLbpluss; 
	  
	  // TFile in which the V values of the 'data' pseudoexps are stored; i.e. in the case of setting exclusion limits on mSugra models, this is TTJets (= "SM only")
  	  TFile *_treeAltHypoFile = TFile::Open(rootFileNameAltHypo.c_str());
  	/*  if (_treeSMPsExpsFile ->IsZombie()) {cout << " FILE "<<_treeSMPsExpsFile << " is a Zombie()....continuing with next point"<<endl; break;break;}
	*/
	  TTree *T = (TTree*)_treeAltHypoFile->Get("Vtree_fr"+TString(Fractionstrstream.str()));
  	
	  float V=0;
  	  T->SetBranchAddress("V",&V);
  	  for(int i=0;i<T->GetEntriesFast();i++){
 		T->GetEntry(i);
 		AltHypo_Vvect.push_back(V);
  	  }
  	  _treeAltHypoFile->Close();
	
	  //cout<<"fraction x = "<<Fractionstrstream.str()<<"  and m0 "<<ll<<"  m12  "<<kk<<"  "<<v_mean/dataVvect.size()<<endl;
	  //reading Valfa, and doing calculations...
	  float Valfa = 0, poweroftest = 0;
	  alfastrstream << alfa;
	  TFile *_treeNullHypoFile = TFile::Open(rootFileNameNullHypo.c_str());

  	  float V_np=0;
  	 /* if (_treeNullHypoFile ->IsZombie()) {cout << " FILE "<<_treeNullHypoFile << " is a Zombie()....continuing with next point"<<endl; break;break;}*/
	  TTree *ValfatreeNullHypo = (TTree*)_treeNullHypoFile->Get("Valfatree_fr"+TString(Fractionstrstream.str()));
	  TTree *T_np = (TTree*)_treeNullHypoFile->Get("Vtree_fr"+TString(Fractionstrstream.str()));
  	  T_np->SetBranchAddress("V",&V_np);
	  for(int i=0;i<T_np->GetEntriesFast();i++){
 		T_np->GetEntry(i);
 		NullHypo_Vvect.push_back(V_np);
	  }

	  sort(NullHypo_Vvect.begin(),NullHypo_Vvect.end()); //moet dit?
	  Valfa = myWeightProbaCalculator.GetProbaPercentile(NullHypo_Vvect, alfa*100); //check if this is correct!!
	  
	  //ValfatreeNullHypo->SetBranchAddress("Valfa"+TString(alfastrstream.str()),&Valfa);
	//for(int i=0;i<ValfatreeNullHypo->GetEntriesFast();i++){
	  //		ValfatreeNullHypo->GetEntry(i);
	//	cout<<"    Valfa = "<<Valfa<<endl;}



        //using power of test
	  poweroftest = myWeightProbaCalculator.PowerOfTest(AltHypo_Vvect,Valfa);
	  PowerOfTestArray[k] = poweroftest;
	  hpoweroftest_tprime[k]->Fill(massestprime[miter],crosssectionstprime[ixsect],poweroftest);

/////	  cout<<  " Filled the exclusion band plot (x = "<<FractionHWEvts[k]<<")  ("<<massestprime[miter]<<" GeV,"<<crosssectionstprime[ixsect]<<" pb), poweroftest = "<<PowerOfTestArray[k]<<endl;

          if(applyondata){
	   //workflow not so much different from determination of 1-beta (see above)
	    TFile *_treeDataFile = TFile::Open(TString(rootFileNameData));
            float V_data=0;
	    TTree *ValfatreeData = (TTree*)_treeDataFile->Get("Valfatree_fr"+TString(Fractionstrstream.str()));
	    TTree *T_data = (TTree*)_treeDataFile->Get("Vtree_fr"+TString(Fractionstrstream.str()));
  	    T_data->SetBranchAddress("V",&V_data);
	    if(T_data->GetEntriesFast() != 1) cout<<"WARNING: DATA FILE MUST HAVE ONLY 1 V-VALUE!!"<<endl;
	    for(int i=0;i<1;i++){
 		T_data->GetEntry(i);		
	    }
	    float pvaluedata = 0;
	    myWeightProbaCalculator.CalculateProba(NullHypo_Vvect, V_data);
	    pvaluedata = 1 - myWeightProbaCalculator.GetProba();//should fix this in the future: extra argument 'left'/'right'.
	    PValueDataArray[k] = pvaluedata;
	    hpvaluedata_tprime[k]->Fill(massestprime[miter],crosssectionstprime[ixsect],pvaluedata);
          
	    //regarding CLs method on data
	    //myWeightProbaCalculatorNew.CalculateProba(NullHypo_Vvect, V_SMonlyMedian);
	    //CLbpluss = 1 - myWeightProbaCalculator.GetProba();
	    CLbpluss = pvaluedata;
	    myWeightProbaCalculator.CalculateProba(AltHypo_Vvect, V_data);
	    CLb = 1 - myWeightProbaCalculator.GetProba();
	    if(likelihood){
	      CLbpluss = 1 - CLbpluss;
	      CLb = 1 - CLb;
	    }
	    CLs = CLbpluss / CLb;
	    cout<<"       *****************>>>>>>>>>>>>>>>>>>>      CLS data = "<<CLs<<endl;
	    CLsArray_Data[k] = CLs;	  
	    hCLs_Data_tprime[k]->Fill(massestprime[miter],crosssectionstprime[ixsect],CLs);
	  }


        //now attempt to use CLs method; to be integrated later in WeightProbaCalculator class (?)	  
	  
	  //median
	  float V_SMonlyMedian = -999;
	  V_SMonlyMedian = myWeightProbaCalculator.GetProbaPercentile(AltHypo_Vvect,50);
	  myWeightProbaCalculator.CalculateProba(NullHypo_Vvect, V_SMonlyMedian);
	  CLbpluss = 1 - myWeightProbaCalculator.GetProba();
	  myWeightProbaCalculator.CalculateProba(AltHypo_Vvect, V_SMonlyMedian);
	  CLb = 1 - myWeightProbaCalculator.GetProba();
	  if(likelihood){
	    CLbpluss = 1 - CLbpluss;
	    CLb = 1 - CLb;
	  }
	  CLs = CLbpluss / CLb;
	  CLsArray_ExpMedian[k] = CLs;
	  hCLs_ExpMedian_tprime[k]->Fill(massestprime[miter],crosssectionstprime[ixsect],CLs);
	  cout<<  " Filled the CLs plot for the median of the SM-only Vdistribution (x = "<<FractionHWEvts[k]<<")  ("<<massestprime[miter]<<" GeV,"<<crosssectionstprime[ixsect]<<" pb), CLbpluss = "<<CLbpluss<<", CLb = "<<CLb<<" ----> CLs = "<<CLs<<endl;

	  //16th percentile
	  float V_SMonly16thPerc = -999;
	  V_SMonly16thPerc = myWeightProbaCalculator.GetProbaPercentile(AltHypo_Vvect,16);
	  myWeightProbaCalculator.CalculateProba(NullHypo_Vvect, V_SMonly16thPerc);
	  CLbpluss = 1 - myWeightProbaCalculator.GetProba();
	  myWeightProbaCalculator.CalculateProba(AltHypo_Vvect, V_SMonly16thPerc);
	  CLb = 1 - myWeightProbaCalculator.GetProba();
	  if(likelihood){
	    CLbpluss = 1 - CLbpluss;
	    CLb = 1 - CLb;
	  }
	  CLs = CLbpluss / CLb;
	  CLsArray_Exp16thPerc[k] = CLs;	  
	  hCLs_Exp16thPerc_tprime[k]->Fill(massestprime[miter],crosssectionstprime[ixsect],CLs);
	  cout<<  " Filled the CLs plot for the 16th percentile of the SM-only Vdistribution (x = "<<FractionHWEvts[k]<<")  ("<<massestprime[miter]<<" GeV,"<<crosssectionstprime[ixsect]<<" pb), CLbpluss = "<<CLbpluss<<", CLb = "<<CLb<<" ----> CLs = "<<CLs<<endl;

	  //84th percentile
	  float V_SMonly84thPerc = -999;
	  V_SMonly84thPerc = myWeightProbaCalculator.GetProbaPercentile(AltHypo_Vvect,84);
	  myWeightProbaCalculator.CalculateProba(NullHypo_Vvect, V_SMonly84thPerc);
	  CLbpluss = 1 - myWeightProbaCalculator.GetProba();
	  myWeightProbaCalculator.CalculateProba(AltHypo_Vvect, V_SMonly84thPerc);
	  CLb = 1 - myWeightProbaCalculator.GetProba();
	  if(likelihood){
	    CLbpluss = 1 - CLbpluss;
	    CLb = 1 - CLb;
	  }
	  CLs = CLbpluss / CLb;
	  CLsArray_Exp84thPerc[k] = CLs;	  
	  hCLs_Exp84thPerc_tprime[k]->Fill(massestprime[miter],crosssectionstprime[ixsect],CLs);          

	 
  	  _treeNullHypoFile->Close();
        
	}////////loop for all fraction
	
	
	
	//cout<<"PowerOfTestArray[0] = "<<PowerOfTestArray[0]<<endl;	
	for(unsigned int k=0;k<nofFractionsHWEvts;k++){
	   PowerOfTestVector.push_back(PowerOfTestArray[k]);
	   CLsVector_ExpMedian.push_back(CLsArray_ExpMedian[k]);
	   CLsVector_Exp16thPerc.push_back(CLsArray_Exp16thPerc[k]);
	   CLsVector_Exp84thPerc.push_back(CLsArray_Exp84thPerc[k]);
	}
	myPowerOfTestValues.push_back(PowerOfTestVector);
	myCLsValues_ExpMedian.push_back(CLsVector_ExpMedian);
	myCLsValues_Exp16thPerc.push_back(CLsVector_Exp16thPerc);
	myCLsValues_Exp84thPerc.push_back(CLsVector_Exp84thPerc);
	if(applyondata){
	  for(unsigned int k=0;k<nofFractionsHWEvts;k++){
	     PValueDataVector.push_back(PValueDataArray[k]);
	     CLsVector_Data.push_back(CLsArray_Data[k]);
	  }
	  myPValuesData.push_back(PValueDataVector);
	  myCLsValues_Data.push_back(CLsVector_Data);
	}
	
	/*for(unsigned int i=0;i<PowerOfTestVector.size();i++){
	   cout<<"      PowerOfTestVector["<<i<<"] = "<<PowerOfTestVector[i]<<endl;
	}
	for(int k=0;k<nofFractionsHWEvts;k++){
	   cout<<" myPowerOfTestValues["<<ixsect<<"]["<<k<<"] = "<<myPowerOfTestValues[ixsect][k]<<endl;
        }*/
	
        /*TCanvas* c1 = new TCanvas("c1","PowerOfTest_vs_Fraction "+TString(rootFileNameNullHypo),500,500);
        TGraph* graph_PowerOfTestvsFraction = new TGraph(nofFractionsHWEvts,FractionHWEvtsArray,PowerOfTestArray);
        graph_PowerOfTestvsFraction->SetMarkerColor(kBlue);
        graph_PowerOfTestvsFraction->SetMarkerStyle(21);
        graph_PowerOfTestvsFraction->GetXaxis()->SetTitle("x");
        graph_PowerOfTestvsFraction->GetYaxis()->SetTitle("1-#beta");
        graph_PowerOfTestvsFraction->Draw("AP");
        c1->Write("PowerOfTest_vs_Fraction",TObject::kOverwrite); //temporary: always overwritten...
        delete c1;
        cout<<endl<<"Output written to "<<_output->GetName()<<endl;*/
	
 
      
     }/////////////////// end of loop for tprime xsections
     
     cout<<"Interpolating for 1-beta = 0.5..."<<endl;
     for(int k=0;k<nofFractionsHWEvts;k++){
       bool InterpolationSucces = false;    
       for (unsigned int ixsect=1;ixsect<crosssectionstprime.size();ixsect++){
           //cout<<" crosssectionstprime["<<ixsect<<"] = "<<crosssectionstprime[ixsect]<<endl;
	   //cout<<"    myPowerOfTestValues["<<ixsect<<"]["<<k<<"] = "<<myPowerOfTestValues[ixsect][k]<<endl;
	   //cout<<"    myPowerOfTestValues["<<ixsect-1<<"]["<<k<<"] = "<<myPowerOfTestValues[ixsect-1][k]<<endl;
	   if (crosssectionstprime[ixsect]>0 && myPowerOfTestValues[ixsect][k]>0.5 && myPowerOfTestValues[ixsect-1][k]<0.5){ 	      
	      xsection_power50[k][miter] = crosssectionstprime[ixsect-1] + ((0.5 - myPowerOfTestValues[ixsect-1][k]) / (myPowerOfTestValues[ixsect][k] - myPowerOfTestValues[ixsect-1][k])) * (crosssectionstprime[ixsect] - crosssectionstprime[ixsect-1]);	
	      //grin->SetPoint(m0point,ll,point_new);
//	      cout<<"*********************************** miter = "<<miter<<", mtprime = "<<massestprime[miter]<<", xsection_power50["<<k<<"]["<<miter<<"] = "<<xsection_power50[k][miter]<<endl;
////	      xaxisarray[miter] = massestprime[miter];
	      InterpolationSucces = true;
	      break;
           }
       }////end of looping interpolation
       if(!InterpolationSucces) cout<<" !! WARNING: interpolation failed, probably extrapolation needed !!"<<endl;
     }

     cout<<"Interpolating for 1-beta = 0.16..."<<endl;
     for(int k=0;k<nofFractionsHWEvts;k++){
       bool InterpolationSucces = false;   
       for (unsigned int ixsect=1;ixsect<crosssectionstprime.size();ixsect++){
//	   cout<<"   mtprime = "<<massestprime[miter]<<", myPowerOfTestValues["<<ixsect<<"]["<<k<<"] = "<<myPowerOfTestValues[ixsect][k]<<",  myPowerOfTestValues["<<ixsect-1<<"]["<<k<<"] = "<<myPowerOfTestValues[ixsect-1][k]<<endl;
	   if (crosssectionstprime[ixsect]>0 && myPowerOfTestValues[ixsect][k]>0.158 && myPowerOfTestValues[ixsect-1][k]<0.158){ 	      
	      xsection_power16[k][miter] = crosssectionstprime[ixsect-1] + ((0.158 - myPowerOfTestValues[ixsect-1][k]) / (myPowerOfTestValues[ixsect][k] - myPowerOfTestValues[ixsect-1][k])) * (crosssectionstprime[ixsect] - crosssectionstprime[ixsect-1]);	
//	      cout<<"      xsection_power16["<<k<<"]["<<miter<<"] = "<<xsection_power16[k][miter]<<endl;	      	      
	      InterpolationSucces = true;
	      break;
           }
       }////end of looping interpolation
       if(!InterpolationSucces) cout<<" !! WARNING: interpolation failed, probably extrapolation needed !!"<<endl;
     }

     cout<<"Interpolating for 1-beta = 0.84..."<<endl;
     for(int k=0;k<nofFractionsHWEvts;k++){   
       bool InterpolationSucces = false; 
       for (unsigned int ixsect=1;ixsect<crosssectionstprime.size();ixsect++){
	   if (crosssectionstprime[ixsect]>0 && myPowerOfTestValues[ixsect][k]>0.842 && myPowerOfTestValues[ixsect-1][k]<0.842){ 	      
	      xsection_power84[k][miter] = crosssectionstprime[ixsect-1] + ((0.842 - myPowerOfTestValues[ixsect-1][k]) / (myPowerOfTestValues[ixsect][k] - myPowerOfTestValues[ixsect-1][k])) * (crosssectionstprime[ixsect] - crosssectionstprime[ixsect-1]);	
	      InterpolationSucces = true;
	      break;
           }
       }////end of looping interpolation
       if(!InterpolationSucces) cout<<" !! WARNING: interpolation failed, probably extrapolation needed !!"<<endl;
     }
     
     if(applyondata){ //be careful! the p-value go from high to low when cross section rises, different compared to 1-beta...
        cout<<"Interpolating for p-value = "<<alfa<<"..."<<endl;
        bool InterpolationSucces = false; 
	for(int k=0;k<nofFractionsHWEvts;k++){   
          for (unsigned int ixsect=1;ixsect<crosssectionstprime.size();ixsect++){
	      //cout<<"myPValuesData["<<ixsect<<"]["<<k<<"] = "<<myPValuesData[ixsect][k]<<endl;
	      if (crosssectionstprime[ixsect]>0 && myPValuesData[ixsect][k]<alfa && myPValuesData[ixsect-1][k]>alfa){ //so here the < and > are switched w.r.t. 1-beta interpolation case... 	      
	         xsection_alfadata[k][miter] = crosssectionstprime[ixsect-1] + ((alfa - myPValuesData[ixsect-1][k]) / (myPValuesData[ixsect][k] - myPValuesData[ixsect-1][k])) * (crosssectionstprime[ixsect] - crosssectionstprime[ixsect-1]);	
	         //cout<<"xsection_alfadata["<<k<<"]["<<miter<<"] = "<<xsection_alfadata[k][miter]<<endl;
		 InterpolationSucces = true;
		 break;
              }
          }////end of looping interpolation
	  if(!InterpolationSucces) cout<<" !! WARNING: interpolation failed, probably extrapolation needed !!"<<endl;
        }	
	cout<<"Interpolating for (data) CLs = "<<alfa<<"..."<<endl;
        InterpolationSucces = false; 
	for(int k=0;k<nofFractionsHWEvts;k++){   
          for (unsigned int ixsect=1;ixsect<crosssectionstprime.size();ixsect++){
	      //cout<<"myPValuesData["<<ixsect<<"]["<<k<<"] = "<<myPValuesData[ixsect][k]<<endl;
	      if (crosssectionstprime[ixsect]>0 && myCLsValues_Data[ixsect][k]<alfa && myCLsValues_Data[ixsect-1][k]>alfa){ //so here the < and > are switched w.r.t. 1-beta interpolation case... 	      
	         xsection_Data_ClsAlpha[k][miter] = crosssectionstprime[ixsect-1] + ((alfa - myCLsValues_Data[ixsect-1][k]) / (myCLsValues_Data[ixsect][k] - myCLsValues_Data[ixsect-1][k])) * (crosssectionstprime[ixsect] - crosssectionstprime[ixsect-1]);	
	         cout<<"      xsection_Data_ClsAlpha["<<k<<"]["<<miter<<"] = "<<xsection_Data_ClsAlpha[k][miter]<<endl;
		 InterpolationSucces = true;
		 break;
              }
          }////end of looping interpolation
	  if(!InterpolationSucces) cout<<" !! WARNING: interpolation failed, probably extrapolation needed !!"<<endl;
        }
     }
     
     //tbc!!
     cout<<"Interpolating for sigma = 2 (counting experiment)..."<<endl;
     //bool InterpolationSucces = false;    
     for (unsigned int ixsect=1;ixsect<crosssectionstprime.size();ixsect++){	   
	   if (crosssectionstprime[ixsect]>0 && mySignSelEvtsValues[ixsect]>2 && mySignSelEvtsValues[ixsect-1]<2){ 	      
	      xsection_counting2sigma[miter] = crosssectionstprime[ixsect-1] + ((2 - mySignSelEvtsValues[ixsect-1]) / (mySignSelEvtsValues[ixsect] - mySignSelEvtsValues[ixsect-1])) * (crosssectionstprime[ixsect] - crosssectionstprime[ixsect-1]);	
	      //InterpolationSucces = true;
	      break;
           }
     }////end of looping interpolation
     //if(!InterpolationSucces) cout<<" !! WARNING: interpolation failed, probably extrapolation needed !!"<<endl;    
        
	
     //for CLs method	  
     cout<<"Interpolating (median) for CLs = "<<alfa<<"..."<<endl;
     for(int k=0;k<nofFractionsHWEvts;k++){
       bool InterpolationSucces = false;   
       for (unsigned int ixsect=1;ixsect<crosssectionstprime.size();ixsect++){
//	   cout<<"   mtprime = "<<massestprime[miter]<<", myCLsValues_ExpMedian["<<ixsect<<"]["<<k<<"] = "<<myCLsValues_ExpMedian[ixsect][k]<<",  myCLsValues_ExpMedian["<<ixsect-1<<"]["<<k<<"] = "<<myCLsValues_ExpMedian[ixsect-1][k]<<endl;
	   if (crosssectionstprime[ixsect]>0 && myCLsValues_ExpMedian[ixsect][k]<alfa && myCLsValues_ExpMedian[ixsect-1][k]>alfa){ 	      
	      xsection_Median_ClsAlpha[k][miter] = crosssectionstprime[ixsect-1] + ((alfa - myCLsValues_ExpMedian[ixsect-1][k]) / (myCLsValues_ExpMedian[ixsect][k] - myCLsValues_ExpMedian[ixsect-1][k])) * (crosssectionstprime[ixsect] - crosssectionstprime[ixsect-1]);	
	      cout<<"      xsection_Median_ClsAlpha["<<k<<"]["<<miter<<"] = "<<xsection_Median_ClsAlpha[k][miter]<<endl;	      	      
/////	      xaxisarray[miter] = massestprime[miter];
	      InterpolationSucces = true;
	      break;
           }
       }////end of looping interpolation
       if(!InterpolationSucces) cout<<" !! WARNING: interpolation failed, probably extrapolation needed !!"<<endl;
     } 
     	  
     cout<<"Interpolating (16th Percentile) for CLs = "<<alfa<<"..."<<endl;
     for(int k=0;k<nofFractionsHWEvts;k++){
       bool InterpolationSucces = false;   
       for (unsigned int ixsect=1;ixsect<crosssectionstprime.size();ixsect++){
	   if (crosssectionstprime[ixsect]>0 && myCLsValues_Exp16thPerc[ixsect][k]<alfa && myCLsValues_Exp16thPerc[ixsect-1][k]>alfa){ 	      
	      xsection_16thPerc_ClsAlpha[k][miter] = crosssectionstprime[ixsect-1] + ((alfa - myCLsValues_Exp16thPerc[ixsect-1][k]) / (myCLsValues_Exp16thPerc[ixsect][k] - myCLsValues_Exp16thPerc[ixsect-1][k])) * (crosssectionstprime[ixsect] - crosssectionstprime[ixsect-1]);	
	      cout<<"      xsection_16thPerc_ClsAlpha["<<k<<"]["<<miter<<"] = "<<xsection_16thPerc_ClsAlpha[k][miter]<<endl;	      	      
	      InterpolationSucces = true;
	      break;
           }
       }////end of looping interpolation
       if(!InterpolationSucces) cout<<" !! WARNING: interpolation failed, probably extrapolation needed !!"<<endl;
     } 

     cout<<"Interpolating (84th Percentile) for CLs = "<<alfa<<"..."<<endl;
     for(int k=0;k<nofFractionsHWEvts;k++){
       bool InterpolationSucces = false;   
       for (unsigned int ixsect=1;ixsect<crosssectionstprime.size();ixsect++){
	   if (crosssectionstprime[ixsect]>0 && myCLsValues_Exp84thPerc[ixsect][k]<alfa && myCLsValues_Exp84thPerc[ixsect-1][k]>alfa){ 	      
	      xsection_84thPerc_ClsAlpha[k][miter] = crosssectionstprime[ixsect-1] + ((alfa - myCLsValues_Exp84thPerc[ixsect-1][k]) / (myCLsValues_Exp84thPerc[ixsect][k] - myCLsValues_Exp84thPerc[ixsect-1][k])) * (crosssectionstprime[ixsect] - crosssectionstprime[ixsect-1]);	
	      cout<<"      xsection_84thPerc_ClsAlpha["<<k<<"]["<<miter<<"] = "<<xsection_84thPerc_ClsAlpha[k][miter]<<endl;	      	      
	      InterpolationSucces = true;
	      break;
           }
       }////end of looping interpolation
       if(!InterpolationSucces) cout<<" !! WARNING: interpolation failed, probably extrapolation needed !!"<<endl;
     } 	  
	  
	  
    }/////////////////// end of loop for tprime masses

    _output = TFile::Open(TString(SaveFile),"RECREATE");
    _output->cd();
    //hNofSelEvts_tprime->SetDrawOption("COLZ");
    hNofSelEvts_tprime->Write();
    hSignSelEvts_tprime->Write();
    for(int k=0;k<nofFractionsHWEvts;k++){    
	   //hpoweroftest_tprime[k]->SetDrawOption("COLZ");
	   hpoweroftest_tprime[k]->Write();
	   //hpvaluedata_tprime[k]->SetDrawOption("COLZ");
	   hpvaluedata_tprime[k]->Write();
	   	   
	   ostringstream Fractionstrstream, lumi_sstream, clstream;
     	   Fractionstrstream << FractionHWEvts[k];
	   lumi_sstream << int(anaEnv.Luminosity + 1); //maybe change this!!!
	   TCanvas *ccc = new TCanvas("cinterpol_x"+TString(Fractionstrstream.str()),"Interpolated Exclusion limit",0,0,700,700);

   /*   TH1F *frame = new TH1F("frame","",numberofmasses,270,470);
      frame->SetMinimum(0); //(1) //0 also for lumi = 36.1 //xsectrange[0]?
      frame->SetMaximum(xsectrange[1]); //cfr 11 for lumi = 250/pb, 26 for lumi = 36.1/pb
      frame->SetDirectory(0);
      frame->SetStats(0);
      frame->GetXaxis()->SetTitle("mass t' (GeV/c^{2})");
      frame->GetYaxis()->SetTitle("#sigma t' pair production (pb)");
      frame->GetXaxis()->SetTickLength(0.02);
      frame->GetXaxis()->SetLabelSize(0.03);
      frame->GetYaxis()->SetTitle("#sigma t' pair production (pb)");
      frame->GetYaxis()->SetMoreLogLabels();
      frame->GetYaxis()->SetLabelSize(0.03);*/
///           gPad->SetLogy();//doesn't work properly??
/////	   ccc->SetLogy();//doesn't work properly?? screws up the axes, axis labels, lines disappear etc...	   
	   ccc->SetTickx();
           ccc->SetTicky();
           ccc->SetGridx();
           ccc->SetGridy();

           //ccc->cd();
/////	   grin[k] = new TGraph(numberofmasses,xaxisarray,xsection_power50[k]);
	   float xassymerror_left[numberofmasses];
	   float xassymerror_right[numberofmasses];
	   float yassymerror_down[numberofmasses];
	   float yassymerror_up[numberofmasses];
	   for(unsigned int u=0;u<numberofmasses;u++){
	       xassymerror_left[u] = 0;
	       xassymerror_right[u] = 0;
	       yassymerror_down[u] = xsection_power50[k][u] - xsection_power16[k][u];
	       yassymerror_up[u] = xsection_power84[k][u] - xsection_power50[k][u];
	   }
           grin[k] = new TGraphAsymmErrors(numberofmasses,xaxisarray,xsection_power50[k],xassymerror_left,xassymerror_right,yassymerror_down,yassymerror_up); //see http://root.cern.ch/root/html/TGraphPainter.html
	   //grin[k]->GetXaxis()->SetMoreLogLabels();
	   grin[k]->SetName("thisTGraph");
	   string title = "CMS Preliminary,     #sqrt{s} = 7 TeV,     L_{int} = " + lumi_sstream.str() + " pb^{-1}"; 	   
	   grin[k]->SetTitle(title.c_str());
	   grin[k]->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
	   //grin[k]->GetYaxis()->SetMoreLogLabels();
	   grin[k]->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");
	   grin[k]->GetYaxis()->SetRangeUser(xsectrange[0],50); //cfr (0,11) for lumi = 250/pb, (0,26) for lumi = 36.1/pb
  //    frame->Draw(" ");
	   grin[k]->SetFillColor(4);
       /////    grin[k]->SetFillStyle(5005);
           grin[k]->SetLineColor(4);
           grin[k]->SetLineWidth(3);
	   grin[k]->SetMarkerColor(4);
           grin[k]->SetMarkerSize(1.);
           grin[k]->SetMarkerStyle(kFullSquare); //kFullSquare = 21  
	   grin[k]->SetFillStyle(3005);
	/////   grin[k]->Draw("ACP4");
	   grin[k]->Draw("ACP3");
	   
	/*   grin_down[k] = new TGraph(numberofmasses,xaxisarray,xsection_power16[k]);
	   grin_down[k]->GetXaxis()->SetTitle("mass t' (GeV/c^{2})");
	   grin_down[k]->GetYaxis()->SetTitle("#sigma t' pair production (pb)");
	   grin_down[k]->GetYaxis()->SetRangeUser(0.5,50); //cfr (0,11) for lumi = 250/pb, (0,26) for lumi = 36.1/pb
	   grin_down[k]->SetFillColor(4);
           grin_down[k]->SetFillStyle(5005);
           grin_down[k]->SetLineColor(4);
           grin_down[k]->SetLineWidth(3);
	   grin_down[k]->SetLineStyle(kDashed);
	   grin_down[k]->Draw("C");
	   
	   grin_up[k] = new TGraph(numberofmasses,xaxisarray,xsection_power84[k]);
	   grin_up[k]->GetXaxis()->SetTitle("mass t' (GeV/c^{2})");
	   grin_up[k]->GetYaxis()->SetTitle("#sigma t' pair production (pb)");
	   grin_up[k]->GetYaxis()->SetRangeUser(0.5,50); //cfr (0,11) for lumi = 250/pb, (0,26) for lumi = 36.1/pb
	   grin_up[k]->SetFillColor(4);
           grin_up[k]->SetFillStyle(5005);
           grin_up[k]->SetLineColor(4);
           grin_up[k]->SetLineWidth(3);
	   grin_up[k]->SetLineStyle(kDashed);
	   grin_up[k]->Draw("C");*/
	   
	   if(applyondata){
	      grin_data[k] = new TGraph(numberofmasses,xaxisarray,xsection_alfadata[k]);
	      grin_data[k]->SetName("thisTGraph_data");
	      grin_data[k]->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
	      grin_data[k]->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");
	      grin_data[k]->GetYaxis()->SetRangeUser(xsectrange[0],50); //cfr (0,11) for lumi = 250/pb, (0,26) for lumi = 36.1/pb
	      grin_data[k]->SetFillColor(3);
              grin_data[k]->SetFillStyle(5005);
              grin_data[k]->SetLineColor(3);
              grin_data[k]->SetLineWidth(3);
	      grin_data[k]->SetMarkerColor(3);
              grin_data[k]->SetMarkerSize(1.);
              grin_data[k]->SetMarkerStyle(kFullSquare);
	      grin_data[k]->Draw("CP"); //if you do this ACP, and grin[k] was also ACP, the grin[k] is not drawn for some reason... Without the 'A' here, it works
	   }
	   
	   TGraph* grin_theo = new TGraph(numberofmasses,xaxisarray,xsection_theory);
	   cout<<"numberofmasses"<<numberofmasses<<endl;
	   for(int ii=0;ii<numberofmasses;ii++){
	     cout<<"      xaxisarray["<<ii<<"] = "<<xaxisarray[ii]<<endl;
	     cout<<"      xsection_theory["<<ii<<"] = "<<xsection_theory[ii]<<endl;
	   }
	   grin_theo->SetName("thisTGraph_theory");
	   grin_theo->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
	   grin_theo->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");
	   grin_theo->GetYaxis()->SetRangeUser(xsectrange[0],50); //cfr (0,11) for lumi = 250/pb, (0,26) for lumi = 36.1/pb
           grin_theo->SetFillStyle(5005);
           grin_theo->SetLineColor(2);
           grin_theo->SetLineWidth(3);
	   grin_theo->SetMarkerColor(2);
           grin_theo->SetMarkerSize(1.);
           grin_theo->SetMarkerStyle(kFullTriangleUp); //kFullTriangleUp = 22
	   grin_theo->Draw("CP"); //"C" worked, then you get only the curve
	   
	   /*TGraph* grin_counting = new TGraph(numberofmasses,xaxisarray,xsection_counting2sigma);
	   grin_counting->GetXaxis()->SetTitle("mass t' (GeV/c^{2})");
	   grin_counting->GetYaxis()->SetTitle("#sigma t' pair production (pb)");
	   grin_counting->GetYaxis()->SetRangeUser(0,xsectrange[1]); //cfr (0,11) for lumi = 250/pb, (0,26) for lumi = 36.1/pb
	   grin_counting->SetFillColor(kMagenta);
           grin_counting->SetFillStyle(5005);
           grin_counting->SetLineColor(kMagenta);
           grin_counting->SetLineWidth(3);
	   grin_counting->SetLineStyle(kDashDotted);
	   grin_counting->Draw("C");*/
	   	   
	   clstream << cl;
	   //string texttex1 = "Expected limit at " + clstream.str() + "% C.L. for L_{int} = " + lumi_sstream.str() + " pb^{-1}";	   
	   
	   TLegend* leg = new TLegend(0.5,0.7,0.9,0.9); //TLegend(0.1,0.7,0.48,0.9)
	   string textExpectedLine = "Expected limit at " + clstream.str() + "% C.L.";
	   leg->AddEntry("thisTGraph",textExpectedLine.c_str(),"lpf"); //"lep"?
	   if(applyondata){
	      string textObservedLine = "Observed limit at " + clstream.str() + "% C.L.";
	      leg->AddEntry("thisTGraph_data",textObservedLine.c_str(),"lp");
	   }
	   leg->AddEntry("thisTGraph_theory","Theoretical prediction (NLO)","lp");
	   leg->Draw();
	   
	   /*string texttex1 = "Expected limit at " + clstream.str() + "% C.L.";
	   TLatex *tex1 = new TLatex (massestprimeHardCoded[0],xsectrange[1]-1,texttex1.c_str());
	   tex1->SetTextColor(4);
	   //tex1->SetTextSize(1);
           tex1->Draw();
	   if(applyondata){
	      string texttex2 = "Observed limit at " + clstream.str() + "% C.L.";
	      TLatex *tex2 = new TLatex (massestprimeHardCoded[0],xsectrange[1]-2,texttex2.c_str());
	      tex2->SetTextColor(3);
	      tex2->Draw();
	   }	   
	   TLatex *tex3 = new TLatex (massestprimeHardCoded[0],xsectrange[1]-3,"Theoretical prediction (NLO)");
	   tex3->SetTextColor(2);
	   tex3->Draw();*/
	   /*TLatex *tex4 = new TLatex (massestprimeHardCoded[0],xsectrange[1]-4,"2 sigma line for counting experiment");
	   tex4->SetTextColor(kMagenta);
	   tex4->Draw();*/
	   	   
	   ccc->Write("InterpolatedExclusion_x"+TString(Fractionstrstream.str()),TObject::kOverwrite);
	   
	   
	   //for CLs method
	   hCLs_ExpMedian_tprime[k]->Write();
	   hCLs_Exp16thPerc_tprime[k]->Write();
	   hCLs_Exp84thPerc_tprime[k]->Write();
	   hCLs_Data_tprime[k]->Write();
	   TCanvas *cCls = new TCanvas("cCls_x"+TString(Fractionstrstream.str()),"Interpolated Exclusion limit with CLs method",0,0,700,700);
           cCls->cd();
	   cCls->SetTickx();
           cCls->SetTicky();
           cCls->SetGridx();
           cCls->SetGridy();
	   //overwriting...
	   for(unsigned int u=0;u<numberofmasses;u++){
	       xassymerror_left[u] = 0;
	       xassymerror_right[u] = 0;
	       yassymerror_down[u] = xsection_Median_ClsAlpha[k][u] - xsection_16thPerc_ClsAlpha[k][u];
	       yassymerror_up[u] = xsection_84thPerc_ClsAlpha[k][u] - xsection_Median_ClsAlpha[k][u];
	   }
	   //yassymerror_up[0] = 3.8 - xsection_Median_ClsAlpha[k][0]; //remember, this is cheating!! only to make plots of ROUND24,25,26 prettier
	   grin_ExpExclusionCls[k] = new TGraphAsymmErrors(numberofmasses,xaxisarray,xsection_Median_ClsAlpha[k],xassymerror_left,xassymerror_right,yassymerror_down,yassymerror_up); //see http://root.cern.ch/root/html/TGraphPainter.html
	   grin_ExpExclusionCls[k]->SetName("thisTGraph_CLs");
	   title = "CMS Preliminary,     #sqrt{s} = 7 TeV,     L_{int} = " + lumi_sstream.str() + " pb^{-1}"; 	   
	   grin_ExpExclusionCls[k]->SetTitle(title.c_str());
	   grin_ExpExclusionCls[k]->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
	   grin_ExpExclusionCls[k]->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");
	   grin_ExpExclusionCls[k]->GetYaxis()->SetRangeUser(xsectrange[0],50); //cfr (0,11) for lumi = 250/pb, (0,26) for lumi = 36.1/pb
	   grin_ExpExclusionCls[k]->SetFillColor(4);
           grin_ExpExclusionCls[k]->SetLineColor(4);
           grin_ExpExclusionCls[k]->SetLineWidth(3);
	   grin_ExpExclusionCls[k]->SetMarkerColor(4);
           grin_ExpExclusionCls[k]->SetMarkerSize(1.);
           grin_ExpExclusionCls[k]->SetMarkerStyle(kFullSquare); //kFullSquare = 21  
	   grin_ExpExclusionCls[k]->SetFillStyle(3005);
	/////   grin[k]->Draw("ACP4");
	   grin_ExpExclusionCls[k]->Draw("ACP3");
	   
	   if(applyondata){
	      grin_dataCls[k] = new TGraph(numberofmasses,xaxisarray,xsection_Data_ClsAlpha[k]);
	      grin_dataCls[k]->SetName("thisTGraph_dataCLs");
	      grin_dataCls[k]->GetXaxis()->SetTitle("t' mass (GeV/c^{2})");
	      grin_dataCls[k]->GetYaxis()->SetTitle("#sigma(pp #rightarrow t'#bar{t}') (pb)");
	      grin_dataCls[k]->GetYaxis()->SetRangeUser(xsectrange[0],50); //cfr (0,11) for lumi = 250/pb, (0,26) for lumi = 36.1/pb
	      grin_dataCls[k]->SetFillColor(3);
              grin_dataCls[k]->SetFillStyle(5005);
              grin_dataCls[k]->SetLineColor(3);
              grin_dataCls[k]->SetLineWidth(3);
	      grin_dataCls[k]->SetMarkerColor(3);
              grin_dataCls[k]->SetMarkerSize(1.);
              grin_dataCls[k]->SetMarkerStyle(kFullSquare);
	      grin_dataCls[k]->Draw("CP"); //if you do this ACP, and grin[k] was also ACP, the grin[k] is not drawn for some reason... Without the 'A' here, it works
	   }
	   
	   grin_theo->Draw("CP"); //"CP"
	   leg->Draw();
	   cCls->Write("InterpolatedExclusion_Cls_x"+TString(Fractionstrstream.str()),TObject::kOverwrite);
	   
    }
    cout<<"Output written to "<<TString(SaveFile)<<endl;
    _output->Close();

/*
   // Kernel Smoother
   // // create new kernel smoother and smooth data with bandwidth = 2.0
    cc->Write(" Exlcusion band2", TObject::kOverwrite);


//c->cd(1);
//////// c->SetTickx();
//////// c->SetTicky();
//////// c->SetGridx();
//////// c->SetGridy();



////frame->Draw(" ");

////grin->SetFillColor(4);
////grin->SetFillStyle(5005);
////grin->SetLineColor(4);
////grin->SetLineWidth(3);

////grin->Draw("C");


    // Kernel Smoother
    // // create new kernel smoother and smooth data with bandwidth = 2.0
    c->cd();
    c->Divide(2,3);
    TGraphSmooth *gs = new TGraphSmooth("normal");
    grout = gs->SmoothKern(grin,"normal",2.0);
    DrawSmooth(1,"Kernel Smoother: bandwidth = 2.0","times","accel");
    // redraw ksmooth with bandwidth = 5.0
    grout = gs->SmoothKern(grin,"normal",5.0);
    DrawSmooth(2,"Kernel Smoother: bandwidth = 5.0","","");
    // Lowess Smoother
    // create new lowess smoother and smooth data with fraction f = 2/3
    grout = gs->SmoothLowess(grin,"",0.67);
    DrawSmooth(3,"Lowess: f = 2/3","","");

    // redraw lowess with fraction f = 0.2
    grout = gs->SmoothLowess(grin,"",0.2);
    DrawSmooth(4,"Lowess: f = 0.2","","");
    //

    // create new super smoother and smooth data with default bass = 0 and span = 0
    grout = gs->SmoothSuper(grin,"",0,0);
    DrawSmooth(5,"Super Smoother: bass = 0","box","");
    // redraw supsmu with bass = 3 (smoother curve)
    grout = gs->SmoothSuper(grin,"",3);
    DrawSmooth(6,"Super Smoother: bass = 3","","");
*/

   cout << "	" << "========================================================" << endl;
}

/*void DrawSmooth(Int_t pad, const char *title, const char *xt, const char *yt)
{
	   c->cd(pad);
	   TH1F *vFrame = gPad->DrawFrame(100,100,500,500);
	   vFrame->SetTitle(title);
	   vFrame->SetTitleSize(0.2);
	   vFrame->SetXTitle(xt);
	   vFrame->SetYTitle(yt);
	   grin->Draw("P");
	   grout->DrawClone("LPX");
}*/



