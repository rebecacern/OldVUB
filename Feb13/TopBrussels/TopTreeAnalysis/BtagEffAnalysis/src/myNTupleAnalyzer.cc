//file to define the function of my analyzer which is containing the b-tag measurement
#include "../interface/myNTupleAnalyzer.h"
//#include "../interface/TDRStyle.h"
#include "../interface/PtDistrRadCase.h"
#include "TMath.h"

double getPtWeight(TFile* f, std::string plot, double pt_jet) {

    
    TH1D* hist = (TH1D*) f->Get(plot.c_str());    
	
    if (hist) {
        
        for (int i=1; i<hist->GetNbinsX(); i++) {
            if (hist->GetBinLowEdge(i) <= pt_jet && hist->GetBinLowEdge(i+1) > pt_jet)
                //if (hist->GetBinContent(i) < 20)
                    return hist->GetBinContent(i);
        }
    }
    
    return 1;

}

//con/de-structor

float round(float depth, float n)
{
    float d;
    int    i;
	
    /* rescale 123.45678 to 12345.678 */ 
    d = fabs(n) * depth;
    /* round off: 12345.678 + 0.5 = 12346.178 -> 12346 */ 
    i = d + 0.5;
    /* restore to its original scale: 12346 -> 123.46 */
    d = (float)i / depth;
	
    if (n>=0)
        return d;   
    else
        return -d;
    return 0;
}

float getNomVal(std::string data_postfix, int fitM,std::string prop,string niwp) {
    
    stringstream fitMode; fitMode << fitM;
    string filename="systematics/"+prop+"_WP"+niwp+"_channel"+data_postfix+"_fitMode"+fitMode.str()+".nominal";
    //cout << filename << endl;
	string line;
	ifstream myfile (filename.c_str());
	if (myfile.is_open()) {
		while ( myfile.good() )
			getline (myfile,line);
		myfile.close();
	} 
    
    //cout << line << endl;
	
	return atof(line.c_str());
	
	
}

float getNomVal(std::string data_postfix, int fitM ,std::string prop,int niwp) {
    stringstream p; p<<niwp;
    return getNomVal(data_postfix,fitM,prop,p.str());
}

void setNomVal(std::string data_postfix, int fitM ,std::string prop,int iwp,float val) {
    
    stringstream niwp; niwp << iwp;
    
    stringstream fitMode; fitMode << fitM;

    string filename="systematics/"+prop+"_WP"+niwp.str()+"_channel"+data_postfix+"_fitMode"+fitMode.str()+".nominal";
	
    ofstream myfile (filename.c_str());
    
    stringstream s; s << val;
    myfile << s.str();
    
    myfile.close();
		
	
}

std::pair<float,float> calcBias(float a, float sa, float b, float sb) {

	float bias = (a-b)/b;
	
	float term1 = pow(1/b,2)*pow(sa,2);
	float term2 = pow((-1/b)-((a-b)/pow(b,2)),2)*pow(sb,2);
	
	float uncBias = sqrt(term1+term2);
	
	return std::pair<float,float>(bias,uncBias);
	
}

vector<double> calcStatUnc (bool silent, vector<double> values) {
	
	vector<double> res;
	
	if (values.size() < 11) return res;
	
	double ntt = values[0]; float sntt = values[1];
	double L = values[2]; float sL = values[3];
	double ebtag = values[4]; float sebtag = values[5];
	double emistag = values[6]; float semistag = values[7];
	double FracBoverTotal = values[8];
	double echi2 = values[9];
	double esel = values[10];
	
	double ebtagcut = ( FracBoverTotal * ebtag ) + ( ( 1 - FracBoverTotal ) * emistag );
	
	double sebtagcut = FracBoverTotal*sebtag;
	
	if (ebtag == 1) {
		
		ebtagcut = 1; 
		sebtagcut = 0; 
	
	}
	
	//cout << ntt << " " << sntt << " " << L << " " << ebtagcut << " " << sebtagcut << " " << emistag << " " << semistag << " " << esel << " " << echi2 << endl;
	if (!silent) cout << "calcStatUnc::UncBtagCut - " << sebtagcut << endl;
	
	res.push_back(sebtagcut);
	
	double exs_term1 = (pow(ntt*sebtagcut,2))/(pow(ebtagcut,4)*pow(echi2,2)*pow(esel,2)*pow(L,2));
	double exs_term2 = (pow(ntt*sL,2))/(pow(ebtagcut,2)*pow(echi2,2)*pow(esel,2)*pow(L,4));
	double exs_term3 = (pow(sntt,2))/(pow(ebtagcut,2)*pow(echi2,2)*pow(esel,2)*pow(L,2));
	
	//cout << exs_term1 << " " << exs_term2 << " " << exs_term3 << endl;
	
	double exs = sqrt(exs_term1+exs_term2+exs_term3);
	
	if (!silent) cout << "calcStatUnc::UncXS - " << exs << endl;
	
	res.push_back(exs);
	
	return res;
}

double calcBtag(double XS, double effCHiSq, double effsel, vector<double> FitForXSResults, float desiredIntLum_,int offset=0) {
	
	//float effsel = 0.039842537;
	
	//float effsel = 0.02647067;
	
	//float effCHiSq = 0.868904;
	
	double misTagRate = FitForXSResults[offset+4];
	
	double f = FitForXSResults[FitForXSResults.size()-1];
		
	double btagCutEff = ((FitForXSResults[offset]/(effsel*effCHiSq))*(1/desiredIntLum_))/XS;
	
	double beff = (btagCutEff-((1-f)*misTagRate))/f;
	
	/*
	bceff = f * eb + (1-f) * eq
	
	-> f * eb = bceff - (1-f) * eq
	 
	 -> eb = [ bceff - (1-f) * eq ] / f
	 
	 */
	
	return beff;
}

vector<double> calcXS(bool silent, float btagEff, float uncbtagEff, double effCHiSq, double effmlb, double effsel, vector<double> FitForXSResults, float desiredIntLum_,int offset=0) {
    
    //btagEff = btagEff/0.926;
	
	//float effsel = 0.039842537;
	
	//float effsel = 0.02647067;
	
	//float effCHiSq = 0.868904;
	
	double misTagRate = FitForXSResults[offset+4];
	
	double FracBoverTotal = FitForXSResults[FitForXSResults.size()-1];
	
	//btagEff=0.81;
	
	if (!silent) cout << "Lumi used for XS calc " << desiredIntLum_ << endl;
	
	vector<double> uncvals;
	
	uncvals.push_back(FitForXSResults[offset]);
	uncvals.push_back(FitForXSResults[offset+1]);
	uncvals.push_back(desiredIntLum_);
	
	//uncvals.push_back(desiredIntLum_*0.04);
	uncvals.push_back(0);
	
	uncvals.push_back(btagEff);
	uncvals.push_back(uncbtagEff);
	uncvals.push_back(misTagRate);
	uncvals.push_back(0.);
	uncvals.push_back(FracBoverTotal);
	uncvals.push_back(effCHiSq*effmlb);
	uncvals.push_back(effsel);
	
	vector<double> UncForXSResults = calcStatUnc(silent,uncvals);
	
	//exit(1);
	
	/*float nttsemimu = ((FitForXSResults[8]*0.8723)/0.25176162);
	 float nttsemimuMC = ((FitForXSResults[0]*0.8723)/0.25176162);
	 
	 float eff = nttsemimu*(27./4.)*(1/desiredIntLum_);
	 float effMC = nttsemimuMC*(27./4.)*(1/desiredIntLum_);
	 
	 vector<float> returnVals;
	 
	 returnVals.push_back(nttsemimu);
	 returnVals.push_back(eff);
	 returnVals.push_back(nttsemimuMC);
	 returnVals.push_back(effMC);*/
    
	
	double btagCutEff = ( FracBoverTotal * btagEff ) + ( ( 1 - FracBoverTotal ) * misTagRate );
	
    //btagCutEff=btagCutEff*1.08;

	if (btagEff == 1) btagCutEff=1; // no btag cut
	
	//float ntt = FitForXSResults[offset+8]/(effsel*btagCutEff*effCHiSq);
	double nttMC = FitForXSResults[offset]/(effsel*btagCutEff*effCHiSq*effmlb);	 
	
	//float exp = ntt*(1/desiredIntLum_);
	double expMC = nttMC*(1/desiredIntLum_);
	
	double expMCOld = (FitForXSResults[offset]/(effsel*btagEff*effCHiSq*effmlb))*(1/desiredIntLum_);
	
	vector<double> returnVals;
	
	//returnVals.push_back(ntt);
	//returnVals.push_back(exp);
	returnVals.push_back(nttMC);
	returnVals.push_back(expMC);
	returnVals.push_back(UncForXSResults[1]);
	
	returnVals.push_back(FracBoverTotal);
	
	returnVals.push_back(desiredIntLum_*157.5*effsel*btagEff*effCHiSq*effmlb);
	
	if (!silent) {
		cout << "+> What should be the btag cut eff :" << (expMCOld*btagEff) << " / 157.5 = " << (expMCOld*btagEff)/157.5 << endl;
		cout << "+> What should be the refsel eff :" << (FitForXSResults[offset]/(157.5*btagEff*effCHiSq))*(1/desiredIntLum_) << endl;
		
		cout << "+> eff(sel) = " << effsel << endl;
		cout << "+> eff(ChiSq) = " << effCHiSq << endl;
		cout << "+> eff(mlb) = " << effmlb << endl;
		cout << "+> eff(btag_cut) = " << btagCutEff << endl;
		cout << "+> eff(btag_est) = " << btagEff << " +- " << uncbtagEff << endl;
		cout << "+> ttbar eff(non-b-btag) = " << misTagRate << endl;
		cout << "+> ttbar b/(nb+b) = " << FracBoverTotal << endl;
		cout << "+> Ntt = " << nttMC << endl;
		//cout << "+> NV_fitted = " << nV << endl;
		cout << " => Using MC - Ratio NttObs/NttExp: " << FitForXSResults[offset] << "/" << returnVals[returnVals.size()-1] << " = " << FitForXSResults[offset]/returnVals[returnVals.size()-1] << endl;
		//cout << " => Using controlSample - Ratio NttObs/NttExp: " << FitForXSResults[offset+8] << "/" << returnVals[4] << " = " << FitForXSResults[offset+8]/returnVals[4] << endl;
		//cout << " => Ntt = " << ntt << " NttMC = " << nttMC << endl;
	}
	
	//exit(1);
	return returnVals;
	
}

double getIsoMu20Weight(double mu_pt, double npv) {
    
    const int nbinsx = 6;
    const int nbinsy = 18;
    
    double bin_max_sf_x[nbinsx] = {5, 10, 15, 20, 25, 30};
    double bin_max_sf_y[nbinsy] = {28, 30, 32, 34, 36, 38, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 400};
    double sf[6][18] = {{1.07171, 0.98514, 0.973092, 1.00493, 0.983205, 0.986829, 0.976192, 0.992651, 0.981572, 0.990406, 1.01543, 0.95207, 1.04386, 0.963883, 1.04895, 1.01818, 1.12045, 1.06852}
        , {1.005, 1.0024, 0.999279, 0.99167, 0.994748, 0.98629, 0.989731, 0.986309, 0.987683, 0.973065, 0.98069, 0.985211, 0.992013, 0.969794, 0.985008, 0.98731, 1.02774, 0.983929}
        , {1.00718, 1.0026, 0.998392, 0.996852, 0.993917, 0.991107, 0.995707, 0.988904, 0.985307, 0.981525, 0.983819, 0.996765, 0.977974, 0.990506, 0.974848, 0.982312, 0.991822, 0.982757}
        , {1.01889, 1.02443, 1.02119, 1.0046, 1.00966, 1.00435, 0.997786, 0.991913, 0.987858, 0.988509, 0.979026, 0.97901, 0.961866, 0.983777, 0.972252, 0.958463, 0.942218, 0.962447}
        , {1.04899, 1.04333, 1.03203, 1.04072, 1.01777, 1.01847, 1.02314, 1.01049, 0.995805, 0.999198, 1.00302, 0.995462, 1.00582, 0.964965, 1.00744, 0.96393, 1.01284, 0.965128}
        , {1.13154, 0.956478, 1.03685, 1.04706, 1.0592, 1.02815, 1.00913, 1.00919, 1.01609, 1.03097, 1.02851, 1.04257, 0.963051, 1.07103, 0.861795, 1.04595, 0.769442, 0.946066}
    };
    
    int bin_x = -1; 
    int bin_y = -1;
    
    for(int i=0; i < nbinsy; ++i){
        if(mu_pt < bin_max_sf_y[i]){
            bin_y = i;
            break;
        }
    }
    if(bin_y == -1){ bin_y = nbinsy - 1; }
    
    for(int i=0; i < nbinsx; ++i){
        if(npv < bin_max_sf_x[i]){
            bin_x = i;
            break;
        }
    }
    if(bin_x == -1){ bin_x = nbinsx - 1; }
    
    return sf[bin_x][bin_y];
    
}

myNTupleAnalyzer::myNTupleAnalyzer(TString *inDir, TString *outDir, int *runSamples, int nRunSamples, double chisqCut, int inBin, int outBin, double desiredIntLum, int nPseudoExp, bool doPseudoExp){
	
	nTTbarBeforeChiSq = 0;
	nTTbarAfterChiSq = 0;
	
    nTTbarBeforeRefSel = 0;
    
	//nTTbarBeforeRefSel = 3701936; // summer11 pt 30
    
    //nTTbarBeforeRefSel = 3698723; // fall11 pt 35
	
	//nTTbarBeforeRefSel = 2947570;
	
    //nTTbarBeforeRefSel = 1203372; // summer12 temp sample
    
    //nTTbarBeforeRefSel = 6442045; // summer12 final sample
    
    //nTTbarBeforeRefSel = 6784910; // summer12_DR53X final sample
    
    //nTTbarBeforeRefSel = 6883735; // Summer12_DR52X toptreeprod V2

    nTTbarAfterRefSel = 0;
    
    nTTbarBeforeChiSq = 0;
    nTTbarAfterChiSq = 0;
    
    nTTbarBeforeMLBCUT = 0;
    nTTbarAfterMLBCUT = 0;
    
    for (unsigned int i=0; i<100; i++) nTTbarAfterRefSel_nPV[i]=0;
	
	inDir_=*inDir;
	outDir_=*outDir;
	runSamples_ = runSamples;
	nRunSamples_ = nRunSamples;
	matchChiSquareCut_=chisqCut;
	
	//for(int i=0; i<nRunSamples_; i++){
	//	cout << runSamples_[i] << endl;
	//}
	
    nTaggers=8;
    //nTaggers=3;
	
    if (doPseudoExp) nTaggers=1;
    
	nWP=nTaggers*3;
	
	for(int i=0; i<nWP; i++){
		eff[i]=-999.;
		effVal[i]=-999.;
		effMCVal[i]=-999.;
		effRRVal[i]=-999.;
		effRRMCVal[i]=-999.;
	}
	
	FMCBias=1.;
	FMCVal=-999.;
	FVal=-999.;
	FMCerrVal=-999.;
	FerrVal=-999.;
	
	changeFbias=false;
	backgroundFraction=1;
	
	inBin_=inBin;
	outBin_=outBin;
	
	nIncreaseSample_=1;
	//int increaseSampletmp_[1]={26};
	//int increaseSampletmp_[1]={28};
	int increaseSampletmp_[1]={-99};
	increaseSample_ = new int[1];
	
	//nIncreaseSample_=7;
	//int increaseSampletmp_[7]={37,38,40,41,42,43,44};
	//increaseSample_ = new int[7];
	
	for(int i=0; i<nIncreaseSample_; i++){
		increaseSample_[i]=increaseSampletmp_[i];
	}
	increaseWeightFraction_=1.5;
	
	desiredIntLum_ = desiredIntLum;
	nPseudoExp_ = nPseudoExp;
	doPseudoExp_ = doPseudoExp;    
	percentiles_ = new int[3];
    
};
myNTupleAnalyzer::~myNTupleAnalyzer(){ 
	cout << "Thanks for using me, byebye." << endl;
	//cout << "I'm going with a bang" << endl;
};


void myNTupleAnalyzer::getEff(double* effReturn){
	for(int i=0; i<nWP; i++){
		effReturn[i]=eff[i];
		//cout <<" myNTupleAnalyzer::getEff : " << eff[i] << endl;
	}
}

void myNTupleAnalyzer::getEffVal(double* effValReturn){
	for(int i=0; i<nWP; i++){
		effValReturn[i]=effVal[i];
	}
}

void myNTupleAnalyzer::getEffMCVal(double* effMCValReturn){
	for(int i=0; i<nWP; i++){
		effMCValReturn[i]=effMCVal[i];
	}
}

void myNTupleAnalyzer::getEffRRVal(double* effRRValReturn){
	for(int i=0; i<nWP; i++){
		effRRValReturn[i]=effRRVal[i]; 
		//cout <<" myNTupleAnalyzer::getEffRRVal : " << effRRVal[i] << endl;
	}
}

void myNTupleAnalyzer::getEffRRMCVal(double* effRRMCValReturn){
	for(int i=0; i<nWP; i++){
		effRRMCValReturn[i]=effRRMCVal[i]; 
		//cout <<" myNTupleAnalyzer::getEffRRMCVal : " << effRRMCVal[i] << endl;
	}
}

void myNTupleAnalyzer::getPercentiles(int* returnPerc){
	for(int i=0; i<3; i++){
		returnPerc[i]=percentiles_[i]; 
	}
}

/*void myNTupleAnalyzer::getFMCVal(double *FMCValReturn){
 FMCValReturn[0]=FMCVal;
 }*/
double myNTupleAnalyzer::getFMCVal(){return FMCVal;}
double myNTupleAnalyzer::getFVal(){return FVal;}

double myNTupleAnalyzer::getFMCerrVal(){return FMCerrVal;}
double myNTupleAnalyzer::getFerrVal(){return FerrVal;}

double myNTupleAnalyzer::getmljMeanVal(){return mljMeanVal;}
double myNTupleAnalyzer::getmljMeannoRWVal(){return mljMeannoRWVal;}
double myNTupleAnalyzer::getmljMeanMCVal(){return mljMeanMCVal;}
double myNTupleAnalyzer::getmljSigmaVal(){return mljSigmaVal;}
double myNTupleAnalyzer::getmljSigmanoRWVal(){return mljSigmanoRWVal;}
double myNTupleAnalyzer::getmljSigmaMCVal(){return mljSigmaMCVal;}

double myNTupleAnalyzer::getmlj_W_MeanVal(){return mlj_W_MeanVal;}
double myNTupleAnalyzer::getmlj_W_MeannoRWVal(){return mlj_W_MeannoRWVal;}
double myNTupleAnalyzer::getmlj_W_MeanMCVal(){return mlj_W_MeanMCVal;}
double myNTupleAnalyzer::getmlj_W_SigmaVal(){return mlj_W_SigmaVal;}
double myNTupleAnalyzer::getmlj_W_SigmanoRWVal(){return mlj_W_SigmanoRWVal;}
double myNTupleAnalyzer::getmlj_W_SigmaMCVal(){return mlj_W_SigmaMCVal;}

double myNTupleAnalyzer::getmlj_R_MeanVal(){return mlj_R_MeanVal;}
double myNTupleAnalyzer::getmlj_R_MeannoRWVal(){return mlj_R_MeannoRWVal;}
double myNTupleAnalyzer::getmlj_R_MeanMCVal(){return mlj_R_MeanMCVal;}
double myNTupleAnalyzer::getmlj_R_SigmaVal(){return mlj_R_SigmaVal;}
double myNTupleAnalyzer::getmlj_R_SigmanoRWVal(){return mlj_R_SigmanoRWVal;}
double myNTupleAnalyzer::getmlj_R_SigmaMCVal(){return mlj_R_SigmaMCVal;}

double myNTupleAnalyzer::get_lepb_b_Counter(){return lepb_b_Counter;}
double myNTupleAnalyzer::get_lepb_nonb_Counter(){return lepb_nonb_Counter;}
double myNTupleAnalyzer::get_lepb_hadqq_Counter(){return lepb_hadqq_Counter;}
double myNTupleAnalyzer::get_lepb_radq_Counter(){return lepb_radq_Counter;}
double myNTupleAnalyzer::get_lepbtag_b_Counter(){return lepbtag_b_Counter;}
double myNTupleAnalyzer::get_lepbtag_nonb_Counter(){return lepbtag_nonb_Counter;}
double myNTupleAnalyzer::get_lepbtag_hadqq_Counter(){return lepbtag_hadqq_Counter;}
double myNTupleAnalyzer::get_lepbtag_radq_Counter(){return lepbtag_radq_Counter;}
double myNTupleAnalyzer::get_q1q2_b_Counter(){return q1q2_b_Counter;}
double myNTupleAnalyzer::get_q1q2_nonb_Counter(){return q1q2_nonb_Counter;}
double myNTupleAnalyzer::get_q1q2_hadqq_Counter(){return q1q2_hadqq_Counter;}
double myNTupleAnalyzer::get_q1q2_radq_Counter(){return q1q2_radq_Counter;}


void myNTupleAnalyzer::setFMCBias(double bias){
	cout << "WARNING    myNTupleAnalyzer::setFBias - you are changing the F ratio with " << bias << " (know what you are doing)" << endl;
	FMCBias=bias;
	changeFbias=true;
	cout << "myNTupleAnalyzer::setFMCBias changeFbias " << changeFbias << endl;
}

void myNTupleAnalyzer::setBackgroundFraction(double fraction){
	cout << "WARNING   myNTupleAnalyzer::setBackgroundFraction - you are changing the background fraction to " << fraction << " (know what you are doing)" << endl;
	backgroundFraction=fraction;
	changeBackgroundFraction=true;
}

//void setTDRStyle();

//change in and out names to accept from function

int nrFiles = 20;

TString inRootFile[20];

TString outRootPath[2]; 
TString outPath;

TString outRootName = "NTupleAnalyzed.root";
TString outRootNamePseudoExp = "NTupleAnalyzed_pseudoExp.root";

//obsolete
/*const double desiredIntLum=1000; //TO CHECK: should this be the same as the value I used for calculating the weights?
 int nPseudoExp=1;//= is there to run on all data
 bool doPseudoExp=false; //to do the pseudo-exps
 */

/*const double desiredIntLum=1000;
 int nPseudoExp=300;
 bool doPseudoExp=true; //to do the pseudo-exps
 */
//end obsolete

// binning for left-right pt-reweighting
const int nBinsVarLR[2]={50,20};
const double lowRangVarLR[2]={0,0};
const double upRangeVarLR[2]={500,2.4};

// binning for SC reweighting (fixed) but obsolete because of variable binning in PTEtaBin.cc

//const int nBinsVarSC[2]={20,10}; //for 1D, then the 10 is not taken into account
//const int nBinsVarSC[2]={51,20}; //for 1D, then the 10 is not taken into account
const int nBinsVarSC[2]={18,20}; //for 1D, then the 10 is not taken into account <-- THIS IS STILL USED
//const int nBinsVarSC[2]={20,10};
const double lowRangVarSC[2]={0,0};
const double upRangeVarSC[2]={500,2.4};
//const double upRangeVarSC[2]={500,3.5};

//const double lowRangVarSC[2]={0,-1};
//const double upRangeVarSC[2]={500,1};

/*const int nBinsVarSC[2]={10,50};
 const double lowRangVarSC[2]={0,0};
 const double upRangeVarSC[2]={500,5};
 */

bool doVarBins=true; // for btaggers
//bool doVarBins=true;

// shift btag distribution instead of reweighting, works only for one tagger so OBSOLETE

bool doShift=false; //warning: when switching this to true make sure the correct shift is applied!
//bool doShift=true;

bool doPtEtaBin=false;
//bool doPtEtaBin=true;

// MLB binning
/*const int nBinsVar0=50;
//const int nBinsVar0=500;
const double lowRangeVar0=0;
const double upRangeVar0=500;*/

int nBinsVar0=50;
//const int nBinsVar0=500;
double lowRangeVar0=0;
double upRangeVar0=500;

// binning btaggers (fix binning)


// OLD TAGGER ORDERING
/*int nBinsBtag[8]={25,25,25,25,25,25,25,25};
double lowRangeBtag[8]={-10,-10,0,0,0,0,0,0};
double upRangeBtag[8]={30,30,8,8,1,1,8,8};
//double wpArray[24]={1.7,3.3,10.2,1.19,1.93,3.41,0.00,1.74,3.05,0.00,0.00,2.00,0.244,0.679,0.898,0.244,0.679,0.898,0.275,0.545,0.790,1.33,2.55,3.74};
double wpArray[24]={1.7,3.3,10.2,1.19,1.93,3.41,0.00,1.74,3.05,0.00,0.00,2.00,0.244,0.679,0.898,0,0,0,1.33,2.55,3.74,0.275,0.545,0.790}; // old ordering
int taggerArray[24]={0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7};
*/

// NEW BTV-11-003 ORDERING

int nBinsBtag[8]={25,25,25,25,25,25,25,25};
double lowRangeBtag[8]={-10,-10,0,0,0,0,0,0};
double upRangeBtag[8]={30,30,3,8,8,8,1,1};

double wpArray[24]={1.7,3.3,10.2,1.19,1.93,3.41,0.275,0.545,0.790,1.33,2.55,3.74,0.00,1.74,3.05,0.00,0.00,2.0,0.244,0.679,0.898,0.455,0.820,0.940}; // old ordering
int taggerArray[24]={0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7};
 

vector< double > bValueBdisc[11]; //needed for the calculation of the variable bin size
bool dobValueBdisc=true; // calculate variable bins??

vector< double > ptValues; // for running method in pt,eta bins
int nPtValues = 5;
int nEtaValues = 5;
double ptBinsSet[6]={30.,45.,60.,80.,110.,999999.};
double etaBinsSet[6]={0,0.48,0.96,1.44,1.92,2.4};
/*
 int nPtValues = 1;
 int nEtaValues = 1;
 double ptBinsSet[1]={0};
 double etaBinsSet[1]={0};
 */
/*int nPtValues = 8;
 int nEtaValues = 8;
 //double ptBinsSet[9]={30.,40.,50.,60.,75.,90.,110.,150.,999999.};
 double ptBinsSet[9]={30.,40.,50.,60.,75.,90.,110.,150.,300.};
 double etaBinsSet[9]={0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4};*/

/*int nPtValues = 1;
 int nEtaValues = 20;
 double etaBinsSet[21]={0,0.12,0.24,0.36,0.48,0.6,0.72,0.84,0.96,1.08,1.2,1.32,1.44,1.56,1.68,1.8,1.92,2.04,2.16,2.28,2.4};
 double ptBinsSet[2]={0,999999};*/


/*int nPtValues = 20;
 int nEtaValues = 1;
 double ptBinsSet[21]={30,34,38,42,46,49,53,58,62,66,71,76,81,86,92,100,108,118,132,155,999999};
 double etaBinsSet[2]={0,999.};*/


vector< double > ptValuesSC; //this are the pt values for the SC reweighing (Don't remember exactly, but I think it is for the variable bins)
int nPtValuesSC = 20;

vector< double > ptValuesBin;
vector< double > ptValuesBinWeight;
vector< double > etaValuesBin; //filled with |absolute| value
vector< double > etaValuesBinWeight; //is here for sanity check

bool doCoutBinning = true;

vector<int> skipFlavours; // obsolete!!!

//vector<pair<int,int> > skipEvents;
bool skipEvents[100][500000]; //there should be not more than 50000 events in one sample
int skipEventsCounter[100]; //there should be not more than 50000 events in one sample
int skipEventsCounter10[100]; //there should be not more than 50000 events in one sample
//bool skipEventsEventSel[100][50000]; //there should be not more than 50000 events in one sample


// new WPArray

//histograms for chiSquare distributions
/*TH1D * ChiSq[10]; 
 TH2D * ChiSq_vs_Mass[2];
 TH2D * E_vs_Mass[2];
 TH2D * Pt_vs_Mass[2];
 TH2D * Pt_vs_E[1];
 TH1D * ChiSq_R_inAll[2]; 
 TH1D * ChiSq_R_inHad[2]; 
 TH1D * ChiSq_R_or_lepb_inHad[2]; 
 TH1D * E_R_inAll[2]; 
 TH1D * E_R_inHad[2]; */

std::map<std::string, TH1D*> MljDistributions; 

//double chisqCutForMass[28]={1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,15.,20.,30.,50.,9999999999.,100.,500.,1000.,5000.,};
//double chisqCutForMass[28]={1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,15.,20.,30.,50.,9999999999.,100.,500.,1000.,5000.,};

//double chisqCutForMass[24]={1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.,20.,25.,30.,40.,50.,60.,70.,85.,100.,500.,1000.,5000.,9999999999.}; // thesis p117 6.22
double chisqCutForMass[24]={1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.,20.,25.,30.,40.,50.,60.,70.,85.,100.,500.,1000.,5000.,10000.};

TH1D * TopMass[24];
TH1D * WMass[24];
double TopMassmax[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double WMassmax[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double TopMassmin[24]={99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.};
double WMassmin[24]={99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.};

//TH1D * cosDelSA[8]; //angle between lepton and leptonic b jet candidate
//TH1D * delSA[8];

TH1D * ptBin;

/*TH1D *delOmegalj[2];
 TH1D *pthadtop[2];
 TH1D *delPhitt[2];
 TH1D *etafifth[2];
 TH1D *ptfifth[2];
 TH1D *btag_trackCountingHighEffBJetTags_hadbPlot[2];
 TH1D *nJetsPlots[2];
 
 //plot fo b-tag efficiency to find bug
 TH1D *btagbug;
 TH1D *btagbugEff;*/

//TH1D * Pt_R_inAll[2]; //replace by PtDistrRadCase
//TH1D * Pt_R_inHad[2]; 

//int main(int argc, char* argv[]){

void myNTupleAnalyzer::run(int verbosity, int leftlimit, int centerleftlimit, int centerrightlimit, int rightlimit, bool doSCreweigh, bool doTwoLights, bool useFit, bool do2D, bool do2Dcontrol, bool doPtEtaBin, bool doJESchange, double JESfactor, bool doNewF, double leftlimitperc, double centerlimitperc, double rightlimitperc, double runNb, bool doFfromMC, int nSystematic,int decay,int fitMode,int nBtag_){ //decay 0: semimu 1: semiel
	
    /*float pom = 0.0024126;
    
    cout << round(10,pom*100) << endl;
    cout << round(100,pom*100) << endl;
    cout << round(10,pom*1000) << endl;
    cout << round(100,pom*1000) << endl;
    
    exit(1);*/

    useTTJetsExcl_=false;

	bool doOnlyMSPlot = false;
    
    bool btagonly = false;
    
    bool doBiasStudy=true;

    bool doFBiasStudy=false;
    
    if (!doPseudoExp_) doBiasStudy=false;
    if (!doPseudoExp_) doFBiasStudy=false;
    
    if (doBiasStudy && doFBiasStudy) {
        cout << "doBiasStudy and doFBiasStudy cannot both be true" << endl;
        exit(1);
    }

	double left_min=leftlimit;
    double left_max=leftlimit+1;
	double mid_min=120;
    double mid_max=200;
    double right_max=400;
    double bias_step=5;
    
    double F_min=2;
    double F_max=2002;
    double F_step=25;
    
    // reweighing of enriched pt to data
    
    bool reweigh_to_data=false; // removes the <pT> difference between data enr and MC enr
    bool reweigh_left_to_all=false; // does not really work, do not use
    
    bool doCorrectForBias=false; // should not be turned on to obtain results
    bool doBtagSFonData=false; // should not be turned on to obtain results
    
    int nBtag = nBtag_;

    if (doPseudoExp_ && nBtag != -1) {        
        
        nTaggers=1;
        
        nWP=3;
        
        stringstream np; np << nBtag;
        
        if (decay == 0)
            outRootNamePseudoExp = "NTupleAnalyzed_pseudoExp_tagger_"+np.str()+"_channel_Mu.root";
        else if (decay == 1)
            outRootNamePseudoExp = "NTupleAnalyzed_pseudoExp_tagger_"+np.str()+"_channel_El.root";

        cout << "wpArray before: " << wpArray[0] << " " << wpArray[1] << " " << wpArray[2] << endl;
        cout << "taggerArray before: " << taggerArray[0] << " " << taggerArray[1] << " " << taggerArray[2] << endl;
        cout << "btag ranges before: " << nBinsBtag[0] << " " << lowRangeBtag[0] << " " << upRangeBtag[0] << endl;
        wpArray[0]=wpArray[3*nBtag];
        wpArray[1]=wpArray[(3*nBtag)+1];
        wpArray[2]=wpArray[(3*nBtag)+2];

        //taggerArray[0]=taggerArray[3*nBtag];
        //taggerArray[1]=taggerArray[(3*nBtag)+1];
        //taggerArray[2]=taggerArray[(3*nBtag)+2];
        
        nBinsBtag[0]=nBinsBtag[nBtag];
        upRangeBtag[0]=upRangeBtag[nBtag];
        lowRangeBtag[0]=lowRangeBtag[nBtag];
        
        cout << "wpArray after: " << wpArray[0] << " " << wpArray[1] << " " << wpArray[2] << endl;
        cout << "taggerArray after: " << taggerArray[0] << " " << taggerArray[1] << " " << taggerArray[2] << endl;
        cout << "btag ranges before: " << nBinsBtag[0] << " " << lowRangeBtag[0] << " " << upRangeBtag[0] << endl;

        //exit(1);
        
    }
		
	cout << "++nSystematic: " << nSystematic << endl; // started as systematic sample switcher, but is now more generally used to switch between sample setups.
	string sampleType = "";
	string sampleSel = "ALL";
	
    cout << "++Decay channel selected: ";
    string data_postfix="";
    if (decay == 0) {
        cout << " 0 (Muon+jets)" << endl;
        data_postfix="_Mu";
    }
    else if (decay == 1) {
        cout << " 1 (Electron+jets)" << endl;
        data_postfix="_El";
    }
    else {
        cout << " " << decay << " (Unknown, ERROR, Exiting) " << endl;
        exit(-1);
    }
    
	float WJets_scale=1;
	float ZJets_scale=1;
    float ST_scale=1;

	if (doPseudoExp_) nSystematic = 0;
	
	/*if (nSystematic != 0 && nSystematic != -2 && nSystematic != -6 && nSystematic != -7) {
		nTaggers = 1;
		nWP = nTaggers*3;
	}*/
	
	int nPDFWeight = -1;
    
	if (nSystematic > 100 && nSystematic < 145) {
		nPDFWeight = nSystematic-100;
		nSystematic = 0;
	}
    
    int pvBinA = -1;
    int pvBinB = -1;
    
    if (nSystematic > 200 && nSystematic < 235) {
		pvBinA = (nSystematic-201)*5;
		pvBinB = pvBinA+5;
		nSystematic = 0;
	}

    int nBinsCSVfix=15;
    
    int i=-1;
    
    
    // NOMINAL SET OF SAMPLES
    nTaggers=1;
    nWP=3;
    
    if (useTTJetsExcl_) {
        nTTbarBeforeRefSel = 24460323/0.438048; //(24460323*0.438048)+(11753993*0.104976)+(31131582*0.456976); // Summer12_DR53X prod V4 3 ttjets samples
        nRunSamples_=15;
        inRootFile[0] = "BtagTree_TTbarJets_SemiLepton";
        inRootFile[1] = "BtagTree_TTbarJets_DiLepton";
        inRootFile[2] = "BtagTree_TTbarJets_FullHadronic";
        inRootFile[3] = "BtagTree_ST_tW_t";
        inRootFile[4] = "BtagTree_ST_tW_tbar";
        inRootFile[5] = "BtagTree_ST_t_t";
        inRootFile[6] = "BtagTree_ST_t_tbar";
        inRootFile[7] = "BtagTree_WJets_3jets";
        inRootFile[8] = "BtagTree_WJets_4jets";
        inRootFile[9] = "BtagTree_ZJets_3jets";
        inRootFile[10] = "BtagTree_ZJets_4jets";
        inRootFile[11] = "BtagTree_WJets_1jets";
        inRootFile[12] = "BtagTree_WJets_2jets";
        inRootFile[13] = "BtagTree_ZJets_1jets";
        inRootFile[14] = "BtagTree_ZJets_2jets";
    } else {
        nTTbarBeforeRefSel = 6830443; // Summer12_DR53X prod V4 incl ttjets sample
        nRunSamples_=13;
        inRootFile[0] = "BtagTree_TTbarJets";
        inRootFile[1] = "BtagTree_ST_tW_t";
        inRootFile[2] = "BtagTree_ST_tW_tbar";
        inRootFile[3] = "BtagTree_ST_t_t";
        inRootFile[4] = "BtagTree_ST_t_tbar";
        inRootFile[5] = "BtagTree_WJets_3jets";
        inRootFile[6] = "BtagTree_WJets_4jets";
        inRootFile[7] = "BtagTree_ZJets_3jets";
        inRootFile[8] = "BtagTree_ZJets_4jets";
        inRootFile[9] = "BtagTree_WJets_1jets";
        inRootFile[10] = "BtagTree_WJets_2jets";
        inRootFile[11] = "BtagTree_ZJets_1jets";
        inRootFile[12] = "BtagTree_ZJets_2jets";
    }
    
    if (doPseudoExp_) {
        
        nRunSamples_=nRunSamples_-4;
            
    }

    // tmp disabled Z+Xjets because one of the syst files failed
    /*nRunSamples_=10;
    inRootFile[0] = "BtagTree_TTbarJets";
    inRootFile[1] = "BtagTree_ST_tW_t";
    inRootFile[2] = "BtagTree_ST_tW_tbar";
    inRootFile[3] = "BtagTree_ST_t_t";
    inRootFile[4] = "BtagTree_ST_t_tbar";
    inRootFile[5] = "BtagTree_ZJets";
    inRootFile[6] = "BtagTree_WJets_1jets";
    inRootFile[7] = "BtagTree_WJets_2jets";
    inRootFile[8] = "BtagTree_WJets_3jets";
    inRootFile[9] = "BtagTree_WJets_4jets";*/

    /*nRunSamples_=3;
    inRootFile[0] = "BtagTree_TTbarJets";
    inRootFile[1] = "BtagTree_ZJets";
    inRootFile[2] = "BtagTree_WJets";
*/
    //inRootFile[0] = "BtagTree_WJets_JESMinus";
    
    // RE-ENABLE
    
    /*
    //if (decay==1) {
        inRootFile[nRunSamples_] = "BtagTree_Data"+data_postfix+"_InvIso";
        //inRootFile[nRunSamples_] = "BtagTree_Data_Mu_InvIso";
        nRunSamples_++;
    //}
    */
    
    /*if (nSystematic>0) {
        nTaggers=3;
        nWP=nTaggers*3;
    }*/
    
    bool doQCD=false;

	switch (nSystematic) {
		case -2:
			sampleType="Data";
			nRunSamples_=1;
			inRootFile[0] = "BtagTree_Data"+data_postfix;
			//inRootFile[0] = "BtagTree_Data1fb";
            doCorrectForBias=false;
            doBtagSFonData=false;
            reweigh_to_data=false;
			break;
        case -9:
			sampleType="Data";
			nRunSamples_=1;
			inRootFile[0] = "BtagTree_Data"+data_postfix;
			//inRootFile[0] = "BtagTree_Data1fb";
            doCorrectForBias=true;
            doBtagSFonData=false;
			break;
		case -3:
			doOnlyMSPlot=true;
			
            nRunSamples_=14;
            if (useTTJetsExcl_)
                nRunSamples_=16;
                        
			sampleType="MSPLotProduction";
            
            //if (decay==1) { // for now no data-driven qcd for muon
                        
            if (doQCD) {
                i=0;
                inRootFile[i] = "BtagTree_Data"+data_postfix+"_InvIso";
                nRunSamples_++;
			}
            
            //nRunSamples_=8;
            //i=0;
            //i++;
			//inRootFile[i] = "BtagTree_ZJets";
			//i++;inRootFile[i] = "BtagTree_WJets";
			i++;inRootFile[i] = "BtagTree_ZJets_1jets";
			i++;inRootFile[i] = "BtagTree_ZJets_2jets";
            i++;inRootFile[i] = "BtagTree_ZJets_3jets";
			i++;inRootFile[i] = "BtagTree_ZJets_4jets";
            
            i++;inRootFile[i] = "BtagTree_WJets_1jets";
			i++;inRootFile[i] = "BtagTree_WJets_2jets";
            i++;inRootFile[i] = "BtagTree_WJets_3jets";
			i++;inRootFile[i] = "BtagTree_WJets_4jets";
            //nRunSamples_--;
            
            i++;inRootFile[i] = "BtagTree_ST_t_t";
            i++;inRootFile[i] = "BtagTree_ST_t_tbar";
            i++;inRootFile[i] = "BtagTree_ST_tW_t";
            i++;inRootFile[i] = "BtagTree_ST_tW_tbar";
            
            if (useTTJetsExcl_) {
                i++;inRootFile[i] = "BtagTree_TTbarJets_SemiLepton";
                i++;inRootFile[i] = "BtagTree_TTbarJets_DiLepton";
                i++;inRootFile[i] = "BtagTree_TTbarJets_FullHadronic";
            } else
                i++;inRootFile[i] = "BtagTree_TTbarJets";

			i++;inRootFile[i] = "BtagTree_Data"+data_postfix;
            
            //cout << i << " " << doQCD << " " << nRunSamples_ << endl; exit(1);
			
            //if (decay==1) 
            if (doQCD)datasets_color.push_back(kYellow); // multijet
			datasets_color.push_back(kAzure-2); // ZJets
			datasets_color.push_back(kAzure-2); // ZJets
			datasets_color.push_back(kAzure-2); // ZJets
			datasets_color.push_back(kAzure-2); // ZJets

			datasets_color.push_back(kGreen-3); // WJets
			datasets_color.push_back(kGreen-3); // WJets
			datasets_color.push_back(kGreen-3); // WJets
			datasets_color.push_back(kGreen-3); // WJets
			
			datasets_color.push_back(kMagenta); // ST t
			datasets_color.push_back(kMagenta); // ST t
			datasets_color.push_back(kMagenta); // ST tW
			datasets_color.push_back(kMagenta); // ST tW
			
			/*datasets_color.push_back(kMagenta+2); // ST t
			datasets_color.push_back(kMagenta+2); // ST t
			datasets_color.push_back(kMagenta); // ST tW
			datasets_color.push_back(kMagenta); // ST tW
             */
			
			datasets_color.push_back(kRed-7); // TT
			datasets_color.push_back(kRed+1); // TT
			
            datasets_color.push_back(kBlack); // Data
			datasets_color.push_back(0); // Data

			//datasets_title.push_back("QCD");
            //if (decay==1) 
            if (doQCD)datasets_title.push_back("multijet");
			datasets_title.push_back("Z/#gamma*#rightarrowl^{+}l^{-}");
			datasets_title.push_back("Z/#gamma*#rightarrowl^{+}l^{-}");
			datasets_title.push_back("Z/#gamma*#rightarrowl^{+}l^{-}");
			datasets_title.push_back("Z/#gamma*#rightarrowl^{+}l^{-}");
			datasets_title.push_back("W#rightarrowl#nu");
			datasets_title.push_back("W#rightarrowl#nu");
			datasets_title.push_back("W#rightarrowl#nu");
			datasets_title.push_back("W#rightarrowl#nu");
			datasets_title.push_back("SingleTop t");
			datasets_title.push_back("SingleTop t");
			datasets_title.push_back("SingleTop tW");
			datasets_title.push_back("SingleTop tW");
			
            datasets_title.push_back("t#bar{t} other");
            datasets_title.push_back("t#bar{t} signal");
    
			datasets_title.push_back("Data");
			
			break;
            
		case 0:
			sampleType="nominal";
            reweigh_to_data=false;
			/*nRunSamples_=7;
			inRootFile[0] = "BtagTree_TTbarJets";
			inRootFile[1] = "BtagTree_ST_SingleTop_tChannel_t";
			inRootFile[2] = "BtagTree_ST_SingleTop_tChannel_tbar";
			inRootFile[3] = "BtagTree_ST_SingleTop_tWChannel_t";
			inRootFile[4] = "BtagTree_ST_SingleTop_tWChannel_tbar";
			inRootFile[5] = "BtagTree_ZJets";
			inRootFile[6] = "BtagTree_WJets";*/
            
			if (nPDFWeight != -1) {
				stringstream s; s<<nPDFWeight;
				sampleType="PDF-Uncertainty-Weight-"+s.str();
			}
            
            if (pvBinA != -1) {
				stringstream s; s<<pvBinA;
				stringstream t; t<<pvBinB;
				sampleType=s.str()+"<nPV<="+t.str();
                btagonly=true;
			}
            
            //doCorrectForBias=true;
			
			break;
            
        case -4:
            sampleType="fixbins";
            reweigh_to_data=false;
			doVarBins=false;
			break;

  
        case -8: // apply pt-rew to data
			sampleType="nominal-datarew";
            reweigh_to_data=false;
            doCorrectForBias=true;
			
			break;

		case 1:
			sampleType="JES-";
//doVarBins=false;
			for (unsigned int i=0; i<nRunSamples_; i++) {
                if (((string)inRootFile[i]).find("InvIso") == -1)
                if (((string)inRootFile[i]).find("InvIso") == -1) inRootFile[i] = inRootFile[i]+"_JESMinus";
            }
            break;
		case 2:
			sampleType="JES+";
			for (unsigned int i=0; i<nRunSamples_; i++) {
                if (((string)inRootFile[i]).find("InvIso") == -1) inRootFile[i] = inRootFile[i]+"_JESPlus";
            }
			break;

		case 3:
			sampleType="JER-";
			for (unsigned int i=0; i<nRunSamples_; i++) {
                if (((string)inRootFile[i]).find("InvIso") == -1) inRootFile[i] = inRootFile[i]+"_JERMinus";
            }
			break;
		case 4:
			sampleType="JER+";
			for (unsigned int i=0; i<nRunSamples_; i++) {
                if (((string)inRootFile[i]).find("InvIso") == -1) inRootFile[i] = inRootFile[i]+"_JERPlus";
            }
			break;
            
		case 5:
			sampleType="PileUp-";
			for (unsigned int i=0; i<nRunSamples_; i++) {
                if (((string)inRootFile[i]).find("InvIso") == -1) inRootFile[i] = inRootFile[i]+"_LessPU";
            }
			break;
		case 6:
			sampleType="PileUp+";
			for (unsigned int i=0; i<nRunSamples_; i++) {
                if (((string)inRootFile[i]).find("InvIso") == -1) inRootFile[i] = inRootFile[i]+"_MorePU";
            }
			break;
			
		case 7:
			sampleType="TTJets-ScaleDown";
            sampleSel="TTJets";
			nTTbarBeforeRefSel=5346767;

            for (int s=0;s<nRunSamples_;s++){
                if (((string)inRootFile[s]).find("TTbarJets") != string::npos) {
                    string rename = ((string)inRootFile[s])+"_ScaleDown";
                    inRootFile[s]=(TString) rename;
                }
			}
            
			break;
			
		case 8:
			sampleType="TTJets-ScaleUp";
            sampleSel="TTJets";
			nTTbarBeforeRefSel=4619133;
			
			for (int s=0;s<nRunSamples_;s++){
                if (((string)inRootFile[s]).find("TTbarJets") != string::npos) {
                    string rename = ((string)inRootFile[s])+"_ScaleUp";
                    inRootFile[s]=(TString) rename;
                }
			}
            
            break;

		case 9:
			sampleType="TTJets-MatchingUp";
            sampleSel="TTJets";
			nTTbarBeforeRefSel=5415003;
                               
			for (int s=0;s<nRunSamples_;s++){
                if (((string)inRootFile[s]).find("TTbarJets") != string::npos) {
                    string rename = ((string)inRootFile[s])+"_MatchingUp";
                    inRootFile[s]=(TString) rename;
                }
			}

            break;
			
		case 10:
			sampleType="TTJets-MatchingDown";
            sampleSel="TTJets";
			nTTbarBeforeRefSel=5476715;
			
			for (int s=0;s<nRunSamples_;s++){
                if (((string)inRootFile[s]).find("TTbarJets") != string::npos) {
                    string rename = ((string)inRootFile[s])+"_MatchingDown";
                    inRootFile[s]=(TString) rename;
                }
			}
            
            break;

		case 11:
			sampleType="WJets-ScaleUp";
            sampleSel="WJets";
            nRunSamples_=8;
			//desiredIntLum_=300;
			inRootFile[0] = "BtagTree_TTbarJets_SemiMuon";
			inRootFile[1] = "BtagTree_TTbarJets_Other";
			inRootFile[2] = "BtagTree_ST_SingleTop_tChannel_t";
			inRootFile[3] = "BtagTree_ST_SingleTop_tChannel_tbar";
			inRootFile[4] = "BtagTree_ST_SingleTop_tWChannel_t";
			inRootFile[5] = "BtagTree_ST_SingleTop_tWChannel_tbar";
			inRootFile[6] = "BtagTree_ZJets";
			inRootFile[7] = "BtagTree_WJets_TuneZ2_scaleup";
			inRootFile[8] = "BtagTree_QCD";
			break;
			
		case 12:
			sampleType="WJets-ScaleDown";
            sampleSel="WJets";
			nRunSamples_=8;
			desiredIntLum_=300;
			inRootFile[0] = "BtagTree_TTbarJets_SemiMuon";
			inRootFile[1] = "BtagTree_TTbarJets_Other";
			inRootFile[2] = "BtagTree_ST_SingleTop_tChannel_t";
			inRootFile[3] = "BtagTree_ST_SingleTop_tChannel_tbar";
			inRootFile[4] = "BtagTree_ST_SingleTop_tWChannel_t";
			inRootFile[5] = "BtagTree_ST_SingleTop_tWChannel_tbar";
			inRootFile[6] = "BtagTree_ZJets";
			inRootFile[7] = "BtagTree_WJets_TuneZ2_scaledown";
			inRootFile[8] = "BtagTree_QCD";
			break;
			
		case 13:
			sampleType="WJets-MatchingUp";
            sampleSel="WJets";
            nRunSamples_=8;
			desiredIntLum_=300;
			inRootFile[0] = "BtagTree_TTbarJets_SemiMuon";
			inRootFile[1] = "BtagTree_TTbarJets_Other";
			inRootFile[2] = "BtagTree_ST_SingleTop_tChannel_t";
			inRootFile[3] = "BtagTree_ST_SingleTop_tChannel_tbar";
			inRootFile[4] = "BtagTree_ST_SingleTop_tWChannel_t";
			inRootFile[5] = "BtagTree_ST_SingleTop_tWChannel_tbar";
			inRootFile[6] = "BtagTree_ZJets";
			inRootFile[7] = "BtagTree_WJets_TuneZ2_matchingup";
			inRootFile[8] = "BtagTree_QCD";
			break;
			
		case 14:
			sampleType="WJets-MatchingDown";
            sampleSel="WJets";
            nRunSamples_=8;
			desiredIntLum_=300;
			inRootFile[0] = "BtagTree_TTbarJets_SemiMuon";
			inRootFile[1] = "BtagTree_TTbarJets_Other";
			inRootFile[2] = "BtagTree_ST_SingleTop_tChannel_t";
			inRootFile[3] = "BtagTree_ST_SingleTop_tChannel_tbar";
			inRootFile[4] = "BtagTree_ST_SingleTop_tWChannel_t";
			inRootFile[5] = "BtagTree_ST_SingleTop_tWChannel_tbar";
			inRootFile[6] = "BtagTree_ZJets";
			inRootFile[7] = "BtagTree_WJets_TuneZ2_matchingdown";
			inRootFile[8] = "BtagTree_QCD";
			break;	
			
		case 17:
			sampleType="WJetsXSDown-50perc";
            sampleSel="WJets";
			WJets_scale=0.5;	
            break;
			
		case 18:
			sampleType="WJetsXSUp-50perc";
            sampleSel="WJets";
			WJets_scale=1.5;
        
			break;

		case 19:
			sampleType="ZJetsXSDown-50perc";
            sampleSel="ZJets";
			ZJets_scale=0.5;	
				
			break;			
		case 20:
			sampleType="ZJetsXSUp-50perc";
            sampleSel="ZJets";
			ZJets_scale=1.5;
			
			break;
			
		case 21:
			sampleType="RightReg-nominal";
			
			leftlimit=70;
			centerleftlimit=170;
			rightlimit=300;
			
			break;
		
		case 22:
			sampleType="RightReg-smaller";
						
			leftlimit=70;
			centerleftlimit=170;
			rightlimit=250;
			
			break;
			
		case 23:
			sampleType="RightReg-bigger";
						
			leftlimit=70;
			centerleftlimit=170;
			rightlimit=350;
			
			break;
			
        case 24:
			sampleType="TopMass-163.5";
            sampleSel="TTJets";
			nRunSamples_=7;
			nTTbarBeforeRefSel=1633191;
			inRootFile[0] = "BtagTree_TTbarJets_Mass_163_5_GeV";
			inRootFile[1] = "BtagTree_ST_SingleTop_tChannel_t";
			inRootFile[2] = "BtagTree_ST_SingleTop_tChannel_tbar";
			inRootFile[3] = "BtagTree_ST_SingleTop_tWChannel_t";
			inRootFile[4] = "BtagTree_ST_SingleTop_tWChannel_tbar";
			inRootFile[5] = "BtagTree_ZJets";
			inRootFile[6] = "BtagTree_WJets";
			
            break;
            
        case 25:
			sampleType="TopMass-181.5";
            sampleSel="TTJets";
			nRunSamples_=7;
			nTTbarBeforeRefSel=1665344;
			inRootFile[0] = "BtagTree_TTbarJets_Mass_181_5_GeV";
			inRootFile[1] = "BtagTree_ST_SingleTop_tChannel_t";
			inRootFile[2] = "BtagTree_ST_SingleTop_tChannel_tbar";
			inRootFile[3] = "BtagTree_ST_SingleTop_tWChannel_t";
			inRootFile[4] = "BtagTree_ST_SingleTop_tWChannel_tbar";
			inRootFile[5] = "BtagTree_ZJets";
			inRootFile[6] = "BtagTree_WJets";
			
			break;
            
        default:
			break;
			
	}
    
	//int nTTbarFilled = 0;
	
	int nGoodChi2Combi = 0;
	int nChi2Combi = 0;	
	
	bool doChi2Check = true;
	
	float ptcutextra = 30;
	
	time_t curr=time(0);
	cout << "++current start time is: " << ctime(&curr) <<endl;
	cout << "++setting TDRStyle()" << endl;
	setMyStyle();
	
	cout << "++defining verbosity" << endl;
	if (doPseudoExp_) verbosity=0;
	
	int verbosity_=verbosity;
	//verbosity_=atoi(argv[1]);
	if(verbosity_==0)  cout << "silent mode" << endl;
	if(verbosity_==1)  cout << "verbosity level: " << verbosity_ << " - all output mode" << endl;
	
	int leftlimit_=0;
	int centerleftlimit_=0;
	int centerrightlimit_=0;
	int rightlimit_=0;
	
	leftlimit_=leftlimit; 
	centerleftlimit_=centerleftlimit;
	centerrightlimit_=centerrightlimit;
	rightlimit_=rightlimit;
	/*  leftlimit_=atoi(argv[2]);
	 centerlimit_=atoi(argv[3]);
	 rightlimit_=atoi(argv[4]);*/
	
	outPath = inDir_ + outDir_;
	
	outRootPath[0] = outPath + outRootName;
	outRootPath[1] = outPath + outRootNamePseudoExp;
	
	cout << outPath << outRootName << " " << endl;
	cout << outPath << outRootNamePseudoExp << " " << endl;
	
	if (verbosity_>-1) {
		cout << "++start" << endl;
		cout << "nTaggers    : " << nTaggers << endl;
		cout << "doPtEtaBin is    : " << doPtEtaBin << endl;
		cout << "doPseudoExp_ is   : " << doPseudoExp_ << endl; 
		cout << "doBiasStudy is   : " << doBiasStudy << endl; 
		cout << "doFBiasStudy is   : " << doFBiasStudy << endl; 
        cout << "nBtag    : " << nBtag << endl;
		cout << "doSCreweigh is   : " << doSCreweigh << endl;
		cout << "doTwoLights is   : " << doTwoLights << endl;
		cout << "useFit is        : " << useFit << endl;
		cout << "doVarBins is     : " << doVarBins << endl;
		cout << "doShift is       : " << doShift << endl;;
		cout << "do2D is          : " << do2D << endl;;
		cout << "do2Dcontrol is   : " << do2Dcontrol << endl;;
		cout << "Reweigh to data is   : " << reweigh_to_data << endl;;
		cout << "Reweigh left to all is   : " << reweigh_left_to_all << endl;
        cout << "doCorrectForBias is: " << doCorrectForBias << endl;
        cout << "useTTJetsExcl is : " << useTTJetsExcl_ << endl;

	}
	
	TRandom3 *rndm = new TRandom3();
	rndm->SetSeed(0);
	double random = 1;
	
	TFile* fout;
	fout = new TFile(outRootPath[0],"RECREATE");
	if(verbosity_>1) cout << "+--> Created output root file: " << outRootPath[0] << endl; 
	fout->cd();
	
	TFile* foutPseudoExp;
	foutPseudoExp = new TFile(outRootPath[1],"RECREATE");
	if(verbosity_>1) cout << "+--> Created output root file for Pseudo Exp: " << outRootPath[1] << endl; 
	foutPseudoExp->cd();
			
	//define storage vectors
	
	vector<string> v_dataSetName;
	
	vector<float> v_weight;
	vector<float> v_weight_nonrew;
    
    vector<int> v_npv;
    
    vector<int> v_runId;

    vector<float> v_ptweight;
    vector<float> v_ptweightControl;
    vector<float> v_ptweightControl2;

	vector<float> v_scaleFactor;
	vector<float> v_njets;
	vector<double> v_met;
	vector<double> v_mvaTrigId;
	vector<double*> v_bTag;
	vector<double*> v_bTagCS1;
	vector<double*> v_bTagCS2;
	vector<double> v_matchChiSquare;
	
	vector<double> v_pt;
	vector<double> v_eta;
	vector<double> v_var;
	vector<double> v_varb;
	vector<double> v_partonFlavour;
    
    vector<double> v_ptMuon;
    vector<double> v_etaMuon;

	
	vector<double> v_ptControl;
	vector<double> v_etaControl;
	vector<double> v_varControl;
	vector<double> v_partonFlavourControl;
	vector<double> v_bTagControl;
	
	vector<double> v_ptControl2;
	vector<double> v_etaControl2;
	vector<double> v_varControl2;
	vector<double> v_partonFlavourControl2;
	vector<double> v_bTagControl2;
	
	vector<double> v_controlVar;
	vector<double> v_controlVar2;
	
	vector<bool*> v_lepb_is;
	vector<bool*> v_q1_is;
	vector<bool*> v_q2_is;
	
	vector<bool> v_skipEvents;
	
	if(verbosity_>-1) cout << "+--> Reading TTrees and storing data" << endl; 
	
	int nCol=0;
	
	// for PDF unc stuff
	float sum_weights_ttbar = 0;
	float n_ttbar = 0;
	
	bool rescalednTTbarBeforeRefsel=false;
    
    TFile* ptrewfile;
    if (decay == 0)
        ptrewfile = new TFile("datamc-ptrew_mu.root","READ");
    else if (decay == 1)
        ptrewfile = new TFile("datamc-ptrew_el.root","READ");

    TFile* tempkle = new TFile("PVREW.root","READ");
    
    string fname = "AllBadHCALLaser.txt";
    string fnameS = "AllBadHCALLaser"+data_postfix+".txt";
    
    ifstream testF(fnameS.c_str());
    
    if (!testF) {

        cout << "BadHCALLaserEvents:: using large file "+fname+" to remove bad events, creating subset "+fnameS << endl;

        testF.close();
            
    } else {
     
        fname = fnameS;
        
        cout << "BadHCALLaserEvents:: using smaller file "+fname+" to remove bad events" << endl;
        
    }
    
    ifstream badev(fname.c_str());

    vector<string> bad_event;
    while (!badev.eof()) {
        
        string ev_;
        
        badev >> ev_;
        
        bad_event.push_back(ev_);
        
        //cout << ev_ << endl;
        
    }
    badev.close();
    
    //cout << bad_event.size() << endl;

    
	for (int iRootFile=0; iRootFile<nrFiles; iRootFile++){
				
		//exit(1);
		
		bool doSkip=true;
		for(int skip=0; skip<nRunSamples_; skip++){
			if(iRootFile==runSamples_[skip]) doSkip=false;
		}
		if(doSkip) continue;
		
		TString inRootPath;
		
		TString PDFWeights;
		
		
		if(!doJESchange) {
			inRootPath = inDir_;
			inRootPath += inRootFile[iRootFile];			
			PDFWeights =  inRootPath+"_PDFWeights.txt";
			inRootPath += ".root";
		} else if(doJESchange) {
			inRootPath = inDir_;
			inRootPath += inRootFile[iRootFile]; 
			inRootPath+="JES_";
			if((JESfactor-1)>0) inRootPath+="+";
			inRootPath+=(int) round((JESfactor-1)*100);
			inRootPath+="%.root";
		}
		
		if(verbosity_>0) cout << "+> reading file: " << inRootPath << endl;
		TFile *f = new TFile(inRootPath);
		
		TRootNTuple *NTuple;
		TTree *tree_ = (TTree*)f->Get("tree");
		NTuple = new TRootNTuple();
		
		TBranch *branch = tree_->GetBranch("TheNTuple");
		branch->SetAddress(&NTuple);
		
		Int_t nEvent = tree_->GetEntries();
		if(verbosity_>3) cout << "+----> nEvent " << nEvent << endl;
        
                                
		//////////////////////////////
		// read PDF stuff if needed //
		//////////////////////////////
		
		
		vector< vector < float > > pdfInfo;

		if (nPDFWeight != -1) {
			if(verbosity_>0) cout << "+> reading PDF weights file: " << PDFWeights << endl;
			
			for(unsigned int iEvt=0; iEvt<nEvent; iEvt++) {
				pdfInfo.push_back(vector<float>(47,0));
			}
			
			fstream file_to_read(PDFWeights, ios::in);
			
			int p=0;
			while (!file_to_read.eof()) {
				
				file_to_read >> pdfInfo[p][0] >> pdfInfo[p][1] >> pdfInfo[p][2] >> pdfInfo[p][3] >> pdfInfo[p][4] >> pdfInfo[p][5] >> pdfInfo[p][6] >> pdfInfo[p][7] >> pdfInfo[p][8] >> pdfInfo[p][9] >> pdfInfo[p][10] >> pdfInfo[p][11] >> pdfInfo[p][12] >> pdfInfo[p][13] >> pdfInfo[p][14] >> pdfInfo[p][15] >> pdfInfo[p][16] >> pdfInfo[p][17] >> pdfInfo[p][18] >> pdfInfo[p][19] >> pdfInfo[p][20] >> pdfInfo[p][21] >> pdfInfo[p][22] >> pdfInfo[p][23] >> pdfInfo[p][24] >> pdfInfo[p][25] >> pdfInfo[p][26] >> pdfInfo[p][27] >> pdfInfo[p][28] >> pdfInfo[p][29] >> pdfInfo[p][30] >> pdfInfo[p][31] >> pdfInfo[p][32] >> pdfInfo[p][33] >> pdfInfo[p][34] >> pdfInfo[p][35] >> pdfInfo[p][36] >> pdfInfo[p][37] >> pdfInfo[p][38] >> pdfInfo[p][39] >> pdfInfo[p][40] >> pdfInfo[p][41] >> pdfInfo[p][42] >> pdfInfo[p][43] >> pdfInfo[p][44] >> pdfInfo[p][45] >> pdfInfo[p][46];
				
				//if (p%1000 == 0)
				//	cout << pdfInfo[p][0] << " " << p << endl;
				p++;
			}
			
			file_to_read.close();
			
			if (((string)PDFWeights).find("TTbar") != string::npos) {
				
				for(unsigned int iEvt=0; iEvt<nEvent; iEvt++) {
					//cout << pdfInfo[iEvt][nPDFWeight+2] << " " << pdfInfo[iEvt][2] << " " << pdfInfo[iEvt][nPDFWeight+2]/pdfInfo[iEvt][2] << endl;
					sum_weights_ttbar+=(pdfInfo[iEvt][nPDFWeight+2]/pdfInfo[iEvt][2]);
				}
				
				n_ttbar+=nEvent;
				
				if(verbosity_>0) cout << " ++> summed ttbar event weights -> " << sum_weights_ttbar << " (so far total of " << n_ttbar << " ttbar events read)" << endl;
				
				//exit(1);
			}			
		}
    
		for(int ev=0; ev < nEvent; ev++){  
            
			float weight = 0;
			float weightfactor=1;
            
			tree_->GetEvent(ev);
            
            //if (fabs(NTuple->partonFlavour()) == 5 && NTuple->mlj() > 170) continue;
            //if (fabs(NTuple->partonFlavour()) != 5) continue;
            
            //if (fabs(NTuple->eta()) < 0.5 && NTuple->mlj() > 140) continue;
            
            //if (decay == 1)
            //    if (NTuple->MET() < 40)
            //        continue;   
            
            //if (decay == 1 && NTuple->mvaTrigID() < 0.8) continue;
            //cout << "ok1" << endl;
            if (decay == 0 && NTuple->MET() <= 30) continue;
            if (decay == 1 && NTuple->MET() <= 40) continue;
            
            //if (NTuple->ptMuon() < 30) continue;
            
            //if (decay == 1 && NTuple->mvaTrigID() < 0.5) continue;
                        
            //cout << "ok1a" << endl;

            // get the muon+jets events
            if (decay==0 && !NTuple->semiMu()) continue;
            //get the electron+jets events
            if (decay==1 && !NTuple->semiEl()) continue;
            //cout << "ok1b" << endl;

            if (NTuple->dataSetName().find("InvIso") != string::npos)
                NTuple->setDataSetName("multijet");
            
            if (NTuple->dataSetName().find("WJets_4jets_") != string::npos)
                NTuple->setDataSetName("WJets_4jets");
            if (NTuple->dataSetName().find("WJets_3jets_") != string::npos)
                NTuple->setDataSetName("WJets_3jets");
            if (NTuple->dataSetName().find("WJets_2jets_") != string::npos)
                NTuple->setDataSetName("WJets_2jets");
            if (NTuple->dataSetName().find("WJets_1jets_") != string::npos)
                NTuple->setDataSetName("WJets_1jets");
            if (doOnlyMSPlot && NTuple->dataSetName().find("WJets_") != string::npos)
                NTuple->setDataSetName("W_"+NTuple->dataSetName());

            if (NTuple->dataSetName().find("ZJets_4jets_") != string::npos)
                NTuple->setDataSetName("ZJets_4jets");
            if (NTuple->dataSetName().find("ZJets_3jets_") != string::npos)
                NTuple->setDataSetName("ZJets_3jets");
            if (NTuple->dataSetName().find("ZJets_2jets_") != string::npos)
                NTuple->setDataSetName("ZJets_2jets");
            if (NTuple->dataSetName().find("ZJets_1jets_") != string::npos)
                NTuple->setDataSetName("ZJets_1jets");
            if (doOnlyMSPlot && NTuple->dataSetName().find("ZJets_") != string::npos)
                NTuple->setDataSetName("Z_"+NTuple->dataSetName());

            if (NTuple->dataSetName().find("Data") == 0 || NTuple->dataSetName().find("data") == 0 || NTuple->dataSetName().find("DATA") == 0) {
             
                NTuple->setDataSetName("Data");
                
                stringstream r; r<<NTuple->runId();
                stringstream t; t<<NTuple->lumiBlockId();
                stringstream u; u<<NTuple->eventId();
                
                string ev = r.str()+":"+t.str()+":"+u.str();
                
                //cout << ev << endl;
                
                int pos = std::find(bad_event.begin(),bad_event.end(),ev)-bad_event.begin();
                
                if (pos < bad_event.size()) {
                    cout << "BadHCALLaserEvents:: Skipping " << ev << " " << pos << endl;
                    if (fnameS != fname) {
                        ofstream fout;
                        fout.open(fnameS.c_str(),ios::out | ios::app);
                        fout << ev << "\n";
                        fout.close();
                    }
                    
                    continue;
                }
            }
            
            // split ttbar signal and ttbar other
            
            if (NTuple->dataSetName().find("TTbarJets") == 0) {
                if ( (decay == 0 && NTuple->truesemiMu()) || (decay == 1 && NTuple->truesemiEl()) ) {
                    NTuple->setDataSetName("TTbarJets_Signal");
                } else { 
                    NTuple->setDataSetName("TTbarJets_Other");
                }
            }
            
            //if (NTuple->ptMuon() <= 35) continue; // for TOP-11-003
              
            //if (NTuple->nPV() > 10) continue;
            //if (NTuple->nJets() != 4) continue;
            
            //cout << NTuple->scalefactor() << " " << NTuple->scalefactor3D() << endl;
            
            //NTuple->setScaleFactor(NTuple->scalefactor3D()); // enable 3D reweighing
            
            //NTuple->setScaleFactor(1); // disable PU reweighing
        
            //if (fabs(NTuple->partonFlavour()) == 5) NTuple->setScaleFactor(1);
            //cout << NTuple->scalefactor() << " ";
            //if (NTuple->dataSetName()!="Data" && NTuple->dataSetName()!="data")
            //    NTuple->setScaleFactor(NTuple->scalefactor()*getPtWeight(tempkle,"pV_weight",NTuple->nPV()));
            //else
              //  cout << NTuple->dataSetName() << " " << NTuple->scalefactor() << endl;
            //cout << NTuple->scalefactor() << " " << endl;
            
            //if (pvBinA != -1) {
            //    if (NTuple->nPV() <= pvBinA || NTuple->nPV() > pvBinB) continue;
            //    NTuple->setScaleFactor(1);
                //cout << pvBinA << " " << pvBinB << endl; exit(1);
            //}
            
            //if (NTuple->pt() < 50 || NTuple->pt() > 80) continue;
			
			//if (NTuple->mlj() > centerleftlimit_ && NTuple->mlj() < rightlimit_ && NTuple->Btag(0) > 1.7) continue;
            
            /*if (NTuple->chiSq()<matchChiSquareCut_ && NTuple->mlj() < 20) {
                cout << "Weird event # " << ev << ", let's dump some info" << endl;
                cout << "mlb = " << NTuple->mlj() << endl;
                cout << "pt bcand = " << NTuple->pt()<< endl;
                cout << "pt Wcand1 = " << NTuple->ptControl() << endl;
                cout << "pt Wcand2 = " << NTuple->ptControl2() << endl;
                cout << "pt muon = " << NTuple->ptMuon() << endl;
                cout << endl;
                exit(1);
            }*/
            
			// define Datasets containers for MultiSamplePlots
            //cout << NTuple->dataSetName() << endl;
            
			if (doOnlyMSPlot) {
				//if (iRootFile == 7 || iRootFile == 8)
					//NTuple->setDataSetName("TT_bla");
				bool present=false;
                bool presentO = false;

				for (unsigned int d=0;d < datasets.size(); d++) {
					if (datasets[d]->Name() == NTuple->dataSetName())
						present=true;
                    if (datasets[d]->Name().find("_Other") != string::npos)
                        presentO=true;
                }
                
                if (NTuple->dataSetName().find("_Signal") != string::npos && !presentO) // first we need to enter other then signal
                    present=true;
                
				if (!present) {
                    
					datasets.push_back(new Dataset(NTuple->dataSetName(),datasets_title[nCol],false,datasets_color[nCol],1,1,1,-1));
                    
                    if (useTTJetsExcl_ && NTuple->dataSetName().find("TTbarJets") != string::npos)
                        datasets[datasets.size()-1]->SetEquivalentLuminosity(250000);
                    else
                        datasets[datasets.size()-1]->SetEquivalentLuminosity((1./NTuple->weight()));
                    
                    if (doOnlyMSPlot && decay == 0 && NTuple->dataSetName().find("multijet") != string::npos) {
                        
                        double lum=(1./NTuple->weight());
                        datasets[datasets.size()-1]->SetEquivalentLuminosity(lum/0.293183);

                        
                    }

                    
                    //continue;

					nCol++;
				}

			}

			if (nSelected_.find(NTuple->dataSetName()) == nSelected_.end()) {
				
				//cout << "POM " << NTuple->dataSetName() << endl;
				//nSelected_[NTuple->dataSetName()].push_back(3140); // data lumi
				//nSelected_[NTuple->dataSetName()].push_back(2140); // data lumi
				//nSelected_[NTuple->dataSetName()].push_back(4568.68); // data lumi
				//nSelected_[NTuple->dataSetName()].push_back(4917.61687); // NEW pixel data lumi
				
                /*if (decay == 0)
                    nSelected_[NTuple->dataSetName()].push_back(905.5); // 2012 lumi mu
                else 
                    nSelected_[NTuple->dataSetName()].push_back(891.415); // 2012 lumi el*/
                nSelected_[NTuple->dataSetName()].push_back(18000.0);
                
                nSelected_[NTuple->dataSetName()].push_back(1/NTuple->weight()); // mc lumi
				nSelected_[NTuple->dataSetName()].push_back(0); // after refsel
				nSelected_[NTuple->dataSetName()].push_back(0); // after refsel + btag TCHEM
				nSelected_[NTuple->dataSetName()].push_back(0); // after refsel (no SF)
				nSelected_[NTuple->dataSetName()].push_back(0); // after refsel + btag TCHEM (no SF)
				
			}
			
			datasetName = NTuple->dataSetName();
            
            // temp reset of data lumi to pixelLumiCalc value
            //if (NTuple->dataSetName()=="Data")
            //    NTuple->setWeight(1/4917.66387);
            
			origLumi = (1./NTuple->weight());
			
			//origLumi = origLumi+(origLumi*0.06);
			
            //cout << "ok2" << endl;

			if(verbosity_>0 || doPseudoExp_){
				if(ev % 10000 == 0)  {
					cout << "+>>> (Sample = " << datasetName << ", L = " << desiredIntLum_ << ", origL = " <<  origLumi << ") analyzing event " << ev << " of " << nEvent <<  flush << "\r" ;//<< endl;
				}
			}    
			
            //NTuple->setScaleFactor(1.);
            
            //cout << NTuple->scalefactor() << endl;
            
			if(doPseudoExp_){ 
				random = rndm->Uniform(1); 
				//if(ev<10) cout << "random "  << random << endl;  
			}

            //WJets_scale=1.;
			if (doPseudoExp_)
				weight = (1/origLumi)*desiredIntLum_*weightfactor;
			else if (NTuple->dataSetName()=="WJets")
				weight = (1/origLumi)*desiredIntLum_*weightfactor*NTuple->scalefactor()*WJets_scale;
			else if (NTuple->dataSetName()=="ZJets")
				weight = (1/origLumi)*desiredIntLum_*weightfactor*NTuple->scalefactor()*ZJets_scale;
			else
				weight = (1/origLumi)*desiredIntLum_*weightfactor*NTuple->scalefactor();
            
            //if (doOnlyMSPlot && NTuple->dataSetName().find("WJets") == 0) {
                //cout << weight << " " << NTuple->scalefactor() <<  endl;
                //weight = weight * (172.2/157.5);
                //cout << weight << endl << endl;
                //exit(1);
            //}
            
            /*if (NTuple->dataSetName().find("Data") != 0 && NTuple->dataSetName().find("data") != 0 && NTuple->dataSetName().find("DATA") != 0) {
                if (decay==0) {
                    weight=weight*getIsoMu20Weight(NTuple->ptMuon(),NTuple->nPV())*1.00*0.999; // Muon trig leg , jet leg, muon ID
                    //weight_nonrew=weight_nonrew*getIsoMu20Weight(NTuple->ptMuon(),NTuple->nPV())*1.00*0.999; // Muon trig leg , jet leg, muon ID
                    //weight=weight-0.032;
                } else if (decay==1) {
                    //weight=weight*0.988*1.000*0.998;
                    weight=weight*0.988*1.000*1.005;
                    //weight_nonrew=weight_nonrew*0.987*1.000*0.998;
                    
                }
            }*/
            
            //cout << NTuple->dataSetName() << " " << weight << " " << NTuple->scalefactor() <<  endl;

            //cout << weightfactor << endl; exit(1);
            
            //cout << ptWeight << endl;
                    
			//weight = NTuple->weight()*desiredIntLum_;
			
			//weight = weight;
			
			if (nPDFWeight != -1) {
				
				if ((pdfInfo[ev][1]-NTuple->mlj())/NTuple->mlj() > 0.00001) {
					cout << pdfInfo[ev][0] << " " << ev << endl;
					cout << pdfInfo[ev][1] << " " << NTuple->mlj() << " " << (pdfInfo[ev][1]-NTuple->mlj())/NTuple->mlj() << endl;
					cout << "ERROR inconsistency between event info from TTree and the info from the PDF uncertainties txt file" << endl;
					exit(1);
				}
				
				if (datasetName.find("TTbar") != string::npos) {
			
					float pdfweight = ((pdfInfo[ev][nPDFWeight+2]/pdfInfo[ev][2])/sum_weights_ttbar)*n_ttbar;
					
					weight = weight*pdfweight;
					
					//weight = weight*(pdfInfo[ev][nPDFWeight+2]/pdfInfo[ev][2]);
					
				} else {
					
					weight = weight*(pdfInfo[ev][nPDFWeight+2]/pdfInfo[ev][2]);
					
				}
				
			}
            
            float weight_nonrew=weight;
            
            //if (ptWeight != -999)
            //    weight = weight*ptWeight;
            
            if (NTuple->dataSetName().find("Data") == 0 || NTuple->dataSetName().find("data") == 0 || NTuple->dataSetName().find("DATA") == 0) {
				weight=1;
                weight_nonrew=1;
				//origLumi = origLumi+(origLumi*0.06);
				//origLumi = origLumi-(origLumi*0.06);
                
                //if (NTuple->mlj() > leftlimit_ && NTuple->mlj() < centerleftlimit_) {
                //    weight = getPtWeight(ptrewfile,"pT_weight_enrdata_alldata",NTuple->pt());
                //    if (weight == -999.) weight =1;
				//}
                
				desiredIntLum_ = origLumi;
				dataLum_ = 	desiredIntLum_;
                                
			}

            
            // reweigh InvIso sample for Data-Driven QCD (for now only e+jets is implemented)
            
                // reweigh for the trigger bins and to the MC yield
			
            if (datasetName.find("multijet") != string::npos) {
                
                weight=1;
                                
                /*if (NTuple->runId() <= 165633) // noniso triggers
                    weight=0.98*0.126; // trigger bin scale factor * mc yield factor
                else if (NTuple->runId() >= 165970 && NTuple->runId() < 175860 ) // iso trigger
                    weight=3.73*0.126;
                else if (NTuple->runId() >= 175860) // 2011B -> high PU
                    weight=3.38*0.126;  */
                
                /*if (doOnlyMSPlot) {
                    NTuple->setScaleFactor(weight);
                    weight=1;
                } else {
                    weight=0.2;
                    //if (decay == 0) weight=0.334;
                    //if (decay == 1) weight=2.524;
                    weight=weight*(1/origLumi)*desiredIntLum_;
                }*/
                
                //cout << weight << endl;
                
                //weight=weight*(1/origLumi)*desiredIntLum_;
                
            }
            
			//cout << "still here?" << endl;
			
			// CHECKING WHICH CHI2 CUT TO TAKE
			if (!doPseudoExp_ && doChi2Check) {
				fout->cd();
				
				string title = "MTopDistribution_NOChiSqCut";
				
				if (chi2CutHistos.find(title) == chi2CutHistos.end()) {
					chi2CutHistos[title] = new TH1D(title.c_str(),title.c_str(),250,0,500);
				}
				chi2CutHistos[title]->Fill(NTuple->HadTopCandMass(),weight);
				
				title = "Chi2";
				
				if (chi2CutHistos.find(title) == chi2CutHistos.end()) {
					chi2CutHistos[title] = new TH1D(title.c_str(),title.c_str(),400,0,4000);
				}
				chi2CutHistos[title]->Fill(NTuple->chiSq(),weight);
				
				
				for (int i=0; i<sizeof(chisqCutForMass)/sizeof(chisqCutForMass[0]); i++) {
					
					stringstream c; c << chisqCutForMass[i];
					string title = "MTopDistribution_ChiSqCut_"+c.str();
					
					if (chi2CutHistos.find(title) == chi2CutHistos.end()) {
						chi2CutHistos[title] = new TH1D(title.c_str(),title.c_str(),250,0,500);
					}
					
					if (NTuple->chiSq() < chisqCutForMass[i])
						chi2CutHistos[title]->Fill(NTuple->HadTopCandMass(),weight);
					
					title = "Chi2_cut_"+c.str();
					
					if (chi2CutHistos.find(title) == chi2CutHistos.end()) {
						chi2CutHistos[title] = new TH1D(title.c_str(),title.c_str(),400,0,4000);
					}
					
					if (NTuple->chiSq() < chisqCutForMass[i])
						chi2CutHistos[title]->Fill(NTuple->chiSq(),weight);
					
					//cout << chi2CutHistos[title]->GetEntries() << " " << chi2CutHistos[title]->GetBinContent(chi2CutHistos[title]->GetMaximumBin()) <<  endl;
					//cout << chi2CutHistos[title]->GetTitle() << endl;
				}
				
			}	


			if (NTuple->dataSetName().find("TTbar") == 0) {
				
				nTTbarBeforeChiSq+=weight_nonrew;	
			
				nTTbarAfterRefSel+=weight_nonrew;
                								
			}
			
			nSelected_[NTuple->dataSetName()][2]+=NTuple->scalefactor();
			nSelected_[NTuple->dataSetName()][4]++;
			
			if (NTuple->Btag(4) > 0.679) {
				nSelected_[NTuple->dataSetName()][3]+=NTuple->scalefactor();
				nSelected_[NTuple->dataSetName()][5]++;
			}
			
			if(!doOnlyMSPlot && !doFBiasStudy && NTuple->chiSq()>matchChiSquareCut_) continue; // TURN ME ON
            
            //if(NTuple->chiSq()>matchChiSquareCut_) continue;

			if (NTuple->dataSetName().find("TTbarJets") == 0 && !rescalednTTbarBeforeRefsel) {
				nTTbarBeforeRefSel = nTTbarBeforeRefSel*(1/origLumi)*desiredIntLum_;
				rescalednTTbarBeforeRefsel=true;
				//cout << nTTbarBeforeRefSel << endl;
				//exit(1);
			}

			if (NTuple->dataSetName().find("TTbarJets_Semi") == 0) {
				if (NTuple->partonFlavour() == 5)
					nGoodChi2Combi++;
				nChi2Combi++;
			}
			
			if (NTuple->dataSetName().find("TTbar") == 0) {
				
				nTTbarAfterChiSq+=weight_nonrew;		
				
				nTTbarBeforeMLBCUT+=weight_nonrew;
				
			}
			
			//extra cut on mlb to remove overflow bin from XS templates
			
			//if (!doOnlyMSPlot && (NTuple->mlj() > upRangeVar0 || NTuple->m3() > upRangeVar0)) continue;
			
            if (!doOnlyMSPlot) {
            
                //if (fitMode == 0 && NTuple->mlj() > upRangeVar0) continue; // TURN ME ON
                
                if (fitMode == 1 && (NTuple->m3() > 500 || NTuple->m3() < 100) ) continue;
                
                if (NTuple->mlj() > upRangeVar0) continue;
                
                //if (NTuple->mlj() > upRangeVar0 || NTuple->m3() > 800) continue;
                //else if (fitMode == 2 && (NTuple->mlj() > upRangeVar0 || NTuple->m3() > upRangeVar0)) continue;
                
            }
            
			if (NTuple->dataSetName().find("TTbar") == 0) {
								
				nTTbarAfterMLBCUT+=weight_nonrew;
				
			}
            
            v_weight_nonrew.push_back(weight_nonrew);
            // reweight MC pt distribution do data 
            
            if ( !doPseudoExp_ && !doOnlyMSPlot) {
                
                float ptw = 1;
                float ptwC1 = 1;
                float ptwC2 = 1;
                
                if (NTuple->dataSetName()=="Data") {
                    
                    reweigh_to_data=false;
                    
                }
                
                if (reweigh_left_to_all && !reweigh_to_data) {
                    
                    if (NTuple->dataSetName()=="Data") {
                        if (NTuple->mlj() > leftlimit_ && NTuple->mlj() < centerleftlimit_)
                            ptw=getPtWeight(ptrewfile,"pT_weight_enrdata_alldata",NTuple->pt());
                        
                        //else if (NTuple->mlj() > centerrightlimit_ && NTuple->mlj() < rightlimit_)
                        //    ptw=getPtWeight(ptrewfile,"pT_weight_deplmc_alldata",NTuple->pt());
                        
                        if (NTuple->mljControl() > leftlimit_ && NTuple->mljControl() < centerleftlimit_)
                            ptwC1=getPtWeight(ptrewfile,"pT_weight_CS_enrdata_alldata",NTuple->ptControl());
                        
                        //else if (NTuple->mljControl() > centerrightlimit_ && NTuple->mljControl() < rightlimit_)
                        //    ptwC1=getPtWeight(ptrewfile,"pT_weight_CS_deplmc_alldata",NTuple->ptControl());
                        
                        if (NTuple->mljControl2() > leftlimit_ && NTuple->mljControl2() < centerleftlimit_)
                            ptwC2=getPtWeight(ptrewfile,"pT_weight_CS_enrdata_alldata",NTuple->ptControl2());
                        
                        //else if (NTuple->mljControl2() > centerrightlimit_ && NTuple->mljControl2() < rightlimit_)
                        //    ptwC2=getPtWeight(ptrewfile,"pT_weight_CS_deplmc_alldata",NTuple->ptControl2());
                    
                    } else {
                        
                        if (NTuple->mlj() > leftlimit_ && NTuple->mlj() < centerleftlimit_)
                            ptw=getPtWeight(ptrewfile,"pT_weight_enrmc_alldata",NTuple->pt());
                        
                        //else if (NTuple->mlj() > centerrightlimit_ && NTuple->mlj() < rightlimit_)
                        //    ptw=getPtWeight(ptrewfile,"pT_weight_deplmc_alldata",NTuple->pt());
                        
                        if (NTuple->mljControl() > leftlimit_ && NTuple->mljControl() < centerleftlimit_)
                            ptwC1=getPtWeight(ptrewfile,"pT_weight_CS_enrmc_alldata",NTuple->ptControl());
                        
                        //else if (NTuple->mljControl() > centerrightlimit_ && NTuple->mljControl() < rightlimit_)
                        //    ptwC1=getPtWeight(ptrewfile,"pT_weight_CS_deplmc_alldata",NTuple->ptControl());
                        
                        if (NTuple->mljControl2() > leftlimit_ && NTuple->mljControl2() < centerleftlimit_)
                            ptwC2=getPtWeight(ptrewfile,"pT_weight_CS_enrmc_alldata",NTuple->ptControl2());
                        
                        //else if (NTuple->mljControl2() > centerrightlimit_ && NTuple->mljControl2() < rightlimit_)
                        //    ptwC2=getPtWeight(ptrewfile,"pT_weight_CS_deplmc_alldata",NTuple->ptControl2());
                    }
                }
                
                if (!reweigh_left_to_all && reweigh_to_data && NTuple->dataSetName()!="Data") {
                    //cout << "broooooooooooooooooooool" << endl; exit(1);
                    //if ( (NTuple->mlj() > leftlimit_ && NTuple->mlj() < centerleftlimit_) || (NTuple->mlj() > centerrightlimit_ && NTuple->mlj() < rightlimit_) )
                        ptw=getPtWeight(ptrewfile,"pT_weight_allmc_alldata",NTuple->pt());
                    
                    
                    /*if ( (NTuple->mljControl() > leftlimit_ && NTuple->mljControl() < centerleftlimit_) || (NTuple->mljControl() > centerrightlimit_ && NTuple->mljControl() < rightlimit_) )
                        ptwC1=getPtWeight(ptrewfile,"pT_weight_CS_allmc_alldata",NTuple->ptControl());
                    if ( (NTuple->mljControl2() > leftlimit_ && NTuple->mljControl2() < centerleftlimit_) || (NTuple->mljControl2() > centerrightlimit_ && NTuple->mljControl2() < rightlimit_) )                      
                        ptwC2=getPtWeight(ptrewfile,"pT_weight_CS_allmc_alldata",NTuple->ptControl2());
                */
                    ptwC1=ptw;
                    ptwC2=ptw;
                } 
                
                //cout << ptw << " " << ptwC1 << endl;
                
                //if (NTuple->mlj() > centerleftlimit_ && NTuple->mlj() < rightlimit_)
                //    ptWeight=getPtWeight(ptrewfile,"pT_weight_enrmc_alldata",NTuple->pt());

                v_ptweight.push_back(ptw);
                v_ptweightControl.push_back(ptwC1);
                v_ptweightControl2.push_back(ptwC2);

            } else {
                
                v_ptweight.push_back(1);
                v_ptweightControl.push_back(1);
                v_ptweightControl2.push_back(1);
            
            }

            // for M3 MLB correlation check
            
            /*string filename = "events_sel_data_channel"+data_postfix+".txt";
            
             //if (NTuple->Btag(7) > 0.545) { // JPM
             if (NTuple->Btag(4) > 0.679) { // CSVM
                
                ofstream ev(filename.c_str(), ios::out | ios::app);
            
                ev << NTuple->runId() << ":" << NTuple->eventId() << ":" << NTuple->lumiBlockId() << ":" << NTuple->mlj() << ":" << NTuple->m3() << endl;
            
                ev.close();
                
            }*/
            
			// now fill vector
			
			v_dataSetName.push_back(NTuple->dataSetName());
			
			v_weight.push_back(weight);
            
            v_runId.push_back(NTuple->runId());
            
            v_npv.push_back(NTuple->nPV());
            
            //v_npv.push_back(1);
            
			//v_weight.push_back(round(1000,weight));
			
			//cout << "*******************" << weight << " " << round(1000,weight) << "*******************" << endl; exit(1);
			
			if (doOnlyMSPlot || doPseudoExp_) {
                
                if (doOnlyMSPlot && useTTJetsExcl_ && NTuple->dataSetName().find("TTbarJets") == 0) {
                    //cout << 250000/(1./NTuple->weight()) << endl;
                    v_scaleFactor.push_back(NTuple->scalefactor()*(250000/(1./NTuple->weight())));
                }
                else
                    v_scaleFactor.push_back(NTuple->scalefactor());
            }
            
            if (doOnlyMSPlot) v_njets.push_back(NTuple->nJets());
            if (doOnlyMSPlot) v_met.push_back(NTuple->MET());
            
            //cout << NTuple->MET() << endl;
            
            if (doOnlyMSPlot) v_mvaTrigId.push_back(NTuple->mvaTrigID());
			
            //cout << NTuple->MET() << endl;
            
			double* bTag = new double[11];
            bTag[0] = NTuple->Btag(0); // TCHE -10 30
            bTag[1] = NTuple->Btag(1); // TCHP -10 30
            bTag[2] = NTuple->Btag(7); // JP
            bTag[3] = NTuple->Btag(6); // JBP
            bTag[4] = NTuple->Btag(2); // SSVHE 0 8
            bTag[5] = NTuple->Btag(3); // SSVHP 0 8
            bTag[6] = NTuple->Btag(4); // CSV
            bTag[7] = NTuple->Btag(8); // CSV retrained
            bTag[8] = NTuple->Btag(0);
            bTag[9] = NTuple->Btag(0);
            bTag[10] = NTuple->Btag(0);
            
            double* bTagC1= new double[11];
            bTagC1[0] = NTuple->BtagCS1(0); // TCHE -10 30
            bTagC1[1] = NTuple->BtagCS1(1); // TCHP -10 30
            bTagC1[2] = NTuple->BtagCS1(7); // JP
            bTagC1[3] = NTuple->BtagCS1(6); // JBP
            bTagC1[4] = NTuple->BtagCS1(2); // SSVHE 0 8
            bTagC1[5] = NTuple->BtagCS1(3); // SSVHP 0 8
            bTagC1[6] = NTuple->BtagCS1(4); // CSV
            bTagC1[7] = NTuple->BtagCS1(8); // CSV retrained
            bTagC1[8] = NTuple->BtagCS1(0);
            bTagC1[9] = NTuple->BtagCS1(0);
            bTagC1[10] = NTuple->BtagCS1(0);
            
            double* bTagC2= new double[11];
            bTagC2[0] = NTuple->BtagCS2(0); // TCHE -10 30
            bTagC2[1] = NTuple->BtagCS2(1); // TCHP -10 30
            bTagC2[2] = NTuple->BtagCS2(7); // JP
            bTagC2[3] = NTuple->BtagCS2(6); // JBP
            bTagC2[4] = NTuple->BtagCS2(2); // SSVHE 0 8
            bTagC2[5] = NTuple->BtagCS2(3); // SSVHP 0 8
            bTagC2[6] = NTuple->BtagCS2(4); // CSV
            bTagC2[7] = NTuple->BtagCS2(8); // CSV retrained
            bTagC2[8] = NTuple->BtagCS2(0);
            bTagC2[9] = NTuple->BtagCS2(0);
            bTagC2[10] = NTuple->BtagCS2(0);
            
            if (doPseudoExp_ && nBtag != -1) {
                bTag[0] = bTag[nBtag];
                bTagC1[0] = bTag[nBtag];
                bTagC2[0] = bTag[nBtag];
            }
            
			v_bTag.push_back(bTag);
			v_bTagCS1.push_back(bTagC1);
			v_bTagCS2.push_back(bTagC2);
			
			v_partonFlavour.push_back(NTuple->partonFlavour());
			v_matchChiSquare.push_back(NTuple->chiSq());
			
			v_var.push_back(NTuple->mlj());
			
            v_varb.push_back(NTuple->m3());
            
			v_pt.push_back(NTuple->pt());
			v_eta.push_back(fabs(NTuple->eta()));
            
            if (doOnlyMSPlot) {
                v_ptMuon.push_back(NTuple->ptMuon());
                v_etaMuon.push_back(NTuple->etaMuon());
                
            }

			v_etaControl.push_back(fabs(NTuple->etaControl()));
			v_ptControl.push_back(NTuple->ptControl());
			v_partonFlavourControl.push_back(NTuple->partonFlavourControl());
			v_bTagControl.push_back(fabs(NTuple->bTagControl()));

			v_etaControl2.push_back(fabs(NTuple->etaControl2()));
			v_ptControl2.push_back(NTuple->ptControl2());
			v_partonFlavourControl2.push_back(NTuple->partonFlavourControl2());
			v_bTagControl2.push_back(fabs(NTuple->bTagControl2()));

			v_controlVar.push_back(NTuple->mljControl());
			v_controlVar2.push_back(NTuple->mljControl2());
			
			bool* lepb_is = new bool[4];
			bool* q1_is = new bool[4];
			bool* q2_is = new bool[4];
			
			lepb_is[0] = NTuple->jet_is_b();	
			lepb_is[1] = NTuple->jet_is_nonb();	
			lepb_is[2] = NTuple->jet_is_hadqq();
			lepb_is[3] = NTuple->jet_is_radq();
			q1_is[0] = NTuple->jetControl_is_b();
			q1_is[1] = NTuple->jetControl_is_nonb();
			q1_is[2] = NTuple->jetControl_is_hadqq();
			q1_is[3] = NTuple->jetControl_is_radq();
			q2_is[0] = NTuple->jetControl2_is_b();
			q2_is[1] = NTuple->jetControl2_is_nonb();
			q2_is[2] = NTuple->jetControl2_is_hadqq();
			q2_is[3] = NTuple->jetControl2_is_radq();
			
			v_lepb_is.push_back(lepb_is);
			v_q1_is.push_back(q1_is);
			v_q2_is.push_back(q2_is);
			
		} // end of event-loop
		//cout << endl;
		f->Close();
	} // end of root-file loop
	//cout << v_weight.size() << " " << v_ptweight.size() << endl; exit(1);

	if (doOnlyMSPlot) {
		
		vector<TH2F*> corr_histos;
		
		double nData = 0;
		double nDataL = 0;
		double nDataR = 0;
		
		double nB = 0;
		double nBL = 0;
		double nBR = 0;
		
		double nNonB= 0;
		double nNonBL= 0;
		double nNonBR= 0;
		
		double nBCSA = 0;
		double nBCS = 0;
		double nBLCS = 0;
		double nBRCS = 0;
		
		double nNonBCSA = 0;
		double nNonBCS= 0;
		double nNonBLCS= 0;
		double nNonBRCS= 0;
		
		double all_b=0;
		double all_c=0;
		double all_dusg=0;
		double tag_b=0;
		double tag_c=0;
		double tag_dusg=0;

		vector<Dataset*> datasets_MC;
		for (unsigned int p=0;p < datasets.size(); p++)
			if (datasets[p]->Name().find("data") == string::npos && datasets[p]->Name().find("Data") == string::npos && datasets[p]->Name().find("DATA") == string::npos)
				datasets_MC.push_back(datasets[p]);
		
		
		cout << "+> FINISHED EVENT LOOP, DUMPING RECORDED DATASETS" << endl << endl; 
		for (unsigned int d=0;d < datasets.size(); d++)
			cout << "\tName: " << datasets[d]->Name() << " XS: " << datasets[d]->Xsection() << " EqLumi: " << datasets[d]->EquivalentLumi() << endl;

		cout << "+> Defining MultiSamplePlots" << endl;
		
		stringstream lum; lum<<round(1,dataLum_);
		//datasets[datasets.size()-1]->SetTitle("Data (L="+lum.str()+"/pb)"); //we allways put data as last dataset
		datasets[datasets.size()-1]->SetTitle("Data"); //we allways put data as last dataset
		
		std::map<string,MultiSamplePlot*> MSPlot;
		std::map<string,TH2D*> profile_PU;
		std::map<string,TH1D*> hist_PU;
		
		MSPlot["nPV"] = new MultiSamplePlot(datasets, "nPV", 36, -0.5, 35.5, "#PV");
		MSPlot["nPV_nonRew"] = new MultiSamplePlot(datasets, "nPV_nonRew", 36, -0.5, 35.5, "#PV");
        
        if (decay==0) {
            MSPlot["pT_lepton"] = new MultiSamplePlot(datasets, "pT_lepton", 100, 0, 500, "p_{T}^{#mu}");
            MSPlot["pT_lepton_btag"] = new MultiSamplePlot(datasets, "pT_lepton_btag", 100, 0, 500, "p_{T}^{#mu}");
            
            MSPlot["Eta_lepton"] = new MultiSamplePlot(datasets, "Eta_lepton", 50, -5, 5, "#eta^{#mu}");
            MSPlot["Eta_lepton_btag"] = new MultiSamplePlot(datasets, "Eta_lepton_btag", 50, -5, 5, "#eta^{#mu}");
        } else if (decay==1) {
            MSPlot["pT_lepton"] = new MultiSamplePlot(datasets, "pT_lepton", 100, 0, 500, "p_{T}^{e}");
            MSPlot["pT_lepton_btag"] = new MultiSamplePlot(datasets, "pT_lepton_btag", 100, 0, 500, "p_{T}^{e}");
            
            MSPlot["Eta_lepton"] = new MultiSamplePlot(datasets, "Eta_lepton", 50, -5, 5, "#eta^{e}");
            MSPlot["Eta_lepton_btag"] = new MultiSamplePlot(datasets, "Eta_lepton_btag", 50, -5, 5, "#eta^{e}");
        }
        
        MSPlot["pT_lepbCand"] = new MultiSamplePlot(datasets, "pT_lepbCand", 100, 0, 500, "p_{T}^{b cand} (GeV)");
        MSPlot["pT_lepbCand_btag"] = new MultiSamplePlot(datasets, "pT_lepbCand_btag", 100, 0, 500, "p_{T}^{b cand} (GeV)");
        MSPlot["pT_WjetCand1"] = new MultiSamplePlot(datasets, "pT_WjetCand1", 100, 0, 500, "p_{T}^{lightjetcand1} (GeV)");
        MSPlot["pT_WjetCand1_btag"] = new MultiSamplePlot(datasets, "pT_WjetCand1_btag", 100, 0, 500, "p_{T}^{lightjetcand1} (GeV)");
        MSPlot["pT_WjetCand2"] = new MultiSamplePlot(datasets, "pT_WjetCand2", 100, 0, 500, "p_{T}^{lightjetcand2} (GeV)");
        MSPlot["pT_WjetCand2_btag"] = new MultiSamplePlot(datasets, "pT_WjetCand2_btag", 100, 0, 500, "p_{T}^{lightjetcand2} (GeV)");
        MSPlot["Eta_lepbCand"] = new MultiSamplePlot(datasets, "Eta_lepbCand", 50, -5, 5, "#eta^{b cand} (GeV)");
        MSPlot["Eta_lepbCand_btag"] = new MultiSamplePlot(datasets, "Eta_lepbCand_btag", 50, -5, 5, "#eta^{b cand} (GeV)");
        MSPlot["Eta_WjetCand1"] = new MultiSamplePlot(datasets, "Eta_WjetCand1", 50, -5, 5, "#eta^{lightjetcand1} (GeV)");
        MSPlot["Eta_WjetCand1_btag"] = new MultiSamplePlot(datasets, "Eta_WjetCand1_btag", 50, -5, 5, "#eta^{lightjetcand1} (GeV)");
        MSPlot["Eta_WjetCand2"] = new MultiSamplePlot(datasets, "Eta_WjetCand2", 50, -5, 5, "#eta^{lightjetcand2} (GeV)");
        MSPlot["Eta_WjetCand2_btag"] = new MultiSamplePlot(datasets, "Eta_WjetCand2_btag", 50, -5, 5, "#eta^{lightjetcand2} (GeV)");

        
        MSPlot["BestJetCombChi2"] = new MultiSamplePlot(datasets, "BestJetCombChi2", 25, 0, 1000, "#chi^{2}");
        
		MSPlot["selectedEventsJetsEta"] = new MultiSamplePlot(datasets, "selectedEventsJetsEta", 13, -2.6, 2.6, "Jet #eta");
        
		MSPlot["nJets"] = new MultiSamplePlot(datasets, "nJets", 8, 2.5, 10.5, "#selected jets");
		MSPlot["nJets_btag"] = new MultiSamplePlot(datasets, "nJets_btag", 8, 2.5, 10.5, "#selected jets");
        
        if (decay==0) {
            MSPlot["MLB_BTV"] = new MultiSamplePlot(datasets, "MLB_BTV", 25, 0, 500, "M_{#mub} (GeV)","Events / 20 GeV");
            MSPlot["MLB_BTV_btag"] = new MultiSamplePlot(datasets, "MLB_BTV_btag", 25, 0, 500, "M_{#mub} (GeV)","Events / 20 GeV");
            MSPlot["MLB_FourJets"] = new MultiSamplePlot(datasets, "MLB_FourJets", 50, 0, 500, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB_FourJets_btag"] = new MultiSamplePlot(datasets, "MLB_FourJets_btag", 50, 0, 500, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB"] = new MultiSamplePlot(datasets, "MLB", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB_btag"] = new MultiSamplePlot(datasets, "MLB_btag", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB_MET"] = new MultiSamplePlot(datasets, "MLB_MET", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB_btag_MET"] = new MultiSamplePlot(datasets, "MLB_btag_MET", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB-bquarks"] = new MultiSamplePlot(datasets_MC, "MLB-bquarks", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB-nonbquarks"] = new MultiSamplePlot(datasets_MC, "MLB-nonbquarks", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB-bquarks_btag"] = new MultiSamplePlot(datasets_MC, "MLB-bquarks_btag", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB-nonbquarks_btag"] = new MultiSamplePlot(datasets_MC, "MLB-nonbquarks_btag", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB_ControlSample"] = new MultiSamplePlot(datasets, "MLB_ControlSample", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB_ControlSample_q1"] = new MultiSamplePlot(datasets, "MLB_ControlSample_q1", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");
            MSPlot["MLB_ControlSample_q2"] = new MultiSamplePlot(datasets, "MLB_ControlSample_q2", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");

        } else {
            MSPlot["MLB_BTV"] = new MultiSamplePlot(datasets, "MLB_BTV", 25, 0, 500, "M_{eb} (GeV)","Events / 20 GeV");
            MSPlot["MLB_BTV_btag"] = new MultiSamplePlot(datasets, "MLB_BTV_btag", 25, 0, 500, "M_{eb} (GeV)","Events / 20 GeV");
            MSPlot["MLB_FourJets"] = new MultiSamplePlot(datasets, "MLB_FourJets", 50, 0, 500, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB_FourJets_btag"] = new MultiSamplePlot(datasets, "MLB_FourJets_btag", 50, 0, 500, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB"] = new MultiSamplePlot(datasets, "MLB", 75, 0, 1200, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB_btag"] = new MultiSamplePlot(datasets, "MLB_btag", 75, 0, 1200, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB_MET"] = new MultiSamplePlot(datasets, "MLB_MET", 75, 0, 1200, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB_btag_MET"] = new MultiSamplePlot(datasets, "MLB_btag_MET", 75, 0, 1200, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB-bquarks"] = new MultiSamplePlot(datasets_MC, "MLB-bquarks", 75, 0, 1200, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB-nonbquarks"] = new MultiSamplePlot(datasets_MC, "MLB-nonbquarks", 75, 0, 1200, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB-bquarks_btag"] = new MultiSamplePlot(datasets_MC, "MLB-bquarks_btag", 75, 0, 1200, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB-nonbquarks_btag"] = new MultiSamplePlot(datasets_MC, "MLB-nonbquarks_btag", 75, 0, 1200, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB_ControlSample"] = new MultiSamplePlot(datasets, "MLB_ControlSample", 75, 0, 1200, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB_ControlSample_q1"] = new MultiSamplePlot(datasets, "MLB_ControlSample_q1", 75, 0, 1200, "M_{eb} (GeV)","Events / 10 GeV");
            MSPlot["MLB_ControlSample_q2"] = new MultiSamplePlot(datasets, "MLB_ControlSample_q2", 75, 0, 1200, "M_{eb} (GeV)","Events / 10 GeV");
        }
        
        MSPlot["M3_BTV"] = new MultiSamplePlot(datasets, "M3_BTV", 25, 0, 500, "m3 (GeV)","Events / 20 GeV");
        MSPlot["M3_BTV_btag"] = new MultiSamplePlot(datasets, "M3_BTV_btag", 25, 0, 500, "m3 (GeV)","Events / 20 GeV");
        
		MSPlot["MLB2"] = new MultiSamplePlot(datasets, "MLB2", 75, 0, 1200, "M_{#mub} (GeV)","Events / 10 GeV");
		MSPlot["MLBL"] = new MultiSamplePlot(datasets, "MLBL", 10 , leftlimit_, centerleftlimit_, "M_{#mub} (GeV)");
		MSPlot["MLBR"] = new MultiSamplePlot(datasets, "MLBR", 10, centerrightlimit_, rightlimit_, "M_{#mub} (GeV)");
		        
        MSPlot["PUWeight"] = new MultiSamplePlot(datasets, "PUWeight", 50, 0, 5, "PUWeight");

        MSPlot["MET"] = new MultiSamplePlot(datasets, "MET", 100, 0, 300, "MET");
        MSPlot["MET_btag"] = new MultiSamplePlot(datasets, "MET_btag", 100, 0, 300, "MET_btag");

        MSPlot["MVATrigID"] = new MultiSamplePlot(datasets, "MVATrigID", 15, 0, 1.5, "MVATrigID");
        
        TH2D* TCHEtrueb_vs_PUWeight= new TH2D("TCHEtrueb_vs_PUWeight","TCHE_vs_PUWeight;Disc;PU weight",50,-10,30,50,0,5);
        
        TH2D* TCHEtrueb_vs_nPV= new TH2D("TCHEtrueb_vs_nPV","TCHEtrueb_vs_nPV;Disc;nPV",50,-10,30,36,-0.5,35);
        TH2D* TCHEtruenonb_vs_nPV= new TH2D("TCHEtruenonb_vs_nPV","TCHEtruenonb_vs_nPV;Disc;nPV",50,-10,30,36,-0.5,35);
        
        TH2D* PTtrueb_vs_nPV= new TH2D("PTtrueb_vs_nPV","PTtrueb_vs_nPV;pT;nPV",100,0,500,36,-0.5,35);
        TH2D* PTtruenonb_vs_nPV= new TH2D("PTtruenonb_vs_nPV","PTtruenonb_vs_nPV;pT;nPV",500,0,500,36,-0.5,35);
        
        TH2D* MLB_VS_PUWeight= new TH2D("MLB_VS_PUWeight","MLB_VS_PUWeight;mlb;PU weight",25,0,500,50,0,5);
        
		corr_histos.push_back(new TH2F("data_mlj_correlation_tche",";m_{#mu j} (GeV); TCHE discriminator",250,0,500,50,-10,30));
		corr_histos.push_back(new TH2F("mlj_correlation_tche",";m_{#mu j} (GeV); TCHE discriminator",250,0,500,50,-10,30));
		corr_histos.push_back(new TH2F("b_mlj_correlation_tche",";m_{#mu j} (GeV); TCHE discriminator",250,0,500,50,-10,30));
		corr_histos.push_back(new TH2F("nonb_mlj_correlation_tche",";m_{#mu j} (GeV); TCHE discriminator",250,0,500,50,-10,30));
		
		TH1D* MLB_ttbar = new TH1D("MLB_ttbar",";M_{#mub} (GeV);a.u.",75,0,1200);
		TH1D* MLB_bkg = new TH1D("MLB_bkg",";M_{#mub} (GeV);a.u.",75,0,1200);
		TH1D* MLB_wjets = new TH1D("MLB_wjets",";M_{#mub} (GeV);a.u.",75,0,1200);
		TH1D* MLB_wjets_noweight = new TH1D("MLB_wjets_noweight",";M_{#mub} (GeV);a.u.",75,0,1200);
		
		TH1D* MLB_MC = new TH1D("MLB_MC",";M_{#mub} (GeV);a.u.",25,0,500);
		TH1D* MLB_DATA = new TH1D("MLB_DATA",";M_{#mub} (GeV);a.u.",25,0,500);
		TH1D* CHI2_MC = new TH1D("CHI2_MC",";#chi^{2};a.u.",10,0,100);
		TH1D* CHI2_DATA = new TH1D("CHI2_DATA",";#chi^{2};a.u.",10,0,100);

		TH1D* MLBCS_MC = new TH1D("MLBCS_MC",";M_{#mub} (GeV);a.u.",25,0,500);
		TH1D* MLBCS_DATA = new TH1D("MLBCS_DATA",";M_{#mub} (GeV);a.u.",25,0,500);
		
		MLB_MC->Sumw2();
		MLBCS_MC->Sumw2();
		CHI2_MC->Sumw2();
		
		TH1D* nTTbar_double = new TH1D("nTTbar_double","M_{#mub} (GeV);a.u.",50,0,500);
		TH1D* nTTbar_float = new TH1D("nTTbar_float","M_{#mub} (GeV);a.u.",50,0,500);
        
        TH2D* MLB_VS_pT_MC = new TH2D("MLB_VS_pT_MC","MLB_VS_pT_MC;M_{#mub} (GeV);p_{T} (GeV/c)",(centerleftlimit_-leftlimit_)/10,leftlimit_,centerleftlimit_,50,0,500);
        TH2D* MLB_VS_pT_DATA = new TH2D("MLB_VS_pT_DATA","MLB_VS_pT_DATA;M_{#mub} (GeV);p_{T} (GeV/c)",(centerleftlimit_-leftlimit_)/10,leftlimit_,centerleftlimit_,50,0,500);
        
        TH2D* MLB_VS_pTbCand = new TH2D("MLB_VS_pTbCand","MLB_VS_pTbCand;M_{#mub} (GeV);p_{T}^{bCand} (GeV/c)",75,0,500,80,0,800);
        
        TH2D* MLB_VS_ETAbCand = new TH2D("MLB_VS_ETAbCand","MLB_VS_ETAbCand;M_{#mub} (GeV);#eta^{bCand}",75,0,500,40,-4,4);   
        
        TH2D* MLB_VS_pTlCand = new TH2D("MLB_VS_pTlCand","MLB_VS_pTlCand;M_{#mub} (GeV);p_{T}^{lCand} (GeV/c)",75,0,500,80,0,800);
        
        TH2D* MLB_VS_ETAlCand = new TH2D("MLB_VS_ETAlCand","MLB_VS_ETAlCand;M_{#mub} (GeV);#eta^{lCand}",75,0,500,40,-4,4);   
        
        TH2D* MLB_VS_chi2bCand = new TH2D("MLB_VS_chi2bCand","MLB_VS_chi2bCand;M_{#mub} (GeV);#chi^{2}_{bCand}",75,0,500,50,0,200);
                
        TH2D* MLB_VS_chi2lCand = new TH2D("MLB_VS_chi2lCand","MLB_VS_chi2lCand;M_{#mub} (GeV);#chi^{2}_{lCand}",75,0,500,50,0,200);

        TH2D* MLB_VS_TCHE = new TH2D("MLB_VS_TCHE","MLB_VS_TCHE;M_{#mub} (GeV);bDisc",75,0,500,25,-10,30);
        TH2D* MLB_VS_TCHE_l = new TH2D("MLB_VS_TCHE_l","MLB_VS_TCHE;M_{#mub} (GeV);bDisc",75,0,500,25,-10,30);
        TH2D* MLB_VS_TCHE_b = new TH2D("MLB_VS_TCHE_b","MLB_VS_TCHE;M_{#mub} (GeV);bDisc",75,0,500,25,-10,30);
        TH2D* MLB_VS_JP = new TH2D("MLB_VS_JP","MLB_VS_CSV;M_{#mub} (GeV);bDisc",75,0,500,25,0,3);
        TH2D* MLB_VS_JP_l = new TH2D("MLB_VS_JP_l","MLB_VS_CSV;M_{#mub} (GeV);bDisc",75,0,500,25,0,3);
        TH2D* MLB_VS_JP_b = new TH2D("MLB_VS_JP_b","MLB_VS_CSV;M_{#mub} (GeV);bDisc",75,0,500,25,0,3);
        TH2D* MLB_VS_CSV = new TH2D("MLB_VS_CSV","MLB_VS_CSV;M_{#mub} (GeV);bDisc",75,0,500,25,-1,1);
        TH2D* MLB_VS_CSV_l = new TH2D("MLB_VS_CSV_l","MLB_VS_CSV;M_{#mub} (GeV);bDisc",75,0,500,25,-1,1);
        TH2D* MLB_VS_CSV_b = new TH2D("MLB_VS_CSV_b","MLB_VS_CSV;M_{#mub} (GeV);bDisc",75,0,500,25,-1,1);

        
        TH2D* TCHE_VS_CSV_b = new TH2D("TCHE_VS_CSV_b","TCHE_VS_CSV_b;TCHE bDisc;CSV bDisc",50,-10,30,50,-1,1);
        TH2D* TCHE_VS_CSV_l = new TH2D("TCHE_VS_CSV_l","TCHE_VS_CSV_l;TCHE bDisc;CSV bDisc",50,-10,30,50,-1,1);
        TH2D* TCHE_VS_JP_b = new TH2D("TCHE_VS_JP_b","TCHE_VS_JP_b;TCHE bDisc;JP bDisc",50,-10,30,50,0,3);
        TH2D* TCHE_VS_JP_l = new TH2D("TCHE_VS_JP_l","TCHE_VS_JP_l;TCHE bDisc;JP bDisc",50,-10,30,50,0,3);
        
        TH1D* CSV_SHAPE_b = new TH1D("CSV_SHAPE_b","CSV_SHAPE_b;bDisc;#jets",50,-1,1);
        TH1D* CSV_SHAPE_l = new TH1D("CSV_SHAPE_l","CSV_SHAPE_b;bDisc;#jets",50,-1,1);
        TH1D* CSV_SHAPE_atag_b = new TH1D("CSV_SHAPE_atag_b","CSV_SHAPE_atag_b;bDisc;#jets",50,-1,1);
        TH1D* CSV_SHAPE_atag_l = new TH1D("CSV_SHAPE_atag_l","CSV_SHAPE_atag_l;bDisc;#jets",50,-1,1);


        TH1D* MET_right_bCand = new TH1D("MET_right_bCand","MET_right_bCand;MET;#events",80,0,800);
        TH1D* MET_right_lCand = new TH1D("MET_right_lCand","MET_right_lCand;MET;#events",80,0,800);
        
        int nBinsPT=50;
        double firstPT=0;
        double lastPT=500;
        
        TH1D* pT_MC = new TH1D("pT_MC","pT_MC;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pT_MC_Rew = new TH1D("pT_MC_Rew","pT_MC_Rew;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pT_DATA = new TH1D("pT_DATA","pT_DATA;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pT_MC_CS = new TH1D("pT_MC_CS","pT_MC;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pT_DATA_CS = new TH1D("pT_DATA_CS","pT_DATA;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        
        TH1D* pTEnr_MC = new TH1D("pTEnr_MC","pT_MC;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTEnr_MCRew = new TH1D("pTEnr_MCRew","pT_MC;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTEnr_MCRew_enrdata= new TH1D("pTEnr_MCRew_enrdata","pT_MC;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTEnr_DATA = new TH1D("pTEnr_DATA","pT_DATA;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTEnr_DATARew = new TH1D("pTEnr_DATARew","pT_DATA;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        
        TH1D* pTEnr_MC_CS = new TH1D("pTEnr_MC_CS","pT_MC;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTEnr_MCRew_CS = new TH1D("pTEnr_MCRew_CS","pT_MC;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTEnr_DATA_CS = new TH1D("pTEnr_DATA_CS","pT_DATA;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTEnr_DATARew_CS = new TH1D("pTEnr_DATARew_CS","pT_DATA;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        
        TH1D* pTDepl_MC = new TH1D("pTDepl_MC","pT_MC;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTDepl_MCRew = new TH1D("pTDepl_MCRew","pT_MC;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTDepl_DATA = new TH1D("pTDepl_DATA","pT_DATA;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTDepl_DATARew = new TH1D("pTDepl_DATARew","pT_DATA;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        
        TH1D* pTDepl_MC_CS = new TH1D("pTDepl_MC_CS","pT_MC;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTDepl_MCRew_CS = new TH1D("pTDepl_MCRew_CS","pT_MC;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTDepl_DATA_CS = new TH1D("pTDepl_DATA_CS","pT_DATA;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        TH1D* pTDepl_DATARew_CS = new TH1D("pTDepl_DATARew_CS","pT_DATA;p_{T} (GeV/c)",nBinsPT,firstPT,lastPT);
        
        pT_MC->Sumw2();
        pT_MC_Rew->Sumw2();
        pTEnr_MC->Sumw2();
        pTEnr_MCRew->Sumw2();
        pTEnr_MCRew_enrdata->Sumw2();
        
        pT_MC_CS->Sumw2();
        pTEnr_MC_CS->Sumw2();
        pTEnr_MCRew_CS->Sumw2();

        pTDepl_MC->Sumw2();
        pTDepl_MCRew->Sumw2();
        
        pTDepl_MC_CS->Sumw2();
        pTDepl_MCRew_CS->Sumw2();
        
        int nBinsPV=36;
        double firstPV=-0.5;
        double lastPV=35.5;
        
        TH1D* pV_MC = new TH1D("pV_MC","#PV",nBinsPV,firstPV,lastPV);
        TH1D* pV_MC_rew = new TH1D("pV_MC_rew","#PV",nBinsPV,firstPV,lastPV);
        TH1D* pV_DATA = new TH1D("pV_DATA","#PV",nBinsPV,firstPV,lastPV);
        
        pV_MC->Sumw2();
        pV_MC_rew->Sumw2();
        pV_DATA->Sumw2(); 

        TH1D* MLB_weight = new TH1D("MLB_weight","MLB_weight;M_{#mub} (GeV);a.u.",25,0,500);
        
		cout << "+> Filling MultiSamplePlots for L="<<dataLum_ << endl;
		
        double scaleFactor_btag=0;
		for (unsigned int n=0;n<v_weight.size();n++) {
            
            // REENABLE IF UPDATED SF's are available
            
            /*if (decay==0) {
                v_scaleFactor[n]=v_scaleFactor[n]*getIsoMu20Weight(v_ptMuon[n],v_npv[n])*1.00*0.999; // Muon trig leg , jet leg, muon ID
                //weight_nonrew=weight_nonrew*getIsoMu20Weight(NTuple->ptMuon(),NTuple->nPV())*1.00*0.999; // Muon trig leg , jet leg, muon ID
                //weight=weight-0.032;
            } else if (decay==1) {
                //weight=weight*0.988*1.000*0.998;
                v_scaleFactor[n]=v_scaleFactor[n]*0.988*1.000*1.005;
                //weight_nonrew=weight_nonrew*0.987*1.000*0.998;
                
            }*/

            /*double sftt=1;
            double sfbkg=1;
            if (decay==0) {
                sftt=0.998027;
                sfbkg=0.515276;
            } else if (decay == 1) {
                sftt=0.882789;
                sfbkg=1.17455;
            }
                
            if (v_dataSetName[n].find("TTbar") != string::npos)
                scaleFactor_btag=v_scaleFactor[n]*sftt;
            else if (v_dataSetName[n] != "data" && v_dataSetName[n] != "Data" && v_dataSetName[n] != "DATA")
                scaleFactor_btag=v_scaleFactor[n]*sfbkg;
            else*/
            scaleFactor_btag=v_scaleFactor[n];


            //cout << v_weight[n] << endl;

			int d=-1;
			
			for (unsigned int p=0;p < datasets.size(); p++)
				if (datasets[p]->Name() == v_dataSetName[n])
					d=p;
			
			if (d == -1) {
			
				cout << "ERROR FINDING DATASET FOR " << v_dataSetName[n] << endl;
				exit(1);
				
			}
				
			double MCWeight = 1;
            double MCWeight2 = 1;
            double MCWeight3 = 1;
            
            double MCWeight2CS1 = 1;
            double MCWeight2CS2 = 1;
            
            double dataWeight = 1;
            double dataWeightCS1 = 1;
            double dataWeightCS2 = 1;
            
            double MCWeight2_dep = 1;
            double MCWeight2CS1_dep = 1;
            double MCWeight2CS2_dep = 1;

            double dataWeight_dep = 1;
            double dataWeightCS1_dep = 1;
            double dataWeightCS2_dep = 1;
            
            if (hist_PU.find(v_dataSetName[n]+"_runs") == hist_PU.end()) {
                hist_PU[v_dataSetName[n]+"_runs"]=new TH1D((v_dataSetName[n]+"_runs").c_str(),"runs;runID;#ev",500,150000,200000);
            }
            hist_PU[v_dataSetName[n]+"_runs"]->Fill(v_runId[n],v_scaleFactor[n]);
            
            //cout << v_dataSetName[n] << " " << v_weight[n] << endl;
            
            if (v_dataSetName[n] != "data" && v_dataSetName[n] != "Data" && v_dataSetName[n] != "DATA") {
                float w = (dataLum_/datasets[d]->EquivalentLumi());
                if (profile_PU.find(v_dataSetName[n]) == profile_PU.end()) {
                    profile_PU[v_dataSetName[n]]=new TH2D(("profile_btagTCHE_PUweight_"+v_dataSetName[n]).c_str(),("profile_btagTCHE_PUweight_"+v_dataSetName[n]+";btag discr;PU weight").c_str(),30,-10,30,20,0,2);
                }
                if (profile_PU.find("All") == profile_PU.end()) {
                    profile_PU["All"]=new TH2D("profile_btagTCHE_PUweight_ALL","profile_btagTCHE_PUweight_ALL;btag discr;PU weight",30,-10,30,20,0,2);
                }
                
                if (profile_PU.find("MLBAll") == profile_PU.end()) {
                    profile_PU["MLBAll"]=new TH2D("profile_MLB_PUweight_ALL","profile_MLB_PUweight_ALL;MLB;PU weight",50,0,500,20,0,2);
                }
                if (profile_PU.find("MLBAll_TAGT") == profile_PU.end()) {
                    profile_PU["MLBAll_TAGT"]=new TH2D("profile_MLB_PUweight_ALL_BTAGT","profile_MLB_PUweight_ALL;MLB;PU weight",50,0,500,20,0,2);
                }
                if (profile_PU.find("MLBAll_trueb") == profile_PU.end()) {
                    profile_PU["MLBAll_trueb"]=new TH2D("profile_MLB_PUweight_ALL_trueb","profile_MLB_PUweight_ALL;MLB;PU weight",50,0,500,20,0,2);
                }
                if (profile_PU.find("MLBAll_CS") == profile_PU.end()) {
                    profile_PU["MLBAll_CS"]=new TH2D("profile_MLB_PUweight_ALL_CS","profile_MLB_PUweight_ALL_CS;MLB;PU weight",50,0,500,20,0,2);
                }
                
                if (fabs(v_partonFlavour[n]) == 5)
                    profile_PU[v_dataSetName[n]]->Fill(v_bTag[n][0],v_scaleFactor[n]);
                
                if (fabs(v_partonFlavour[n]) == 5)
                    profile_PU["All"]->Fill(v_bTag[n][0],v_scaleFactor[n],w);
            
                profile_PU["MLBAll"]->Fill(v_var[n],v_scaleFactor[n],w);
                if (v_bTag[n][0] > 10.2)
                    profile_PU["MLBAll_TAGT"]->Fill(v_var[n],v_scaleFactor[n],w);
                if (fabs(v_partonFlavour[n]) == 5)
                    profile_PU["MLBAll_trueb"]->Fill(v_var[n],v_scaleFactor[n],w);
                if (v_ptControl[n] >= ptcutextra && v_bTagControl[n] < 3.3)
                    profile_PU["MLBAll_CS"]->Fill(v_controlVar2[n],v_scaleFactor[n],w);
                if (v_ptControl2[n] >= ptcutextra && v_bTagControl2[n] < 3.3)
                    profile_PU["MLBAll_CS"]->Fill(v_controlVar2[n],v_scaleFactor[n],w);
                
                if (profile_PU.find("PT_Lepton_VS_MET") == profile_PU.end()) {
                    profile_PU["PT_Lepton_VS_MET"]=new TH2D("profile_PT_Lepton_VS_MET","profilePT_Lepton_VS_MET;p_{T}^{lepton};MET",50,0,500,80,0,800);
                }
                profile_PU["PT_Lepton_VS_MET"]->Fill(v_ptMuon[n],v_met[n],w);

                
                if (hist_PU.find("MLBL") == hist_PU.end()) {
                    hist_PU["MLBL"]=new TH1D("MLB_LEFT","MLB_LEFT;MLB",50,0,500);
                    hist_PU["MLBL_NoPURew"]=new TH1D("MLBL_NoPURew","MLBL_NoPURew;MLB",50,0,500);
                    hist_PU["MLBR"]=new TH1D("MLB_RIGHT","MLB_RIGHT;MLB",50,0,500);
                    hist_PU["MLBR_NoPURew"]=new TH1D("MLBR_NoPURew","MLBR_NoPURew;MLB",50,0,500);
                }
                if (v_var[n] > leftlimit_ && v_var[n] < centerleftlimit_) {
                    hist_PU["MLBL"]->Fill(v_var[n],w*v_scaleFactor[n]);
                    hist_PU["MLBL_NoPURew"]->Fill(v_var[n],w);
                }
                if (v_var[n] > centerleftlimit_ && v_var[n] < rightlimit_) {
                    hist_PU["MLBR"]->Fill(v_var[n],w*v_scaleFactor[n]);
                    hist_PU["MLBR_NoPURew"]->Fill(v_var[n],w);
                }
                if (hist_PU.find("MLBLCS") == hist_PU.end()) {
                    hist_PU["MLBLCS"]=new TH1D("MLBCS_LEFT","MLB_LEFT;MLB",50,0,500);
                    hist_PU["MLBLCS_NoPURew"]=new TH1D("MLBLCS_NoPURew","MLBL_NoPURew;MLB",50,0,500);
                    hist_PU["MLBRCS"]=new TH1D("MLBCS_RIGHT","MLB_RIGHT;MLB",50,0,500);
                    hist_PU["MLBRCS_NoPURew"]=new TH1D("MLBRCS_NoPURew","MLBR_NoPURew;MLB",50,0,500);
                }
                if (v_ptControl[n] >= ptcutextra && v_bTagControl[n] < 3.3) {
                    if (v_controlVar[n] > leftlimit_ && v_controlVar[n] < centerleftlimit_) {
                        hist_PU["MLBLCS"]->Fill(v_controlVar[n],w*v_scaleFactor[n]);
                        hist_PU["MLBLCS_NoPURew"]->Fill(v_controlVar[n],w);
                    }
                    if (v_controlVar[n] > centerleftlimit_ && v_controlVar[n] < rightlimit_) {
                        hist_PU["MLBRCS"]->Fill(v_controlVar[n],w*v_scaleFactor[n]);
                        hist_PU["MLBRCS_NoPURew"]->Fill(v_controlVar[n],w);
                    }
                }
                if (v_ptControl2[n] >= ptcutextra && v_bTagControl2[n] < 3.3) {
                    if (v_controlVar2[n] > leftlimit_ && v_controlVar2[n] < centerleftlimit_) {
                        hist_PU["MLBLCS"]->Fill(v_controlVar2[n],w*v_scaleFactor[n]);
                        hist_PU["MLBLCS_NoPURew"]->Fill(v_controlVar2[n],w);
                    }
                    if (v_controlVar2[n] > centerleftlimit_ && v_controlVar2[n] < rightlimit_) {
                        hist_PU["MLBRCS"]->Fill(v_controlVar2[n],w*v_scaleFactor[n]);
                        hist_PU["MLBRCS_NoPURew"]->Fill(v_controlVar2[n],w);
                    }
                }
            
                
                //if (v_dataSetName[n].find("TTbarJets") == -1) {
                    
                    string title="MLB_Template_"+v_dataSetName[n];
                    string title2="MLB_Template_BackGround";
                
                if (v_dataSetName[n].find("TTbarJets") == 0) title="MLB_Template_TTbar";
                if (v_dataSetName[n].find("ST_t_") == 0) title="MLB_Template_ST_t";
                if (v_dataSetName[n].find("ST_tW_") == 0) title="MLB_Template_ST_tW";
                    
                    if (hist_PU.find(title) == hist_PU.end()) {
                        hist_PU[title]=new TH1D(title.c_str(),(title+";MLB").c_str(),50,0,500);
                        hist_PU[title+"_JPM"]=new TH1D((title+"_JPM").c_str(),(title+"_JPM;MLB").c_str(),50,0,500);
                        hist_PU[title+"_NonRew"]=new TH1D((title+"_NonRew").c_str(),(title+"_NonRew;MLB").c_str(),50,0,500);
                        hist_PU[title+"_JPM_NonRew"]=new TH1D((title+"_JPM_NonRew").c_str(),(title+"_JPM_NonRew;MLB").c_str(),50,0,500);
                    }
                    if (hist_PU.find(title2) == hist_PU.end()) {
                        hist_PU[title2]=new TH1D(title2.c_str(),(title2+";MLB").c_str(),50,0,500);
                        hist_PU[title2+"_JPM"]=new TH1D((title2+"_JPM").c_str(),(title2+"_JPM;MLB").c_str(),50,0,500);
                        hist_PU[title2+"_NonRew"]=new TH1D((title2+"_NonRew").c_str(),(title2+"_NonRew;MLB").c_str(),50,0,500);
                        hist_PU[title2+"_JPM_NonRew"]=new TH1D((title2+"_JPM_NonRew").c_str(),(title2+"_JPM_NonRew;MLB").c_str(),50,0,500);
                    }
                    //if (v_dataSetName[n] == "WJets") {
                    
                if (v_var[n] < upRangeVar0 &&  v_matchChiSquare[n]<matchChiSquareCut_ ) {

                    float w = (400./datasets[d]->EquivalentLumi());
                    
                    if (v_dataSetName[n].find("multijet") == 0) w=1;
                    
                    //cout << dataLum_ << endl; exit(2);
                    hist_PU[title]->Fill(v_var[n],w*v_scaleFactor[n]);
                    hist_PU[title+"_NonRew"]->Fill(v_var[n],w);
                    
                    if (v_dataSetName[n].find("TTbarJets") == string::npos) {
                        hist_PU[title2]->Fill(v_var[n],w*v_scaleFactor[n]);
                        hist_PU[title2+"_NonRew"]->Fill(v_var[n],w);
                    }
                    
                    //if (v_bTag[n][6] > 0.679) { //CSVM
                    if (v_bTag[n][2] > 0.545) { // JPM
                        hist_PU[title+"_JPM"]->Fill(v_var[n],w*v_scaleFactor[n]);
                        hist_PU[title+"_JPM_NonRew"]->Fill(v_var[n],w);
                        if (v_dataSetName[n].find("TTbarJets") == string::npos) {
                            hist_PU[title2+"_JPM"]->Fill(v_var[n],w*v_scaleFactor[n]);
                            hist_PU[title2+"_JPM_NonRew"]->Fill(v_var[n],w);
                        }
                    }
                }
                //}
            }
            if (v_dataSetName[n] != "data" && v_dataSetName[n] != "Data" && v_dataSetName[n] != "DATA") {
                MCWeight = (dataLum_/datasets[d]->EquivalentLumi())*v_scaleFactor[n];
                
                // LEFT
                double ptWeight=getPtWeight(ptrewfile,"pT_weight_enrmc_alldata",v_pt[n]);
                if (ptWeight != -999)
                    MCWeight2 = MCWeight*ptWeight;
                double ptWeightb=getPtWeight(ptrewfile,"pT_weight_allmc_alldata",v_pt[n]);
                if (ptWeightb != -999)
                    MCWeight3 = MCWeight*ptWeightb;
                
                double ptWeightCS1=getPtWeight(ptrewfile,"pT_weight_CS_enrmc_alldata",v_ptControl[n]);
                if (ptWeightCS1 != -999)
                    MCWeight2CS1 = MCWeight*ptWeightCS1;
                
                double ptWeightCS2=getPtWeight(ptrewfile,"pT_weight_CS_enrmc_alldata",v_ptControl2[n]);
                if (ptWeightCS2 != -999)
                    MCWeight2CS2 = MCWeight*ptWeightCS2;
                
                // RIGHT
                ptWeight=getPtWeight(ptrewfile,"pT_weight_deplmc_alldata",v_pt[n]);
                if (ptWeight != -999)
                    MCWeight2_dep = MCWeight*ptWeight;
                ptWeightCS1=getPtWeight(ptrewfile,"pT_weight_CS_deplmc_alldata",v_ptControl[n]);
                if (ptWeightCS1 != -999)
                    MCWeight2CS1_dep = MCWeight*ptWeightCS1;
                
                ptWeightCS2=getPtWeight(ptrewfile,"pT_weight_CS_deplmc_alldata",v_ptControl2[n]);
                if (ptWeightCS2 != -999)
                    MCWeight2CS2_dep = MCWeight*ptWeightCS2;

			} else {
                // LEFT
                dataWeight=getPtWeight(ptrewfile,"pT_weight_enrdata_alldata",v_pt[n]);
                if (dataWeight == -999)
                    dataWeight=1;
                dataWeightCS1=getPtWeight(ptrewfile,"pT_weight_CS_enrdata_alldata",v_ptControl[n]);
                if (dataWeightCS1 == -999)
                    dataWeightCS1=1;
                dataWeightCS2=getPtWeight(ptrewfile,"pT_weight_CS_enrdata_alldata",v_ptControl2[n]);
                if (dataWeightCS2 == -999)
                    dataWeightCS2=1;
                
                // RIGHT
                dataWeight_dep=getPtWeight(ptrewfile,"pT_weight_depldata_alldata",v_pt[n]);
                if (dataWeight_dep == -999)
                    dataWeight_dep=1;
                dataWeightCS1_dep=getPtWeight(ptrewfile,"pT_weight_CS_depldata_alldata",v_ptControl[n]);
                if (dataWeightCS1_dep == -999)
                    dataWeightCS1_dep=1;
                dataWeightCS2_dep=getPtWeight(ptrewfile,"pT_weight_CS_depldata_alldata",v_ptControl2[n]);
                if (dataWeightCS2_dep == -999)
                    dataWeightCS2_dep=1;

            }
            
			if (v_matchChiSquare[n]<matchChiSquareCut_ ) {//&& v_var[n] > upRangeVar0) {
				nTTbar_double->Fill(v_var[n],v_weight[n]);
				nTTbar_float->Fill(v_var[n],v_weight[n]);
			}
			
            if (v_matchChiSquare[n]<matchChiSquareCut_ ) {
                if (v_dataSetName[n] != "data" && v_dataSetName[n] != "Data" && v_dataSetName[n] != "DATA") {
                    MSPlot["PUWeight"]->Fill(v_scaleFactor[n], datasets[d], true, dataLum_);
                    if (fabs(v_partonFlavour[n] == 5)) {
                        TCHEtrueb_vs_PUWeight->Fill(v_bTag[n][0],v_scaleFactor[n],v_weight[n]/v_scaleFactor[n]);
                        TCHEtrueb_vs_nPV->Fill(v_bTag[n][0],v_npv[n],v_weight[n]/v_scaleFactor[n]);
                        PTtrueb_vs_nPV->Fill(v_pt[n],v_npv[n],v_weight[n]/v_scaleFactor[n]);
                    } else {
                        TCHEtruenonb_vs_nPV->Fill(v_bTag[n][0],v_npv[n],v_weight[n]/v_scaleFactor[n]);
                        PTtruenonb_vs_nPV->Fill(v_pt[n],v_npv[n],v_weight[n]/v_scaleFactor[n]);
                    }
                }
            }

            if (v_dataSetName[n] != "data" && v_dataSetName[n] != "Data" && v_dataSetName[n] != "DATA")

                MLB_VS_PUWeight->Fill(v_var[n],v_scaleFactor[n],v_weight[n]/v_scaleFactor[n]);

            MSPlot["nPV"]->Fill(v_npv[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
            MSPlot["nPV_nonRew"]->Fill(v_npv[n], datasets[d], true, dataLum_);
            
            if (v_matchChiSquare[n]<matchChiSquareCut_ ) {
                MSPlot["pT_lepton"]->Fill(v_ptMuon[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_bTag[n][2] > 0.545) MSPlot["pT_lepton_btag"]->Fill(v_ptMuon[n], datasets[d], true, dataLum_*scaleFactor_btag);
                
                MSPlot["pT_lepbCand"]->Fill(v_pt[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_bTag[n][2] > 0.545) MSPlot["pT_lepbCand_btag"]->Fill(v_pt[n], datasets[d], true, dataLum_*scaleFactor_btag);
                
                MSPlot["pT_WjetCand1"]->Fill(v_ptControl[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_bTag[n][2] > 0.545) MSPlot["pT_WjetCand1_btag"]->Fill(v_ptControl[n], datasets[d], true, dataLum_*scaleFactor_btag);
                
                MSPlot["pT_WjetCand2"]->Fill(v_ptControl2[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_bTag[n][2] > 0.545) MSPlot["pT_WjetCand2_btag"]->Fill(v_ptControl2[n], datasets[d], true, dataLum_*scaleFactor_btag);
                
                MSPlot["Eta_lepton"]->Fill(v_etaMuon[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_bTag[n][2] > 0.545) MSPlot["Eta_lepton_btag"]->Fill(v_etaMuon[n], datasets[d], true, dataLum_*scaleFactor_btag);
                
                MSPlot["Eta_lepbCand"]->Fill(v_eta[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_bTag[n][2] > 0.545) MSPlot["Eta_lepbCand_btag"]->Fill(v_eta[n], datasets[d], true, dataLum_*scaleFactor_btag);
                
                MSPlot["Eta_WjetCand1"]->Fill(v_etaControl[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_bTag[n][2] > 0.545) MSPlot["Eta_WjetCand1_btag"]->Fill(v_etaControl[n], datasets[d], true, dataLum_*scaleFactor_btag);
                
                MSPlot["Eta_WjetCand2"]->Fill(v_etaControl2[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_bTag[n][2] > 0.545) MSPlot["Eta_WjetCand2_btag"]->Fill(v_etaControl2[n], datasets[d], true, dataLum_*scaleFactor_btag);

                MSPlot["nJets"]->Fill(v_njets[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_bTag[n][2] > 0.545) MSPlot["nJets_btag"]->Fill(v_njets[n], datasets[d], true, dataLum_*scaleFactor_btag);
                
                MSPlot["MET"]->Fill(v_met[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_bTag[n][2] > 0.545) MSPlot["MET_btag"]->Fill(v_met[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                
                //cout << v_met[n] << endl;
                
                MSPlot["MVATrigID"]->Fill(v_mvaTrigId[n], datasets[d], true, dataLum_*v_scaleFactor[n]);


            }
            
            MSPlot["BestJetCombChi2"]->Fill(v_matchChiSquare[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
            
			MSPlot["selectedEventsJetsEta"]->Fill(v_eta[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
            
            //** Plot some taggers **//
        
            if (MSPlot.find("bTagger_TCHE") == MSPlot.end()) {
                MSPlot["bTagger_TCHE"] = new MultiSamplePlot(datasets, "bTagger_TCHE", 50, -10, 30, "TCHE Discriminator","Number of jets"); 
                MSPlot["bTagger_TCHP"] = new MultiSamplePlot(datasets, "bTagger_TCHP", 50, -10, 30, "TCHP Discriminator","Number of jets");
                MSPlot["bTagger_JP"] = new MultiSamplePlot(datasets, "bTagger_JP", 50, 0, 3, "JP Discriminator","Number of jets"); 
                MSPlot["bTagger_JBP"] = new MultiSamplePlot(datasets, "bTagger_JBP", 50, 0, 8, "JBP Discriminator","Number of jets");
                MSPlot["bTagger_SSVHE"] = new MultiSamplePlot(datasets, "bTagger_SSVHE", 50, 0, 8, "SSVHE Discriminator","Number of jets"); 
                MSPlot["bTagger_SSVHP"] = new MultiSamplePlot(datasets, "bTagger_SSVHP", 50, 0, 8, "SSVHP Discriminator","Number of jets");
                MSPlot["bTagger_CSV"] = new MultiSamplePlot(datasets, "bTagger_CSV", 50, 0, 1, "CSV Discriminator","Number of jets"); 
                MSPlot["bTagger_CSVR"] = new MultiSamplePlot(datasets, "bTagger_CSVR", 50, 0, 1, "CSVR Discriminator","Number of jets"); 
            }
            
            MSPlot["bTagger_TCHE"]->Fill(v_bTag[n][0], datasets[d], true, dataLum_*v_scaleFactor[n]);   
            MSPlot["bTagger_TCHP"]->Fill(v_bTag[n][1], datasets[d], true, dataLum_*v_scaleFactor[n]);   
            MSPlot["bTagger_JP"]->Fill(v_bTag[n][2], datasets[d], true, dataLum_*v_scaleFactor[n]);   
            MSPlot["bTagger_JBP"]->Fill(v_bTag[n][3], datasets[d], true, dataLum_*v_scaleFactor[n]); 
            MSPlot["bTagger_SSVHE"]->Fill(v_bTag[n][4], datasets[d], true, dataLum_*v_scaleFactor[n]);   
            MSPlot["bTagger_SSVHP"]->Fill(v_bTag[n][5], datasets[d], true, dataLum_*v_scaleFactor[n]); 
            MSPlot["bTagger_CSV"]->Fill(v_bTag[n][6], datasets[d], true, dataLum_*v_scaleFactor[n]);   
            MSPlot["bTagger_CSVR"]->Fill(v_bTag[n][7], datasets[d], true, dataLum_*v_scaleFactor[n]);   



			if (v_dataSetName[n] != "data" && v_dataSetName[n] != "Data" && v_dataSetName[n] != "DATA") {
				CHI2_MC->Fill(v_matchChiSquare[n],MCWeight);

				if (v_bTagControl[n] < 3 && v_bTagControl2[n] < 3) {
					MLBCS_MC->Fill(v_controlVar[n],v_weight[n]);					
					MLBCS_MC->Fill(v_controlVar2[n],v_weight[n]);
				}

			} else {
				MLB_DATA->Fill(v_var[n]);
				CHI2_DATA->Fill(v_matchChiSquare[n]);
				
				if (v_bTagControl[n] < 3 && v_bTagControl2[n] < 3) {
					MLBCS_DATA->Fill(v_controlVar[n]);					
					MLBCS_DATA->Fill(v_controlVar2[n]);					
				}
				
			}
			
            if (v_matchChiSquare[n]<matchChiSquareCut_ ) {
                MSPlot["MLB2"]->Fill(v_var[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
            MSPlot["MLB"]->Fill(v_var[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
			MSPlot["MLB_BTV"]->Fill(v_var[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_met[n]>40) {
                    MSPlot["MLB_MET"]->Fill(v_var[n], datasets[d], true, dataLum_*scaleFactor_btag);
                }

            if (v_njets[n] == 4)
                MSPlot["MLB_FourJets"]->Fill(v_var[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
			MSPlot["M3_BTV"]->Fill(v_varb[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
			}
            
            if (v_matchChiSquare[n]<matchChiSquareCut_ ) {
                if (v_var[n] > leftlimit_ && v_var[n] <= centerleftlimit_)
                    MSPlot["MLBL"]->Fill(v_var[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
                if (v_var[n] > centerrightlimit_ && v_var[n] < rightlimit_)
                    MSPlot["MLBR"]->Fill(v_var[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
			}

			if (v_ptControl[n] > ptcutextra && v_ptControl2[n] > ptcutextra) {
				MSPlot["MLB_ControlSample"]->Fill(v_controlVar[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
				MSPlot["MLB_ControlSample"]->Fill(v_controlVar2[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
				MSPlot["MLB_ControlSample_q1"]->Fill(v_controlVar[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
				MSPlot["MLB_ControlSample_q2"]->Fill(v_controlVar2[n], datasets[d], true, dataLum_*v_scaleFactor[n]);
			}
			
            if (v_matchChiSquare[n]<matchChiSquareCut_) {

			//if (v_bTag[n][0] > 3.3) { // TCHEM
            //if (v_bTag[n][6] > 0.679) { // CSVM
                if (v_bTag[n][2] > 0.545) { // JPM
                    MSPlot["MLB_btag"]->Fill(v_var[n], datasets[d], true, dataLum_*scaleFactor_btag);
                    if (v_met[n]>40) {
                        MSPlot["MLB_btag_MET"]->Fill(v_var[n], datasets[d], true, dataLum_*scaleFactor_btag);
                    }
                    MSPlot["MLB_BTV_btag"]->Fill(v_var[n], datasets[d], true, dataLum_*scaleFactor_btag);
                    if (v_njets[n] == 4)
                        MSPlot["MLB_FourJets_btag"]->Fill(v_var[n], datasets[d], true, dataLum_*scaleFactor_btag);
                    MSPlot["M3_BTV_btag"]->Fill(v_varb[n], datasets[d], true, dataLum_*scaleFactor_btag);
                }
            }
            
            if (v_dataSetName[n] != "data" && v_dataSetName[n] != "Data" && v_dataSetName[n] != "DATA") {
                
                MLB_VS_pT_MC->Fill(v_var[n],v_pt[n],MCWeight);
                
                if (v_matchChiSquare[n]<matchChiSquareCut_) {
                    MLB_VS_TCHE->Fill(v_var[n],v_bTag[n][0],MCWeight);
                    MLB_VS_CSV->Fill(v_var[n],v_bTag[n][6],MCWeight);
                    MLB_VS_JP->Fill(v_var[n],v_bTag[n][2],MCWeight);
                    
                    if (fabs(v_partonFlavour[n]) == 5) {
                        MLB_VS_pTbCand->Fill(v_var[n],v_pt[n],MCWeight);
                        MLB_VS_ETAbCand->Fill(v_var[n],v_eta[n],MCWeight);
                        MLB_VS_chi2bCand->Fill(v_var[n],v_matchChiSquare[n],MCWeight);
                        MLB_VS_TCHE_b->Fill(v_var[n],v_bTag[n][0],MCWeight);
                        MLB_VS_CSV_b->Fill(v_var[n],v_bTag[n][6],MCWeight);
                        MLB_VS_JP_b->Fill(v_var[n],v_bTag[n][2],MCWeight);
                        
                        TCHE_VS_CSV_b->Fill(v_bTag[n][0],v_bTag[n][6],MCWeight);
                        TCHE_VS_JP_b->Fill(v_bTag[n][0],v_bTag[n][2],MCWeight);
                        
                        CSV_SHAPE_b->Fill(v_bTag[n][6],MCWeight);
                        if (v_bTag[n][0] < 3)
                            CSV_SHAPE_atag_b->Fill(v_bTag[n][6],MCWeight);

                    } else {
                        MLB_VS_pTlCand->Fill(v_var[n],v_pt[n],MCWeight);
                        MLB_VS_ETAlCand->Fill(v_var[n],v_eta[n],MCWeight);
                        MLB_VS_chi2lCand->Fill(v_var[n],v_matchChiSquare[n],MCWeight);
                        MLB_VS_TCHE_l->Fill(v_var[n],v_bTag[n][0],MCWeight);
                        MLB_VS_CSV_l->Fill(v_var[n],v_bTag[n][6],MCWeight);
                        MLB_VS_JP_l->Fill(v_var[n],v_bTag[n][2],MCWeight);
                    
                        TCHE_VS_CSV_l->Fill(v_bTag[n][0],v_bTag[n][6],MCWeight);
                        TCHE_VS_JP_l->Fill(v_bTag[n][0],v_bTag[n][2],MCWeight);
                        
                        CSV_SHAPE_l->Fill(v_bTag[n][6],MCWeight);
                        if (v_bTag[n][0] < 3)
                            CSV_SHAPE_atag_l->Fill(v_bTag[n][6],MCWeight);

                    }
                    
                    if (v_var[n] > centerrightlimit_ && v_var[n] < rightlimit_) {
                        if (fabs(v_partonFlavour[n]) == 5) {
                            MET_right_bCand->Fill(v_met[n],MCWeight);
                        } else {
                            MET_right_lCand->Fill(v_met[n],MCWeight);
                        }
                    }
                }
                
                pT_MC->Fill(v_pt[n],MCWeight);
                pT_MC_Rew->Fill(v_pt[n],MCWeight*getPtWeight(ptrewfile,"pT_weight_allmc_alldata",v_pt[n]));
                
                pV_MC->Fill(v_npv[n],MCWeight);
                pV_MC_rew->Fill(v_npv[n],MCWeight*getPtWeight(ptrewfile,"pV_weight",v_npv[n]));
                

                if (v_ptControl[n] >= 30) pT_MC_CS->Fill(v_ptControl[n],MCWeight);
                if (v_ptControl2[n] >= 30) pT_MC_CS->Fill(v_ptControl2[n],MCWeight);
                
                //LEFT
                
                if (v_var[n] > leftlimit_ && v_var[n] < centerleftlimit_) {
                    pTEnr_MCRew->Fill(v_pt[n],MCWeight2);
                    pTEnr_MCRew_enrdata->Fill(v_pt[n],MCWeight3);
                    pTEnr_MC->Fill(v_pt[n],MCWeight);
                }
                
                if (v_controlVar[n] > leftlimit_ && v_controlVar[n] < centerleftlimit_) {
                    if (v_ptControl[n] >= 30) pTEnr_MC_CS->Fill(v_ptControl[n],MCWeight);
                    if (v_ptControl[n] >= 30) pTEnr_MCRew_CS->Fill(v_ptControl[n],MCWeight2CS1);
                }
                if (v_controlVar2[n] > leftlimit_ && v_controlVar2[n] < centerleftlimit_) {
                    if (v_ptControl2[n] >= 30) pTEnr_MC_CS->Fill(v_ptControl2[n],MCWeight);
                    if (v_ptControl2[n] >= 30) pTEnr_MCRew_CS->Fill(v_ptControl2[n],MCWeight2CS2);
                }
                
                // RIGHT
                
                if (v_var[n] > centerrightlimit_ && v_var[n] < rightlimit_) {
                    pTDepl_MCRew->Fill(v_pt[n],MCWeight2_dep);
                    pTDepl_MC->Fill(v_pt[n],MCWeight);
                }
                
                if (v_controlVar[n] > centerrightlimit_ && v_controlVar[n] < rightlimit_) {
                    if (v_ptControl[n] >= 30) pTDepl_MC_CS->Fill(v_ptControl[n],MCWeight);
                    if (v_ptControl[n] >= 30) pTDepl_MCRew_CS->Fill(v_ptControl[n],MCWeight2CS1_dep);
                }
                if (v_controlVar2[n] > centerrightlimit_ && v_controlVar2[n] < rightlimit_) {
                    if (v_ptControl2[n] >= 30) pTDepl_MC_CS->Fill(v_ptControl2[n],MCWeight);
                    if (v_ptControl2[n] >= 30) pTDepl_MCRew_CS->Fill(v_ptControl2[n],MCWeight2CS2_dep);
                }

                
            } else {
                
                MLB_VS_pT_DATA->Fill(v_var[n],v_pt[n]);

                pT_DATA->Fill(v_pt[n]);
                pV_DATA->Fill(v_npv[n]);

                if (v_ptControl[n] >= 30) pT_DATA_CS->Fill(v_ptControl[n]);
                if (v_ptControl2[n] >= 30) pT_DATA_CS->Fill(v_ptControl2[n]);
                
                // LEFT
                if (v_var[n] > leftlimit_ && v_var[n] <= centerleftlimit_) {
                    pTEnr_DATA->Fill(v_pt[n]);
                    pTEnr_DATARew->Fill(v_pt[n],dataWeight);
                }
                
                if (v_controlVar[n] > leftlimit_ && v_controlVar[n] < centerleftlimit_) {
                    if (v_ptControl[n] >= 30) pTEnr_DATA_CS->Fill(v_ptControl[n]);
                    if (v_ptControl[n] >= 30) pTEnr_DATARew_CS->Fill(v_ptControl[n],dataWeightCS1);
                }
                if (v_controlVar2[n] > leftlimit_ && v_controlVar2[n] < centerleftlimit_) {
                    if (v_ptControl2[n] >= 30) pTEnr_DATA_CS->Fill(v_ptControl2[n]);
                    if (v_ptControl2[n] >= 30) pTEnr_DATARew_CS->Fill(v_ptControl2[n],dataWeightCS2);
                }
                
                // RIGHT
                if (v_var[n] > centerrightlimit_ && v_var[n] < rightlimit_) {
                    pTDepl_DATA->Fill(v_pt[n]);
                    pTDepl_DATARew->Fill(v_pt[n],dataWeight_dep);
                }
                
                if (v_controlVar[n] > centerrightlimit_ && v_controlVar[n] < rightlimit_) {
                    if (v_ptControl[n] >= 30) pTDepl_DATA_CS->Fill(v_ptControl[n]);
                    if (v_ptControl[n] >= 30) pTDepl_DATARew_CS->Fill(v_ptControl[n],dataWeightCS1_dep);
                }
                if (v_controlVar2[n] > centerrightlimit_ && v_controlVar2[n] < rightlimit_) {
                    if (v_ptControl2[n] >= 30) pTDepl_DATA_CS->Fill(v_ptControl2[n]);
                    if (v_ptControl2[n] >= 30) pTDepl_DATARew_CS->Fill(v_ptControl2[n],dataWeightCS2_dep);
                }

            }

			// testing mistag rate
			
			if (v_dataSetName[n] != "data" && v_dataSetName[n] != "Data" && v_dataSetName[n] != "DATA") {

				if (v_matchChiSquare[n]<matchChiSquareCut_ && v_var[n] < 500) {
                    
                    MLB_MC->Fill(v_var[n],MCWeight);

					if (fabs(v_partonFlavour[n]) == 5) {
						all_b+=v_weight[n];
						if (v_bTag[n][0] > 3.3)
							tag_b+=v_weight[n];
					}
					else if (fabs(v_partonFlavour[n]) == 4) {
						all_c+=v_weight[n];
						if (v_bTag[n][0] > 3.3)
							tag_c+=v_weight[n];
					}
					else {
						all_dusg+=v_weight[n];
						if (v_bTag[n][0] > 3.3)
							tag_dusg+=v_weight[n];
					}
					
				}
				
			}
			
			if (d < datasets_MC.size()) {
				if (fabs(v_partonFlavour[n]) == 5) {
                    MSPlot["MLB-bquarks"]->Fill(v_var[n], datasets_MC[d], true, dataLum_*v_scaleFactor[n]);
                    if (v_bTag[n][0] > 1.7)
                        MSPlot["MLB-bquarks_btag"]->Fill(v_var[n], datasets_MC[d], true, dataLum_*v_scaleFactor[n]);
                }
				if (fabs(v_partonFlavour[n]) != 5) {
                    MSPlot["MLB-nonbquarks"]->Fill(v_var[n], datasets_MC[d], true, dataLum_*v_scaleFactor[n]);
                    if (v_bTag[n][0] > 1.7)
                        MSPlot["MLB-nonbquarks_btag"]->Fill(v_var[n], datasets_MC[d], true, dataLum_*v_scaleFactor[n]);
                }
				if (v_dataSetName[n].find("TTbar") == 0) {
					MLB_ttbar->Fill(v_var[n],v_weight[n]);
				} else if (v_dataSetName[n].find("WJets") == 0) {
					MLB_wjets->Fill(v_var[n],v_weight[n]);
					MLB_wjets_noweight->Fill(v_var[n],1);
				} else {
					MLB_bkg->Fill(v_var[n],v_weight[n]);
				}
				
				if (v_matchChiSquare[n]<matchChiSquareCut_) {
					if (fabs(v_partonFlavour[n]) == 5) nB+=v_weight[n];
					if (fabs(v_partonFlavour[n]) != 5) nNonB+=v_weight[n];
					
					nData+=v_weight[n];
					
					if (v_var[n] > leftlimit_ && v_var[n] < centerleftlimit_) {		
						nDataL+=v_weight[n];
						if (fabs(v_partonFlavour[n]) == 5) nBL+=v_weight[n];
						if (fabs(v_partonFlavour[n]) != 5) nNonBL+=v_weight[n];
					}
					if (v_var[n] >= centerrightlimit_ && v_var[n] < rightlimit) {		
						nDataR+=v_weight[n];
						if (fabs(v_partonFlavour[n]) == 5) nBR+=v_weight[n];
						if (fabs(v_partonFlavour[n]) != 5) nNonBR+=v_weight[n];
					}
					
					if (fabs(v_partonFlavourControl[n]) == 5) nBCSA+=v_weight[n];
					if (fabs(v_partonFlavourControl[n]) != 5) nNonBCSA+=v_weight[n];
					if (fabs(v_partonFlavourControl2[n]) == 5) nBCSA+=v_weight[n];
					if (fabs(v_partonFlavourControl2[n]) != 5) nNonBCSA+=v_weight[n];
					
					//if (v_bTagControl[n] < 3 && v_bTagControl2[n] < 3) { 
                    if (v_ptControl[n] > ptcutextra && v_ptControl2[n] > ptcutextra) {
						// CS JET 1
						if (fabs(v_partonFlavourControl[n]) == 5) nBCS+=v_weight[n];
						if (fabs(v_partonFlavourControl[n]) != 5) nNonBCS+=v_weight[n];
						
						if (v_controlVar[n] > leftlimit && v_controlVar[n] < centerleftlimit) {			
							if (fabs(v_partonFlavourControl[n]) == 5) nBLCS+=v_weight[n];
							if (fabs(v_partonFlavourControl[n]) != 5) nNonBLCS+=v_weight[n];
						}
						if (v_controlVar[n] > centerrightlimit_ && v_controlVar[n] < rightlimit) {			
							if (fabs(v_partonFlavourControl[n]) == 5) nBRCS+=v_weight[n];
							if (fabs(v_partonFlavourControl[n]) != 5) nNonBRCS+=v_weight[n];
						}
						
						// CS JET 2
						if (fabs(v_partonFlavourControl2[n]) == 5) nBCS+=v_weight[n];
						if (fabs(v_partonFlavourControl2[n]) != 5) nNonBCS+=v_weight[n];
						
						if (v_controlVar2[n] > leftlimit && v_controlVar2[n] < centerleftlimit) {			
							if (fabs(v_partonFlavourControl2[n]) == 5) nBLCS+=v_weight[n];
							if (fabs(v_partonFlavourControl2[n]) != 5) nNonBLCS+=v_weight[n];
						}
						if (v_controlVar2[n] > centerrightlimit_ && v_controlVar2[n] < rightlimit) {			
							if (fabs(v_partonFlavourControl2[n]) == 5) nBRCS+=v_weight[n];
							if (fabs(v_partonFlavourControl2[n]) != 5) nNonBRCS+=v_weight[n];
						}
					}
					
				}
			}
			
			if (v_dataSetName[n] != "data" && v_dataSetName[n] != "Data" && v_dataSetName[n] != "DATA") {
				corr_histos[1]->Fill(v_var[n],v_bTag[n][0]);
				if (fabs(v_partonFlavour[n]) == 5) corr_histos[2]->Fill(v_var[n],v_bTag[n][0],v_weight[n]);
				if (fabs(v_partonFlavour[n]) != 5) corr_histos[3]->Fill(v_var[n],v_bTag[n][0],v_weight[n]);
			} else {
				corr_histos[0]->Fill(v_var[n],v_bTag[n][0]);
			}

		}

		cout << "+> Saving MultiSamplePlots" << endl;
			
		mkdir("MSPlot",0777);
		mkdir(("MSPlot/channel"+data_postfix).c_str(),0777);
		//TFile* fMS = new TFile("JetMuonMass_StackPlot.root","RECREATE");
		TFile* fMS = new TFile(("StackPlots"+data_postfix+".root").c_str(),"RECREATE");
		
		for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
		{
			MultiSamplePlot *temp = it->second;
			string name = it->first;
			temp->showNumberEntries(false);
			temp->setDataLumi(round(1,dataLum_));
			temp->Draw(false, name, false, true, true, true, true);
			temp->Write(fout, name, true, ("MSPlot/channel"+data_postfix+"/").c_str());
			//if (name == "MLB") 
            temp->Write(fMS, name, true, ("MSPlot/channel"+data_postfix+"/").c_str());

		}
		fMS->Close();
		
		TFile* fout = new TFile("corrHistos.root","RECREATE");
		for (unsigned int p=0;p<corr_histos.size();p++) {
			
			fout->cd();
			
			corr_histos[p]->Write();
			
		}
        
        for(map<string,TH2D*>::const_iterator it = profile_PU.begin(); it != profile_PU.end(); it++) {
            it->second->Write();
            TProfile* profile = (TProfile*) it->second->ProfileX(("TPRofile_btag_vs_PUWeight_"+it->first).c_str());
            profile->Write();
            delete profile;

        }

        for(map<string,TH1D*>::const_iterator it = hist_PU.begin(); it != hist_PU.end(); it++)
            it->second->Write();

		MLB_ttbar->Scale(1./MLB_ttbar->Integral());
		MLB_bkg->Scale(1./MLB_bkg->Integral());

		MLB_wjets->Scale(1./MLB_wjets->Integral());
		MLB_wjets_noweight->Scale(1./MLB_wjets_noweight->Integral());

		MLB_ttbar->Write();
		MLB_bkg->Write();

		MLB_wjets->Write();
		MLB_wjets_noweight->Write();
		
		MLB_DATA->Write();
		MLBCS_DATA->Write();
		CHI2_DATA->Write();
		
		MLB_MC->Write();
		MLBCS_MC->Write();
		CHI2_MC->Write();
		
		TH1D* MLB_RATIO = (TH1D*) MLB_DATA->Clone();
		MLB_RATIO->SetName("MLB_RATIO");
		MLB_RATIO->Divide(MLB_MC);
		MLB_RATIO->Write();
		
		TH1D* MLBCS_RATIO = (TH1D*) MLBCS_DATA->Clone();
		MLBCS_RATIO->SetName("MLBCS_RATIO");
		MLBCS_RATIO->Divide(MLBCS_MC);
		MLBCS_RATIO->Write();
		
		TH1D* CHI2_RATIO = (TH1D*) CHI2_DATA->Clone();
		CHI2_RATIO->SetName("CHI2_RATIO");
		CHI2_RATIO->Divide(CHI2_MC);
		CHI2_RATIO->Write();
		
		nTTbar_float->Write();
		nTTbar_double->Write();
        
        MLB_VS_pT_MC->Write();
        MLB_VS_pT_DATA->Write();
        
        MLB_VS_pTbCand->Write();
        MLB_VS_ETAbCand->Write();
        MLB_VS_pTlCand->Write();
        MLB_VS_ETAlCand->Write();
        MLB_VS_chi2bCand->Write();
        MLB_VS_chi2lCand->Write();
        
        MLB_VS_TCHE->Write();
        MLB_VS_TCHE_l->Write();
        MLB_VS_TCHE_b->Write();
        MLB_VS_CSV->Write();
        MLB_VS_CSV_l->Write();
        MLB_VS_CSV_b->Write();
        MLB_VS_JP->Write();
        MLB_VS_JP_l->Write();
        MLB_VS_JP_b->Write();
        
        TCHE_VS_CSV_b->Write();
        TCHE_VS_JP_b->Write();
        TCHE_VS_CSV_l->Write();
        TCHE_VS_JP_l->Write();
        
        CSV_SHAPE_b->Write();
        CSV_SHAPE_atag_b->Write();
        CSV_SHAPE_l->Write();
        CSV_SHAPE_atag_l->Write();
        
        TProfile* MLB_VS_TCHE_profile = (TProfile*) MLB_VS_TCHE->ProfileX("MLB_VS_TCHE_profile");
        TProfile* MLB_VS_TCHE_l_profile = (TProfile*) MLB_VS_TCHE_l->ProfileX("MLB_VS_TCHE_l_profile");
        TProfile* MLB_VS_TCHE_b_profile = (TProfile*) MLB_VS_TCHE_b->ProfileX("MLB_VS_TCHE_b_profile");
        TProfile* MLB_VS_CSV_profile = (TProfile*) MLB_VS_CSV->ProfileX("MLB_VS_CSV_profile");
        TProfile* MLB_VS_CSV_l_profile = (TProfile*) MLB_VS_CSV_l->ProfileX("MLB_VS_CSV_l_profile");
        TProfile* MLB_VS_CSV_b_profile = (TProfile*) MLB_VS_CSV_b->ProfileX("MLB_VS_CSV_b_profile");
        TProfile* MLB_VS_JP_profile = (TProfile*) MLB_VS_JP->ProfileX("MLB_VS_JP_profile");
        TProfile* MLB_VS_JP_l_profile = (TProfile*) MLB_VS_JP_l->ProfileX("MLB_VS_JP_l_profile");
        TProfile* MLB_VS_JP_b_profile = (TProfile*) MLB_VS_JP_b->ProfileX("MLB_VS_JP_b_profile");

        MLB_VS_TCHE_profile->Write();
        MLB_VS_TCHE_l_profile->Write();
        MLB_VS_TCHE_b_profile->Write();
        MLB_VS_CSV_profile->Write();
        MLB_VS_CSV_l_profile->Write();
        MLB_VS_CSV_b_profile->Write();
        MLB_VS_JP_profile->Write();
        MLB_VS_JP_l_profile->Write();
        MLB_VS_JP_b_profile->Write();
        
        TProfile* profile_bcandpt = (TProfile*) MLB_VS_pTbCand->ProfileX("TPRofile_MLB_VS_pTbCand");
        TProfile* profile_bcandeta = (TProfile*) MLB_VS_ETAbCand->ProfileX("TPRofile_MLB_VS_ETAbCand");
        TProfile* profile_lcandpt = (TProfile*) MLB_VS_pTlCand->ProfileX("TPRofile_MLB_VS_pTlCand");
        TProfile* profile_lcandeta = (TProfile*) MLB_VS_ETAlCand->ProfileX("TPRofile_MLB_VS_ETAlCand");

        TProfile* profile_bcandchi2 = (TProfile*) MLB_VS_chi2bCand->ProfileX("TPRofile_MLB_VS_chi2bCand");
        TProfile* profile_lcandchi2 = (TProfile*) MLB_VS_chi2lCand->ProfileX("TPRofile_MLB_VS_chi2lCand");

        profile_bcandpt->Write();
        profile_bcandeta->Write();
        profile_lcandpt->Write();
        profile_lcandeta->Write();
        profile_bcandchi2->Write();
        profile_lcandchi2->Write();
        
        MET_right_bCand->Write();
        MET_right_lCand->Write();
        
        TCHEtrueb_vs_PUWeight->Write();
        TCHEtrueb_vs_nPV->Write();
        TCHEtruenonb_vs_nPV->Write();
        PTtrueb_vs_nPV->Write();
        PTtruenonb_vs_nPV->Write();
        MLB_VS_PUWeight->Write();
        
 		
		fout->Close();
        
		TFile* fout2 = new TFile("datamc-ptrew.root","RECREATE");
        
        TDirectory* dir = fout2->mkdir("ControlPlots");
        
        dir->cd();
        
        pTEnr_MC->Write();
        pTEnr_MCRew->Write();
        pTEnr_MCRew_enrdata->Write();
        pTEnr_DATA->Write();
        pTEnr_DATARew->Write();
        
        pTEnr_MC_CS->Write();
        pTEnr_MCRew_CS->Write();
        pTEnr_DATA_CS->Write();
        pTEnr_DATARew_CS->Write();

        pTDepl_MC->Write();
        pTDepl_MCRew->Write();
        pTDepl_DATA->Write();
        pTDepl_DATARew->Write();
        
        pTDepl_MC_CS->Write();
        pTDepl_MCRew_CS->Write();
        pTDepl_DATA_CS->Write();
        pTDepl_DATARew_CS->Write();

        pV_MC->Scale(1./pV_MC->Integral());
        pV_MC_rew->Scale(1./pV_MC_rew->Integral());
        pV_DATA->Scale(1./pV_DATA->Integral());

        pV_MC->Write();
        pV_MC_rew->Write();
        pV_DATA->Write();
        
        pT_MC->Write();
        pT_MC_Rew->Write();
        pT_DATA->Write();
        
        pT_MC_CS->Write();
        pT_DATA_CS->Write();
        
        MLB_MC->Write();
        MLB_DATA->Write();
        
        MLB_weight->Add(MLB_DATA);
        MLB_weight->Divide(MLB_MC);
        
        MLB_weight->Write();
        
        fout2->cd();
        
        TH1D* pT_weight_data = new TH1D("pT_weight_enrdata_alldata","pT_weight;p_{T} (GeV/c);weight",nBinsPT,firstPT,lastPT);
        
        TH1D* pT_weight_mc = new TH1D("pT_weight_enrmc_alldata","pT_weight;p_{T} (GeV/c);weight",nBinsPT,firstPT,lastPT);
        
        TH1D* pT_weight = new TH1D("pT_weight_allmc_alldata","pT_weight;p_{T} (GeV/c);weight",nBinsPT,firstPT,lastPT);

        TH1D* pT_weight_mc2 = new TH1D("pT_weight_enrmc_enrdata","pT_weight;p_{T} (GeV/c);weight",nBinsPT,firstPT,lastPT);

        TH1D* pT_weight_dataCS = new TH1D("pT_weight_CS_enrdata_alldata","pT_weight;p_{T} (GeV/c);weight",nBinsPT,firstPT,lastPT);
        
        TH1D* pT_weight_mcCS = new TH1D("pT_weight_CS_enrmc_alldata","pT_weight;p_{T} (GeV/c);weight",nBinsPT,firstPT,lastPT);
        
        TH1D* pT_weightCS = new TH1D("pT_weight_CS_allmc_alldata","pT_weight;p_{T} (GeV/c);weight",nBinsPT,firstPT,lastPT);
        
        TH1D* pT_weight_data_depl = new TH1D("pT_weight_depldata_alldata","pT_weight;p_{T} (GeV/c);weight",nBinsPT,firstPT,lastPT);
        
        TH1D* pT_weight_mc_depl = new TH1D("pT_weight_deplmc_alldata","pT_weight;p_{T} (GeV/c);weight",nBinsPT,firstPT,lastPT);
        
        TH1D* pT_weight_data_deplCS = new TH1D("pT_weight_CS_depldata_alldata","pT_weight;p_{T} (GeV/c);weight",nBinsPT,firstPT,lastPT);
        
        TH1D* pT_weight_mc_deplCS = new TH1D("pT_weight_CS_deplmc_alldata","pT_weight;p_{T} (GeV/c);weight",nBinsPT,firstPT,lastPT);

        pT_weight_mc->Add(pT_DATA);
        pT_weight_mc->Divide(pTEnr_MC);
        pT_weight_mc->Write();
        
        pT_weight_mc2->Add(pTEnr_DATA);
        pT_weight_mc2->Divide(pTEnr_MC);
        pT_weight_mc2->Write();
        
        pT_weight_data->Add(pT_DATA);
        pT_weight_data->Divide(pTEnr_DATA);
        pT_weight_data->Write();
        
        pT_weight->Add(pT_DATA);
        pT_weight->Divide(pT_MC);
        pT_weight->Write();

        pT_weight_mcCS->Add(pT_DATA_CS);
        pT_weight_mcCS->Divide(pTEnr_MC_CS);
        pT_weight_mcCS->Write();
        
        pT_weight_dataCS->Add(pT_DATA_CS);
        pT_weight_dataCS->Divide(pTEnr_DATA_CS);
        pT_weight_dataCS->Write();
        
        pT_weightCS->Add(pT_DATA_CS);
        pT_weightCS->Divide(pT_MC_CS);
        pT_weightCS->Write();
        
        pT_weight_mc_depl->Add(pT_DATA);
        pT_weight_mc_depl->Divide(pTDepl_MC);
        pT_weight_mc_depl->Write();
        
        pT_weight_data_depl->Add(pT_DATA);
        pT_weight_data_depl->Divide(pTDepl_DATA);
        pT_weight_data_depl->Write();
        
        pT_weight_mc_deplCS->Add(pT_DATA_CS);
        pT_weight_mc_deplCS->Divide(pTDepl_MC_CS);
        pT_weight_mc_deplCS->Write();
        
        pT_weight_data_deplCS->Add(pT_DATA_CS);
        pT_weight_data_deplCS->Divide(pTDepl_DATA_CS);
        pT_weight_data_deplCS->Write();
        
        TH1D* pV_weight = new TH1D("pV_weight","pV_weight;#PV;weight",nBinsPV,firstPV,lastPV);

        pV_weight->Add(pV_DATA);
        pV_weight->Divide(pV_MC);
        pV_weight->Write();

        fout2->Close();
		
		cout << "# jets = " << nData << " # L = " << nDataL << " # R = " << nDataR << endl;
		cout << "B = " << nB << " NonB = " << nNonB << " B/nonB = " << (double)nB/(double)(nB+nNonB) << endl;
		cout << "BL = " << nBL << " NonBL = " << nNonBL << " BL/nonBL = " << (double)nBL/(double)(nBL+nNonBL) << " nonBL/BL = " << (double)nNonBL/(double)(nBL+nNonBL) << endl;
		cout << "BR = " << nBR << " NonBR = " << nNonBR << " BR/nonBE = " << (double)nBR/(double)(nBR+nNonBR) << " nonBR/BR = " << (double)nNonBR/(double)(nBR+nNonBR) << endl;
		
		cout << "CSA - B = " << nBCSA << " NonB = " << nNonBCSA << " nonB/B = " << (double)nNonBCSA/(double)(nBCSA+nNonBCSA) << endl;
		cout << "CS - B = " << nBCS << " NonB = " << nNonBCS << " nonB/B = " << (double)nNonBCS/(double)(nBCS+nNonBCS) << endl;
		cout << "CS - BL = " << nBLCS << " NonBL = " << nNonBLCS << " BL/nonBL = " << (double)nBLCS/(double)(nBLCS+nNonBLCS) << " nonBL/BL = " << (double)nNonBLCS/(double)(nBLCS+nNonBLCS) << endl;
		cout << "CS - BR = " << nBRCS << " NonBR = " << nNonBRCS << " BR/nonBE = " << (double)nBRCS/(double)(nBRCS+nNonBRCS) << " nonBR/BR = " << (double)nNonBRCS/(double)(nBRCS+nNonBRCS) << endl;
		
		cout << "F = " << nNonBL/nNonBR << endl;
		cout << "FCS = " << nNonBLCS/nNonBRCS << endl;
		cout << "#L = " << nBL+nNonBL << endl;
		cout << "#LCS = " << nBLCS+nNonBLCS << endl;
		cout << "#R = " << nBR+nNonBR << endl;
		cout << "#RCS = " << nBRCS+nNonBRCS << endl;
		
		cout << endl << "**** Mistag rate test (TCHEM) ****" << endl;
		
		cout << "eff b: " << tag_b/all_b << endl; 
		cout << "eff c: " << tag_c/all_c << endl; 
		cout << "eff dusg: " << tag_dusg/all_dusg << endl; 
        
        cout << endl << "**** Efficiencies ****" << endl;

        cout << "+> RefSel cut efficiency (MC based) " << endl;
        cout << " ++> nttbar before = " << nTTbarBeforeRefSel << endl;
        cout << " ++> nttbar after = " << nTTbarAfterRefSel << endl;
        cout << " ++> nttbar after2 = " << nTTbarAfterRefSel2 << endl;
        cout << " ++> refsel eff = " << (double)nTTbarAfterRefSel/(double)nTTbarBeforeRefSel  << endl;	
        
		
		exit(1);
	
	}
	
	std::map<int, vector<float> > rangesbTag;

	if (doVarBins) {		
		
		if (sampleType.find("nominal") != string::npos && !doPseudoExp_) {
		//if  (doVarBins) {
            
            for (unsigned int nt=0; nt < nTaggers; nt++) {
				
				double max = upRangeBtag[nt];
				double min = lowRangeBtag[nt];
				
				if (verbosity_ > 0) cout << "+--> Defining variable bins for tagger " << nt << " (nBins: " << nBinsBtag[nt] << " min: " << min << " max: "<< max << ") "<< endl;
				
                // OLD METHOD 
                
				/*std::vector<float> bTagVals;
				
				for (unsigned int s=0; s< v_bTag.size(); s++) 
					if (v_bTag[s][nt] > min && v_bTag[s][nt] <= max)
						bTagVals.push_back(v_bTag[s][nt]);
				
				//for (unsigned int s=0; s< v_bTag.size(); s++) 
				//	cout << v_bTag[s][nt] << " ";
				//c//out << endl;
				
				sort(bTagVals.begin(),bTagVals.end());
				
				int len = bTagVals.size();
				
				int nPerBin = len/nBinsBtag[nt];
				
				for (int i=0; i<bTagVals.size()-1; i++) {
					
					if (bTagVals[i] > bTagVals[i+1])
						cout << "ERROR " << i << " " << i+1 << endl;
					
				}
                				
				for (int i=0; i<nBinsBtag[nt]; i++) {
					
					float local_min = 0;
					
					if (i*nPerBin < len) {
						//cout << i*nPerBin << endl;
						
						local_min = bTagVals[i*nPerBin];
					}
					
					if (rangesbTag[nt].size() == 0 && local_min < min )
						rangesbTag[nt].push_back(min);
					else 
						rangesbTag[nt].push_back(local_min);
					
					
				}*/
                
                // NEW METHOD
                
                float binSize = (max-min)/nBinsBtag[nt];
                
                cout << binSize << endl;
                
                for (int i=0; i<nBinsBtag[nt]; i++) {
                
                    if (min+(i*binSize) < max)
                        rangesbTag[nt].push_back(min+(i*binSize));
                    else
                        rangesbTag[nt].push_back(max);
                }
                
                // PUT IN WPS
				
				vector<float> wps;
				
				for (unsigned int i=0; i<sizeof(wpArray)/sizeof(wpArray[0]); i++)
					if (taggerArray[i] == nt)
						wps.push_back(wpArray[i]);
                
                map<int,bool> listIndexClosest;
				
				for (unsigned int i=0;i<wps.size(); i++) {
					
					int indexClosest = 0;
					
					for (unsigned int j=0;j<rangesbTag[nt].size(); j++) {
						
						if (fabs(rangesbTag[nt][j]-wps[i]) < fabs(rangesbTag[nt][indexClosest]-wps[i])) {
							
                            if (listIndexClosest.find(j) == listIndexClosest.end()) {
                                indexClosest = j;   
                            } else {
                                //cout << j << endl;
                                if (wps[i] > rangesbTag[nt][j])
                                    indexClosest=j+1;
                                else if (wps[i] < rangesbTag[nt][j])
                                    indexClosest=j-1;
                            }
                        }
					}
					//cout << wps[i] <<  " " << indexClosest << " "<< binLowEdges[indexClosest] << endl;
					
                    //cout << indexClosest << " " << rangesbTag[nt][indexClosest] << endl;

                    if (indexClosest == nBinsBtag[nt]) {
                        
                        //cout<< "lol" << endl;
                        
                        nBinsBtag[nt] = nBinsBtag[nt]+1;
                        
                        //rangesbTag[nt][indexClosest+1] = rangesbTag[nt][indexClosest];
                        
                        rangesbTag[nt].push_back(wps[i]);

                        //rangesbTag[nt].push_back(wps[i]+(upRangeBtag[nt]-wps[i])/2);

                        listIndexClosest[indexClosest]=true;
                        
                    } else {
                        rangesbTag[nt][indexClosest] = wps[i];
                        
                        listIndexClosest[indexClosest]=true;                        
                    }
                    
                    //cout << indexClosest << " " << rangesbTag[nt][indexClosest] << endl;
					
					//cout << "BIS " << wps[i] <<  " " << indexClosest << " "<< binLowEdges[indexClosest] << endl;
					
				}
				
				rangesbTag[nt].push_back(max);
				
				stringstream t; t << nt;
				fstream bins(("binning/tagger_"+t.str()+"_channel"+data_postfix+".bins").c_str(),ios::out | ios::trunc);
				
				for (unsigned int i=0;i<rangesbTag[nt].size(); i++) {
					
					if (verbosity_ > 0) cout << rangesbTag[nt][i] << " ";
					
					bins << rangesbTag[nt][i] << "\n";
					
				}if (verbosity_ > 0) cout << endl;
				
				bins.close();
                				
			}
		} else {
			
			for (unsigned int nt=0; nt < nTaggers; nt++) {
				
				rangesbTag[nt].clear();
				
				stringstream t; 
                if (doPseudoExp_ && nBtag != -1)
                    t << nBtag;
                else
                    t << nt;
                
				string filename = "binning/tagger_"+t.str()+"_channel"+data_postfix+".bins";
                
				cout << "+--> Reading variable bins for tagger " << nt << " from " << filename << endl;
				//exit(1);
				fstream bins(filename.c_str(), ios::in);
				
				while(!bins.eof()) {
				
					double boundary = 0;
					
					bins >> boundary;
					
					if (rangesbTag[nt].size() <= nBinsBtag[nt])
						rangesbTag[nt].push_back(boundary);
					
				}
				
				bins.close();
				
				for (unsigned int i=0;i<rangesbTag[nt].size(); i++) {
					
					cout << rangesbTag[nt][i] << " ";
										
				}cout << endl;
				
				
			}				
		
			//exit(1);
			
		}
	}
    
    //exit(1);
        
    vector<double> pexp_left;
    vector<double> pexp_mid;
    vector<double> pexp_right;
    vector<double> pexp_f;
    
    if (doPseudoExp_ && doBiasStudy) {
        
        double step = bias_step;
        
        //step=10;
        
        for (double c=left_min; c<left_max; c+=step) {
            for (double a=mid_min; a<mid_max; a+=step) {
                for (double b=a+10; b<right_max; b+=step) {
                    cout << c << " " << a << " " << b << endl;
                    //for (int p=0; p<5;p++) {
                    pexp_left.push_back(c);
                    pexp_mid.push_back(a);
                    pexp_right.push_back(b);
                //}
                }
            }
        }
        
        if (step == 5) {
            
            nBinsVar0=100;
            lowRangeVar0=0;
            upRangeVar0=500;
            
        }
        
        cout << endl << pexp_f.size() << " -- " << pexp_left.size() << " " << pexp_mid.size() << " " << pexp_right.size() << endl;
        
        nPseudoExp_ = pexp_mid.size();
        
        //nPseudoExp_ = 10;
        
    }
    
    if (doPseudoExp_ && doFBiasStudy) {
        for (double a=F_min; a<F_max; a+=F_step) {
            pexp_f.push_back(a);
            cout << a << endl;
        }
        cout << endl << pexp_f.size() << " -- " << pexp_left.size() << " " << pexp_mid.size() << " " << pexp_right.size() << endl;
        
        nPseudoExp_ = pexp_f.size();

    }
    //exit(1);
    
	for(int iPE=0; iPE<nPseudoExp_; iPE++){
        
        //cout << "here" << endl;
        
        if (doPseudoExp_ && doBiasStudy) {
            leftlimit_ = pexp_left[iPE];
            centerleftlimit_ = pexp_mid[iPE];
            centerrightlimit_ = pexp_mid[iPE];
            rightlimit_ = pexp_right[iPE];
		}
        if (doPseudoExp_ && doFBiasStudy) {
            matchChiSquareCut_=pexp_f[iPE];
            //matchChiSquareCut_=100.0;
        }
        
		std::map<int,vector<double> > FitForXSResults;
		
		//initialize the skipEvent array
		/*for(int i=0; i<nrFiles; i++){
			skipEventsCounter[i]=0;
			skipEventsCounter10[i]=0;
			
			for(int j=0; j<500000; j++){
				skipEvents[i][j]=false;
			}
		}*/  
		
		fout->cd();
		
		//define a dynamic array
		PtEtaBinContainer *aContainer=NULL;
		
		aContainer = new PtEtaBinContainer(verbosity_,1,0,nTaggers,doVarBins,doShift,false,inBin_,outBin_,0);  
        
        aContainer->setLumi(desiredIntLum_);
		
		if(verbosity_>1) cout << "+--> Pseudo-Experiment " << iPE << " Created PtEtaBinContainer " << endl; 

		if (doVarBins) {
			aContainer->SetVarBins(rangesbTag);
			if (verbosity_>1) cout << "+-->Pseudo-Experiment " << iPE << " Defined variable bin-ranges for all btaggers" << endl;
		}
			
		std::map<int,vector<double> > WPMap;
		for (unsigned int iwp=0; iwp<nWP; iwp++) {
			
			WPMap[taggerArray[iwp]].push_back(wpArray[iwp]);						
		}
		

		TH1D* template_fit_data = new TH1D("data","data",nBinsVar0,lowRangeVar0,upRangeVar0);
		TH1D* template_fit_cs_data = new TH1D("data_cs","data_cs",nBinsVar0,lowRangeVar0,upRangeVar0);
		
		double global_shift = 0;
		
		for (int loop=-1; loop<2; loop++) {
			
			if (loop == -1) {
				
				//exit(1);
			
				aContainer->DefineSignalSamplePlots(nBinsVarLR[0],lowRangVarLR[0],upRangeVarLR[0],nBinsVarLR[1],lowRangVarLR[1],upRangeVarLR[1],nBinsVarSC[0],lowRangVarSC[0],upRangeVarSC[0],nBinsVarSC[1],lowRangVarSC[1],upRangeVarSC[1],nBinsVar0,lowRangeVar0,upRangeVar0,nBinsBtag,lowRangeBtag,upRangeBtag);
				aContainer->DefineControlSamplePlots(nBinsVarSC[0],lowRangVarSC[0],upRangeVarSC[0],nBinsVarSC[1],lowRangVarSC[1],upRangeVarSC[1],nBinsVar0,lowRangeVar0,upRangeVar0,nBinsBtag,lowRangeBtag,upRangeBtag);
				
				if(verbosity_>1) cout << "+--> Pseudo-Experiment " << iPE << " Defined Signal and Control sample plots" << endl;
				
				if(verbosity_>-1) cout << "+--> Pseudo-Experiment " << iPE << " Processing stored events" << endl; 
				
			}
			
			if (loop == 0) {
				
				
				if(verbosity_>-1) cout << "  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Filling Signal Sample, Control Sample and XS templates " << endl; 
				
				
				for (unsigned int n=0;n<v_weight.size();n++) {
					
					if (iPE == 0) {
						if (n==0)cout << "  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Setting up v_skipEvents vector" << endl;
						v_skipEvents.push_back(false);
					}

					//if(doPseudoExp_ && iPE == 0){ // FIXME REMOVE iPE == 0
                    
					if((doPseudoExp_ && !doBiasStudy && !doFBiasStudy) || (doPseudoExp_ && (doBiasStudy || doFBiasStudy) && iPE == 0)){ 
						
						random = rndm->Uniform(1); 
						
						if(v_weight[n]>1) { 
							cout << n << " (dataset: " << v_dataSetName[n] << ") " << v_weight[n] << endl;
							cout << "warning, the weight is higher than 1, the selection for the random sample will not be correct" << endl;
							exit(1);
						}
						v_skipEvents[n]=false;
						if(random>v_weight[n]) {
							v_skipEvents[n]=true;
						}
                        
					}
					
					/*random = rndm->Uniform(1);
					
					if (v_weight[n] > 1)
					cout << random << " " << v_weight[n] << endl;
					
					for (unsigned int a=0;a<v_skipEvents.size();a++)
						if (a < 500)
							v_skipEvents[a]=false;
						else 
							v_skipEvents[a]=true;*/

					if (v_skipEvents[n]) continue;
                    
                    string title = "PEXP_check_exp";
                    
                    if (doFBiasStudy && v_matchChiSquare[n] > matchChiSquareCut_) continue;
				
					//cout << "n: " << n << endl;
					
					/*ofstream chiSqfile;
					chiSqfile.open ("tmpcheck1.txt", ios::out | ios::app );
					
					// (n < 50000)
					chiSqfile << iPE << " " << n  << "\n";
					
					chiSqfile.close();*/
					
					
					double weight = v_weight[n]*v_ptweight[n];
					double weightC1 = v_weight[n]*v_ptweightControl[n];
					double weightC2 = v_weight[n]*v_ptweightControl2[n];
                    
                    //cout << v_weight[n] << " " << weight << " " << weightC1 << " " << weightC2 << endl;

					double weight_nonrew = v_weight[n];
					
					if (doPseudoExp_) {
                        
                        weight=1;
                        weight_nonrew=1;
                        //weight=v_scaleFactor[n];
                        //weight_nonrew=v_scaleFactor[n];
                        
                        weightC1=1;
                        weightC2=1;
                        
                    }
					// fill signal sample
                    
					aContainer->FillSignalSamplePlots(weight,weight_nonrew,v_partonFlavour[n],v_lepb_is[n][2],v_lepb_is[n][3],v_matchChiSquare[n],v_bTag[n],WPMap,v_pt[n],v_eta[n],v_var[n]+global_shift,v_varb[n],leftlimit_,centerleftlimit_,centerrightlimit_,rightlimit_);
					
					// fill control sample
					if (v_ptControl[n] > ptcutextra && v_ptControl2[n] > ptcutextra) {
							aContainer->FillControlSamplePlots(weightC1,v_partonFlavourControl[n],v_q1_is[n][2],v_q1_is[n][3],v_matchChiSquare[n],v_bTagCS1[n],v_ptControl[n],v_etaControl[n],v_controlVar[n]+global_shift,leftlimit_,centerleftlimit_,centerrightlimit_,rightlimit_);
							aContainer->FillControlSamplePlots(weightC2,v_partonFlavourControl2[n],v_q2_is[n][2],v_q2_is[n][3],v_matchChiSquare[n],v_bTagCS2[n],v_ptControl2[n],v_etaControl2[n],v_controlVar2[n]+global_shift,leftlimit_,centerleftlimit_,centerrightlimit_,rightlimit_);
					}
					//}
					
					/*	if (v_ptControl[n] > ptcutextra && v_bTagControl[n] < 100)
							aContainer->FillControlSamplePlots(weight,v_partonFlavourControl[n],v_q1_is[n][2],v_q1_is[n][3],v_matchChiSquare[n],v_ptControl[n],v_etaControl[n],v_controlVar[n],leftlimit_,centerleftlimit_,centerrightlimit_,rightlimit_);
						if (v_ptControl2[n] > ptcutextra && v_bTagControl2[n] < 100)
							aContainer->FillControlSamplePlots(weight,v_partonFlavourControl2[n],v_q2_is[n][2],v_q2_is[n][3],v_matchChiSquare[n],v_ptControl2[n],v_etaControl2[n],v_controlVar2[n],leftlimit_,centerleftlimit_,centerrightlimit_,rightlimit_);
					*/
					// fill templates for XS template fit
					
                    // JES TEST
                    
                    float extraw=1;
                    //if (v_dataSetName[n].find("TTbar") == string::npos) {
                    //    extraw = 1;
                    //}
                    
					aContainer->FillXStemplates(weight_nonrew*extraw,v_dataSetName[n],v_partonFlavour[n],v_bTag[n],WPMap,v_var[n],v_varb[n],leftlimit_,centerleftlimit_,centerrightlimit_,rightlimit_);
				
					if (v_dataSetName[n].find("TTbar") == 0)
						
						nTTbarAfterRefSel2+=weight;
					
				}
				
				
				aContainer->SetErrorsSignalSamples();
				aContainer->SetErrorsControlSamples();
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Set errors for 1D Var and 1D b-tag plots "  << endl;
				
				aContainer->MakeSoverSBPlots(); // plot 6.4 p 102
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Made S over S plus B plots "  << endl;
				
				aContainer->MakeProfileXplots(); 
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Made ProfileX plots "  << endl;
				
				aContainer->Make1DXplots(); 
				aContainer->SetError1DXplots();
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Made 1DX plots "  << endl;
				
				aContainer->Make1DYplots(); 
				aContainer->SetError1DYplots();
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Made 1DY plots "  << endl;
				
				aContainer->MakeMCEffPlots();
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Made MC eff plots "  << endl;
				
				aContainer->Make1DYVar2plots(); // obsolete
				aContainer->SetError1DYVar2plots(); // obsolete
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Made 1DY Var2 plots "  << endl;
				
				
				//in this function I fit all the bins:
				aContainer->MakeXRatioPlot(true); // pt-reweighting between left and right region plot 6.11 p 108
				
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Done fitting, now changing the parameters "  << endl;
				//in this function I change the fit parameters of the bins to that one of the global bin
				//aContainer->ChangeLeftRightPars();
				
				aContainer->Make2DRatioPlot(); // obsolete, doesn't work well to reweight both in pt and eta for left-right
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Made LR pt ratio and 2D ratio plots " << endl;
				
				aContainer->MakeSCprojectionPlots();
				aContainer->MakeSCVar12RatioPlot();
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Made SC pt ratio and 2D ratio plots " << endl;
				
				aContainer->ReweighRight(); // define plots to draw reweighted btag,... distributions for substraction in the enriched region. Call it defineReweightPlots? Also the SC reweight plots are defined
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Made Right Reweighted plots "  << endl;
                
                aContainer->ReweighLeft(); // define plots to draw reweighted btag,... distributions for substraction in the enriched region. Call it defineReweightPlots? Also the SC reweight plots are defined
				if(verbosity_>1) cout<<"  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Made Left Reweighted plots "  << endl;
			}
			
			else if (loop == 1) {
				
				if(verbosity_>-1) cout << "  +--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Running analysis" << endl; 
				
				//cout << "weight? " << weight << endl;
				
				for (unsigned int n=0;n<v_weight.size();n++) {

					if (v_skipEvents[n]) continue;
					
                    double weight = v_weight[n]*v_ptweight[n];
					double weightC1 = v_weight[n]*v_ptweightControl[n];
					double weightC2 = v_weight[n]*v_ptweightControl2[n];
                    
					//double weight_nonrew = v_weight[n];
					
                    if (doPseudoExp_) {
                        
                        weight=1;
                        weightC1=1;
                        weightC2=1;
                        
                        //weight=v_scaleFactor[n];
                    }					
					//cout << "2nd loop n: " << n << " " << weight<< endl;

					aContainer->FillReweighRight(do2D,useFit,weight,v_partonFlavour[n],v_bTag[n],v_pt[n],v_eta[n],v_var[n],leftlimit_,centerleftlimit_,centerrightlimit_,rightlimit_);//original
					
					//aContainer->FillReweighControl(do2Dcontrol, weight,partonFlavour,lepb_is[2],lepb_is[3],pt,eta,var,leftlimit_,centerleftlimit_,centerrightlimit_,rightlimit_,matchChiSquare); //this one is for the reweighing of the non-b mlj distribution in signal sample between different chi-square cuts
					double* bTagCuts = new double[3];
					
					bTagCuts[0]=1.7; bTagCuts[1]=3.3; bTagCuts[2]=10.2; // VERY BAD!!!!!!
					if (v_ptControl[n] > ptcutextra && v_ptControl2[n] > ptcutextra) {
						aContainer->FillReweighControl(v_bTagCS1[n],bTagCuts,do2Dcontrol, weightC1,v_partonFlavourControl[n],v_q1_is[n][2],v_q1_is[n][3],v_ptControl[n],v_etaControl[n],v_controlVar[n],leftlimit_,centerleftlimit_,centerrightlimit_,rightlimit_,v_matchChiSquare[n]);
						aContainer->FillReweighControl(v_bTagCS2[n],bTagCuts,do2Dcontrol, weightC2,v_partonFlavourControl2[n],v_q2_is[n][2],v_q2_is[n][3],v_ptControl2[n],v_etaControl2[n],v_controlVar2[n],leftlimit_,centerleftlimit_,centerrightlimit_,rightlimit_,v_matchChiSquare[n]);
						
					}
				}
				
				aContainer->GetLRratio(FMCBias,doNewF,newFpt,newFeta);
				if(verbosity_>1) cout<<"+--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Calculated the LR ratio "  << endl;
				aContainer->MeasureEff(doSCreweigh); // not reweighted
				if(verbosity_>1) cout<<"+--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Measured the efficiency "  << endl;
				aContainer->MakeReweighRatio(); 
				if(verbosity_>1) cout<<"+--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Made reweigh ratio"  << endl;
				aContainer->MeasureEffRR(doSCreweigh); // RR stands for Right Reweighted
				if(verbosity_>1) cout<<"+--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Measured the efficiency Right reweighted"  << endl;
                aContainer->MeasureMistagEffRR(doSCreweigh); // RR stands for Right Reweighted
				if(verbosity_>1) cout<<"+--> (L=" << desiredIntLum_ << ") Pseudo-Experiment " << iPE << " Measured the mistag efficiency Right reweighted"  << endl;
				
				stringstream matchChiSquareCutSTR; matchChiSquareCutSTR << matchChiSquareCut_;
				
				if (doVarBins && !btagonly)
					FitForXSResults = aContainer->doMLJTemplateFit(matchChiSquareCutSTR.str(),fitMode,data_postfix,nSystematic);
				else {
				
					for (unsigned int m=0;m<20;m++)
						for (unsigned int n=0;n<20;n++)
							FitForXSResults[m].push_back(n);	
				
				}
								
			}
			
		} // end of 'loop' loop
		
		//exit(1);
		
		stringstream matchChiSquareCutSTR; matchChiSquareCutSTR << matchChiSquareCut_;
		
		double chiSQCutEff;
		double refSelEff;
		double MLBCutEff;
		
		double chiSQCutEff_sample = (double)nTTbarAfterChiSq/(double)nTTbarBeforeChiSq;
		double MLBCutEff_sample = (double)nTTbarAfterMLBCUT/(double)nTTbarBeforeMLBCUT;		
		double refSelEff_sample = (double)nTTbarAfterRefSel/(double)nTTbarBeforeRefSel;
        
        /*cout << nTTbarBeforeRefSel << endl;
        cout << nTTbarAfterRefSel << endl;
        cout << nTTbarBeforeChiSq << endl;
        cout << nTTbarAfterChiSq << endl;
        cout << nTTbarAfterMLBCUT << endl;
        cout << nTTbarBeforeMLBCUT << endl;
        
        
        exit(1);
		*/
        
		if (sampleType.find("nominal") != string::npos) {
			
			ofstream effs; effs.open(("efficiencies_channel"+data_postfix+".nominal").c_str(), ios::out | ios::trunc );
			
			effs << refSelEff_sample << " " << chiSQCutEff_sample << " " << MLBCutEff_sample;
			
			effs.close();
			
		}			
		
		if (sampleType.find("nominal") == string::npos) {
			
			double chiSQCutEff_nom = 0;
			double MLBCutEff_nom = 0;		
			double refSelEff_nom = 0;
			
			ifstream effs; effs.open(("efficiencies_channel"+data_postfix+".nominal").c_str(), ios::in );
			
			while (!effs.eof())
				effs >> refSelEff_nom >> chiSQCutEff_nom >> MLBCutEff_nom;
			
			cout << refSelEff_nom << " " << chiSQCutEff_nom << " " << MLBCutEff_nom << endl;
												
			effs.close();
			
			chiSQCutEff = chiSQCutEff_nom; 
			MLBCutEff = MLBCutEff_nom;
			refSelEff = refSelEff_nom;
            
            //chiSQCutEff = chiSQCutEff_sample;
			//MLBCutEff = MLBCutEff_sample;
			//refSelEff = refSelEff_sample;

			
			//if (datasetName == "Data") // FIXME : get new value
                //refSelEff=refSelEff*0.98544968;// lepton ID and muon eff SF
                //refSelEff=refSelEff*0.98;// lepton ID and muon eff SF
			
			
		} else {
			
			chiSQCutEff = chiSQCutEff_sample;
			MLBCutEff = MLBCutEff_sample;
			refSelEff = refSelEff_sample;
			
		}
        
        /*if (datasetName == "Data") {
			// old
			
			//spring11
			//chiSQCutEff = 0.894470;
			//refSelEff = 0.0178175;
			
			//summer11
			//chiSQCutEff = 0.9113*0.994836; // mupt > 20
			//refSelEff = 0.0253855;
			//chiSQCutEff = 0.909642*0.992768; // mupt > 35
			//refSelEff = 0.0181441;
			
			//refSelEff=refSelEff*0.98544968; // lepton ID and muon eff SF
			
			// TEMPORARY!
			if (nSystematic >= 7 && nSystematic <=13) { // Fall10 values for ttjets systematic samples
             chiSQCutEff = 0.930668;
             refSelEff = 0.0261852;
             }
             
             if (nSystematic >= 19 && nSystematic <=20) { // Spring11 values for ttjets systematic samples
             chiSQCutEff = 0.91274;
             refSelEff = 0.0235167;
             }
            
            //refSelEff=0.0272082+(0.0272082*1.0*0.01);
            //chiSQCutEff=1;
            //MLBCutEff=1;
            
		}*/
        
        //if (decay == 1 && sampleType.find("Data") != string::npos)
        //    refSelEff=refSelEff*0.9341*1.004;
		
		if (!doPseudoExp_) {
			
			std::map<int,int> taggerCounter; // first int nTagger, second int nTimesTagger
			
			std::map<int,std::vector<float> > infoForPDF,InfoForCombinedFit;
			
			int lastTagger=0;
			int nTimesWP=1;
			
			for(int iwp=0; iwp < nWP; iwp++){

				if (verbosity_ > 1) cout << "+--------------------------------------------------------------------------" << endl;
				if (verbosity_ > 1) cout << "+---- MLJ template fit results (Tagger: " << taggerArray[iwp] << " WP: " << wpArray[iwp] << ")" << endl;
				if (verbosity_ > 1) cout << "+--------------------------------------------------------------------------" << endl;
				
				double* results = new double[30];
				
				taggerCounter[taggerArray[iwp]]++;
				
				//cout << taggerCounter[taggerArray[iwp]] << endl;
                				
				aContainer->GetWPEff(true,doSCreweigh,wpArray[iwp],taggerArray[iwp],results,doPseudoExp_,true,0,0,false,runNb);
  
                //results[8]=0.605;
                //results[9]=0.017;
                
                /* SINCE BTAG IS DATA_DRIVEN we want to fix the beff for theory systematics to the nominal value */
                
                /*if (sampleType.find("Data") != string::npos) {
                
                    setNomVal(data_postfix,fitMode,"DATA_eff_meas_fdata",iwp,results[8]);
                    setNomVal(data_postfix,fitMode,"DATA_ueff_meas_fdata",iwp,results[9]);
                    
                    if (data_postfix=="_El") {
                        results[8] = getNomVal("_Mu","DATA_eff_meas_fdata",iwp);
                        results[9] = getNomVal("_Mu","DATA_ueff_meas_fdata",iwp);
                    }
                    
                    //results[8] = getNomVal(data_postfix,fitMode,"eff_true",iwp);
                    //results[9] = getNomVal(data_postfix,fitMode,"ueff_true",iwp);
                
                }*/
                
                if (sampleType.find("nominal") != string::npos) {
                 
                    setNomVal(data_postfix,fitMode,"eff_true",iwp,results[0]);
                    setNomVal(data_postfix,fitMode,"ueff_true",iwp,results[1]);
                    
                    setNomVal(data_postfix,fitMode,"eff_meas_fmc",iwp,results[6]);
                    setNomVal(data_postfix,fitMode,"ueff_meas_fmc",iwp,results[7]);
                    
                    setNomVal(data_postfix,fitMode,"eff_meas_fdata",iwp,results[8]);
                    setNomVal(data_postfix,fitMode,"ueff_meas_fdata",iwp,results[9]);
                    
                } else if (( nSystematic > 6 && nSystematic < 21 ) || nPDFWeight > -1) {
                //} else if (nSystematic > 0) {
                    
                    //cout << "i'm here" << endl; exit(1);
                
                    results[0] = getNomVal(data_postfix,fitMode,"eff_true",iwp);
                    results[1] = getNomVal(data_postfix,fitMode,"ueff_true",iwp);
                    
                    results[6] = getNomVal(data_postfix,fitMode,"eff_meas_fmc",iwp);
                    results[7] = getNomVal(data_postfix,fitMode,"ueff_meas_fmc",iwp);
                    
                    results[8] = getNomVal(data_postfix,fitMode,"eff_meas_fdata",iwp);
                    results[9] = getNomVal(data_postfix,fitMode,"ueff_meas_fdata",iwp);
                    
                
                } else if (datasetName == "Data" && doBtagSFonData) {
                    
                    cout << "Never set doBtagSFonData to true" << endl;
                    exit(1);
                    double sf = results[8]/getNomVal(data_postfix,fitMode,"eff_meas_fdata",iwp);
                    
                    setNomVal(data_postfix,fitMode,"sf_for_data",iwp,sf);
                    
                    cout << "Applying measured btag SF for data of " << sf << endl;
                    
                    results[8] = results[8]/sf;

                    // after the calcXS calls below, the results[8] value is put back to it's original
                }
                
                /* for bias correction in BTV-11-003 SF's */
                std::pair<float,float> biasforcorr = calcBias(results[0],results[1],results[8],results[9]);
                                   
                if (doCorrectForBias && sampleType.find("nominal") != string::npos) {
                
                    if (verbosity_ > 0) cout << "Storing measured btag eff bias of " << round(100,biasforcorr.first*100) << "+-" << round(100,biasforcorr.second*100) << endl;
                    
                    setNomVal(data_postfix,fitMode,"biascorr",iwp,biasforcorr.first);
                    
                }
                
                if (datasetName == "Data" && doCorrectForBias) {
                    
                    stringstream p; p<<iwp;
                    double bias = getNomVal(data_postfix,fitMode,"biascorr",p.str());
                    cout << "Correcting data measured btag eff for known bias of " << bias << endl;
                    cout << "E_b = " << results[8] << endl;
                    if (bias > 0)
                        results[8] = results[8]+fabs((results[8]*bias));
                    else
                        results[8] = results[8]-fabs((results[8]*bias));
                    
                    cout << "E_bCorr = " << results[8] << endl;
                    
				}
                //exit(1);
                
				//for (unsigned int i=0; i<FitForXSResults[taggerArray[iwp]].size(); i+=1)
				//	cout << "FitForXSResults[" << taggerArray[iwp] <<"][" << i << "] = " << FitForXSResults[taggerArray[iwp]][i] << endl;
				
				vector<double> XSresults;

				float lumiToCalc = 0;
				if (datasetName =="Data")
					lumiToCalc = origLumi;
				else
					lumiToCalc = desiredIntLum_;
				
				if (iwp==0 || taggerArray[iwp] != lastTagger) {
					cout << "++> CALCULATING XS FOR NOMINAL FIT RESULTS (no btag cut)" << endl;
					XSresults = calcXS(doPseudoExp_,1,0,chiSQCutEff,MLBCutEff,refSelEff,FitForXSResults[taggerArray[iwp]],lumiToCalc); 
					
					cout << " => SigmaTTbar = " << XSresults[1] << " +- " << XSresults[2] << "(stat) pb" <<  endl << endl << endl;
					//cout << "+> TTbar Cross Section with TTbar MC and V data templates" << endl;
					//cout << "SigmaTTbar = " << XSresults[1] << endl;cout << endl;
					
					infoForPDF[iwp].push_back(round(10,XSresults[1]));
					infoForPDF[iwp].push_back(round(10,XSresults[2]));
					
				} else if (iwp != 0) {
					
					infoForPDF[iwp].push_back(infoForPDF[iwp-1][0]);
					infoForPDF[iwp].push_back(infoForPDF[iwp-1][1]);
				}
				
				if (taggerArray[iwp] != lastTagger) {
					lastTagger=taggerArray[iwp];
					nTimesWP=1;
				}
				
                FitForXSResults[taggerArray[iwp]][((nTimesWP)*5)+4]=results[24];
				cout << "++> CALCULATING XS WITH BTAG CUT AND Btag EFF Fdata (using array offset: "<< (nTimesWP)*5 << ")" << endl;
				//XSresults = calcXS(doPseudoExp_,results[8],FitForXSResults[taggerArray[iwp]],lumiToCalc,(nTimesWP)*5); // 16 entries for each cut
				XSresults = calcXS(doPseudoExp_,results[8],results[9],chiSQCutEff,MLBCutEff,refSelEff,FitForXSResults[taggerArray[iwp]],lumiToCalc,(nTimesWP)*5); // 16 entries for each cut
                
                // if applying the SF for data now revert to original value
                
                /*if (datasetName == "Data" && doBtagSFonData) {
                                        
                    double sf = getNomVal(data_postfix,fitMode,"sf_for_data",iwp);
                    
                    cout << "Removing measured btag SF for data of " << sf << endl;
                    
                    results[8] = results[8]*sf;
                    
                    // after the calcXS calls below, the results[8] value is put back to it's original
                }*/ // DISABLED FOR SECURITY


				stringstream s; s << iwp;
			
				string title = "2DMeasurement_FDATA_XS_btag_IWP_"+s.str();
				
				string title2 = "ALT2DMeasurement_FDATA_XS_btag_IWP_"+s.str();
				if (histos1D.find(title2) == histos1D.end()) {
					histos1D[title2] = new TH1D(title2.c_str(),title2.c_str(),200,125,325);
				}
				
				double contours[3];
				contours[0] = 1; //1 sigma
				contours[1] = 4; //2 sigma
				
                int b=50;
                double bmin=0.5;
                double bmax=1;
                

                string wpstr = "";
                if (iwp%3==0) {
                    b=50;
                    bmin=0.5;
                    bmax=1;
                }else if ( (iwp-1) %3==0) {
                    b=60;
                    bmin=0.4;
                    bmax=1;
                    //exit(0);
                }else if ( (iwp-2) %3==0){
                    b=100;
                    bmin=0.0;
                    bmax=1;
                }
                                
                cout << "AAA " << results[8] << " " << bmin << endl;
				if (histos2D.find(title) == histos2D.end()) {
					histos2D[title] = new TH2D(("ellipse_"+title).c_str(),(title+";#hat{#epsilon}_{b};#hat{#sigma}_{t#bar{t}} (pb)").c_str(),b,bmin,bmax,200,125,325);
					histos2D["contour_"+title] = new TH2D(("contour_ellipse_"+title).c_str(),(title+";#hat{#epsilon}_{b};#hat{#sigma}_{t#bar{t}} (pb)").c_str(),b,bmin,bmax,200,125,325);
					histos2D["contour_"+title]->SetContour(2,contours);
				}
				
				if (histos2D.find(title2) == histos2D.end()) {
					histos2D[title2] = new TH2D(("ellipse_"+title2).c_str(),(title2+";#hat{#sigma}_{t#bar{t}} (pb);#hat{#epsilon}_{b}").c_str(),200,125,325,100,0,1);
					histos2D["contour_"+title2] = new TH2D(("contour_ellipse_"+title2).c_str(),(title2+";#hat{#sigma}_{t#bar{t}} (pb);#hat{#epsilon}_{b}").c_str(),200,125,325,100,0,1);
					histos2D["contour_"+title2]->SetContour(2,contours);
				}
				
				//results[9]=results[9]*2; // just to test the alt method
				
				double gx[1]; double egx[1];
				double gy[1]; double egy[1];
				
				gx[0] = results[8]; egx[0] = results[9];
				gy[0] = XSresults[1]; egy[0] = XSresults[2];
				
				cout << gx[0] << "+-" << egx[0] << endl;
				cout << gy[0] << "+-" << egy[0] << endl;
				
				///exit(1);
				
				if (graphs.find(title) == graphs.end()) {
					
					graphs[title] = new TGraphErrors(1,gx,gy,egx,egy);
					graphs[title]->SetNameTitle(("dataresult"+title).c_str(),(title+";#epsilon_{b};#sigma_{t#bar{t}} (pb)").c_str());
					
				}
				
				if (graphs.find(title2) == graphs.end()) {
					
					graphs[title2] = new TGraphErrors(1,gy,gx,egy,egx);
					graphs[title2]->SetNameTitle(("dataresult"+title2).c_str(),(title2+";#sigma_{t#bar{t}} (pb);#epsilon_{b}").c_str());
					
				}
				
				if (histos1D.find(title) == histos1D.end()) {
					histos1D[title] = new TH1D(title.c_str(),title.c_str(),(int)((2*egx[0])/0.005),gx[0]-egx[0],gx[0]+egx[0]);
				}
				
				for (unsigned int i=0; i<histos1D[title]->GetNbinsX()+1; i++) {
					
					// first we calc the horizontal band
					
					//if (histos1D[title]->GetBinLowEdge(i) > results[8]+results[9] || histos1D[title]->GetBinLowEdge(i-1) < results[8]-results[9]) continue;
					
					vector<double> tmpRes = calcXS(true,histos1D[title]->GetXaxis()->GetBinCenter(i),results[9],chiSQCutEff,MLBCutEff,refSelEff,FitForXSResults[taggerArray[iwp]],lumiToCalc,(nTimesWP)*5);
					
					histos1D[title]->SetBinContent(i,tmpRes[1]);
					
					//cout << histos1D[title]->GetXaxis()->GetBinCenter(i) << " +- " << results[9] << " ==> " << tmpRes[1] << " +- " << tmpRes[3] << endl;
					
				}
				
				for (unsigned int i=1; i<histos1D[title2]->GetNbinsX()+1; i++) {
					
					// first we calc the horizontal band
					
					if (histos1D[title2]->GetBinLowEdge(i) > XSresults[1]+XSresults[2] || histos1D[title2]->GetBinLowEdge(i-1) < XSresults[1]-XSresults[2]) continue;
					
					double tmpRes = calcBtag(histos1D[title2]->GetXaxis()->GetBinCenter(i),chiSQCutEff,refSelEff,FitForXSResults[taggerArray[iwp]],lumiToCalc,(nTimesWP)*5);
					
					histos1D[title2]->SetBinContent(i,tmpRes);
					
					//cout << tmpRes[1] << " +- " << tmpRes[3] << endl;
					
				}
				
				TF1* fit = new TF1("fit","pol1",histos1D[title]->GetBinLowEdge(1),histos1D[title]->GetBinLowEdge(histos1D[title]->GetNbinsX()+1));
				TF1* fit2 = new TF1("fit2","pol1",histos1D[title2]->GetBinLowEdge(1),histos1D[title2]->GetBinLowEdge(histos1D[title2]->GetNbinsX()+1));
				
				histos1D[title]->Fit(fit,"RQ");
				histos1D[title2]->Fit(fit2,"RQ");
				
				//cout << "IWP " << iwp << " fit pars: " << fit->GetParameter(0) << " " << fit->GetParameter(1) << endl;
				
				for (unsigned int i=1; i<histos2D[title]->GetXaxis()->GetNbins()+1; i++) {
					
					//cout << histos1D[title]->GetXaxis()->GetBinCenter(i) << " " << fit->Eval(histos1D[title]->GetXaxis()->GetBinCenter(i)) << endl;
					
					for (unsigned int j=1; j<histos2D[title]->GetYaxis()->GetNbins()+1; j++) {

						float XS = histos2D[title]->GetYaxis()->GetBinCenter(j);
						float btag = histos2D[title]->GetXaxis()->GetBinCenter(i);
						
						float sigma_xs = pow((XS-fit->Eval(btag))/XSresults[2],2);
						float sigma_btag = pow((btag-results[8])/results[9],2);
						
						float sigma = sigma_btag+sigma_xs;
						
						if (sigma <= 4) {
							histos2D["contour_"+title]->Fill(btag,XS,sigma);
                            //cout << btag << " " << XS << " " << sigma << endl;
                        }
						histos2D[title]->Fill(btag,XS,sigma);
						
					}
					
				}
				
				for (unsigned int i=1; i<histos2D[title2]->GetXaxis()->GetNbins()+1; i++) {
					for (unsigned int j=1; j<histos2D[title2]->GetYaxis()->GetNbins()+1; j++) {
						
						float XS = histos2D[title2]->GetXaxis()->GetBinCenter(i);
						float btag = histos2D[title2]->GetYaxis()->GetBinCenter(j);
												
						float sigma_xs = pow((XS-XSresults[1])/XSresults[2],2);
						float sigma_btag = pow((btag-fit2->Eval(XS))/results[9],2);
						
						//cout << XS << " -> " << fit2->Eval(XS) << endl;
						
						float sigma = sigma_btag+sigma_xs;
						
						if (sigma <= 4)
							histos2D["contour_"+title2]->Fill(XS,btag,sigma);
						histos2D[title2]->Fill(XS,btag,sigma);
						
					}
										
				}

				// make contour plot 
				
				if (doVarBins) {
					
					TCanvas *c1 = new TCanvas("c1","results");
					TPad    *z1 = new TPad("z1","z1", 0., 0., 1, 1); z1->Draw();
					TPad    *z2 = new TPad("z2","z2", 0., 0., 1, 1); z2->Draw();
					z2->SetFillStyle(4000); // z2 in transparent
					
					histos2D[title]->SetTitle("");
					histos2D["contour_"+title]->SetTitle("");
					
					z1->cd();
					histos2D[title]->Draw("cont4");
                    histos2D[title]->GetYaxis()->SetTitleOffset(1.3);

					
					z2->cd();
					//z2->Range(0.851562,-0.000616667,2.09844,0.0115167);
					histos2D["contour_"+title]->SetLineColor(kWhite);
					histos2D["contour_"+title]->SetLineWidth(4.);
					histos2D["contour_"+title]->Draw("cont3");
					graphs[title]->SetMarkerColor(kRed);
					graphs[title]->SetLineColor(kRed);
					graphs[title]->SetLineWidth(2);
					graphs[title]->Draw("same P");
					
					float XSmin = 99999999;
					float XSmax = 0;
					
					float effmin=0;
					float effmax=0;
					
					for (unsigned int x=0; x <histos2D["contour_"+title]->GetXaxis()->GetNbins(); x++) {
						for (unsigned int y=0; y <histos2D["contour_"+title]->GetYaxis()->GetNbins(); y++) {
							
							if (fabs(histos2D["contour_"+title]->GetBinContent(x,y)-1) < 0.05) {
								
								if (histos2D["contour_"+title]->GetYaxis()->GetBinCenter(y) > XSmax) {
									XSmax = histos2D["contour_"+title]->GetYaxis()->GetBinCenter(y);
									effmax = histos2D["contour_"+title]->GetXaxis()->GetBinCenter(x);
								}
								if (histos2D["contour_"+title]->GetYaxis()->GetBinCenter(y) < XSmin) {
									XSmin = histos2D["contour_"+title]->GetYaxis()->GetBinCenter(y);
									effmin = histos2D["contour_"+title]->GetXaxis()->GetBinCenter(x);
								}
								
								//cout << "XS " << histos2D["contour_"+title]->GetYaxis()->GetBinCenter(y) << " dCHi2 " << histos2D["contour_"+title]->GetBinContent(x,y) << endl;
								
								
							}
						}
					}
					
					cout << "Cross section uncertainty interval: [" << XSmin <<","<< XSmax << "]" << endl;
					
					/*TLine* line = new TLine(0,XSmax,effmax,XSmax);
					line->SetLineColor(kGreen);
					line->SetLineWidth(1);
					line->Draw("same");
					
					TLine* line2 = new TLine(0,XSmin,effmin,XSmin);
					line2->SetLineWidth(1);
					line2->SetLineColor(kGreen);
					line2->Draw("same");
					*/
					
					//histos1D[title]->Draw("same");
					
                    if (nSystematic == -2) {
                        stringstream slum; slum << round(10,desiredIntLum_/1000);
                        /*TLatex* text = new TLatex(0.13,0.96,("CMS Preliminary, #sqrt{s} = 8 TeV, #int Ldt = "+slum.str()+" fb^{-1}").c_str());
                         text->SetNDC();
                         //text->SetTextFont(42);
                         text->SetTextSize(0.04);
                         
                         text->Draw();*/
                        // Draw the CMS line
                        
                        TLatex* latex = new TLatex();
                        latex->SetNDC();
                        latex->SetTextSize(0.04);
                        latex->SetTextAlign(31); // align right
                        latex->DrawLatex(0.45, 0.95, "CMS Preliminary");
                        
                        TLatex* latex2 = new TLatex();
                        latex2->SetNDC();
                        latex2->SetTextSize(0.04);
                        latex2->SetTextAlign(31); // align right
                        latex2->DrawLatex(0.87, 0.95, (slum.str() + " fb^{-1} at #sqrt{s} = 8 TeV").c_str());
					} else {
                        
                        TLatex* text = new TLatex(0.13,0.955,"CMS Simulation");
                        
                        text->SetTextSize(0.05);
                        
                        text->SetNDC();
                        
                        text->Draw(); 
                        
                    }
					//fit->Draw("same");
					
					mkdir("2DMeasurement",0777);
					
					histos2D[title]->GetYaxis()->SetTitleOffset(0.8);
					histos2D[title]->GetYaxis()->SetLabelSize(0.04);
					histos2D["contour_"+title]->GetYaxis()->SetTitleOffset(0.8);
					histos2D["contour_"+title]->GetYaxis()->SetLabelSize(0.04);
					
					histos2D[title]->GetXaxis()->SetTitleOffset(0.8);
					histos2D[title]->GetXaxis()->SetLabelSize(0.04);
					histos2D["contour_"+title]->GetXaxis()->SetTitleOffset(0.8);
					histos2D["contour_"+title]->GetXaxis()->SetLabelSize(0.04);

					c1->Update();
					
					TPaveStats* sb = (TPaveStats*)histos2D[title]->GetListOfFunctions()->FindObject("stats");
					sb->SetTextColor(kWhite);
					sb->SetLineColor(kWhite);
					sb->SetX1NDC(.99);
					sb->SetX2NDC(1.);
					sb->SetY1NDC(0.99);
					sb->SetY2NDC(1.);
					
					c1->Draw();
					
					c1->Update();
					
					sb = (TPaveStats*)histos2D["contour_"+title]->GetListOfFunctions()->FindObject("stats");
					sb->SetTextColor(kWhite);
					sb->SetLineColor(kWhite);
					sb->SetX1NDC(.99);
					sb->SetX2NDC(1.);
					sb->SetY1NDC(0.99);
					sb->SetY2NDC(1.);
					
					c1->Draw();
					
					
					/*if (iwp == 1) {
						TFile* f2D = new TFile("2DMeasurement/2D_Data_TCHEM.root","RECREATE");
						f2D->cd();
						c1->Write();
						f2D->Close();
						
					}*/
					stringstream fitM; fitM<<fitMode;
                    
					c1->SaveAs(("2DMeasurement/"+title+"_channel"+data_postfix+"_fitMode"+fitM.str()+".png").c_str());
					c1->SaveAs(("2DMeasurement/"+title+"_channel"+data_postfix+"_fitMode"+fitM.str()+".pdf").c_str());
					c1->SaveAs(("2DMeasurement/"+title+"_channel"+data_postfix+"_fitMode"+fitM.str()+".root").c_str());
					
					//exit(1);
					
					// make contour plot ALT METHOD
					
					TCanvas *c2 = new TCanvas("c2","results");
					TPad    *z3 = new TPad("z3","z3", 0., 0., 1, 1); z3->Draw();
					TPad    *z4 = new TPad("z4","z4", 0., 0., 1, 1); z4->Draw();
					z4->SetFillStyle(4000); // z3 in transparent
					
					histos2D[title2]->SetTitle("");
					histos2D["contour_"+title2]->SetTitle("");
					
					z3->cd();
					histos2D[title2]->Draw("cont4");
					
					z4->cd();
					//z2->Range(0.851562,-0.000616667,2.09844,0.0115167);
					histos2D["contour_"+title2]->Draw("cont3");
					graphs[title2]->Draw("same P");
					
					//fit->Draw("same");
					
					
					histos2D[title2]->GetYaxis()->SetTitleOffset(0.8);
					histos2D["contour_"+title2]->GetYaxis()->SetTitleOffset(0.8);
					histos2D[title2]->GetXaxis()->SetTitleOffset(0.8);
					histos2D["contour_"+title2]->GetXaxis()->SetTitleOffset(0.8);
					
					//c2->SaveAs(("2DMeasurement/"+title2+".C").c_str());
					//c2->SaveAs(("2DMeasurement/"+title2+".png").c_str());
					
					//delete z1;
					//delete z2;
					delete c1;
					delete z3;
					delete z4;
					delete c2;
					
					delete fit;
					delete fit2;
				}
				
				cout << "+> TTbar Cross Section with TTbar+V MC templates" << endl;
				cout << " => SigmaTTbar = " << XSresults[1] << " +- " << XSresults[2] << "(stat) pb" <<  endl << endl << endl;
				//cout << "+> TTbar Cross Section with TTbar MC and V data templates" << endl;
				//cout << "SigmaTTbar = " << XSresults[1] << endl;cout << endl;
				
				infoForPDF[iwp].push_back(round(10,XSresults[1]));
				infoForPDF[iwp].push_back(round(10,XSresults[2]));
				infoForPDF[iwp].push_back(round(10,results[8]*100));
				infoForPDF[iwp].push_back(round(10,results[9]*100));
				infoForPDF[iwp].push_back(round(1000,results[14]));
				infoForPDF[iwp].push_back(round(1000,results[15]));
                
                InfoForCombinedFit[iwp].push_back(results[8]);
                InfoForCombinedFit[iwp].push_back(results[9]);
				
				cout << "++> CALCULATING XS WITH BTAG CUT AND Btag EFF Fmc (using array offset: "<< (nTimesWP)*5 << ")" << endl;
				//XSresults = calcXS(doPseudoExp_,results[8],FitForXSResults[taggerArray[iwp]],lumiToCalc,(nTimesWP)*5); // 16 entries for each cut
				XSresults = calcXS(doPseudoExp_,results[6],results[7],chiSQCutEff,MLBCutEff,refSelEff,FitForXSResults[taggerArray[iwp]],lumiToCalc,(nTimesWP)*5); // 16 entries for each cut
				
				cout << "+> TTbar Cross Section with TTbar+V MC templates" << endl;
				cout << " => SigmaTTbar = " << XSresults[1] << " +- " << XSresults[2] << "(stat) pb" <<  endl << endl << endl;
				//cout << "+> TTbar Cross Section with TTbar MC and V data templates" << endl;
				//cout << "SigmaTTbar = " << XSresults[1] << endl;cout << endl;
				
				infoForPDF[iwp].push_back(round(10,XSresults[1]));
				infoForPDF[iwp].push_back(round(10,XSresults[2]));
				infoForPDF[iwp].push_back(round(10,results[6]*100));
				infoForPDF[iwp].push_back(round(10,results[7]*100));
				infoForPDF[iwp].push_back(round(1000,results[10]));
				infoForPDF[iwp].push_back(round(1000,results[11]));
				
				
				cout << "++> CALCULATING XS WITH BTAG CUT AND Btag EFF MCTRUTH (using array offset: "<< (nTimesWP)*5 << ")" << endl;
				//XSresults = calcXS(doPseudoExp_,results[8],FitForXSResults[taggerArray[iwp]],lumiToCalc,(nTimesWP)*5); // 16 entries for each cut
				XSresults = calcXS(doPseudoExp_,results[0],results[1],chiSQCutEff,MLBCutEff,refSelEff,FitForXSResults[taggerArray[iwp]],lumiToCalc,(nTimesWP)*5); // 5 entries for each cut
				
				cout << "+> TTbar Cross Section with TTbar+V MC templates" << endl;
				cout << " => SigmaTTbar = " << XSresults[1] << " +- " << XSresults[2] << "(stat) pb" <<  endl << endl << endl;
				//cout << "+> TTbar Cross Section with TTbar MC and V data templates" << endl;
				//cout << "SigmaTTbar = " << XSresults[1] << endl;cout << endl;
				
				infoForPDF[iwp].push_back(round(10,XSresults[1]));
				infoForPDF[iwp].push_back(round(10,XSresults[2]));
				infoForPDF[iwp].push_back(round(10,results[0]*100));
				infoForPDF[iwp].push_back(round(10,results[1]*100));
				
				infoForPDF[iwp].push_back(round(10,FitForXSResults[taggerArray[iwp]][((nTimesWP)*5)+4]*100));
				infoForPDF[iwp].push_back(round(10,FitForXSResults[taggerArray[iwp]][FitForXSResults[taggerArray[iwp]].size()-1]*100));
				
                infoForPDF[iwp].push_back(round(10,results[22]*100)); // true mistag
				infoForPDF[iwp].push_back(round(10,results[23]*100));
                infoForPDF[iwp].push_back(round(10,results[24]*100)); // measured mistag
				infoForPDF[iwp].push_back(round(10,results[25]*100));

                
                InfoForCombinedFit[iwp].push_back(FitForXSResults[taggerArray[iwp]][((nTimesWP)*5)+4]);
				InfoForCombinedFit[iwp].push_back(FitForXSResults[taggerArray[iwp]][FitForXSResults[taggerArray[iwp]].size()-1]);
				
				
				//cout << FitForXSResults[taggerArray[iwp]][((nTimesWP)*5)+4] <<" " << infoForPDF[iwp][14] <<  endl;
				//cout << FitForXSResults[taggerArray[iwp]][FitForXSResults[taggerArray[iwp]].size()-1] << endl;
				
				cout << infoForPDF[iwp].size() << endl;
				nTimesWP++;
				
				delete results;
				
			}
			
			// build latex file
			
			// open latex file to write out the results for btag eff and XS
			
			lastTagger=0;
			
			ofstream resultsFile;
			ofstream btagANtableFile;
			ofstream btagRESFile;
			ofstream mistagRESFile;
			ofstream xsRESFile;
			ofstream xsANtableFile;
			
			stringstream cut; cut << matchChiSquareCut_;
			stringstream fitm; fitm << fitMode;
			
			string resultsFileName = "2DBtagEff_XS_analysisResults_chisqCut"+cut.str()+"_sampleType_"+sampleType+"_fitMode"+fitm.str()+"_channel"+data_postfix+".tex";
			//string btagANtableFileName = "BTAG_ANTables_chisqCut"+cut.str()+"_sampleType_"+sampleType+"_fitMode"+fitm.str()+"_channel"+data_postfix+".tex";
			string btagANtableFileName = "BTAG_ANTables_chisqCut"+cut.str()+"_sampleType_"+sampleType+"_channel"+data_postfix+".tex";
			//string btagRESFileName = "BTAG_RESULTS_chisqCut"+cut.str()+"_sampleType_"+sampleType+"_fitMode"+fitm.str()+"_channel"+data_postfix+".txt";
			string btagRESFileName = "BTAG_RESULTS_chisqCut"+cut.str()+"_sampleType_"+sampleType+"_channel"+data_postfix+".txt";
			string mistagRESFileName = "MISTAG_RESULTS_chisqCut"+cut.str()+"_sampleType_"+sampleType+"_channel"+data_postfix+".txt";
			string xsRESFileName = "XS_RESULTS_chisqCut"+cut.str()+"_sampleType_"+sampleType+"_fitMode"+fitm.str()+"_channel"+data_postfix+".txt";
			string xsANtableFileName = "XS_ANTables_chisqCut"+cut.str()+"_sampleType_"+sampleType+"_fitMode"+fitm.str()+"_channel"+data_postfix+".tex";
			
			resultsFile.open (resultsFileName.c_str(), ios::out | ios::trunc );
			btagANtableFile.open (btagANtableFileName.c_str(), ios::out | ios::trunc );
			btagRESFile.open (btagRESFileName.c_str(), ios::out | ios::trunc );
			mistagRESFile.open (mistagRESFileName.c_str(), ios::out | ios::trunc );
			xsRESFile.open (xsRESFileName.c_str(), ios::out | ios::trunc );
			xsANtableFile.open (xsANtableFileName.c_str(), ios::out | ios::trunc );
			
			resultsFile << "\\documentclass[landscape,14pt]{article}\n";
			resultsFile << "\\usepackage{graphicx}\n";
			resultsFile << "\\usepackage[usenames,dvipsnames]{color}\n";
			resultsFile << "\\usepackage[margin=0.5cm]{geometry}\n";
			resultsFile << "\\title{Results for the 2D ($\\epsilon_{b},\\sigma_{t\\bar{t}}$) measurement}\n";
			resultsFile << "\\author{Michael Maes}\n";
			resultsFile << "\\date{\\today}\n";
			resultsFile << "\\begin{document}\n";
			//resultsFile << "\\maketitle\n";
			//resultsFile << "\\newpage\n";
			
			/*
			 NEW PAGE FOR EACH TAGGER
			 
			 */
			
			lastTagger=0;
			
			string discName = "";
			string sdiscName = "";
			
			// BTAG SUMMARY TABLE FOR AN!
			
			if (datasetName != "Data") {
				btagANtableFile << "& WP & $\\epsilon_{b}^{MCtruth} (\\%)$ & $F^{exp}$ & $\\hat{\\epsilon}_{b}(F^{exp}) (\\%)$ & $F$ & $\\hat{\\epsilon}_{b}(F) (\\%)$  \\\\ \n";
			} else {
				btagANtableFile << "& WP & $F$ & $\\hat{\\epsilon}_{b}(F) (\\%)$  \\\\ \n";					
			}
			
			// XS SUMMARY TABLE FOR AN!

			xsANtableFile << "& WP & Type & $\\epsilon_{b}$ & $\\epsilon_{mistag}$ & $\\frac{\\#b}{\\#b+\\#non-b}$ & $\\sigma_{t\\bar{t}}$ (pb) \\\\ \n";

			// XS SUMMARY TABLE FOR AN
			
			xsANtableFile << "\\hline\n";
			xsANtableFile << " & ";
			xsANtableFile << "No cut & ";
			xsANtableFile << "- & ";
			xsANtableFile << "- & ";
			xsANtableFile << "- & ";
			xsANtableFile << "- & ";
			xsANtableFile << infoForPDF[0][0]<< " $\\pm$ " << infoForPDF[0][1];				
			xsANtableFile << "\\\\\n";
            
            stringstream sys; sys<<nSystematic;
            
            TString fileEff = "FitOutput/Efficiencies_nSystematic_"+sys.str()+data_postfix+".root";
            
            TFile* effs = new TFile(fileEff,"RECREATE");
            			
			for(int iwp=0; iwp < nWP; iwp++){
				
				if (iwp != 0 && taggerArray[iwp] == lastTagger) continue;
				
				lastTagger = taggerArray[iwp];
				
				vector<string> PlotFileNames = aContainer->GetFitPlotPaths(lastTagger);
				
                // OLD ORDERING
				/*if (taggerArray[iwp] == 0) {
                 discName="TrackCountingHighEfficiency tagger";
                 sdiscName="TCHE";
                 } else if (taggerArray[iwp] == 1) {
                 discName="TrackCountingHighPurity tagger";
                 sdiscName="TCHP";
                 } else if (taggerArray[iwp] == 2) {
                 discName="SimpleSecondaryVertexHighEff tagger";
                 sdiscName="SSVHE";
                 } else if (taggerArray[iwp] == 3) {
                 discName="SimpleSecondaryVertexHighPurity tagger";
                 sdiscName="SSVHP";
                 } else if (taggerArray[iwp] == 4) {
                 discName="Combined Secondary Vertex tagger";
                 sdiscName="CSV";
                 } else if (taggerArray[iwp] == 5) {
                 discName="Combined Secondary Vertex MVA tagger";
                 sdiscName="CSVMVA";
                 } else if (taggerArray[iwp] == 7) {
                 discName="Jet Probability tagger";
                 sdiscName="JP";
                 } else if (taggerArray[iwp] == 6) {
                 discName="Jet BProbability tagger";
                 sdiscName="JBP";*/
                
                if (taggerArray[iwp] == 0) {
                    discName="TrackCountingHighEfficiency tagger";
                    sdiscName="TCHE";
                } else if (taggerArray[iwp] == 1) {
                    discName="TrackCountingHighPurity tagger";
                    sdiscName="TCHP";
                } else if (taggerArray[iwp] == 4) {
                    discName="SimpleSecondaryVertexHighEff tagger";
                    sdiscName="SSVHE";
                } else if (taggerArray[iwp] == 5) {
                    discName="SimpleSecondaryVertexHighPurity tagger";
                    sdiscName="SSVHP";
                } else if (taggerArray[iwp] == 6) {
                    discName="Combined Secondary Vertex tagger";
                    sdiscName="CSV";
                } else if (taggerArray[iwp] == 7) {
                    discName="Combined Secondary Vertex tagger (New)";
                    sdiscName="CSVR";
                } else if (taggerArray[iwp] == 2) {
                    discName="Jet Probability tagger";
                    sdiscName="JP";
                } else if (taggerArray[iwp] == 3) {
                    discName="Jet BProbability tagger";
                    sdiscName="JBP";
				} else {
					discName="Unknown tagger";
					sdiscName="N/A";
				}
				
				if (datasetName == "Data")
					discName += " (Results on Collision data) ";
                
                string chan="Unknown";
                
                if (decay==0) chan="$\\mu$+jets";
                else if (decay==1) chan="e+jets";
                
                string var="Unknown";
                if (fitMode==0) var="$m_{lj}$";
                if (fitMode==1) var="$M3$";
                if (fitMode==2) var="2D $(m_{lj},M3)$";
    
				resultsFile << "\\newpage\n"; // start new page for each tagger*/
				//resultsFile << "\\begin{center}\n";
				//resultsFile << "\\section{" << discName << " L=" << desiredIntLum_ << " $pb^{-1}$ ($\\chi^{2}_{cut}=" << cut.str() << "$, $\\epsilon_{\\chi^{2}_{cut}}=" << round(1000,chiSQCutEff)*100 << "\\%$, $\\epsilon_{sel}=" << round(1000,refSelEff)*100 << "\\%$)}\n"; 
				resultsFile << "\\section{" << discName << " L=" << round(10,desiredIntLum_/1000.0) << " $fb^{-1}$ ($\\chi^{2}_{cut}=" << cut.str() << "$, channel: " << chan << ", Var: " << var << ", sample: " << sampleType << ")}\n"; 
				//resultsFile << "\\begin{center}\n";
				
				// TABLE WITH BTAG EFF INFORMATION 
				
				resultsFile << "\\subsection{Results on $\\epsilon_{b}$ and $\\epsilon_{q}$}\n";
				resultsFile << "\\begin{center}\n";
				resultsFile << "\t\\begin{tabular}{|c|c|c|c|c|c||c|c|}\n";
				resultsFile << "\t\t\\hline\n";
				//resultsFile << "\t\tWP & $\\epsilon_{b}^{True} (\\%)$ & $F^{exp}$ & $\\hat{\\epsilon}_{b}(F^{exp}) (\\%)$ & $F_{data}$ & $\\hat{\\epsilon}_{b}(F^{CS}) (\\%)$ & Bias[$F^{CS}$,$F^{exp}$] (\\%)& Bias[$\\hat{\\epsilon}_{b}(F^{exp})$,$\\epsilon_{b}^{true}$] (\\%)& Bias[$\\hat{\\epsilon}_{b}(F^{CS})$,$\\epsilon_{b}^{True}$] (\\%)\\\\\n";
				resultsFile << "\t\tWP & $\\epsilon_{b}^{True} (\\%)$ & $F^{exp}$ & $\\hat{\\epsilon}_{b}(F^{exp}) (\\%)$ & $F_{data}$ & $\\hat{\\epsilon}_{b}(F^{CS}) (\\%)$ & $\\epsilon_{q}^{True} (\\%)$ & $\\hat{\\epsilon_{q}} (\\%)$ \\\\\n";
				
                resultsFile << "\t\t\\hline\n";
				resultsFile << "\t\t\\hline\n"; 
				
				bool putSDiscName=false;
				
				for(int niwp=0; niwp < nWP; niwp++){
					
					if (taggerArray[niwp] != lastTagger) continue;
					
					if (wpArray[niwp] == 0) continue;
					
					resultsFile << "\t\t" << wpArray[niwp] << " & ";
					
					if (datasetName != "Data") {
						resultsFile << "{\\color{ForestGreen} " << infoForPDF[niwp][16] << " $\\pm$ " << infoForPDF[niwp][17] << "} & ";
						resultsFile << infoForPDF[niwp][12] << " $\\pm$ " << infoForPDF[niwp][13] << " & ";
						resultsFile << "{\\color{BurntOrange}" << infoForPDF[niwp][10] << " $\\pm$ " << infoForPDF[niwp][11] << "} & ";
					} else {
						resultsFile << "{\\color{ForestGreen} - }& ";
						resultsFile << " - & ";
						resultsFile << "{\\color{BurntOrange} - }& ";
					}
					
					resultsFile << infoForPDF[niwp][6] << " $\\pm$ " << infoForPDF[niwp][7] << " & ";
					resultsFile << "{\\color{blue}" << infoForPDF[niwp][4] << " $\\pm$ " << infoForPDF[niwp][5] << "} & ";
                    
                    if (datasetName == "Data")
                        resultsFile << " - & ";
                    else    
                        resultsFile << "{\\color{OliveGreen} " << infoForPDF[niwp][20] << " $\\pm$ " << infoForPDF[niwp][21] << "} & ";
                    
                    resultsFile << "{\\color{Purple} " << infoForPDF[niwp][22] << " $\\pm$ " << infoForPDF[niwp][23] << "} ";

                    resultsFile << "\\\\ \n";

					/*if (datasetName != "Data") {

						std::pair<float,float> biasF = calcBias(infoForPDF[niwp][6],infoForPDF[niwp][7],infoForPDF[niwp][12],infoForPDF[niwp][13]);
						std::pair<float,float> biasEpsFEXP = calcBias(infoForPDF[niwp][10],infoForPDF[niwp][11],infoForPDF[niwp][16],infoForPDF[niwp][17]);
						std::pair<float,float> biasEpsF = calcBias(infoForPDF[niwp][4],infoForPDF[niwp][5],infoForPDF[niwp][16],infoForPDF[niwp][17]);
						
						resultsFile << round(10,biasF.first*100) << "$\\pm$"<< round(10,biasF.second*100) << " & ";

						resultsFile << round(10,biasEpsFEXP.first*100) << "$\\pm$"<< round(10,biasEpsFEXP.second*100) << " & ";
                    
						resultsFile << round(10,biasEpsF.first*100) << "$\\pm$"<< round(10,biasEpsF.second*100);

						resultsFile << "\\\\ \n";
					
					} else {
					
					  resultsFile << " - & ";
					  resultsFile << " - & - \\\\ \n ";
						
					}*/
					
					string wpstr = "";
					if (niwp%3==0)
						wpstr="L";
					else if ( (niwp-1) %3==0)
						wpstr="M";
					else if ( (niwp-2) %3==0)
						wpstr="T";
					
					string tstr = "";
					
					if (!putSDiscName) {
						tstr = sdiscName;
						putSDiscName=true;
					}
					
					if (tstr != "") btagANtableFile << "\\hline\n"; 

                    // commented now because we define the SF to the E_b (method)
                    if (datasetName == "Data") {											
                        btagRESFile << sdiscName << " " << wpstr << " " << infoForPDF[niwp][4] << " " << infoForPDF[niwp][5] << "\n"; //method
                        mistagRESFile << sdiscName << " " << wpstr << " " << infoForPDF[niwp][22] << " " << infoForPDF[niwp][23] << "\n"; //method
                    }
                    else {
                        btagRESFile << sdiscName << " " << wpstr << " " << infoForPDF[niwp][16] << " " << infoForPDF[niwp][17]; //truth
                        btagRESFile << " " << infoForPDF[niwp][4] << " " << infoForPDF[niwp][5] << "\n"; //meas
                        mistagRESFile << sdiscName << " " << wpstr << " " << infoForPDF[niwp][20] << " " << infoForPDF[niwp][21]; //truth
                        mistagRESFile << " " << infoForPDF[niwp][22] << " " << infoForPDF[niwp][23] << "\n"; //meas
					}
					if (datasetName != "Data") {												
						btagANtableFile << tstr << " & " << wpstr << " & " << infoForPDF[niwp][16] << " $\\pm$ " << infoForPDF[niwp][17] << " & ";
						btagANtableFile << infoForPDF[niwp][12] << " $\\pm$ " << infoForPDF[niwp][13] << " & ";
						btagANtableFile << infoForPDF[niwp][10] << " $\\pm$ " << infoForPDF[niwp][11] << " & ";
						btagANtableFile << infoForPDF[niwp][6] << " $\\pm$ " << infoForPDF[niwp][7] << " & ";
						btagANtableFile << infoForPDF[niwp][4] << " $\\pm$ " << infoForPDF[niwp][5] << " \\\\ \n";
					} else {
						btagANtableFile << tstr << " & " << wpstr << " & " << infoForPDF[niwp][6] << " $\\pm$ " << infoForPDF[niwp][7] << " & ";
						btagANtableFile << infoForPDF[niwp][4] << " $\\pm$ " << infoForPDF[niwp][5] << " \\\\ \n";						
					}
					
					resultsFile << "\t\t\\hline\n"; 
					
				}
				resultsFile << "\t\\end{tabular}\n";
				resultsFile << "\\end{center}\n";
				
				// TABLE WITH XS MEASUREMENT INFORMATION 

                if (datasetName != "Data") resultsFile << "\\begin{minipage}{0.5\\textwidth}\n";
				resultsFile << "\\subsection{Results on ($\\epsilon_{b},\\sigma_{t\\bar{t}}$)}\n";
                resultsFile << "\\vspace{-7mm}\n";
				resultsFile << "\\begin{center}\n";
				resultsFile << "\t$\\begin{array}{cc}\n";

				resultsFile << "\t\t\\begin{tabular}{|c|c|c|}\n";
				resultsFile << "\t\t\t\\hline\n";
				resultsFile << "\t\t\t Efficiency & Value & Nominal value \\\\ \n";	
				resultsFile << "\t\t\t\\hline\n"; 
				resultsFile << "\t\t\t\\hline\n";
				//$\\epsilon_{\\chi^{2}_{cut}}=" << round(1000,chiSQCutEff)*100 << "\\%$, $\\epsilon_{sel}=" << round(1000,refSelEff)*100 << "\\%$)}
				if (datasetName != "Data") {
					resultsFile << "\t\t\t$\\epsilon_{\\chi^{2}_{cut}}$ & " << round(1000,chiSQCutEff_sample)*100 << " & " << round(1000,chiSQCutEff)*100 << " \\\\ \n";
					resultsFile << "\t\t\t\\hline\n";
					resultsFile << "\t\t\t$\\epsilon_{m_{\\mu j}^{cut}}$ & " << round(1000,MLBCutEff_sample)*100 << " & " << round(1000,MLBCutEff)*100 << " \\\\ \n";
					resultsFile << "\t\t\t\\hline\n";
					resultsFile << "\t\t\t$\\epsilon_{sel}$ & " << round(1000,refSelEff_sample)*100 << " & " << round(1000,refSelEff)*100 << " \\\\ \n";
				} else {
					resultsFile << "\t\t\t$\\epsilon_{\\chi^{2}_{cut}}$ & - & " << round(1000,chiSQCutEff)*100 << " \\\\ \n";
					resultsFile << "\t\t\t\\hline\n";
					resultsFile << "\t\t\t$\\epsilon_{m_{\\mu j}^{cut}}$ & - & " << round(1000,MLBCutEff)*100 << " \\\\ \n";
					resultsFile << "\t\t\t\\hline\n";
					resultsFile << "\t\t\t$\\epsilon_{sel}$ & - & " << round(1000,refSelEff)*100 << " \\\\ \n";
				}
				
				resultsFile << "\t\t\t\\hline\n";

				resultsFile << "\t\t\\end{tabular} &\n";
				
				resultsFile << "\t\t\\begin{tabular}{|c|c|c|c|c|c|}\n";
				resultsFile << "\t\t\t\\hline\n";
				resultsFile << "\t\t\tWP & Type & $\\epsilon_{b}$ & $\\epsilon_{mistag}$ & $\\frac{\\#b}{\\#b+\\#non-b}$ & $\\sigma_{t\\bar{t}}$ (pb) \\\\ \n";	
				resultsFile << "\t\t\t\\hline\n"; 
				resultsFile << "\t\t\t\\hline\n";
				
				resultsFile << "\t\t\t No cut & ";
				resultsFile << "- & ";
				resultsFile << "- & ";
				resultsFile << "- & ";
				resultsFile << "- & ";
				resultsFile << infoForPDF[iwp][0]<< " $\\pm$ " << infoForPDF[iwp][1];				
				resultsFile << "\\\\\n";
				resultsFile << "\t\t\t\\hline\n"; 
				
				 
				putSDiscName = false;
				
				for(int niwp=0; niwp < nWP; niwp++){
					
					if (taggerArray[niwp] != lastTagger) continue;
					
					if (wpArray[niwp] == 0) continue;
					
					resultsFile << "\t\t\t" << wpArray[niwp] << " & ";
					resultsFile << "{\\color{blue}$\\hat{\\epsilon}_{b}(F^{CS})$} & ";
					resultsFile << "{\\color{blue}" << infoForPDF[niwp][4] << " $\\pm$ " << infoForPDF[niwp][5] << "} & ";
					resultsFile << "{\\color{blue}" << infoForPDF[niwp][18] << "} & ";
					resultsFile << "{\\color{blue}" << infoForPDF[niwp][19] << "} & ";
					resultsFile << "{\\color{blue}" << infoForPDF[niwp][2] << " $\\pm$ " << infoForPDF[niwp][3] << "}";				
					resultsFile << "\\\\\n";
					resultsFile << "\t\t\t\\hline\n";
					
					
					
					if (datasetName != "Data") {
						resultsFile << "\t\t\t & ";
						resultsFile << "{\\color{BurntOrange}$\\hat{\\epsilon}_{b}(F^{exp})$} & ";
						resultsFile << "{\\color{BurntOrange}" << infoForPDF[niwp][10] << " $\\pm$ " << infoForPDF[niwp][11] << "} & ";
						resultsFile << "{\\color{BurntOrange}" << infoForPDF[niwp][18] << "} & ";
						resultsFile << "{\\color{BurntOrange}" << infoForPDF[niwp][19] << "} & ";
						resultsFile << "{\\color{BurntOrange}" << infoForPDF[niwp][8] << " $\\pm$ " << infoForPDF[niwp][9] << "}";				
						resultsFile << "\\\\\n";
						resultsFile << "\t\t\t\\hline\n"; 
						
						resultsFile << "\t\t\t & ";
						resultsFile << "{\\color{OliveGreen}$\\epsilon_{b}^{true}$} & ";
						resultsFile << "{\\color{OliveGreen}" << infoForPDF[niwp][16] << " $\\pm$ " << infoForPDF[niwp][17] << "} & ";
						resultsFile << "{\\color{OliveGreen}" << infoForPDF[niwp][18] << "} & ";
						resultsFile << "{\\color{OliveGreen}" << infoForPDF[niwp][19] << "} & ";
						resultsFile << "{\\color{OliveGreen}" << infoForPDF[niwp][14] << " $\\pm$ " << infoForPDF[niwp][15] << "}";				
						resultsFile << "\\\\\n";
						resultsFile << "\t\t\t\\hline\n"; 
					}
					
					//XS SUMMARY TABLE FOR AN
					
					string wpstr = "";
					if (niwp%3==0)
						wpstr="L";
					else if ( (niwp-1) %3==0)
						wpstr="M";
					else if ( (niwp-2) %3==0)
						wpstr="T";
					
					string tstr = "";
					
					if (!putSDiscName) {
						tstr = sdiscName;
						putSDiscName=true;
					}
					
					if (tstr != "") xsANtableFile << "\\hline\n"; 
					
					xsANtableFile << tstr << " & " << wpstr << " & ";
					xsANtableFile << "$\\hat{\\epsilon}_{b}(F^{CS})$& ";
					xsANtableFile << infoForPDF[niwp][4] << " $\\pm$ " << infoForPDF[niwp][5] << " & ";
					xsANtableFile << infoForPDF[niwp][18] << " & ";
					xsANtableFile << infoForPDF[niwp][19] << " & ";
					xsANtableFile << infoForPDF[niwp][2] << " $\\pm$ " << infoForPDF[niwp][3] << "";				
					xsANtableFile << "\\\\\n";					
					
					
					if (datasetName != "Data") {
						
						xsANtableFile << "& & ";
						xsANtableFile << "$\\hat{\\epsilon}_{b}(F^{exp})$& ";
						xsANtableFile << infoForPDF[niwp][10] << " $\\pm$ " << infoForPDF[niwp][11] << " & ";
						xsANtableFile << infoForPDF[niwp][18] << " & ";
						xsANtableFile << infoForPDF[niwp][19] << " & ";
						xsANtableFile << infoForPDF[niwp][8] << " $\\pm$ " << infoForPDF[niwp][9] << "";				
						xsANtableFile << "\\\\\n";
						
						xsANtableFile << "& & ";
						xsANtableFile << "$\\epsilon_{b}^{true}$& ";
						xsANtableFile << infoForPDF[niwp][16] << " $\\pm$ " << infoForPDF[niwp][17] << " & ";
						xsANtableFile << infoForPDF[niwp][18] << " & ";
						xsANtableFile << infoForPDF[niwp][19] << " & ";
						xsANtableFile << infoForPDF[niwp][14] << " $\\pm$ " << infoForPDF[niwp][15] << "";				
						xsANtableFile << "\\\\\n";
						
					}
                    
                    xsRESFile << sdiscName << " " << wpstr << " " << infoForPDF[niwp][16] << " " << infoForPDF[niwp][17] << " " << infoForPDF[niwp][4] << " " << infoForPDF[niwp][5] << " " << infoForPDF[niwp][2] << " " << infoForPDF[niwp][3] << "\n"; //method
					
                    
                    stringstream countwp; countwp<<niwp;
                    TString effname="eff_WP_"+countwp.str();
                    TH1D* effHisto = new TH1D(effname,sampleType.c_str(),8,-0.5,7.5);
                    
                    effHisto->SetBinContent(1,desiredIntLum_/1000);
                    effHisto->SetBinContent(2,refSelEff);
                    effHisto->SetBinContent(3,chiSQCutEff);
                    effHisto->SetBinContent(4,MLBCutEff);
                    effHisto->SetBinContent(5,chiSQCutEff);
                    /*effHisto->SetBinContent(5,infoForPDF[niwp][19]/100);
                     effHisto->SetBinContent(6,infoForPDF[niwp][18]/100);
                     effHisto->SetBinContent(7,infoForPDF[niwp][4]/100);
                     effHisto->SetBinError(7,infoForPDF[niwp][5]/100);*/
                    effHisto->SetBinContent(5,InfoForCombinedFit[niwp][3]);
                    effHisto->SetBinContent(6,InfoForCombinedFit[niwp][2]);
                    effHisto->SetBinContent(7,InfoForCombinedFit[niwp][0]);
                    effHisto->SetBinError(7,InfoForCombinedFit[niwp][1]);

                    effHisto->Write();
					
					
				}
                
				
				resultsFile << "\t\t\\end{tabular}\n";
				
				resultsFile << "\t\\end{array}$\n";
				resultsFile << "\\end{center}\n";
                
                if (datasetName != "Data") {
                    resultsFile << "\\end{minipage}\n";
                    resultsFile << "\\qquad\n";

                    resultsFile << "\\begin{minipage}{0.57\\textwidth}\n";
                    
                    resultsFile << "\\begin{center}\n";
                    resultsFile << "\t\\begin{tabular}{|c|c|c|c|c|}\n";
                    resultsFile << "\t\t\\hline\n";

                    resultsFile << "\t\tWP & $\\frac{F^{CS}-F^{exp}}{F^{exp}}$ & $\\frac{\\epsilon_{b}^{exp}-\\epsilon_{b}^{True}}{\\epsilon_{b}^{True}}$ & $\\frac{\\epsilon_{b}^{CS}-\\epsilon_{b}^{True}}{\\epsilon_{b}^{True}}$ &  $\\frac{\\epsilon_{q}-\\epsilon_{q}^{True}}{\\epsilon_{q}^{True}}$\\\\\n";
                    
                    resultsFile << "\t\t\\hline\n";
                    resultsFile << "\t\t\\hline\n"; 
                    
                    for(int niwp=0; niwp < nWP; niwp++){
                        
                        if (taggerArray[niwp] != lastTagger) continue;
                        
                        if (wpArray[niwp] == 0) continue;
                        
                        resultsFile << "\t\t" << wpArray[niwp] << " & ";
                        
                        std::pair<float,float> biasF = calcBias(infoForPDF[niwp][6],infoForPDF[niwp][7],infoForPDF[niwp][12],infoForPDF[niwp][13]);
                        std::pair<float,float> biasEpsFEXP = calcBias(infoForPDF[niwp][10],infoForPDF[niwp][11],infoForPDF[niwp][16],infoForPDF[niwp][17]);
                        std::pair<float,float> biasEpsF = calcBias(infoForPDF[niwp][4],infoForPDF[niwp][5],infoForPDF[niwp][16],infoForPDF[niwp][17]);
                        std::pair<float,float> biasQEps = calcBias(infoForPDF[niwp][22],infoForPDF[niwp][23],infoForPDF[niwp][20],infoForPDF[niwp][21]);
                        
                        //resultsFile << "{\\color{LimeGreen} " << infoForPDF[niwp][20] << " $\\pm$ " << infoForPDF[niwp][21] << "} & ";

                        resultsFile << "{\\color{Black} " << round(10,biasF.first*100) << "$\\pm$"<< round(10,biasF.second*100) << " } & ";
                        
                        resultsFile << "{\\color{BurntOrange} " << round(10,biasEpsFEXP.first*100) << "$\\pm$"<< round(10,biasEpsFEXP.second*100) << " } & ";
                        
                        resultsFile << "{\\color{blue} " << round(10,biasEpsF.first*100) << "$\\pm$"<< round(10,biasEpsF.second*100) << " } & ";
                        
                        resultsFile << "{\\color{Purple} " << round(10,biasQEps.first*100) << "$\\pm$"<< round(10,biasQEps.second*100) << " } ";
                        
                        resultsFile << "\\\\ \n";
                        
                        resultsFile << "\t\t\\hline\n"; 
                        
                    }
                    
                    
                    resultsFile << "\t\\end{tabular}\n";
                    resultsFile << "\\end{center}\n";

                
                    resultsFile << "\\end{minipage}\n";
                
                }
				
				// TEMPLATE FIT PLOTS
                resultsFile << "\\vspace{-0.2cm}\n";
				
				resultsFile << "\\subsection{Template Fit and 2D correlation plots (with $F_{data}$)}\n";
				
				resultsFile << "\\vspace{-0.5cm}\n";
				resultsFile << "\\begin{figure}[!h]\n";
				resultsFile << "\\begin{center}\n";
				resultsFile << "$\\begin{array}{ccc}\n";
				
				//cout << PlotFileNames.size() << endl; exit(1);
				
				vector<float> discwp;
				
				for(int niwp=0; niwp < nWP; niwp++)
					if (taggerArray[niwp] == lastTagger)
						discwp.push_back(wpArray[niwp]);
				
				//for (int p=0;p<discwp.size();p++)
				//	cout << iwp << " " << p << " " << discwp[p] << endl;
				
				//exit(1);
				
				for (unsigned int file=0;file<PlotFileNames.size(); file++) {
					//cout << PlotFileNames[file] << endl;
					if (PlotFileNames[file].find("NOCUT") != -1) continue;
					
					if (discwp[0] == 0 && PlotFileNames[file].find("LCut") != -1) continue;
					if (discwp[1] == 0 && PlotFileNames[file].find("MCut") != -1) continue;
					if (discwp[2] == 0 && PlotFileNames[file].find("TCut") != -1) continue;
					
					if (file != PlotFileNames.size()-1)
						resultsFile << "\\includegraphics[width=0.18\\textwidth]{"<<PlotFileNames[file]<<"} &\n";
					else 
						resultsFile << "\\includegraphics[width=0.18\\textwidth]{"<<PlotFileNames[file]<<"} \\\\ \n";
					
				}
				
				discwp.clear();
								
				for(int niwp=0; niwp < nWP; niwp++){
					
					if (taggerArray[niwp] != lastTagger) continue;
					
					if (wpArray[niwp] == 0) continue;
					
					stringstream wp; wp << niwp;
                    stringstream fitm; fitm << fitMode;

					//cout << PlotFileNames[file] << endl;"2DMeasurement/"+title+"_channel"+data_postfix+".png"
					string filename = "2DMeasurement/2DMeasurement_FDATA_XS_btag_IWP_"+wp.str()+"_channel"+data_postfix+"_fitMode"+fitm.str()+".png";
					
					if (niwp != nWP-1) {
						if (taggerArray[niwp+1] == lastTagger) {
							resultsFile << "\\includegraphics[width=0.18\\textwidth]{"<<filename<<"} &\n";
						} else {
							resultsFile << "\\includegraphics[width=0.18\\textwidth]{"<<filename<<"} \\\\\n";
						}
					} else {
						resultsFile << "\\includegraphics[width=0.18\\textwidth]{"<<filename<<"} \\\\\n";
					}
				}
				
				for(int niwp=0; niwp < nWP; niwp++){
					
					if (taggerArray[niwp] != lastTagger) continue;
					
					if (wpArray[niwp] == 0) continue;
					
					stringstream wp; wp << wpArray[niwp];
					
					if (niwp != nWP-1) {
						if (taggerArray[niwp+1] == lastTagger) {
							resultsFile << "Figure: WP "<<wp.str()<<"&\n";
						} else {
							resultsFile << "Figure: WP "<<wp.str()<<"\n";
						}
					} else {
						resultsFile << "Figure: WP "<<wp.str()<<"\n";
					}					
				}
				
				resultsFile << "\\end{array}$\n";
				resultsFile << "\\end{center}\n";
				resultsFile << "\\end{figure}\n";
				
			}
			resultsFile << "\\end{document}\n";
                        
            effs->Close();
			
			resultsFile.close();
			btagANtableFile.close();
			btagRESFile.close();
			mistagRESFile.close();

			xsRESFile.close();

            /*if (sampleType.find("Data") != string::npos) {
				ofstream Eb; Eb.open("Eb.data", ios::out );
				Eb << infoForPDF[0][4];
                Eb.close();
            }*/
            
			if (sampleType.find("nominal") != string::npos) {
                
				for(int niwp=0; niwp < nWP; niwp++){
                    
                    stringstream p; p<<niwp;

                    std::pair<float,float> biasEpsF = calcBias(infoForPDF[niwp][4],infoForPDF[niwp][5],infoForPDF[niwp][16],infoForPDF[niwp][17]);
                    biasEpsF.first = round(100,biasEpsF.first*100);
                    biasEpsF.second = round(100,biasEpsF.second*100);
                    
                    std::pair<float,float> biasMistag = calcBias(infoForPDF[niwp][22],infoForPDF[niwp][23],infoForPDF[niwp][20],infoForPDF[niwp][21]);
                    biasMistag.first = round(100,biasMistag.first*100);
                    biasMistag.second = round(100,biasMistag.second*100);
                    
                    //float dataEb = getNomVal(data_postfix,fitMode,"Eb.data");
                    
                    /*ofstream F; F.open(("systematics/F_WP"+p.str()+"_channel"+data_postfix+".nominal").c_str(), ios::out );
                    ofstream Eb; Eb.open(("systematics/Eb_WP"+p.str()+"_channel"+data_postfix+".nominal").c_str(), ios::out );
                    ofstream XS; XS.open(("systematics/XS_WP"+p.str()+"_channel"+data_postfix+".nominal").c_str(), ios::out );
                    ofstream uF; uF.open(("systematics/uF_WP"+p.str()+"_channel"+data_postfix+".nominal").c_str(), ios::out );
                    ofstream uEb; uEb.open(("systematics/uEb_WP"+p.str()+"_channel"+data_postfix+".nominal").c_str(), ios::out );
                    ofstream uXS; uXS.open(("systematics/uXS_WP"+p.str()+"_channel"+data_postfix+".nominal").c_str(), ios::out );
                    ofstream SF; SF.open(("systematics/SF_WP"+p.str()+"_channel"+data_postfix+".nominal").c_str(), ios::out );
                    
                    SF << infoForPDF[niwp][4]/infoForPDF[niwp][16];
                    
                    F << infoForPDF[niwp][6];
                    //Eb << infoForPDF[niwp][4];
                    XS << infoForPDF[niwp][2];
                    
                    uF << infoForPDF[niwp][7];
                    //uEb << infoForPDF[niwp][5];
                    uXS << infoForPDF[niwp][3];
                    
                    Eb << biasEpsF.first;
                    uEb << biasEpsF.second;
                    
                    F.close(); Eb.close(); XS.close();
                    uF.close(); uEb.close(); uXS.close();
                    SF.close();*/
                    
                    double uSF = round(1000,sqrt((pow(infoForPDF[niwp][5],2)/pow(infoForPDF[niwp][16],2) )+( (pow(infoForPDF[niwp][4],2)*pow(infoForPDF[niwp][17],2))/ pow(infoForPDF[niwp][16],4))));
                    double uSFq = round(1000,sqrt((pow(infoForPDF[niwp][23],2)/pow(infoForPDF[niwp][20],2) )+( (pow(infoForPDF[niwp][22],2)*pow(infoForPDF[niwp][21],2))/ pow(infoForPDF[niwp][20],4))));
                    
                    setNomVal(data_postfix,fitMode,"F",niwp,infoForPDF[niwp][6]);
                    setNomVal(data_postfix,fitMode,"Eb",niwp,biasEpsF.first);
                    setNomVal(data_postfix,fitMode,"Eq",niwp,biasMistag.first);
                    setNomVal(data_postfix,fitMode,"XS",niwp,infoForPDF[niwp][2]);
                    setNomVal(data_postfix,fitMode,"uF",niwp,infoForPDF[niwp][7]);
                    setNomVal(data_postfix,fitMode,"uEb",niwp,biasEpsF.second);
                    setNomVal(data_postfix,fitMode,"uEq",niwp,biasMistag.second);
                    setNomVal(data_postfix,fitMode,"uXS",niwp,infoForPDF[niwp][3]);
                    setNomVal(data_postfix,fitMode,"SF",niwp,infoForPDF[niwp][4]/infoForPDF[niwp][16]);
                    setNomVal(data_postfix,fitMode,"SFq",niwp,infoForPDF[niwp][22]/infoForPDF[niwp][20]);
                    setNomVal(data_postfix,fitMode,"uSF",niwp,uSF);
                    setNomVal(data_postfix,fitMode,"uSFq",niwp,uSFq);

                    stringstream fitm; fitm << fitMode;
                    ofstream csv; csv.open(("systematics/systematics_results_WP"+p.str()+"_channel"+data_postfix+"_fitMode"+fitm.str()+".csv").c_str(), ios::out | ios::app );
                    
                    
                    if (nWP > 1) {
                        //F
                        csv << sampleType << " " << sampleSel << " " << infoForPDF[niwp][6] << " " << infoForPDF[niwp][7] << " 0 0 ";

                        // Eb + SFb
                        csv << biasEpsF.first << " " << biasEpsF.second << " 0 0 ";
                        csv << infoForPDF[niwp][4]/infoForPDF[niwp][16]  << " " << uSF << " 0 0 ";
                        
                        // Eq + SFq
                        csv << biasMistag.first << " " << biasMistag.second << " 0 0 ";
                        csv << infoForPDF[niwp][22]/infoForPDF[niwp][20] << " " << uSFq << " 0 0 ";
                        
                        // XS
                        csv	<< infoForPDF[niwp][2] << " " << infoForPDF[niwp][3] << " 0 0 " << "\n";
                    }
                    csv.close();
                }
				
			} else {
                
                for(int niwp=0; niwp < nWP; niwp++){
                    
                    stringstream p; p<<niwp;

                    std::pair<float,float> biasEpsF = calcBias(infoForPDF[niwp][4],infoForPDF[niwp][5],infoForPDF[niwp][16],infoForPDF[niwp][17]);
                    biasEpsF.first = round(100,biasEpsF.first*100);
                    biasEpsF.second = round(100,biasEpsF.second*100);
                    
                    std::pair<float,float> biasMistag = calcBias(infoForPDF[niwp][22],infoForPDF[niwp][23],infoForPDF[niwp][20],infoForPDF[niwp][21]);
                    biasMistag.first = round(100,biasMistag.first*100);
                    biasMistag.second = round(100,biasMistag.second*100);
                    
                    float nomF = getNomVal(data_postfix,fitMode,"F",p.str()); 
                    float nomEb = getNomVal(data_postfix,fitMode,"Eb",p.str()); 
                    float nomEq = getNomVal(data_postfix,fitMode,"Eq",p.str()); 
                    float nomXS = getNomVal(data_postfix,fitMode,"XS",p.str()); 
                    
                    float nomuF = getNomVal(data_postfix,fitMode,"uF",p.str()); 
                    float nomuEb = getNomVal(data_postfix,fitMode,"uEb",p.str()); 
                    float nomuEq = getNomVal(data_postfix,fitMode,"uEq",p.str()); 
                    float nomuXS = getNomVal(data_postfix,fitMode,"uXS",p.str()); 
                    
                    float nomSF = getNomVal(data_postfix,fitMode,"SF",p.str());
                    float nomSFq = getNomVal(data_postfix,fitMode,"SFq",p.str());

                    float nomuSF = getNomVal(data_postfix,fitMode,"uSF",p.str());
                    float nomuSFq = getNomVal(data_postfix,fitMode,"uSFq",p.str());

                    float diffF = infoForPDF[niwp][6]-nomF;
                    //float diffEb = infoForPDF[niwp][4]-nomEb;
                    float diffXS = infoForPDF[niwp][2]-nomXS;
                    
                    float udiffF = sqrt(pow(nomuF,2)+pow(infoForPDF[niwp][7],2));
                    //float udiffEb = sqrt(pow(nomuEb,2)+pow(infoForPDF[niwp][5],2));
                    float udiffXS = sqrt(pow(nomuXS,2)+pow(infoForPDF[niwp][3],2));
                    
                    float diffEb = biasEpsF.first-nomEb;
                    float udiffEb = sqrt(pow(nomuEb,2)+pow(biasEpsF.second,2));

                    float diffEq = biasMistag.first-nomEq;
                    float udiffEq = sqrt(pow(nomuEq,2)+pow(biasMistag.second,2));
                    
                    float SF = infoForPDF[niwp][4]/infoForPDF[niwp][16];
                    float SFq = infoForPDF[niwp][22]/infoForPDF[niwp][20];
                    
                    double uSF = round(1000,sqrt((pow(infoForPDF[niwp][5],2)/pow(infoForPDF[niwp][16],2) )+( (pow(infoForPDF[niwp][4],2)*pow(infoForPDF[niwp][17],2))/ pow(infoForPDF[niwp][16],4))));
                    double uSFq = round(1000,sqrt((pow(infoForPDF[niwp][23],2)/pow(infoForPDF[niwp][20],2) )+( (pow(infoForPDF[niwp][22],2)*pow(infoForPDF[niwp][21],2))/ pow(infoForPDF[niwp][20],4))));
                    
                    float diffSF = SF-nomSF;
                    float diffSFq = SFq-nomSFq;
                    
                    float udiffSF = sqrt(pow(nomuSF,2)+pow(uSF,2));
                    float udiffSFq = sqrt(pow(nomuSFq,2)+pow(uSFq,2));

                    stringstream fitm; fitm << fitMode;
                    ofstream csv; csv.open(("systematics/systematics_results_WP"+p.str()+"_channel"+data_postfix+"_fitMode"+fitm.str()+".csv").c_str(), ios::out | ios::app );
                    
                    if (nWP > 1) {
                        // F
                        csv << sampleType << " " << sampleSel << " " << infoForPDF[niwp][6] << " " << infoForPDF[niwp][7] << " " << diffF << " " << udiffF << " ";

                        // Eb + SFb
                        
                        csv << biasEpsF.first << " " << biasEpsF.second << " " << diffEb << " " << udiffEb << " ";
                        csv	<< SF << " " << uSF << " " << diffSF << " " << udiffSF << " ";

                        // Eq + SFq
                        csv << biasMistag.first << " " << biasMistag.second << " " << diffEq << " " << udiffEq << " ";
                        csv	<< SFq << " " << uSFq << " " << diffSFq << " " << udiffSFq << " ";

                        // XS
                        csv	<< infoForPDF[niwp][2] << " " << infoForPDF[niwp][3] << " " << diffXS << " " << udiffXS << "\n";
                    }
                    
                    csv.close();
                    
                }
				
			}
			
			
			
			//exit(1);
			
			// do some chi2 checks
			/*ofstream chiSqfile;
			chiSqfile.open ("chi2cutcheck.txt", ios::out | ios::app );
			
			chiSqfile << matchChiSquareCut_ << " " << wpArray[3] << " " << infoForPDF[3][2] << " " << infoForPDF[3][8] << " " << (infoForPDF[3][2]-infoForPDF[3][8])/infoForPDF[3][8] << "\n";
			
			chiSqfile.close();
			//exit(1);
			*/
			if(verbosity_>1) cout << "+--> MLJ template fits for XS done " << endl;
			
			fout->cd();   
			
			if(!doPseudoExp_){  
				
				fout->mkdir("2DMeasurementHistos");
				
				fout->cd("2DMeasurementHistos");
				
				for(std::map<std::string,TH1D*>::const_iterator it = histos1D.begin(); it != histos1D.end(); it++) {
					it->second->Write();
					
					//cout << "++> Writing " << it->first << endl;
					
					//exit(1);
					delete histos1D[it->first];
				}
				
				for(std::map<std::string,TH2D*>::const_iterator it = histos2D.begin(); it != histos2D.end(); it++) {
					it->second->Write();
					
					//cout << "++> Writing " << it->first << endl;
					
					//exit(1);
					delete histos2D[it->first];
				}
				
				for(std::map<std::string,TGraphErrors*>::const_iterator it = graphs.begin(); it != graphs.end(); it++) {
					it->second->Write();
					
					//cout << "++> Writing " << it->first << endl;
					
					//exit(1);
					delete graphs[it->first];
				}
				
				for(std::map<std::string,TCanvas*>::const_iterator it = canvas_2Dmeas.begin(); it != canvas_2Dmeas.end(); it++) {
				//	it->second->Write();
					
					//cout << "++> Writing " << it->first << endl;
					
					//exit(1);
					delete canvas_2Dmeas[it->first];
				}
				fout->cd();
				
				if (!doPseudoExp_ && doChi2Check) {
					if(verbosity_>0) cout<<"+> writing the Chi2 cutted MTop histos " << endl;
					fout->mkdir("Chi2CutCheck");
					fout->cd("Chi2CutCheck");
					for(std::map<std::string,TH1D*>::const_iterator it = chi2CutHistos.begin(); it != chi2CutHistos.end(); it++) {
						
						cout << it->first << " "  << it->second->GetTitle() << endl;
						
						it->second->Write();
						
						//delete chi2CutHistos[it->first];
					}
					
					int nCuts = (int)(sizeof(chisqCutForMass)/sizeof(chisqCutForMass[0]));
					
					float xVals[100], yVals[100];
					float xVals2[100], yVals2[100];
					float xVals3[100], yVals3[100];
					
					for (int i=0; i<nCuts; i++) {
						
						stringstream c; c << chisqCutForMass[i];
						string title = "MTopDistribution_ChiSqCut_"+c.str();
						
						// definition 1
						int maxBin = chi2CutHistos[title]->GetMaximumBin();
						
						float peakMtop = chi2CutHistos[title]->GetBinCenter(maxBin);
						
						float maxMtop = 0;
						float minMtop = 99999999;
						
						for (unsigned int b=1;b<chi2CutHistos[title]->GetNbinsX();b++) {
							
							if (chi2CutHistos[title]->GetBinContent(b) > 0) {
								if (chi2CutHistos[title]->GetBinCenter(b) > maxMtop)
									maxMtop = chi2CutHistos[title]->GetBinCenter(b);
								
								if (chi2CutHistos[title]->GetBinCenter(b) < minMtop)
									minMtop = chi2CutHistos[title]->GetBinCenter(b);
							}
						}
						
						xVals[i]=chisqCutForMass[i];
						yVals[i]=fabs(fabs(maxMtop-peakMtop)-fabs(peakMtop-minMtop));
						
						cout << "Def1: " << i << " " << xVals[i] << " " << minMtop << " " << peakMtop << " "  << maxMtop << " " << yVals[i] << endl;
						
						// definition2
						
						Float_t maxMTOP = chi2CutHistos[title]->GetMean();
						Float_t sigmaMTOP = chi2CutHistos[title]->GetRMS();
						
						TF1 *fitfunc;
						
						string func_title = title+"_Fitted";
						
						double rms = 5*chi2CutHistos[title]->GetBinWidth(1);
						double maxbin =  chi2CutHistos[title]->GetBinCenter(chi2CutHistos[title]->GetMaximumBin());
						
						fitfunc = new TF1(func_title.c_str(),"gaus");
						fitfunc->SetRange(maxbin-rms,maxbin+rms);
						chi2CutHistos[title]->Fit(fitfunc,"RQ");
						
						fitfunc->Write();	
						
						Float_t maxMTOP2 = fitfunc->GetParameter(1);
						Float_t sigmaMTOP2 = fitfunc->GetParameter(2);
						
						int BinLow = 0; int BinHigh = 0;
						int BinLow2 = 0; int BinHigh2 = 0;
						
						for (Int_t bin = 1; bin < chi2CutHistos[title]->GetNbinsX(); bin++) {
							
							float up = chi2CutHistos[title]->GetBinCenter(bin)+(chi2CutHistos[title]->GetBinWidth(bin)/2);
							float low = chi2CutHistos[title]->GetBinCenter(bin)-(chi2CutHistos[title]->GetBinWidth(bin)/2);
							
							if (maxMTOP-(2*sigmaMTOP) > low && maxMTOP-(2*sigmaMTOP) < up) BinLow = bin;
							if (maxMTOP+(2*sigmaMTOP) > low && maxMTOP+(2*sigmaMTOP) < up) BinHigh = bin;
							
							if (maxMTOP2-(2*sigmaMTOP2) > low && maxMTOP2-(2*sigmaMTOP2) < up) BinLow2 = bin;
							if (maxMTOP2+(2*sigmaMTOP2) > low && maxMTOP2+(2*sigmaMTOP2) < up) BinHigh2 = bin;
							
						}
						
						//cout << maxMTOP << " " << sigmaMTOP << " " << histo1D_btag["BestJetCombMTop"]->GetBinCenter(BinLow) << " " << histo1D_btag["BestJetCombMTop"]->GetBinCenter(BinHigh) << endl;
						
						Int_t origFraction = chi2CutHistos["MTopDistribution_NOChiSqCut"]->Integral(BinLow,BinHigh);
						Int_t newFraction = chi2CutHistos[title]->Integral(BinLow,BinHigh);
						Int_t newFraction2 = chi2CutHistos[title]->Integral(BinLow2,BinHigh2);
						
						xVals2[i]=chisqCutForMass[i];
						yVals2[i]=(double)newFraction/(double)origFraction;
						
						xVals3[i]=chisqCutForMass[i];
						yVals3[i]=(double)newFraction2/(double)origFraction;
						
						cout << "Def2: " << i << " " << xVals2[i] << " " << yVals2[i]  << " " << newFraction << " " << origFraction << endl;
						
						
					}
					
					cout << nCuts << endl;
					
					//exit(1);
					
					TGraph* MtopSymmetry_trough_diff = new TGraph((int)(sizeof(chisqCutForMass)/sizeof(chisqCutForMass[0])),xVals,yVals);
					
					MtopSymmetry_trough_diff->SetNameTitle("Graph_MtopSymmetry_trough_diff","Symmetry of mtop distribution;#chi^{2} cut; Difference");
					
					MtopSymmetry_trough_diff->Write();
					
					TGraph* NeventsRemoved = new TGraph((int)(sizeof(chisqCutForMass)/sizeof(chisqCutForMass[0])),xVals2,yVals2);
					
					NeventsRemoved->SetNameTitle("Graph_NeventsRemoved_FromTopMassPeak","Fraction of events in the m_{t} mass peak;#chi^{2} cut; Fraction");
					
					NeventsRemoved->Write();
					
					TGraph* NeventsRemoved_fit = new TGraph((int)(sizeof(chisqCutForMass)/sizeof(chisqCutForMass[0])),xVals3,yVals3);
					
					NeventsRemoved_fit->SetNameTitle("Graph_NeventsRemoved_FromTopMassPeak_fit","Fraction of events in the m_{t} mass peak;#chi^{2} cut; Fraction");
					
					NeventsRemoved_fit->Write();
					
				}
				
				
			}
			
		} // end of XS calc
		
		if (!doPseudoExp_) {
			cout << "+> WriteContainerToRootFile: done" << endl;
		
			aContainer->WriteContainerToRootFile(fout,true,true,true,true,true,true,true,true,true,true,false,true); //write the pt/eta bins
		
		}
		
		cout<<"+--> Histograms written to root file " << outRootPath[0] << endl;
		
		cout<<"+--> Doing pseudo-experiment stuff " << endl;
		
		cout << "************************************************"<< endl << endl;;
		//    for(int iwp=0; iwp < nWP; iwp++){
		
		int psLastTagger=0;
		int psnTimesWP=1;
		for(int iwp=0; iwp < nWP; iwp++){
			if(iwp!=-1) {//0 ||iwp==1 ||iwp==2 ) { //||iwp==3 ||iwp==4 ||iwp==15 ||iwp==19){
				double *wpEffArray = new double[30];
				//if(!doPseudoExp_) 
				if(verbosity_>1) cout << "****** Efficiency values at working point "<< iwp <<" (Tagger: " << taggerArray[iwp] << " WP: "<< wpArray[iwp] << ") ******"<< endl;	
				if(verbosity_>1) aContainer->CoutWPEff(true,doSCreweigh,wpArray[iwp],taggerArray[iwp],wpEffArray,doPseudoExp_,true,0,0,true,runNb);
				
				if (doPseudoExp_) {	 
					
					foutPseudoExp->cd();
					
					//if(verbosity_>1) cout << "DOING PSEUDOEXP!!!!!!!!!!!!!!!!!" << endl;
					aContainer->GetWPEff(true,doSCreweigh,wpArray[iwp],taggerArray[iwp],wpEffArray,doPseudoExp_,true,0,0,false,runNb);
					stringstream s; s << wpArray[iwp];
                    
					string title = "Efficiency_PseudoExps_WP_"+s.str();
					
					if (pullHistos.find(title) == pullHistos.end())
						pullHistos[title] = new TH1D(title.c_str(),(title+";#hat{#epsilon}_{b};pseudo experiments").c_str(),50,0,2);
					pullHistos[title]->Fill(wpEffArray[8]);
					
					title = "UncEfficiency_PseudoExps_WP_"+s.str();
					
					if (pullHistos.find(title) == pullHistos.end())
						pullHistos[title] = new TH1D(title.c_str(),(title+";#delta#hat{#epsilon}_{b};pseudo experiments").c_str(),50,0,2);
					pullHistos[title]->Fill(wpEffArray[9]);
					
					title = "Efficiency_Unc_2D_PseudoExps_WP_"+s.str();
					
					if (pullHistos2D.find(title) == pullHistos2D.end())
						pullHistos2D[title] = new TH2D(title.c_str(),(title+";#hat{#epsilon}_{b};#delta#hat{#epsilon}_{b}").c_str(),80,0,8,80,0,8);
					pullHistos2D[title]->Fill(wpEffArray[8],wpEffArray[9]);
                    
                    // bias Eff(Fdata,ShapeData) vs Eff(Fmc,ShapeData)
                    title = "EffBias_VS_centerLimit_rightlimit_WP_"+s.str();
                    
                    if (pullHistos2D.find(title) == pullHistos2D.end())
						pullHistos2D[title] = new TH2D(title.c_str(),(title+";CenterLimit;Right Limit").c_str(),(mid_max-mid_min)/bias_step,mid_min,mid_max,(right_max-mid_min-20)/bias_step,mid_min+20,right_max);
                    pullHistos2D[title]->Fill(centerleftlimit_,rightlimit_,(wpEffArray[8]-wpEffArray[6])/wpEffArray[6]);
                    
                    int minb = 0, maxb=1.5;
                    int minl = 0, maxl=0.5;
                    
                    if (iwp == 0) {
                        minb = 0.6;
                        maxb=1.2;
                        minl = 0.1;
                        maxl=0.3;
                    } else if(iwp == 1) {
                        minb = 0.4;
                        maxb=0.8;
                        minl = 0.0;
                        maxl=0.01;
                    } else if(iwp == 2) {
                        minb = 0.1;
                        maxb=0.5;
                        minl = 0.;
                        maxl=0.05;
                    }  
                    
                    title = "MC_MisTagEfficiency_PseudoExps_WP_"+s.str();
					
					if (pullHistos.find(title) == pullHistos.end())
						pullHistos[title] = new TH1D(title.c_str(),(title+";#epsilon_{q};pseudo experiments").c_str(),75,minl,maxl);
					pullHistos[title]->Fill(wpEffArray[22]);

                    title = "Meas_MisTagEfficiency_PseudoExps_WP_"+s.str();
					
					if (pullHistos.find(title) == pullHistos.end())
						pullHistos[title] = new TH1D(title.c_str(),(title+";#epsilon_{q};pseudo experiments").c_str(),75,minl,maxl);
					pullHistos[title]->Fill(wpEffArray[24]);

                    title = "MC_MisTagEfficiency_VS_BtagEffMeas_WP_"+s.str();
                    
                    if (pullHistos2D.find(title) == pullHistos2D.end())
						pullHistos2D[title] = new TH2D(title.c_str(),(title+";#hat{#epsilon}_{b};#hat{#epsilon}_{#slash{b}}").c_str(),75,minb,maxb,75,minl,maxl);
                    pullHistos2D[title]->Fill(wpEffArray[8],wpEffArray[22]);
                    

                    title = "Meas_MisTagEfficiency_VS_BtagEffMeas_WP_"+s.str();
                    
                    if (pullHistos2D.find(title) == pullHistos2D.end())
						pullHistos2D[title] = new TH2D(title.c_str(),(title+";#hat{#epsilon}_{b};#hat{#epsilon}_{#slash{b}}").c_str(),75,minb,maxb,75,minl,maxl);
                    pullHistos2D[title]->Fill(wpEffArray[8],wpEffArray[24]);
                    

                    /*title = "3DEffBias_VS_centerLimit_rightlimit_WP_"+s.str();
                    
                    if (pullHistos3D.find(title) == pullHistos3D.end())
						pullHistos3D[title] = new TH3D(title.c_str(),(title+";Left Limit;CenterLimit;Right Limit").c_str(),(left_max-left_min)/bias_step,left_min,left_max,(mid_max-mid_min)/bias_step,mid_min,mid_max,(right_max-mid_min-20)/bias_step,mid_min+20,right_max);
                    pullHistos3D[title]->Fill(leftlimit_,centerleftlimit_,rightlimit_,(wpEffArray[8]-wpEffArray[6])/wpEffArray[6]);
                    */
                    
                    // bias Eff(Fdata,ShapeData) vs Eff(True)
                    if (doBiasStudy) {
                        title = "EffBiasTrue_VS_centerLimit_rightlimit_WP_"+s.str();
                        
                        if (pullHistos2D.find(title) == pullHistos2D.end())
                            pullHistos2D[title] = new TH2D(title.c_str(),(title+";CenterLimit;Right Limit").c_str(),(mid_max-mid_min)/bias_step,mid_min,mid_max,(right_max-mid_min-20)/bias_step,mid_min+20,right_max);
                        pullHistos2D[title]->Fill(centerleftlimit_,rightlimit_,(wpEffArray[8]-wpEffArray[0])/wpEffArray[0]);
                        
                        // bias Eff(Fmc,ShapeData) vs Eff(True)
                        title = "EffMCBiasTrue_VS_centerLimit_rightlimit_WP_"+s.str();
                        
                        if (pullHistos2D.find(title) == pullHistos2D.end())
                            pullHistos2D[title] = new TH2D(title.c_str(),(title+";CenterLimit;Right Limit").c_str(),(mid_max-mid_min)/bias_step,mid_min,mid_max,(right_max-mid_min-20)/bias_step,mid_min+20,right_max);
                        pullHistos2D[title]->Fill(centerleftlimit_,rightlimit_,(wpEffArray[6]-wpEffArray[0])/wpEffArray[0]);
                        
                        // bias Eff(Fmc,ShapeMC) vs Eff(True)
                        title = "EffMC_ShapeMC_BiasTrue_VS_centerLimit_rightlimit_WP_"+s.str();
                        
                        if (pullHistos2D.find(title) == pullHistos2D.end())
                            pullHistos2D[title] = new TH2D(title.c_str(),(title+";CenterLimit;Right Limit").c_str(),(mid_max-mid_min)/bias_step,mid_min,mid_max,(right_max-mid_min-20)/bias_step,mid_min+20,right_max);
                        pullHistos2D[title]->Fill(centerleftlimit_,rightlimit_,(wpEffArray[18]-wpEffArray[0])/wpEffArray[0]);
                        
                        // bias Eff(Fdata,ShapeMC) vs Eff(True)
                        title = "Eff_SchapeMC_BiasTrue_VS_centerLimit_rightlimit_WP_"+s.str();
                        
                        if (pullHistos2D.find(title) == pullHistos2D.end())
                            pullHistos2D[title] = new TH2D(title.c_str(),(title+";CenterLimit;Right Limit").c_str(),(mid_max-mid_min)/bias_step,mid_min,mid_max,(right_max-mid_min-20)/bias_step,mid_min+20,right_max);
                        pullHistos2D[title]->Fill(centerleftlimit_,rightlimit_,(wpEffArray[20]-wpEffArray[0])/wpEffArray[0]);
                        
                        title = "FBias_VS_centerLimit_rightlimit_WP_"+s.str();
                        
                        if (pullHistos2D.find(title) == pullHistos2D.end())
                            pullHistos2D[title] = new TH2D(title.c_str(),(title+";CenterLimit;Right Limit").c_str(),(mid_max-mid_min)/bias_step,mid_min,mid_max,(right_max-mid_min-20)/bias_step,mid_min+20,right_max);
                        pullHistos2D[title]->Fill(centerleftlimit_,rightlimit_,(wpEffArray[14]-wpEffArray[10])/wpEffArray[10]);
                        
                        title = "Binning_Xcheck_centerLimit_rightlimit_WP_"+s.str();
                        
                        if (pullHistos2D.find(title) == pullHistos2D.end())
                            pullHistos2D[title] = new TH2D(title.c_str(),(title+";CenterLimit;Right Limit").c_str(),(mid_max-mid_min)/bias_step,mid_min,mid_max,(right_max-mid_min-20)/bias_step,mid_min+20,right_max);
                        pullHistos2D[title]->Fill(centerleftlimit_,rightlimit_);
                    }
                    
                    if (doFBiasStudy) {
                        
                        title = "FBias_VS_chi2cut_WP_"+s.str();
                        
                        if (pullHistos.find(title) == pullHistos.end())
                            pullHistos[title] = new TH1D(title.c_str(),(title+";Chi2 cut;Bias F").c_str(),(F_max-F_min)/F_step,F_min,F_max);
                        pullHistos[title]->Fill(matchChiSquareCut_,(wpEffArray[14]-wpEffArray[10])/wpEffArray[10]);
                        

                        
                    }
                         
					                    
                    //cout << centerleftlimit << " " << rightlimit << endl; exit(1);
					//cout << (wpEffArray[8]-wpEffArray[6])/wpEffArray[6] << endl; exit(1);
                    
					bTagPseudoExpResults[wpArray[iwp]].push_back(std::pair<float,float>(wpEffArray[8],wpEffArray[9]));
					misTagPseudoExpResults[wpArray[iwp]].push_back(std::pair<float,float>(wpEffArray[24],wpEffArray[25]));
					
					if(verbosity_>1) cout << "*********************" << wpEffArray[8] << " +- " << wpEffArray[9] << endl;
					// exit(0);
					
					// correlate Epsilon B and sigma TTbar
					
					/*float exp = desiredIntLum_*23.3333*0.25176162;
					 
					 float nttsemimu = ((FitForXSResults[8]*0.8723)/0.25176162);
					 float nttsemimuMC = ((FitForXSResults[0]*0.8723)/0.25176162);
					 
					 float eff = nttsemimu*(27./4.)*(1/desiredIntLum_);
					 float effMC = nttsemimuMC*(27./4.)*(1/desiredIntLum_);
					 */
					
					
					// no btag!!!
					vector<double> sigmaRes;
					
					if (iwp ==0) {
						
						sigmaRes = calcXS(doPseudoExp_,1,0,chiSQCutEff,MLBCutEff,refSelEff,FitForXSResults[taggerArray[iwp]],desiredIntLum_);
						
						//cout << "Control sample - Expected: " << sigmaRes[4] << " Obtained: " << FitForXSResults[8] << " RATIO: " << FitForXSResults[8]/sigmaRes[4] << endl;
						//cout << "MC - Expected: " << sigmaRes[4] << " Obtained: " << FitForXSResults[0] << " RATIO: " << FitForXSResults[0]/sigmaRes[4] << endl;
						
						title = "MC_XS_VS_BtagEfficiency_PseudoExps";//_WP_"+s.str();
						
						if (XSHistos.find(title) == XSHistos.end())
							XSHistos[title] = new TH2D(title.c_str(),(title+";#epsilon_{b}^{estimated};#sigma_{t#bar{t}} (pb)").c_str(),40,0,2,100,100,300);
						
						XSHistos[title]->Fill(wpEffArray[8],sigmaRes[1]);
						
						bTagPseudoExpResultsForXS["XS"].push_back(sigmaRes[1]);
						bTagPseudoExpResultsForXS["eXS"].push_back(sigmaRes[2]);
						
					}
					
					if (taggerArray[iwp] != psLastTagger) {
						psLastTagger=taggerArray[iwp];
						psnTimesWP=1;
					}
					
					// btag cut
					
					if (verbosity_ > 0) cout << "++> Using btag cut to calc XS using array offset " << (psnTimesWP*5) << endl;
					
					stringstream nd, nwp; nd << taggerArray[iwp]; nwp << wpArray[iwp];
					
					// epsilon_b estimated FData
					
					sigmaRes.clear();
					
					sigmaRes = calcXS(doPseudoExp_,wpEffArray[8],wpEffArray[9],chiSQCutEff,MLBCutEff,refSelEff,FitForXSResults[taggerArray[iwp]],desiredIntLum_,(psnTimesWP*5));
					
					title = "MC_XS_VS_BtagEfficiency_PseudoExps_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str();
					
					if (XSHistos.find(title) == XSHistos.end())
						XSHistos[title] = new TH2D(title.c_str(),(title+";#epsilon_{b}^{estimated};#sigma_{t#bar{t}} (pb)").c_str(),40,0,2,100,100,300);
					
					XSHistos[title]->Fill(wpEffArray[8],sigmaRes[1]);
					
					ofstream chiSqfile;
					chiSqfile.open ("tmpcheck.txt", ios::out | ios::app );
					
					if (nwp.str() == "1.7")
					chiSqfile << nwp.str() << " " << wpEffArray[8] << " " << sigmaRes[1] << "\n";
					
					chiSqfile.close();
					
					bTagPseudoExpResultsForXS["XS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str()].push_back(sigmaRes[1]);
					bTagPseudoExpResultsForXS["eXS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str()].push_back(sigmaRes[2]);
					
					// epsilon_b estimated FMC
					
					sigmaRes.clear();
					
					sigmaRes = calcXS(doPseudoExp_,wpEffArray[6],wpEffArray[7],chiSQCutEff,MLBCutEff,refSelEff,FitForXSResults[taggerArray[iwp]],desiredIntLum_,(psnTimesWP*5));
					
					title = "MC_XS_VS_BtagEfficiency_PseudoExps_FMC_nDisc_"+nd.str()+"_WP_"+nwp.str();
					
					if (XSHistos.find(title) == XSHistos.end())
						XSHistos[title] = new TH2D(title.c_str(),(title+";#epsilon_{b}^{estimated};#sigma_{t#bar{t}} (pb)").c_str(),40,0,2,100,100,300);
					
					XSHistos[title]->Fill(wpEffArray[6],sigmaRes[1]);
					
					sigmaRes.clear();
					
					// epsilon_b MC Truth
					
					sigmaRes = calcXS(doPseudoExp_,wpEffArray[0],wpEffArray[1],chiSQCutEff,MLBCutEff,refSelEff,FitForXSResults[taggerArray[iwp]],desiredIntLum_,(psnTimesWP*5));
					
					title = "MC_XS_VS_BtagEfficiency_PseudoExps_BTAGMCTRUTH_nDisc_"+nd.str()+"_WP_"+nwp.str();
					
					if (XSHistos.find(title) == XSHistos.end())
						XSHistos[title] = new TH2D(title.c_str(),(title+";#epsilon_{b}^{estimated};#sigma_{t#bar{t}} (pb)").c_str(),40,0,2,100,100,300);
					
					XSHistos[title]->Fill(wpEffArray[0],sigmaRes[1]);
					
					psnTimesWP++;
					
					//exit(1);
					// plots to compare with e+jets
					
					std::map<string,float> EffCalcDetails_ = aContainer->GetEffCalcDetails(wpArray[iwp],0,0,0,false);
					
					if (iwp == 0) {
						
						bTagPseudoExpResultsForF["F"].push_back(EffCalcDetails_["F"]);
						bTagPseudoExpResultsForF["Fe"].push_back(EffCalcDetails_["Fe"]);
						bTagPseudoExpResultsForF["Fmc"].push_back(EffCalcDetails_["Fmc"]);
						bTagPseudoExpResultsForF["Fmce"].push_back(EffCalcDetails_["Fmce"]);
						
					}
					
					/*for (std::map<string,float>::const_iterator it = EffCalcDetails_.begin(); it != EffCalcDetails_.end(); ++it)
					 cout << "From myNupleAnalyzer.cc " << it->first << " -> " << it->second << endl;
					 
					 title = "Events_Below_Cut_WP_"+s.str();
					 
					 if (pullHistos.find(title) == pullHistos.end())
					 pullHistos[title] = new TH1D(title.c_str(),(title+";#A;pseudo experiments").c_str(),5000,-5000,5000);
					 pullHistos[title]->Fill(EffCalcDetails_["min"]);
					 
					 title = "FMCEvents_Below_Cut_WP_"+s.str();
					 
					 if (pullHistos.find(title) == pullHistos.end())
					 pullHistos[title] = new TH1D(title.c_str(),(title+";#A;pseudo experiments").c_str(),5000,-5000,5000);
					 pullHistos[title]->Fill(EffCalcDetails_["minMC"]);
					 
					 title = "Events_Above_Cut_WP_"+s.str();
					 
					 if (pullHistos.find(title) == pullHistos.end())
					 pullHistos[title] = new TH1D(title.c_str(),(title+";#B;pseudo experiments").c_str(),5000,-5000,5000);
					 pullHistos[title]->Fill(EffCalcDetails_["plus"]);
					 
					 title = "FMCEvents_Above_Cut_WP_"+s.str();
					 
					 if (pullHistos.find(title) == pullHistos.end())
					 pullHistos[title] = new TH1D(title.c_str(),(title+";#B;pseudo experiments").c_str(),5000,-5000,5000);
					 pullHistos[title]->Fill(EffCalcDetails_["plusMC"]);
					 
					 title = "F_DataDriven_WP"+s.str();
					 
					 if (pullHistos.find(title) == pullHistos.end())
					 pullHistos[title] = new TH1D(title.c_str(),(title+";F;pseudo experiments").c_str(),40,0,4);
					 pullHistos[title]->Fill(EffCalcDetails_["F"]);
					 
					 title = "F_MC_WP"+s.str();
					 
					 if (pullHistos.find(title) == pullHistos.end())
					 pullHistos[title] = new TH1D(title.c_str(),(title+";F_{MC};pseudo experiments").c_str(),40,0,4);
					 pullHistos[title]->Fill(EffCalcDetails_["Fmc"]);*/
					
					//cout << "/////------------------------------------------------------------------------------------------------------" << EffCalcDetails_["Fmc"]<<endl;
					
					
					//exit(0);
					
					
				}
			}
		}
		if(verbosity_>1) cout << "************************************************"<< endl << endl;;
		
		if (doPseudoExp_) {
			
			string title;
			
			/*vector<float> sigmaRes = calcXS(doPseudoExp_,1,FitForXSResults,desiredIntLum_);
			 
			 string title = "BIAS_XS";
			 
			 if (pullHistos.find(title) == pullHistos.end())
			 pullHistos[title] = new TH1D(title.c_str(),(title+";nTTbar fitted - nTTbar Input;pseudo experiments").c_str(),500,-1000,1000);
			 
			 pullHistos[title]->Fill(FitForXSResults[8]-nTTbarFilled);
			 
			 title = "BIAS_XS_MC";
			 
			 if (pullHistos.find(title) == pullHistos.end())
			 pullHistos[title] = new TH1D(title.c_str(),(title+";nTTbar fitted - nTTbar Input;pseudo experiments").c_str(),500,-1000,1000);
			 
			 pullHistos[title]->Fill(FitForXSResults[0]-nTTbarFilled);*/
			
			/////////
			
			std::map<string,float> EffCalcDetails_ = aContainer->GetEffCalcDetails(1.7,0,0,0,false);
			
			title = "NL";
			
			if (pullHistos.find(title) == pullHistos.end())
				pullHistos[title] = new TH1D(title.c_str(),(title+";# events left;pseudo experiments").c_str(),500,-1000,1000);
			pullHistos[title]->Fill(EffCalcDetails_["NLeft"]);
			
			title = "NR";
			
			if (pullHistos.find(title) == pullHistos.end())
				pullHistos[title] = new TH1D(title.c_str(),(title+";# events right;pseudo experiments").c_str(),500,-1000,1000);
			pullHistos[title]->Fill(EffCalcDetails_["NRight"]);
			
			//cout << "------------------------------------------------------------------------------------------------------" << EffCalcDetails_["NLeft"]<<endl;
			
			//exit(0);
		}
		
		if(verbosity_>0 || doPseudoExp_) cout<<"+> deleting the pointer "<< endl;
		//RM::delete Container1;
		
		delete aContainer;
		
		/*if (!doPseudoExp_) {
		 cout << "+> Checking the ttbar/non-ttbar ratio in the pseudo-data mlj distribution" << endl;
		 
		 fout->cd();
		 for (int p=0; p<1; p++) {
		 
		 stringstream s; s<<p;
		 
		 cout << endl << "  ====== Tagger " << p << " ====== " << endl;
		 
		 float IntTTbar = MljDistributions["Mlj_SignalSample_TTbar"]->Integral();
		 float IntOther = MljDistributions["Mlj_SignalSample_Other"]->Integral();
		 
		 cout << endl << "    -> || NO Cut || nTTbar " << IntTTbar << " || nOther " << IntOther << " || nOther/(nOther+nTTbar) " << IntOther/(IntOther+IntTTbar) << " || " << endl;
		 
		 IntTTbar = MljDistributions["Mlj_SignalSample_TTbar_bTagger_"+s.str()+"_bTagL"]->Integral();
		 IntOther = MljDistributions["Mlj_SignalSample_Other_bTagger_"+s.str()+"_bTagL"]->Integral();
		 
		 cout << endl << "    -> || Loose Cut || nTTbar " << IntTTbar << " || nOther " << IntOther << " || nOther/(nOther+nTTbar) " << IntOther/(IntOther+IntTTbar) << " || " << endl;
		 
		 IntTTbar = MljDistributions["Mlj_SignalSample_TTbar_bTagger_"+s.str()+"_bTagM"]->Integral();
		 IntOther = MljDistributions["Mlj_SignalSample_Other_bTagger_"+s.str()+"_bTagM"]->Integral();
		 
		 cout <<  endl << "    -> || Medium Cut || nTTbar " << IntTTbar << " || nOther " << IntOther << " || nOther/(nOther+nTTbar) " << IntOther/(IntOther+IntTTbar) << " || " << endl;
		 
		 
		 IntTTbar = MljDistributions["Mlj_SignalSample_TTbar_bTagger_"+s.str()+"_bTagT"]->Integral();
		 IntOther = MljDistributions["Mlj_SignalSample_Other_bTagger_"+s.str()+"_bTagT"]->Integral();
		 
		 cout << endl << "    -> || Tight Cut || nTTbar " << IntTTbar << " || nOther " << IntOther << " || nOther/(nOther+nTTbar) " << IntOther/(IntOther+IntTTbar) << " || " << endl;
		 
		 }cout << endl;
		 }*/
		
		if(verbosity_>0) cout<<"+> closing the root-file " << outRootPath[0] << endl;
		fout->Close();
		
		//cout << " +> nttbar in histogram to do template fit on = " << nTTbarFilled << endl;
		
		//nTTbarFilled = 0;
		
	}	// end of PSEUDOEXP loop
	
	if(verbosity_>0) cout<<"+> Now make the histogram of the pull for each WP" <<endl;
	
	// do pullhistos for btag & btag unc
	if (doPseudoExp_) {
		
		foutPseudoExp->cd();
		
		for(int iwp=0; iwp < 3; iwp++){
			
			int nTot = 0;
			int nNeg = 0;
			
			cout << "++> Making the PullHisto for WP " << wpArray[iwp] << endl;
			
			float meanExp = 0;
			for (unsigned int w=0; w<bTagPseudoExpResults[wpArray[iwp]].size(); w++) {
				
				//if (fabs(bTagPseudoExpResults[wpArray[iwp]][w].first > 1))
				//	continue;
				
				//cout << "  Experiment " << w << " result: " << bTagPseudoExpResults[wpArray[iwp]][w].first << " +- " << bTagPseudoExpResults[wpArray[iwp]][w].second << endl;
				
				meanExp += bTagPseudoExpResults[wpArray[iwp]][w].first;
				
			}
			
			meanExp = meanExp/bTagPseudoExpResults[wpArray[iwp]].size();

			float meanExpMisTag = 0;
			for (unsigned int w=0; w<misTagPseudoExpResults[wpArray[iwp]].size(); w++) {
				
				meanExpMisTag += misTagPseudoExpResults[wpArray[iwp]][w].first;
				
			}
			
			meanExpMisTag = meanExpMisTag/misTagPseudoExpResults[wpArray[iwp]].size();

			cout << "  Mean of the btag experiments: " << meanExp << endl;
			cout << "  Mean of the mistag experiments: " << meanExpMisTag << endl;
			
			stringstream s; s << wpArray[iwp];
			string title = "Pull_WP_"+s.str();
			string title2 = "Pull_MisTag_WP_"+s.str();
			
			if (pullHistos.find(title) == pullHistos.end())
				pullHistos[title] = new TH1D(title.c_str(),title.c_str(),40,-10,10);

			if (pullHistos.find(title2) == pullHistos.end())
				pullHistos[title2] = new TH1D(title2.c_str(),title2.c_str(),40,-10,10);

			for (unsigned int w=0; w<bTagPseudoExpResults[wpArray[iwp]].size(); w++) {
				
				//if (fabs(bTagPseudoExpResults[wpArray[iwp]][w].first > 1))
				//	continue;
				
				nTot++;
				float pull = (bTagPseudoExpResults[wpArray[iwp]][w].first-meanExp)/bTagPseudoExpResults[wpArray[iwp]][w].second ;
				float mpull = (misTagPseudoExpResults[wpArray[iwp]][w].first-meanExpMisTag)/misTagPseudoExpResults[wpArray[iwp]][w].second ;
				
				if (pull < 0) {
					cout << "  ++> Pull experiment " << w << ": (" << bTagPseudoExpResults[wpArray[iwp]][w].first << " - " << meanExp << ") / " << bTagPseudoExpResults[wpArray[iwp]][w].second << " = " << pull << endl;
					nNeg++;
				}
				pullHistos[title]->Fill((bTagPseudoExpResults[wpArray[iwp]][w].first-meanExp)/bTagPseudoExpResults[wpArray[iwp]][w].second);
				pullHistos[title2]->Fill(mpull);
				
			}
			
			cout << "+> Of the " << nTot << " experiments, " << nNeg << " gave a negative pull (" << (double)nNeg/(double)nTot << "%)" << endl;
			
			cout << endl;
            
            /*title = "EffBias_VS_centerLimit_rightlimit_WP_"+s.str();
            
            double avgEffBias=0;
            for (int res=0;res<bTagEffBiasResults.size();res++)
                avgEffBias+=bTagEffBiasResults[res];
            avgEffBias=avgEffBias/bTagEffBiasResults.size();
            
            if (pullHistos2D.find(title) == pullHistos2D.end())
                pullHistos2D[title] = new TH2D(title.c_str(),(title+";CenterLimit;Right Limit").c_str(),16,120,200,52,140,600);
            pullHistos2D[title]->Fill(centerleftlimit_,rightlimit_,avgEffBias);*/

			
		}
		
		cout << endl;
		
	}	
	
	// do pullhistos for F
	
	if (doPseudoExp_) {
		
		foutPseudoExp->cd();
		
		cout << "++> Making the PullHisto for F " << endl;
		
		float meanExpF = 0;
		float meanExpFmc = 0;
		for (unsigned int w=0; w<bTagPseudoExpResultsForF["F"].size(); w++) {
			
			//if (fabs(bTagPseudoExpResults[wpArray[iwp]][w].first > 1))
			//	continue;
			
			//cout << " FOR F  Experiment " << w << " result: " << bTagPseudoExpResultsForF["F"][w] << " +- " << bTagPseudoExpResultsForF["Fe"][w] << endl;
			//cout << " FOR FMC  Experiment " << w << " result: " << bTagPseudoExpResultsForF["Fmc"][w] << " +- " << bTagPseudoExpResultsForF["Fmce"][w] << endl;
			
			meanExpF += bTagPseudoExpResultsForF["F"][w];
			meanExpFmc += bTagPseudoExpResultsForF["Fmc"][w];
			
		}
		
		meanExpF = meanExpF/bTagPseudoExpResultsForF["F"].size();
		meanExpFmc = meanExpFmc/bTagPseudoExpResultsForF["Fmc"].size();
		
		cout << "  Mean of the experiments for F: " << meanExpF << endl;
		cout << "  Mean of the experiments for Fmc: " << meanExpFmc << endl;
		
		string title = "Pull_F_DataDriven";
		string titleMC = "Pull_F_MC";
		
		if (pullHistos.find(title) == pullHistos.end())
			pullHistos[title] = new TH1D(title.c_str(),title.c_str(),160,-10,10);
		
		if (pullHistos.find(titleMC) == pullHistos.end())
			pullHistos[titleMC] = new TH1D(titleMC.c_str(),titleMC.c_str(),160,-10,10);
		
		for (unsigned int w=0; w<bTagPseudoExpResultsForF["F"].size(); w++) {
			
			//if (fabs(bTagPseudoExpResults[wpArray[iwp]][w].first > 1))
			//	continue;
			
			float pull = (bTagPseudoExpResultsForF["F"][w]-meanExpF); ///bTagPseudoExpResultsForF["Fe"][w];
			float pullmc = (bTagPseudoExpResultsForF["Fmc"][w]-meanExpFmc); ///bTagPseudoExpResultsForF["Fmce"][w];
			
			//if (pull < 0) {
			//cout << "  ++> F Pull experiment " << w << ": (" << bTagPseudoExpResultsForF["F"][w] << " - " << meanExpF << ") / " << bTagPseudoExpResultsForF["Fe"][w] << " = " << pull << endl;
			//cout << "  ++> Fmc Pull experiment " << w << ": (" << bTagPseudoExpResultsForF["Fmc"][w] << " - " << meanExpFmc << ") / " << bTagPseudoExpResultsForF["Fmce"][w] << " = " << pull << endl;
			//}
			pullHistos[title]->Fill(pull);
			pullHistos[titleMC]->Fill(pullmc);
			
		}
		
		//cout << "+> Of the " << nTot << " experiments, " << nNeg << " gave a negative pull (" << (double)nNeg/(double)nTot << "%)" << endl;
		
		cout << endl;
		
	}
	
	// do pullhistos for XS
	
	if (doPseudoExp_) {
		
		cout << "++> Making the PullHisto for XS " << endl;
		
		for(int iwp=0; iwp < nWP; iwp++){
			
			stringstream nd, nwp; nd << taggerArray[iwp]; nwp << wpArray[iwp];
			
			float meanExp = 0;
			
			for (unsigned int w=0; w<bTagPseudoExpResultsForXS["XS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str()].size(); w++) {
				
				cout << "XS for exp " << w << " iwp " << iwp << " " << bTagPseudoExpResultsForXS["XS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str()][w] << endl;
				meanExp += bTagPseudoExpResultsForXS["XS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str()][w];
				
			}
			
			meanExp = meanExp/bTagPseudoExpResultsForXS["XS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str()].size();
			
			//cout << "  IWP " << iwp << " -> Mean XS in experiments " << meanExp << endl;
			
			string title = "Pull_XS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str();
			string title2 = "UNC_XS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str();
			string title3 = "XS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str();
			
			if (pullHistos.find(title) == pullHistos.end())
				pullHistos[title] = new TH1D(title.c_str(),title.c_str(),40,-10,10);
			
			if (pullHistos.find(title2) == pullHistos.end())
				pullHistos[title2] = new TH1D(title2.c_str(),(title2+";#delta#hat{#sigma_{t#bar{t}}} (pb)").c_str(),50,0,100);
			
			if (pullHistos.find(title3) == pullHistos.end())
				pullHistos[title3] = new TH1D(title3.c_str(),(title3+";#hat{#sigma_{t#bar{t}}} (pb)").c_str(),50,0,300);
			
			for (unsigned int w=0; w<bTagPseudoExpResultsForXS["XS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str()].size(); w++) {
				
				float pull = (bTagPseudoExpResultsForXS["XS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str()][w]-meanExp)/bTagPseudoExpResultsForXS["eXS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str()][w];
				
				pullHistos[title]->Fill(pull);
				pullHistos[title2]->Fill(bTagPseudoExpResultsForXS["eXS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str()][w]);
				pullHistos[title3]->Fill(bTagPseudoExpResultsForXS["XS_FDATA_nDisc_"+nd.str()+"_WP_"+nwp.str()][w]);
				
			}
			
		}
		
	}
	
	if(verbosity_>0) cout<<"+> writing the Pseudo Exp plots Still need to do this for all taggers/wp's " << outRootPath[1] << endl;
	foutPseudoExp->cd();
	
	for(std::map<std::string,TH1D*>::const_iterator it = pullHistos.begin(); it != pullHistos.end(); it++) {
		
		it->second->Write();
		
		//cout << "++> Writing " << it->first << endl;
		
		//delete pullHistos[it->first];
	}
	
	for(std::map<std::string,TH2D*>::const_iterator it = pullHistos2D.begin(); it != pullHistos2D.end(); it++) {
		
		it->second->Write();
		
		//cout << "++> Writing " << it->first << endl;
		
		//delete pullHistos[it->first];
	}
    
    for(std::map<std::string,TH3D*>::const_iterator it = pullHistos3D.begin(); it != pullHistos3D.end(); it++) {
		
		it->second->Write();
		
		//cout << "++> Writing " << it->first << endl;
		
		//delete pullHistos[it->first];
	}
	
	
	for(std::map<std::string,TH2D*>::const_iterator it = XSHistos.begin(); it != XSHistos.end(); it++) {
		it->second->Write();
		
		//cout << "++> Writing " << it->first << endl;
		
		delete XSHistos[it->first];
	}
	
	if(verbosity_>0) cout<<"+> closing the Pseudo Exp root-file " << outRootPath[1] << endl;
	foutPseudoExp->Close();
	
	cout << "written: " << outRootPath[0] << endl;
	cout << "written: " << outRootPath[1] << endl;
	cout << "finished NTupleAnalyzer.C" << endl;

	cout << "+> RefSel cut efficiency (MC based) " << endl;
	cout << " ++> nttbar before = " << nTTbarBeforeRefSel << endl;
	cout << " ++> nttbar after = " << nTTbarAfterRefSel << endl;
	cout << " ++> nttbar after2 = " << nTTbarAfterRefSel2 << endl;
	cout << " ++> refsel eff = " << (double)nTTbarAfterRefSel/(double)nTTbarBeforeRefSel  << endl;	
	
	cout << "+> CHI2 cut efficiency (MC based) (after refsel)" << endl;
	cout << " ++> nttbar before cut = " << nTTbarBeforeChiSq << endl;
	cout << " ++> nttbar after cut = " << nTTbarAfterChiSq << endl;
	cout << " ++> cut eff = " << (double)nTTbarAfterChiSq/(double)nTTbarBeforeChiSq  << endl;	
	
	cout << "+> MLB cut efficiency (MC based) (after chi2)" << endl;
	cout << " ++> nttbar before cut = " << nTTbarBeforeMLBCUT << endl;
	cout << " ++> nttbar after cut = " << nTTbarAfterMLBCUT << endl;
	cout << " ++> cut eff = " << (double)nTTbarAfterMLBCUT/(double)nTTbarBeforeMLBCUT  << endl;	
	

	cout << endl << " ++> nttbar GOOD combi after cut: " << 	nGoodChi2Combi << " on total " << nChi2Combi << " combis = " << (float)nGoodChi2Combi/(float)nChi2Combi << "%" << endl;

	cout << " " << endl;
	cout << "+> Selection yields" << endl;

	fstream yields("yield_selection.csv", ios::out | ios::trunc);
	yields << "Sample Ntot Nsel Effsel Yield NselBtag EffselBtag Yield\n\n";
	for (std::map<std::string, std::vector<float> >::const_iterator it=nSelected_.begin(); it != nSelected_.end(); ++it) {
		
		cout << it->first << " - # events after refsel: " << it->second[2] << endl;
		cout << it->first << " - # events after refsel (no SF): " << it->second[4] << endl;
		cout << it->first << " - (L="<<it->second[0]<<")# events after refsel: " << it->second[2]*(it->second[0]/it->second[1]) << endl;
		cout << it->first << " - # events after refsel+btag CSVM: " << it->second[3] << endl;
		cout << it->first << " - # events after refsel+btag CSVM (no SF): " << it->second[5] << endl;
		cout << it->first << " - (L="<<it->second[0]<<")# events after refsel+btag CSVM: " << it->second[3]*(it->second[0]/it->second[1]) << endl;
		
		yields << it->first << " - " << it->second[4] << " - " << it->second[2]*(it->second[0]/it->second[1]) << " ";
		yields << it->second[5] << " - " << it->second[3]*(it->second[0]/it->second[1]) << "\n";

		cout << "" << endl;
	} yields.close();
	
	time_t currEnd=time(0);
	cout << "current end time is: " << ctime(&currEnd) <<endl;
	
	time_t currDiff=currEnd-curr;
	cout << "calculation time is: " << currDiff << " seconds = " << (double) currDiff/60 << " minutes" <<endl;
	cout << "time per Pseudo Exp: " << (double) currDiff/ (double) nPseudoExp_ << " seconds" << endl;
	
	
}

