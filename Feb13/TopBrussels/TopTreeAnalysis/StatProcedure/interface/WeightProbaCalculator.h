#ifndef WeightProbaCalculator_h
#define WeightProbaCalculator_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TObject.h"
#include "TMath.h"

#include "RooStats/HybridCalculator.h" /*just include this because otherwise compiling complains about namespace*/

using namespace std;
using namespace RooStats;

/**
	//	Aim: Calculation of p-values for pseudoexperiments or real data, creating corresponding V and p distributions, 
	//	and calculating p95 or P05Data.\n
	//
	//	A p value of 'data' (= (pseudo)experiment) is the integral of a reference V distribution larger than the corresponding V value
	//	of the 'data', divided by the total integral of the reference V distribution. p95 is  the value of p for which 95% of the pseudoexperiments 
	//	result in a smaller p value. P05Data is the probability for 'data' to have a p value smaller than 0.05.\n
	//
	//	For explanation of the V values, see other class SampleCombinedWeightCalculator. 
	//
	//	Example of usage:
	//	- call constructor
	//	- initialize the distributions with the Init() functions
	//	- CalculateProba() in a loop over 'data'  V values, and fill a vector with the resulting GetProba() value
	//	- FillDistributions()
	//	- WriteDistributions()
	//	- get the value of p95 with GetProba95() or tha value of P05Data with GetP05Data()

*/

struct sortDescending {
     bool operator() (float a, float b) { return (a>b);}
};

struct sortAscending {
     bool operator() (float a, float b) { return (a<b);}
};

class WeightProbaCalculator {
 public:
  
     WeightProbaCalculator();
     /** default constructor */
     WeightProbaCalculator(string refmcname, string dataname);
     /** names will appear in the titles of the created histograms */
     ~WeightProbaCalculator();
     /** destructor */
     
     void InitRefVDistr(int nbinsx, double xlow, double xup);  /** initialize reference Monte Carlo V distribution*/   
     void InitVDistr(int nbinsx, double xlow, double xup); /** initialize 'data' (= (pseudo)experiment(s)) V distribution*/
     void InitpDistr(int nbinsx, double xlow, double xup); /** initialize ('data') p distribution*/
     void InitNofSigmaDistr(int nbinsx, double xlow, double xup); /** initialize ('data') sigma distribution*/
      
     void CalculateProba(vector<float> Vvalues_MC, float Vvalue_data); /** main method; calculate p value*/    
     void FillDistributions(vector<float> Vvalues_MC, vector<float> Vvalues_data, vector<float> pvalues);
     /** fill V (both reference and 'data') and p distributions*/
     void WriteDistributions(); /** write distributions*/
     float GetProba95(vector <float> pvalues); /** get p95; the value for which 95% of the pseudoexperiments result in a smaller p value*/
     float GetP05Data(vector <float> pvalues); /** get P05Data; the probability for 'data' to have a p value smaller than 0.05*/
     float GetProbaPercentile(vector <float> pvalues, double percent);
     /** get the percentile (related to the specified second argument) of p-values, i.e. the p-value below which a percentage of all p-values fall.
    The 50th percentile is the median, and in principle, the 95th percentile is p95*/
     
     string GetRefMCName() const {return RefMCName_;}; /** get name of reference (Monte Carlo) dataset(s)*/
     string GetDataName() const {return DataName_;}; /** get name of 'data'set(s)*/
     TH1F* GetRefVDistribution() const {return RefVDistribution_;}; /** get reference (Monte Carlo) V distribution*/
     TH1F* GetVDistribution() const {return VDistribution_;}; /** get 'data' V distribution*/
     TH1F* GetpDistribution() const {return pDistribution_;}; /** get ('data') p distribution*/
     TH2F* GetCorrelationHisto_pV() const {return CorrelationHisto_pV_;}; 
     /** get correlation histogram between p and V of the 'data'. Can be used, but most of the time not interesting*/
     float GetProba() const {return Proba_;}; /** get p-value*/
     float GetSignificance() const {return Significance_;}; /** get the significance, the familiar name for the null p-value in terms of 1-sided Gaussian significance*/
     float ProbaToSignificance(float proba); /** convert manually a p-value to significance*/
     float SignificanceToProba(float significance); /** convert manually a significance to p-value*/
     float PowerOfTest(vector<float> Vvalues_data, float Valfa); /** calculate and get the power of the test, 1-beta, given a V distribution and V_alfa*/
     bool isNewPhysics() {return isNewPhysics_;};  /** not used*/ 
          
 private:
     void InitCorrHisto_pV(); 
     /** private method to initialize correlation histogram between p and V according to previous initializations. Can be used, but most of the time not interesting*/
     void SortVectorDescending(vector<float> &pvalues); /** sorts vector of p values such that the first element is the largest*/
     void SortVectorAscending(vector<float> &pvalues); /** sorts vector of p values such that the first element is the smallest*/
      
     string RefMCName_;
     string DataName_;
     TH1F* RefVDistribution_; 
     TH1F* VDistribution_;
     TH1F* pDistribution_;
     TH1F* NofSigmaDistribution_;
     TH2F* CorrelationHisto_pV_;
     float Proba_;
     float Significance_;
     bool isNewPhysics_;
};

#endif
