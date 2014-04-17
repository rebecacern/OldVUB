#ifndef EventCombinedWeightCalculator_h
#define EventCombinedWeightCalculator_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "TMath.h"
#include "TH1F.h"
#include "TopTreeAnalysis/Content/interface/MCObsExpectation.h"


using namespace std;

/**
	//	Aim: Calculation of combined weight for an event. 
	//	
	//	A bin significance is computed by computing the deviation of bin contents between the data and the estimation histograms (of an observable)
	//	to which the event belongs. Next the significances of all observables are combined for that event, to result in a combined weight S.
	//
	//	Example of usage: inside a loop over events;
	//	- call constructor
	//      - SetObsExp(), this is to set the MCObsExpectation objects that should contain the histograms that are needed
	//	- CalculateWeight()
	//	- fill a vector with the resulting GetCombinedWeight() value (this is the so-called 'S-weight')

*/

class EventCombinedWeightCalculator {
 private:

 public:  
     EventCombinedWeightCalculator(string SignificanceCombination);
     /** default constructor */
     ~EventCombinedWeightCalculator();
     /** destructor */    

     void SetSystematics(float LumiError, vector<MCObsExpectation*> MCObsExp_JESSyst, bool doLumiSyst, bool doXsectionSyst, bool doJESSyst);   
     /*void CalculateWeight(vector<pair < string, float > > & event_string, vector<TH1F*>& histos_Obs_data, vector<TH1F*>& histos_Obs_estim);*/
     void SetObsExp(vector<MCObsExpectation*> MCSMObsExp, vector<MCObsExpectation*> MCNPObsExp);
     void SetObsExp(vector<MCObsExpectation*> MCSMObsExp);
     void CalculateWeight(vector<pair < string, float > > & event_string, vector<TH1F*>& histos_Obs_data, vector < Dataset* > datasets);
     /** main method; calculate the combined weight for the event (event_string = list of values of observables, with the name of the observable)*/
     float GetSignificance() const {return Significance_;};
     /** get the squared significance for the bin to which the event belongs for an observable*/
     float GetCombinedWeight() const {return CombinedWeight_;};
     /** get the combined (w.r.t. the observables) weight for the event*/
     bool isUnderflowWarning() const {return isUnderflowWarning_;}; 
     /** a boolean to indicate if the specified event in the CalculateWeight method belongs to the underflow bin of at least one observable histogram*/

 private:
     void Init(); /** initializes some private members*/
     float Significance_;
     float CombinedWeight_;
     bool isUnderflowWarning_;
     string SignificanceCombination_;
     vector<MCObsExpectation*> MCSMObsExp_; /** should contain the observable histograms obtained by running on the Standard Model MC samples*/
     vector<MCObsExpectation*> MCNPObsExp_; /** should contain the observable histograms obtained by running on the Standard Model and New Physics MC samples*/
     float LumiError_;
     vector<MCObsExpectation*> MCObsExp_JESSyst_; /** should contain the observable histograms obtained by varying the JES...*/
     bool doLumiSyst_;
     bool doXsectionSyst_;
     bool doJESSyst_;    
};

#endif
