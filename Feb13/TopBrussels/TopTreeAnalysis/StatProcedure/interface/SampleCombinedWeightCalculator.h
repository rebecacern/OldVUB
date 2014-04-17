#ifndef SampleCombinedWeightCalculator_h
#define SampleCombinedWeightCalculator_h

#include <iostream>
#include <iomanip>
#include <vector>
#include "TMath.h"
#include "TString.h"

using namespace std;

/**
	//	Aim: Calculation of a combined weight V for a sample (list of events). 
	//
	//	This V value is an estimator which later can be used for hypothesis tests (see other class WeightProbaCalculator) of the consistency
	//	of data with the Standard Model.
	//
	//	Example of usage:
	//	For example, within a loop over pseudoexperiments, and after a loop over events:
	//	- call constructor
	//	- CalculateV() according to a vector of weights S which result from the EventCombinedWeightCalculator class in a loop over
	//	selected events, or according to a vector of 'eventinfos', to take all information of the events (weight, New Physics, values of variables) along
	//	during the ranking of the events
	//	- fill a vector with the resulting GetV() value

*/
   
struct comparestruct {
     bool operator() (float a, float b) { return (a>b);}
};

struct eventinfo {
    float weight;
    unsigned int isnp;
    vector<float>* variableValues;
}; /** this struct groups information of an event: the values of the variables, the weight, and if it is New Physics (isnp = 1) or not (isnp = 0)*/

struct comparestruct_eventinfo {
     bool operator() (eventinfo a, eventinfo b) { return (a.weight>b.weight);}
};

class SampleCombinedWeightCalculator {

 public: 
     SampleCombinedWeightCalculator();
     /** default constructor */
     ~SampleCombinedWeightCalculator();     
     /** destructor */

     void RankEvents(vector<float>& weights_per_event); /** ranks the events such that the first element is the highest weight*/
     void RankEvents(vector<eventinfo>& eventsinfos); /** use to take other event info along during ranking*/
     void CalculateV(vector<float>& weights_per_event, TString combination, float fraction);
     /** main method; calculate the combined weight V for the sample, according to a combination specified by a TString and a fraction of the 
     events with highest weight*/
     void CalculateV(vector<eventinfo>& eventsinfos, TString combination, float fraction); /** use to take other event info along*/
     vector<float> GetWeights_per_event() const {return Weights_per_event_;};
     vector<eventinfo> GetEventsInfos() const {return EventsInfos_;};
     TString GetCombination() const {return Combination_;};
     /** get the string which specifies how to combine the weigths of the events to obtain one weight for the sample*/
     float GetFraction() const {return Fraction_;}; /** get the chosen fraction of the events with highest weight*/
     float GetV() const {return V_;}; /** get the V value (= the combined (w.r.t. the events) weight) of the sample*/
     unsigned int GetnEventsInFinalSample() const {return nEventsInFinalSample_;};
     /** get the total number of events in final sample*/
     unsigned int GetnEventsInHighestWeightSubsample() const {return nEventsInHWSubsample_;};
     /** get the total number of events in highest-weight subsample*/
     unsigned int GetnNPEventsInFinalSample() const {return nNPEventsInFinalSample_;};
     /** get the total number of New Physics events in final sample*/
     unsigned int GetnNPEventsInHWSubsample() const {return nNPEventsInHWSubsample_;};
     /** get the total number of New Physics events in highest-weight subsample*/
     float GetFracNPEventsInFinalSample() const {return FracNPEventsInFinalSample_;};
     /** get the fraction of events in the final sample that are New Physics events*/
     float GetFracNPEventsInHWSubsample() const {return FracNPEventsInHWSubsample_;};
     /** get the fraction of events in the highest-weight subsample that are New Physics events*/
     float GetFracNPEventsEndUpInHWSubsample() const {return FracNPEventsEndUpInHWSubsample_;};
     /** get the fraction of New Physics events in the final sample that end up in the highest-weight subsample*/
     float GetSignalSignificanceInHWSubSample() const {return SignalSignificanceInHWSubSample_;};
     /** get the significance of signal events, s/sqrt(b), in the highest-weight subsample*/
     
 private:
     comparestruct myComparison;
     comparestruct_eventinfo myComparison_eventinfo;
     vector<float> Weights_per_event_;
     vector<eventinfo> EventsInfos_;
     TString Combination_;
     float Fraction_;
     float V_;
     unsigned int nEventsInFinalSample_; /** total number of events in final sample*/
     unsigned int nEventsInHWSubsample_; /** total number of events in highest-weight subsample*/
     unsigned int nNPEventsInFinalSample_; /** total number of New Physics events in final sample*/
     unsigned int nNPEventsInHWSubsample_; /** total number of New Physics events in highest-weight subsample*/
     float FracNPEventsInFinalSample_; /** fraction of events in the final sample that are New Physics events*/
     float FracNPEventsInHWSubsample_; /** fraction of events in the highest-weight subsample that are New Physics events*/
     float FracNPEventsEndUpInHWSubsample_; /** fraction of New Physics events in the final sample that end up in the highest-weight subsample*/
     float SignalSignificanceInHWSubSample_; /** significance of signal events, s/sqrt(b), in highest-weight subsample*/
     void CalculatenEvents(); /** do calculations related to number of events in sample*/
};

#endif
