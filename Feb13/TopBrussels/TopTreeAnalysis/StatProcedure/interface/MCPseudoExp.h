#ifndef MCPseudoExp_h
#define MCPseudoExp_h

#include "TopTreeAnalysis/Content/interface/Dataset.h"
#include "TopTreeAnalysis/Content/interface/AnalysisEnvironment.h"
/*#include "TRandom1.h"*/  /*Ranlux Random number generator class with periodicity > 10**14*/
#include "TRandom3.h"  /*Random number generator class based on M. Matsumoto and T. Nishimura, Mersenne Twistor, period 2**19937-1*/
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "TChain.h"
#include "TTree.h"

using namespace std;

/**
	//	Aim: Selecting events from a Monte Carlo sample to do a pseudoexperiment.
	//
	//	Example of usage:
	//	For example within a loop over pseudoexperiments, and before a loop over events:
	//	- call constructor
	//	- CalculatePseudoExp()
	//	- obtain the resulting vector of integers GetRandomNumbers(), which can be used inside a loop over events as a list of entry
	//	numbers to pick out events e.g.: 
	//	eventTree->GetEvent(myMCPseudoExp.GetRandomNumbers()[ievt]), where eventTree is a TChain* for instance.
	//	Attention: the vector of integers retrieved by GetRandomNumbers() is ranked (from low to high); this seemed necessary due
	//	to the way ROOT runs over the eventTree (otherwise much slower).

*/
struct rankstruct {
     bool operator() (int a, int b) { return (a<b);}
};

class MCPseudoExp {

 public:  
     MCPseudoExp();
     /** default constructor */
     ~MCPseudoExp();
     /** destructor */   

     void CalculatePseudoExp(const Dataset* dataset, float& luminosity);
     /** main method; calculates the amount events to pick out of a dataset according to a luminosity of the dataset, and generates 
     corresponding random integers which can be used in a later stage to effectively pick out the events (doing the pseudoexperiment)*/
     void CalculatePseudoExp(const Dataset* dataset, float& luminosity, int& ievtmin, bool  nPseudoSeries=false);
     /** This constructor is dedicated to the use of continuous events in order to increase speed. 
     The idea is to split the sample in a continuous block starting from the event ievtmin. 
     The number of events is randomly taken from a poisson distribution.
     If there is no event left or if ievtmin>sample.size, this constructor is equivalent to the previous one.*/
     void CalculatePseudoExp(TTree* seleventtree, float& force_luminosity, const Dataset* dataset, AnalysisEnvironment anaEnv);
     /** when running on 'TTrees' with (selected) events (inside a file corresponding to fout) rather than toptrees*/
     vector<int> GetRandomNumbers() const {return RandomNumbers_;}; 
     /** vector of random numbers generated inside the CalculatePseudoExp() method*/
     float GetLuminosity() const {return Luminosity_;}; /** get the specified luminosity*/
     float GetPreSelEfficiency() const {return PreSelEfficiency_;}; /** get the preselection efficiency of the used dataset*/
     float GetXsection() const {return Xsection_;}; /** get the cross section of the used dataset (in 1/pb)*/
     float GetMeanNEvents() const {return MeanNEvents_;}; /** get the mean number of events to be picked in the pseudoexperiment*/
     int GetNEvents() const {return NEvents_;}; /** get the actual number of picked events in the pseudoexperiment*/
     
 private:
     void GenerateRandomNumbers(int nNumbers, int eventtree_nentries); /** generate random integers*/
     void RankVector(vector<int> &randomnumbers); /** rank the randomnumbers vector from low (first element) to high*/
     vector<int> RandomNumbers_;
     float Luminosity_;
     float PreSelEfficiency_;
     float Xsection_; 
     float MeanNEvents_;
     int NEvents_; 
};

#endif
