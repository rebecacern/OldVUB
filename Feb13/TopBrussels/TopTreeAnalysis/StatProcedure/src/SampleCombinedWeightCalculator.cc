#include "../interface/SampleCombinedWeightCalculator.h"
  
SampleCombinedWeightCalculator::SampleCombinedWeightCalculator(){
      //initialize private members
	Combination_ = TString("");
	Fraction_ = 0.1;
	V_ = 1.0;
	nEventsInFinalSample_ = 0;
	nEventsInHWSubsample_ = 0;
	nNPEventsInFinalSample_ = 0;
	nNPEventsInHWSubsample_ = 0;
	FracNPEventsInFinalSample_ = 0;
	FracNPEventsInHWSubsample_ = 0;
	FracNPEventsEndUpInHWSubsample_ = 0;
	SignalSignificanceInHWSubSample_ = 0;
}

SampleCombinedWeightCalculator::~SampleCombinedWeightCalculator(){}

void SampleCombinedWeightCalculator::RankEvents(vector<float>& weights_per_event) {
	sort(weights_per_event.begin(),weights_per_event.end(),myComparison);	
}

void SampleCombinedWeightCalculator::RankEvents(vector<eventinfo>& eventsinfos) {
	sort(eventsinfos.begin(),eventsinfos.end(),myComparison_eventinfo);	
}

void SampleCombinedWeightCalculator::CalculatenEvents(){
	// used to do calculations related to number of events in sample, when you take other event info along
	nEventsInFinalSample_ = EventsInfos_.size();
	nEventsInHWSubsample_ = int(Fraction_*EventsInfos_.size());
	if(nEventsInHWSubsample_ == 0){
		cout << "Warning: number of events in highest-weight subsample = 0!" << endl;
		nEventsInHWSubsample_++;
	}
	nNPEventsInFinalSample_ = 0;
	nNPEventsInHWSubsample_ = 0;
	for(unsigned int i=0;i<nEventsInFinalSample_;i++){
		nNPEventsInFinalSample_ = nNPEventsInFinalSample_ + EventsInfos_[i].isnp;
		if(i<nEventsInHWSubsample_) nNPEventsInHWSubsample_ = nNPEventsInHWSubsample_ + EventsInfos_[i].isnp;
	}
	FracNPEventsInFinalSample_ = float(nNPEventsInFinalSample_)/float(nEventsInFinalSample_);
	FracNPEventsInHWSubsample_ = float(nNPEventsInHWSubsample_)/float(nEventsInHWSubsample_);
	FracNPEventsEndUpInHWSubsample_ = float(nNPEventsInHWSubsample_)/float(nNPEventsInFinalSample_);
	SignalSignificanceInHWSubSample_ = float(nNPEventsInHWSubsample_)/float(sqrt(nEventsInHWSubsample_-nNPEventsInHWSubsample_)); // signal/sqrt(background)
}

void SampleCombinedWeightCalculator::CalculateV(vector<float>& weights_per_event, TString combination, float fraction) {
      ////////////////////////////////////////
      // Calculate the combined (w.r.t. the events) weight V for a sample
      ////////////////////////////////////////	
	Weights_per_event_ = weights_per_event;
	Combination_ = combination;
	Fraction_ = fraction;
	nEventsInFinalSample_ = Weights_per_event_.size();
	nEventsInHWSubsample_ = int(Fraction_*Weights_per_event_.size());
	if(nEventsInHWSubsample_ == 0){
		cout << "Warning: number of events in highest-weight subsample = 0!" << endl;
		nEventsInHWSubsample_++;
	}
	
	if(Combination_ == "sum") V_ = 0.0;
	else if(Combination_ == "product") V_ = 1.0;
	else if(Combination_ == "average") V_ = 0.0;
	else if(Combination_ == "logproduct") V_ = 0.0;
	else cout << "Undefined combination of weights (options: \"sum\", \"product\", \"average\", \"logproduct\")" << endl;
	
	RankEvents(Weights_per_event_);  //ik heb het veranderd; vroeger niet de berekeningen gedaan op de private member...
	nEventsInFinalSample_ = Weights_per_event_.size();
	nEventsInHWSubsample_ = int(Fraction_*Weights_per_event_.size());
	if(nEventsInHWSubsample_ == 0){
		cout << "Warning: number of events in highest-weight subsample = 0!" << endl;
		nEventsInHWSubsample_++;
	}
	for(unsigned int i=0; i<nEventsInHWSubsample_; i++){
		if(Combination_ == "sum") V_ = V_ + Weights_per_event_[i];
		else if(Combination_ == "product") V_ = V_*Weights_per_event_[i];
		else if(Combination_ == "average" && nEventsInHWSubsample_ != 0) V_ = V_ + (Weights_per_event_[i] / nEventsInHWSubsample_);
		else if(Combination_ == "logproduct") V_ = V_ + TMath::Log(EventsInfos_[i].weight);
	}	 	
}

void SampleCombinedWeightCalculator::CalculateV(vector<eventinfo>& eventsinfos, TString combination, float fraction) {
      ////////////////////////////////////////
      // Calculate the combined (w.r.t. the events) weight V for a sample; take event info along
      ////////////////////////////////////////	
	EventsInfos_ = eventsinfos;
	Combination_ = combination;
	Fraction_ = fraction;
	
	if(Combination_ == "sum") V_ = 0.0;
	else if(Combination_ == "product") V_ = 1.0;
	else if(Combination_ == "average") V_ = 0.0;
	else if(Combination_ == "logproduct") V_ = 0.0;
	else cout << "Undefined combination of weights (options: \"sum\", \"product\", \"average\", \"logproduct\")" << endl;
	
	RankEvents(EventsInfos_);
	CalculatenEvents();
	//cout<<"EventsInfos_.size() = "<<EventsInfos_.size()<<endl;
	//cout<<"-----> nEventsInHWSubsample_ = "<<nEventsInHWSubsample_<<endl;
	for(unsigned int i=0; i<nEventsInHWSubsample_; i++){
		if(Combination_ == "sum"){
		  V_ = V_ + EventsInfos_[i].weight;
		  //cout<<" added "<<EventsInfos_[i].weight<<" to V_, updated V_ = "<<V_<<endl;
		  //cout<<"    some info about the event: "<<endl;
		  //for(unsigned int j=0; j<(*(EventsInfos_[i].variableValues)).size(); j++){
		  //   cout<<"      var "<<j<<": "<<(*(EventsInfos_[i].variableValues))[j]<<endl;
		  //}
		}
		else if(Combination_ == "product") V_ = V_*EventsInfos_[i].weight;
		else if(Combination_ == "average" && nEventsInHWSubsample_ != 0) V_ = V_ + (EventsInfos_[i].weight / nEventsInHWSubsample_);
		else if(Combination_ == "logproduct") V_ = V_ + TMath::Log(EventsInfos_[i].weight);
	}		 	
}
