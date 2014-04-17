#include "../interface/EventCombinedWeightCalculator.h"
  
EventCombinedWeightCalculator::EventCombinedWeightCalculator(string SignificanceCombination){
	SignificanceCombination_ = SignificanceCombination;
	doLumiSyst_ = false;
        doXsectionSyst_ = false;
        doJESSyst_ = false;
}

EventCombinedWeightCalculator::~EventCombinedWeightCalculator(){}

void EventCombinedWeightCalculator::Init(){
      // initialize some private members	
   	if(SignificanceCombination_=="PRODUCT"){
	   Significance_ = 1.;
	   CombinedWeight_ = 1.;
	}	   
	if(SignificanceCombination_=="SUM"||SignificanceCombination_=="SUM_LRatioMinded"){
	   Significance_ = 0.;
	   CombinedWeight_ = 0.;
	}
   	isUnderflowWarning_ = false;
}

void EventCombinedWeightCalculator::SetSystematics(float LumiError, vector<MCObsExpectation*> MCObsExp_JESSyst, bool doLumiSyst, bool doXsectionSyst, bool doJESSyst){
      // setting members regarding systematics
	LumiError_ = LumiError;
	MCObsExp_JESSyst_ = MCObsExp_JESSyst;
	doLumiSyst_ = doLumiSyst;
        doXsectionSyst_ = doXsectionSyst;
        doJESSyst_ = doJESSyst;
}

void EventCombinedWeightCalculator::SetObsExp(vector<MCObsExpectation*> MCSMObsExp,vector<MCObsExpectation*> MCNPObsExp){
      // setting the observable histograms obtained by running on the Standard Model MC samples
        SetObsExp(MCSMObsExp);
	if(SignificanceCombination_=="SUM_LRatioMinded"){
	   MCNPObsExp_ = MCNPObsExp;
	}
}

void EventCombinedWeightCalculator::SetObsExp(vector<MCObsExpectation*> MCSMObsExp){
      // setting the observable histograms obtained by running on the Standard Model and New Physics MC samples
        MCSMObsExp_ = MCSMObsExp;
}

void EventCombinedWeightCalculator::CalculateWeight(vector<pair <string,float> > & event_string, vector<TH1F*>& histos_Obs_data, vector < Dataset* > datasets) {
      ////////////////////////////////////////
      // Calculate the combined (w.r.t. the observables) S-weight for an event
      ////////////////////////////////////////	
	int verbosity = 1;
	Init();
	vector<float> significances; //squared bin significances
        //loop over variables
        //check if the event is in an underflow bin of one of the observable histograms

//		if (event_string.size()!=event.size()) cout <<" Possible Problem with the collection of the events......CHECK statistical macro.....!!!!!!!!!!!!!"<<event.size()<<"  "<<event_string.size()<<endl;
		//cout<<" Will try to find event "<<" with  "<<event_string.size()<<" for Obs_Data "<<histos_Obs_data[1]->GetName()<<" and estimation "<<histos_Obs_estim[1]->GetName()<<"  "<<endl;	
	for(unsigned int i = 0; i < event_string.size(); i++) {
	   int bin=0;
	   int bin_est_SM=0, bin_est_NP=0, bin_JESSyst_est=0;
	   float binContent_data=0.;
	   float binContent_estim_SM=0., binContent_estim_NP=0., binContent_JESSyst_estim_SM=0.;
	   float binContent_estim_SMContributions[datasets.size()]; //different contributions to a bin, of the datasamples (SM 'background' processes)
	   float binVar=0.;
	   float binVar_Stat=0., binVar_Syst=0.;
	   TH1F* estim_SM=0;
	   TH1F* estim_NP=0;
	   TH1F* estim_JESSyst_SM=0;
	   vector<TH1F*> estim_SMContributions; //vector of histograms of the different contributions of the datasamples (SM 'background' processes)
	   //copy from MCObsExpectation to the vector
	   for (unsigned int l=0;l<MCSMObsExp_.size();l++){
		if (MCSMObsExp_[l]->GetHistoSMProcesses()->GetName()=="SMProcess_"+event_string[i].first){
	           if(verbosity>1) 
		    cout<<endl<<"  found MCSMObsExp plots....... "<<MCSMObsExp_[l]->GetHistoSMProcesses()->GetName()<<"   "<<l<<endl;           
		   estim_SM = MCSMObsExp_[l]->GetHistoSMProcesses();
		   bin_est_SM = estim_SM->FindBin(event_string[i].second);
		   int overflowbin_estim_SM = estim_SM->GetNbinsX() + 1;
		   if(bin_est_SM == overflowbin_estim_SM) bin_est_SM = overflowbin_estim_SM - 1;//because I added the overflow bin to the last bin in the histo, but the value of the variable for the event can still get in the overflowbin!! Note: I still need to test this more carefully!?	   
		   binContent_estim_SM = estim_SM->GetBinContent(bin_est_SM);	
		   for(unsigned int d=0;d<datasets.size();d++){
		        string dsname;
			dsname = datasets[d]->Name();
                        if(!(dsname.find("SUSY")<=dsname.size() || dsname.find("LM")<=dsname.size() || dsname.find("Zp")<=dsname.size() || dsname.find("NP")<=dsname.size() || dsname.find("prime")<=dsname.size())){
			  TH1F* temphist_SMContribution=0;
			  temphist_SMContribution = MCSMObsExp_[l]->GetHistogram(dsname);
		   	  estim_SMContributions.push_back(temphist_SMContribution);
			  binContent_estim_SMContributions[d] = estim_SMContributions[d]->GetBinContent(bin_est_SM);
		        }
		   }		   		   
	           break;
	        }
           }
	   
	   bin = histos_Obs_data[i]->FindBin(event_string[i].second);
	   int overflowbin_Obs_data = histos_Obs_data[i]->GetNbinsX() + 1;
	   if(bin == overflowbin_Obs_data) bin = overflowbin_Obs_data - 1;//because I added the overflow bin to the last bin in the histo, but the value of the variable for the event can still get in the overflowbin!! Note: I still need to test this more carefully!?	   
	   
	   if(SignificanceCombination_!="SUM_LRatioMinded"){
	      binContent_data = histos_Obs_data[i]->GetBinContent(bin);
	   }
	   else{
	       //in the case SignificanceCombination_=="SUM_LRatioMinded", you don't care at all about 'bin' and 'binContent_data', you work with the bin contents of the SM and SM+NP expectations
	       for (unsigned int l=0;l<MCNPObsExp_.size();l++){
		if (MCNPObsExp_[l]->GetHistoAll()->GetName()=="All_"+event_string[i].first){
	           if(verbosity>1) 
		    cout<<endl<<"  found MCNPObsExp plots....... "<<MCNPObsExp_[l]->GetHistoAll()->GetName()<<"   "<<l<<endl;           
		   estim_NP = MCNPObsExp_[l]->GetHistoAll();
		   bin_est_NP = estim_NP->FindBin(event_string[i].second);
		   int overflowbin_estim_NP = estim_NP->GetNbinsX() + 1;
		   if(bin_est_NP == overflowbin_estim_NP) bin_est_NP = overflowbin_estim_NP - 1;	   
		   binContent_estim_NP = estim_NP->GetBinContent(bin_est_NP);		   		   
		   break;
	        }
              }	   
	   }	

     	   if(bin<1 || bin_est_SM<1){
		isUnderflowWarning_ = true;
	 	cout<<" !!!!!!!!!!!!!!!!!!!!!!      Found event in WRONG BIN --------------------->>>>>>>>>>>>>>>>>"<<i<<" with  "<<event_string[i].second<<" for Obs_Data "<<histos_Obs_data[i]->GetName()<<"in bin  "<<bin<<endl;
		break;
	   }
	   else {
		binVar=0.0;
		
		//cout<< " Sanity check....." <<histos_Obs_data[i]->GetName()<<"   should match  this  "<<event_string[i].first<<"   and value is ... "<<event_string[i].second<<"  bin of data "<<bin<<" bin of estim  "<<bin_est<<"  data bin contents "<<binContent_data<<" binContent est  "<<binContent_estim<<" data bins  "<<histos_Obs_data[i]->GetNbinsX()<<" estim bins  "<<histos_Obs_estim[i]->GetNbinsX()<<endl;
		binVar_Stat=0.0;
		if(SignificanceCombination_!="SUM_LRatioMinded") binVar_Stat = TMath::Power(histos_Obs_data[i]->GetBinError(bin),2) + TMath::Power(estim_SM->GetBinError(bin_est_SM),2);	//sum of variances
		else binVar_Stat = TMath::Power(estim_NP->GetBinError(bin_est_NP),2) + TMath::Power(estim_SM->GetBinError(bin_est_SM),2);
		//binVar = TMath::Power(histos_Obs_data[i]->GetBinError(bin),2); //neglect estim error temporary!!!!
		
		binVar_Syst=0.0; //the lumierror correctly included...??
	        if(doXsectionSyst_){
		  if(verbosity>2) cout<<"    --- Adding Xsection systematics into squared bin significance ---"<<endl; 
		  for(unsigned int d=0;d<datasets.size();d++){
		        if(verbosity>1) cout<<"            cross section (dataset "<<d<<") contribution to binVar_Syst: the square of "<<datasets[d]->XsectionError()<<" * "<<binContent_estim_SMContributions[d]<<endl;
			if(verbosity>1) cout<<"            normfactor error (dataset "<<d<<") contribution to binVar_Syst: the square of "<<LumiError_<<" * "<<binContent_estim_SMContributions[d]<<endl;
			binVar_Syst = binVar_Syst + TMath::Power((datasets[d]->XsectionError())*binContent_estim_SMContributions[d],2);			
		        if(verbosity>1) cout<<"              updated binVar_Syst: "<<binVar_Syst<<endl;
		  }
		}
		if(doLumiSyst_){
	           if(verbosity>2) cout<<"    --- Adding Luminosity systematics into squared bin significance ---"<<endl; 
		   binVar_Syst = binVar_Syst + TMath::Power(LumiError_*binContent_estim_SM,2); //lumierror contribution
		}
		if(doJESSyst_){
	           if(verbosity>2) cout<<"    --- Adding JES systematics into squared bin significance ---"<<endl; 
	   	   for (unsigned int l=0;l<MCObsExp_JESSyst_.size();l++){
			  if (MCObsExp_JESSyst_[l]->GetHistoSMProcesses()->GetName()=="SMProcess_"+event_string[i].first){
	           	    if(verbosity>1) 
		              cout<<endl<<"  found MCObsExp_JESSsyst_ plots....... "<<MCObsExp_JESSyst_[l]->GetHistoSMProcesses()->GetName()<<"   "<<l<<endl;           
		            estim_JESSyst_SM = MCObsExp_JESSyst_[l]->GetHistoSMProcesses();
		            bin_JESSyst_est = estim_JESSyst_SM->FindBin(event_string[i].second);
		            binContent_JESSyst_estim_SM = estim_JESSyst_SM->GetBinContent(bin_est_SM);		   		   
	           	    break;
	        	  }
           	   }
		   if(verbosity>1) cout<<"  binVar_Syst of JES = "<<TMath::Power(binContent_estim_SM-binContent_JESSyst_estim_SM,2)<<endl;
		   binVar_Syst = binVar_Syst + TMath::Power(binContent_estim_SM-binContent_JESSyst_estim_SM,2); //ok?		
		}
		binVar = binVar_Stat + binVar_Syst;
		
		if(verbosity>1){
		   cout<<" Found event "<<i<<" with  "<<event_string[i].second<<" for Obs_Data "<<histos_Obs_data[i]->GetName()<<"in bin  "<<bin<<" event size "<<event_string.size()<<"  binContent  data "<<binContent_data<<" binContent est  "<<binContent_estim_SM<<endl; 
		   cout << "   value of observable " << i << ": " << event_string[i].second << endl;
		   cout << "   bin data " << bin << " bin estim (SM)  "<<bin_est_SM<<endl;
		   cout << "   binVar: " << binVar << endl;
		   cout << "      binVar_Stat = "<<binVar_Stat<<endl;
		   cout << "      binVar_Syst = "<<binVar_Syst<<endl;
		   cout << "   binContent_data: " << binContent_data << endl;
		   cout << "   binContent_estim_SM: " << binContent_estim_SM << endl;
		   cout << "   binContent_estim_NP: " << binContent_estim_NP << endl;	  
		   cout << "   Error histos_Obs_data["<<i<<"]: " << histos_Obs_data[i]->GetBinError(bin) << endl;
		   cout << "   Error estim_SM["<<i<<"]: " << estim_SM->GetBinError(bin_est_SM) << endl;
		   cout << "   Error estim_NP["<<i<<"]: " << estim_NP->GetBinError(bin_est_NP) << endl;
		}
		
		if(binVar!=0) {
		   if(SignificanceCombination_!="SUM_LRatioMinded") Significance_ = TMath::Power((binContent_data - binContent_estim_SM),2) / binVar;  // deviation between data and estimation
	      	   else Significance_ = TMath::Power((binContent_estim_NP - binContent_estim_SM),2) / binVar;
		  //  cout << "Significance_: " << Significance_ << endl;		  
		   significances.push_back(Significance_);
		}
		else {
		   cout << "[Significance_] Warning: trying to divide by 0!!" << endl;
		   Significance_ = 0;  // ignore events with related binVar = 0 (by doing this you ignore the events because in the end you take the highest weights)
		   // in normal operation, this should not be possible, since an event cannot appear in an empty bin.
		}
	    }
	} //end loop over event_string vector
	
	//loop over variables
	for(unsigned int i = 0; i < significances.size(); i++) {
	    if(SignificanceCombination_=="PRODUCT") CombinedWeight_ = CombinedWeight_*significances[i];
	    if(SignificanceCombination_=="SUM"||SignificanceCombination_=="SUM_LRatioMinded"){ 
	      CombinedWeight_ = CombinedWeight_ + significances[i]; 
	      if(verbosity>3) cout<<"   "<<i<<" added "<<significances[i]<<", updated CombinedWeight_ to "<<CombinedWeight_<<endl;
	    }
	}
	     
	if(verbosity>1){
	     cout << "----> CombinedWeight_: " << CombinedWeight_ << endl;
	}
	if(significances.size()>0) significances.clear();
}

/*void EventCombinedWeightCalculator::CalculateWeight(vector<pair <string,float> > & event_string, vector<TH1F*>& histos_Obs_data, vector<TH1F*>& histos_Obs_estim) {
      ////////////////////////////////////////
      // Calculate the combined (w.r.t. the observables) weight for an event
      ////////////////////////////////////////	
	int verbosity = 1;
	Init();
	vector<float> significances; //(squared) bin significances
        //loop over variables
        //check if the event is in an underflow bin of one of the observable histograms

//		if (event_string.size()!=event.size()) cout <<" Possible Problem with the collection of the events......CHECK statistical macro.....!!!!!!!!!!!!!"<<event.size()<<"  "<<event_string.size()<<endl;
		//cout<<" Will try to find event "<<" with  "<<event_string.size()<<" for Obs_Data "<<histos_Obs_data[1]->GetName()<<" and estimation "<<histos_Obs_estim[1]->GetName()<<"  "<<endl;	
	for(unsigned int i = 0; i < event_string.size(); i++) {
	   int bin=0;
	   int bin_est=0;
	   float binContent_data=0.;
	   float binContent_estim=0.;
	   float binVar=0.;
	   //for (unsigned int ll=0;ll<histos_Obs_data.size();ll++){		
	   bin = histos_Obs_data[i]->FindBin(event_string[i].second);
	   bin_est = histos_Obs_estim[i]->FindBin(event_string[i].second);

     	   if(bin<1 || bin_est<1){
		isUnderflowWarning_ = true;
	 	cout<<" !!!!!!!!!!!!!!!!!!!!!!      Found event in WRONG BIN --------------------->>>>>>>>>>>>>>>>>"<<i<<" with  "<<event_string[i].second<<" for Obs_Data "<<histos_Obs_data[i]->GetName()<<"in bin  "<<bin<<endl;
		break;
	   }
	   else {
		binContent_data = histos_Obs_data[i]->GetBinContent(bin);
		binContent_estim = histos_Obs_estim[i]->GetBinContent(bin_est); 
		//cout<< " Sanity check....." <<histos_Obs_data[i]->GetName()<<"   should match  this  "<<event_string[i].first<<"   and value is ... "<<event_string[i].second<<"  bin of data "<<bin<<" bin of estim  "<<bin_est<<"  data bin contents "<<binContent_data<<" binContent est  "<<binContent_estim<<" data bins  "<<histos_Obs_data[i]->GetNbinsX()<<" estim bins  "<<histos_Obs_estim[i]->GetNbinsX()<<endl;
		binVar=0.0;
		binVar = TMath::Power(histos_Obs_data[i]->GetBinError(bin),2) + TMath::Power(histos_Obs_estim[i]->GetBinError(bin_est),2);	//sum of variances
		//binVar = TMath::Power(histos_Obs_data[i]->GetBinError(bin),2); //neglect estim error temporary!!!!
		//cout<<"NEGLECTING UNCERTAINTY ON ESTIMATION"<<endl;
		if(verbosity>1){
			cout<<" Found event "<<i<<" with  "<<event_string[i].second<<" for Obs_Data "<<histos_Obs_data[i]->GetName()<<"in bin  "<<bin<<" event size "<<event_string.size()<<"  binContent  data "<<binContent_data<<" binContent est  "<<binContent_estim<<endl; 
		  	cout << endl;
			cout << "   value of observable " << i << ": " << event_string[i].second << endl;
		  	cout << "   bin data " << bin << " bin estim  "<<bin_est<<endl;
		  	cout << "   binVar: " << binVar << endl;
		  	cout << "   binContent_data: " << binContent_data << endl;
		  	cout << "   binContent_estim: " << binContent_estim << endl;	  
		  	cout << "   Error histos_Obs_data["<<i<<"]: " << histos_Obs_data[i]->GetBinError(bin) << endl;
		  	cout << "   Error histos_Obs_estim["<<i<<"]: " << histos_Obs_estim[i]->GetBinError(bin_est) << endl;
		}
		
		if(binVar!=0) {
			Significance_ = TMath::Power((binContent_data - binContent_estim),2) / binVar;  // deviation between data and estimation
	      		//  cout << "Significance_: " << Significance_ << endl;		  
			significances.push_back(Significance_);
		}
		else {
			if(verbosity>1){
				cout << "[Significance_] Warning: trying to divide by 0" << endl;
				cout << "   value of observable " << i << ": " << event_string[i].second << endl;
		  		cout << "   bin: " << bin << endl;
				cout << "   binVar: " << binVar << endl;
		  		cout << "   binContent_data: " << binContent_data << endl;
		  		cout << "   binContent_estim: " << binContent_estim << endl;	  
		  		cout << "   Error histos_Obs_data["<<i<<"]: " << histos_Obs_data[i]->GetBinError(bin) << endl;
		  		cout << "   Error histos_Obs_estim["<<i<<"]: " << histos_Obs_estim[i]->GetBinError(bin) << endl;
			}
			Significance_ = 0;  // ignore events with related binVar = 0 (by doing this you ignore the events because in the end you take the highest weights)
			// in normal operation, this should not be possible, since an event cannot appear in an empty bin.
		}
		//	cout << "Significance_: " << Significance_ << endl;		  
		//break;
	    }
            //}/////////loop in all the histos to find the event

	} //end loop over event_string vector
	//loop over variables, and check which observables (which entries of the 'weights' vector) are to be considered as sensitive
	for(unsigned int i = 0; i < significances.size(); i++) {
	    if(SignificanceCombination_=="PRODUCT") CombinedWeight_ = CombinedWeight_*significances[i];
	    if(SignificanceCombination_=="SUM"){ 
	      CombinedWeight_ = CombinedWeight_ + significances[i]; 
	      if(verbosity>3) cout<<"   "<<i<<" added "<<significances[i]<<", updated CombinedWeight_ to "<<CombinedWeight_<<endl;
	    }
	}
	     
	if(verbosity>3){
	     cout << "CombinedWeight_: " << CombinedWeight_ << endl;
	}
	if(significances.size()>0) significances.clear();
}
*/
