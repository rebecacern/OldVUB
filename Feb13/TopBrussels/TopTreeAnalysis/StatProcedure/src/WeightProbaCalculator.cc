#include "../interface/WeightProbaCalculator.h"
  
WeightProbaCalculator::WeightProbaCalculator(){
	RefMCName_ = "";
	DataName_ = "";
	InitRefVDistr(50,0,700);
	InitVDistr(50,0,700);
	InitpDistr(10,0,1);
	InitNofSigmaDistr(20,0,5);
}

WeightProbaCalculator::WeightProbaCalculator(string refmcname, string dataname){
      // names will appear in the titles and names of the created histograms
	RefMCName_ = refmcname;
	DataName_ = dataname;
	InitRefVDistr(50,0,700);
	InitVDistr(50,0,700);
	InitpDistr(10,0,1);
	InitNofSigmaDistr(20,0,5);
}

WeightProbaCalculator::~WeightProbaCalculator(){}

void WeightProbaCalculator::InitRefVDistr(int nbinsx, double xlow, double xup){
	RefVDistribution_ = new TH1F("RefVDistribution "+TString(RefMCName_),"Monte Carlo V distribution "+TString(RefMCName_),nbinsx,xlow,xup);
	RefVDistribution_->GetXaxis()->SetTitle("V");
}

void WeightProbaCalculator::InitVDistr(int nbinsx, double xlow, double xup){
	VDistribution_ = new TH1F("VDistribution "+TString(DataName_ ),"V distribution "+TString(DataName_ ),nbinsx,xlow,xup);
	VDistribution_->GetXaxis()->SetTitle("V");
}

void WeightProbaCalculator::InitpDistr(int nbinsx, double xlow, double xup){
	pDistribution_ = new TH1F("pDistribution "+TString(DataName_ ),"p distribution "+TString(DataName_ ),nbinsx,xlow,xup);
	pDistribution_->GetXaxis()->SetTitle("p");
}

void WeightProbaCalculator::InitNofSigmaDistr(int nbinsx, double xlow, double xup){
	NofSigmaDistribution_ = new TH1F("NofSigmaDistribution "+TString(DataName_ ),"number of sigma distribution "+TString(DataName_ ),nbinsx,xlow,xup);
	NofSigmaDistribution_->GetXaxis()->SetTitle("number of sigma");
}

void WeightProbaCalculator::InitCorrHisto_pV(){
      // private method 
      //can be used, but most of the time not interesting
	CorrelationHisto_pV_ = new TH2F("CorrelationHisto_pV "+TString(DataName_ ),"Correlation p and V "+TString(DataName_ ),VDistribution_->GetNbinsX(),VDistribution_->GetXaxis()->GetXmin(),VDistribution_->GetXaxis()->GetXmax(),pDistribution_->GetNbinsX(),pDistribution_->GetXaxis()->GetXmin(),pDistribution_->GetXaxis()->GetXmax());
	CorrelationHisto_pV_->GetXaxis()->SetTitle("V");
	CorrelationHisto_pV_->GetYaxis()->SetTitle("p");
}

void WeightProbaCalculator::CalculateProba(vector<float> Vvalues_MC, float Vvalue_data){
      ////////////////////////////////////////
      // Calculate the p value of the data corresponding to the given V value
      ////////////////////////////////////////
	int verbosity = 0;

	int integral_larger = 0;
	int integral_total = 0;
	for(unsigned int i=0; i<Vvalues_MC.size(); i++){
		if(!(Vvalues_MC[i] != Vvalues_MC[i]) && !(Vvalue_data != Vvalue_data)){  // to exclude nan
			if(Vvalues_MC[i]>Vvalue_data) 
				integral_larger++;
			integral_total++;
		}
		else cout << "Warning: nan V value" << endl;
	}
	if(verbosity>0) 
		cout << "Vvalue_data: " << Vvalue_data << ", integral_larger: " << integral_larger << ", integral_total: " << integral_total << endl;
	if(integral_total!=0){
		Proba_ = float(integral_larger)/float(integral_total);
		//Significance_ = RooStats::PValueToSignificance(Proba_);  // for the moment HERE not automatic, because of issue with p-value or 1 - p-value, which is also not automatic... I'll have to see/fix this; possibly function on Proba_ or 1-Proba_ is the same anyway...!
	}
	else 
		cout << "Warning: integral Monte Carlo V distribution = 0" << endl;
	
}

void WeightProbaCalculator::FillDistributions(vector<float> Vvalues_MC, vector<float> Vvalues_data, vector<float> pvalues){
	for(unsigned int i=0; i<Vvalues_MC.size(); i++){
		if(!(Vvalues_MC[i] != Vvalues_MC[i])){
			RefVDistribution_->Fill(Vvalues_MC[i]);
		}
	}
	//InitCorrHisto_pV(); //can be used, but most of the time not interesting
	for(unsigned int i=0; i<Vvalues_data.size(); i++){
		if(!(Vvalues_data[i] != Vvalues_data[i])){
			VDistribution_->Fill(Vvalues_data[i]);
			pDistribution_->Fill(pvalues[i]);
			//CorrelationHisto_pV_->Fill(Vvalues_data[i],pvalues[i]); //can be used, but most of the time not interesting
			NofSigmaDistribution_->Fill(ProbaToSignificance(pvalues[i]));
		}
	}
}

void WeightProbaCalculator::WriteDistributions(){
	cout << "Writing distributions..." << endl;
	RefVDistribution_->Write(RefVDistribution_->GetName(),TObject::kOverwrite);
	VDistribution_->Write(VDistribution_->GetName(),TObject::kOverwrite);
	pDistribution_->Write(pDistribution_->GetName(),TObject::kOverwrite);
	//CorrelationHisto_pV_->Write(CorrelationHisto_pV_->GetName(),TObject::kOverwrite); //can be used, but most of the time not interesting
	NofSigmaDistribution_->Write(NofSigmaDistribution_->GetName(),TObject::kOverwrite);
}

void WeightProbaCalculator::SortVectorDescending(vector<float> &pvalues){
      // sort such that first element is the largest
	sortDescending myComparison;
	sort(pvalues.begin(),pvalues.end(),myComparison);
}
     
void WeightProbaCalculator::SortVectorAscending(vector<float> &pvalues){
      // sort such that first element is the smallest
	sortAscending myComparison;
	sort(pvalues.begin(),pvalues.end(),myComparison);
}

float WeightProbaCalculator::GetProba95(vector<float> pvalues){
      // calulate p95; the value for which 95% of the pseudoexperiments have a smaller p value
	SortVectorDescending(pvalues);
	return pvalues[int(0.05*pvalues.size())+1];
}

float WeightProbaCalculator::GetP05Data(vector<float> pvalues){
      // calulate P05Data; the probability for 'data' to have a p value smaller than 0.05
      	unsigned int nofpvaluesSmaller = 0;
	for(unsigned int i=0;i<pvalues.size();i++){
		if(pvalues[i]<0.05) nofpvaluesSmaller++;
	}
	return float(nofpvaluesSmaller)/float(pvalues.size());
}

float WeightProbaCalculator::GetProbaPercentile(vector<float> pvalues, double percent){
     // calculate the percentile (related to the specified second argument) of e.g. p-values;
     // the 50th percentile is the median, and in principle, the 95th percentile is p95
	SortVectorAscending(pvalues);
	double pvaluesArray[pvalues.size()];  // needs to be an array of doubles (Double_t) to make use of TMath::Quantiles();
	for(unsigned int i=0;i<pvalues.size();i++){
		pvaluesArray[i] = pvalues[i];
	}
	double probability[1], quantile[1];
	probability[0] = percent/100;
	TMath::Quantiles((int) pvalues.size(),1,pvaluesArray,quantile,probability);  //see ROOT documentation for non-default calculation of quantiles
	return (float) quantile[0];
}

float WeightProbaCalculator::ProbaToSignificance(float proba){
	return RooStats::PValueToSignificance(proba);
}

float WeightProbaCalculator::SignificanceToProba(float significance){
	return RooStats::SignificanceToPValue(significance);
}

float WeightProbaCalculator::PowerOfTest(vector<float> Vvalues_data, float Valfa){
	// calulate the power of the test (1-beta); the probability for 'data' to have a V value smaller than Valfa. This is calculated in the same way as P05Data,
	// so this could be standardized in the future 
	unsigned int nofVvaluesSmaller = 0;
	for(unsigned int i=0;i<Vvalues_data.size();i++){
		if(Vvalues_data[i]<Valfa) nofVvaluesSmaller++;
	}
	return float(nofVvaluesSmaller)/float(Vvalues_data.size());
}


