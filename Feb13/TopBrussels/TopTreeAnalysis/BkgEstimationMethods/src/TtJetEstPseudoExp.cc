#include "../interface/TtJetEstPseudoExp.h"

TtJetEstPseudoExp::TtJetEstPseudoExp(int NbOfPseudoExp, TtJetEstimation *ttjEst){

	NbOfPseudoExp_    = NbOfPseudoExp;
	ttjEst_            = ttjEst;
	rand              = new TRandom3(0);
	hEstMinusExp = 0;
	hSystError = 0;

	if(ttjEst){
		ExpPlotCRBkg = ttjEst->PlotByDatasetsCRBkg();
		ExpPlotSRBkg = ttjEst->PlotByDatasetsSRBkg();
		ExpPlotCRTtJet = ttjEst->PlotByDatasetsCRTtJets();
		ExpPlotSRTtJet = ttjEst->PlotByDatasetsSRTtJets();
		ExpPlotCRNP = ttjEst->PlotByDatasetsCRNP();
		ExpPlotSRNP = ttjEst->PlotByDatasetsSRNP();
		for(unsigned int i=0;i<ExpPlotCRBkg.size();i++)
			PEPlotCRBkg.push_back((TH1F*) ExpPlotCRBkg[i]->Clone());
		for(unsigned int i=0;i<ExpPlotSRBkg.size();i++)
			PEPlotSRBkg.push_back((TH1F*) ExpPlotSRBkg[i]->Clone());
		for(unsigned int i=0;i<ExpPlotCRBkg.size();i++)
			PEPlotBkgEstim.push_back((TH1F*) ExpPlotCRBkg[i]->Clone());
		for(unsigned int i=0;i<ExpPlotCRTtJet.size();i++)
			PEPlotCRTtJet.push_back((TH1F*) ExpPlotCRTtJet[i]->Clone());
		for(unsigned int i=0;i<ExpPlotSRTtJet.size();i++)
			PEPlotSRTtJet.push_back((TH1F*) ExpPlotSRTtJet[i]->Clone());
		for(unsigned int i=0;i<ExpPlotCRNP.size();i++)
			PEPlotCRNP.push_back((TH1F*) ExpPlotCRNP[i]->Clone());
		for(unsigned int i=0;i<ExpPlotSRNP.size();i++)
			PEPlotSRNP.push_back((TH1F*) ExpPlotSRNP[i]->Clone());
		if(ttjEst->PlotEstimated()->GetXaxis()->GetXbins ()->fN>0){
			hEstMinusExp = new TProfile("hEstMinusExp","Estimation - Expectation", ttjEst->PlotEstimated()->GetXaxis()->GetNbins(),ttjEst->PlotEstimated()->GetXaxis()->GetXbins ()->fArray);
			hSystError = new TH1F("hSystError","Mean(Estimation - Expectation) bin per bin",ttjEst->PlotEstimated()->GetXaxis()->GetNbins(),ttjEst->PlotEstimated()->GetXaxis()->GetXbins ()->fArray);
		}
		else{
			hEstMinusExp = new TProfile("hEstMinusExp","Estimation - Expectation", ttjEst->PlotEstimated()->GetNbinsX(), ttjEst->PlotEstimated()->GetXaxis()->GetXmin(), ttjEst->PlotEstimated()->GetXaxis()->GetXmax());
			hSystError = new TH1F("hSystError","Mean(Estimation - Expectation) bin per bin", ttjEst->PlotEstimated()->GetNbinsX(), ttjEst->PlotEstimated()->GetXaxis()->GetXmin(), ttjEst->PlotEstimated()->GetXaxis()->GetXmax());
		}
	}
	hPEPlotEstim = new TH1F("hPEPlotEstim","N_{t#bar{t} expected}",100,0,10000);
	hPEPlotExpect = new TH1F("hPEPlotExpect","N_{t#bar{t} expected}",100,0,10000);
        hPEData = new TH1F("hPEData","Pseudo-data",100,0,10000);
        hNttEst = new TH1F("hNttEst","N_{t#bar{t} estimated}",50,0,5000);
	hNttEstMinusNttExp = new TH1F("hNttEstMinusNttExp","",50,-1000,1000);
	hPull = new TH1F("hPull","Pull",50,-10,10);
}

TtJetEstPseudoExp::~TtJetEstPseudoExp(){	
	delete hPEPlotEstim;
	delete hPEPlotExpect;
	delete hPEData;
	delete hNttEst;
	delete hNttEstMinusNttExp;
	delete hPull;
}

void TtJetEstPseudoExp::PrintResults(){
}

void TtJetEstPseudoExp::RollTheDice(bool verbose, bool addBkgEstimError, vector<float> BkgEstimError){
	//Compute a fake bkg estimation according to BkgEstimError vector which should contain the relative error in % expected for the different bkg
	//prevent case where vector size is different for bkg
	if(addBkgEstimError && BkgEstimError.size()!= PEPlotCRBkg.size()) addBkgEstimError = false;
	
	if(verbose) std::cout<<"// Let's start the RollTheDice method..."<<std::endl;
	for(int i=0; i<NbOfPseudoExp_; i++)
	{
		if(verbose && i%50==0) std::cout<<"RollTheDice method : running the "<<i+1<<"th pseudo-exp"<<std::endl;
			//first reset the histograms
			hPEData->Reset();
			hPEPlotEstim->Reset();
			hPEPlotExpect->Reset();
			//CR
			for(unsigned int d=0;d<PEPlotCRBkg.size();d++) PEPlotCRBkg[d]->Reset();
			for(unsigned int d=0;d<PEPlotCRTtJet.size();d++) PEPlotCRTtJet[d]->Reset();
			for(unsigned int d=0;d<PEPlotCRNP.size();d++) PEPlotCRNP[d]->Reset();
			//SR
			for(unsigned int d=0;d<PEPlotSRBkg.size();d++) PEPlotSRBkg[d]->Reset();
			for(unsigned int d=0;d<PEPlotSRTtJet.size();d++) PEPlotSRTtJet[d]->Reset();
			for(unsigned int d=0;d<PEPlotSRNP.size();d++) PEPlotSRNP[d]->Reset();
			//
			for(unsigned int d=0;d<PEPlotBkgEstim.size();d++) PEPlotBkgEstim[d]->Reset();
			///////////////////////

			///////////////////////////////////////// 
			// Fill histos with randomized numbers
			//CR
			for(unsigned int d=0;d<ExpPlotCRBkg.size();d++)//loop over datasets
				for(int b=1;b<ExpPlotCRBkg[d]->GetNbinsX()+1;b++){
					PEPlotCRBkg[d]->SetBinContent(b,rand->PoissonD(ExpPlotCRBkg[d]->GetBinContent(b)));
					hPEPlotEstim->SetBinContent(b,hPEPlotEstim->GetBinContent(i)+PEPlotCRBkg[d]->GetBinContent(b));
					if(addBkgEstimError) hPEData->SetBinContent(b,hPEData->GetBinContent(i)+PEPlotCRBkg[d]->GetBinContent(b));
				}
			for(unsigned int d=0;d<ExpPlotCRTtJet.size();d++)//loop over datasets
				for(int b=1;b<ExpPlotCRTtJet[d]->GetNbinsX()+1;b++){
					PEPlotCRTtJet[d]->SetBinContent(b,rand->PoissonD(ExpPlotCRTtJet[d]->GetBinContent(b)));
					hPEPlotEstim->SetBinContent(b,hPEPlotEstim->GetBinContent(i)+PEPlotCRTtJet[d]->GetBinContent(b));
					hPEData->SetBinContent(b,hPEData->GetBinContent(i)+PEPlotCRTtJet[d]->GetBinContent(b));
				}
			for(unsigned int d=0;d<ExpPlotCRNP.size();d++)//loop over datasets
				for(int b=1;b<ExpPlotCRNP[d]->GetNbinsX()+1;b++){
					PEPlotCRNP[d]->SetBinContent(b,rand->PoissonD(ExpPlotCRNP[d]->GetBinContent(b)));
					hPEPlotEstim->SetBinContent(b,hPEPlotEstim->GetBinContent(i)+PEPlotCRNP[d]->GetBinContent(b));
					if(addBkgEstimError) hPEData->SetBinContent(b,hPEData->GetBinContent(i)+PEPlotCRNP[d]->GetBinContent(b));
				}
			//SR minus CR
			for(unsigned int d=0;d<ExpPlotSRBkg.size();d++)//loop over datasets
				for(int b=1;b<ExpPlotSRBkg[d]->GetNbinsX()+1;b++){
					PEPlotSRBkg[d]->SetBinContent(b,rand->PoissonD(ExpPlotSRBkg[d]->GetBinContent(b)-ExpPlotCRBkg[d]->GetBinContent(b)));
					if(addBkgEstimError){ 
						hPEData->SetBinContent(b,hPEData->GetBinContent(i)+PEPlotSRBkg[d]->GetBinContent(b));
						//PEPlotBkgEstim[d]->SetBinContent(b,PEPlotBkgEstim[d]->GetBinContent(i)+rand->Gaus((PEPlotCRBkg[d]->GetBinContent(b)+PEPlotSRBkg[d]->GetBinContent(b)),(PEPlotCRBkg[d]->GetBinContent(b)+PEPlotSRBkg[d]->GetBinContent(b))*BkgEstimError[d]));
						PEPlotBkgEstim[d]->SetBinContent(b,PEPlotBkgEstim[d]->GetBinContent(i)+PEPlotCRBkg[d]->GetBinContent(b)+PEPlotSRBkg[d]->GetBinContent(b));
					}
				}
			for(unsigned int d=0;d<ExpPlotSRTtJet.size();d++)//loop over datasets
				for(int b=1;b<ExpPlotSRTtJet[d]->GetNbinsX()+1;b++){
					PEPlotSRTtJet[d]->SetBinContent(b,rand->PoissonD(ExpPlotSRTtJet[d]->GetBinContent(b)-ExpPlotCRTtJet[d]->GetBinContent(b)));
					hPEPlotExpect->SetBinContent(b,hPEPlotExpect->GetBinContent(i)+PEPlotSRTtJet[d]->GetBinContent(b));
					hPEData->SetBinContent(b,hPEData->GetBinContent(i)+PEPlotSRTtJet[d]->GetBinContent(b));
				}
			for(unsigned int d=0;d<ExpPlotSRNP.size();d++)//loop over datasets
				for(int b=1;b<ExpPlotSRNP[d]->GetNbinsX()+1;b++){
					PEPlotSRNP[d]->SetBinContent(b,rand->PoissonD(ExpPlotSRNP[d]->GetBinContent(b)-ExpPlotCRNP[d]->GetBinContent(b)));
					if(addBkgEstimError) hPEData->SetBinContent(b,hPEData->GetBinContent(i)+PEPlotSRNP[d]->GetBinContent(b));
				}
		float NttExp = -999.;
		float NttEst = -999.;
		float NttEstErr = -999.;
		pair<float,float> NRegion = ttjEst_->NRegion();
		//compute numbers
 		float n1 = ttjEst_->NofEvents (hPEPlotEstim, NRegion).first;
       		float n2 = ttjEst_->NofEvents (hPEData, NRegion).first;
	        float nBkgEstim = 0;
		if(addBkgEstimError){ 
			for(unsigned int d=0;d<PEPlotBkgEstim.size();d++)
				nBkgEstim+=ttjEst_->NofEvents (PEPlotBkgEstim[d], NRegion).first;
		}
		float factor = (n2 - nBkgEstim) / n1;
		float factorError = (sqrt(n2)*n1+n1*sqrt(n2))/(n1*n1);
		//if(TMath::IsNaN(factor)) factor = 1;
		NttEstErr = sqrt(hPEPlotEstim->Integral())*factor+hPEPlotEstim->Integral()*factorError;
		hPEPlotEstim->Scale (factor);
		//fill values
		NttExp =  hPEPlotExpect->Integral();
		NttEst =  hPEPlotEstim->Integral();
		//fill plots
		hNttEst->Fill(NttEst);
		hNttEstMinusNttExp->Fill(NttEst-NttExp);
		hPull->Fill((NttEst-NttExp)/NttEstErr);
		if(hEstMinusExp){
			for(int b=1; b<hPEPlotEstim->GetNbinsX()+1; b++){
				hEstMinusExp->Fill(hPEPlotExpect->GetBinCenter(b),hPEPlotEstim->GetBinContent(b)-hPEPlotExpect->GetBinContent(b));
			}
		}
		
		}
	if(hSystError){
			for(int b=1; b<hPEPlotEstim->GetNbinsX()+1; b++)
				hSystError->SetBinContent(b, hEstMinusExp->GetBinContent(b));
	}
}


void TtJetEstPseudoExp::Write(TFile* fout, string label){
	string dirname = "TtJetEstimation" + label;
	string dirname2 = "PseudoExp";
	if(!fout->cd (dirname.c_str ())) {
		fout->mkdir (dirname.c_str ());
	}
	fout->cd(dirname.c_str());
	if(!fout->cd (dirname2.c_str())) {
		fout->GetDirectory(dirname.c_str())->mkdir (dirname2.c_str());
	}
	fout->cd((dirname+string("/")+dirname2).c_str());
	hNttEst->Write();
	hNttEstMinusNttExp->Write();
	hPull->Write();
	if(hEstMinusExp) hEstMinusExp->Write();
	if(hSystError) hSystError->Write();
}
