#include "../interface/WJetShapeEstimation.h"

WJetShapeEstimation::WJetShapeEstimation(TAxis* axis, int NofIterations, string VarLabel){
	histo_inclusif_ = new TH1F (TString ("histo_inclusif_") + VarLabel, "histo_inclusif", axis->GetNbins(), axis->GetXbins()->fArray);
	histo_0bj_ = new TH1F (TString ("histo_0bj_") + VarLabel, "histo_0bj", axis->GetNbins(), axis->GetXbins()->fArray);
	histo_1bj_ = new TH1F (TString ("histo_1bj_") + VarLabel, "histo_1bj", axis->GetNbins(), axis->GetXbins()->fArray);
	histo_Wjets_ = new TH1F (TString ("histo_Wjets_") + VarLabel, "histo_Wjets", axis->GetNbins(), axis->GetXbins()->fArray);
	histo_ttjets_ = new TH1F (TString ("histo_ttjets_") + VarLabel, "histo_ttjets", axis->GetNbins(), axis->GetXbins()->fArray);
	histo_Wlike_ = new TH1F (TString ("histo_Wlike_") + VarLabel, "histo_Wlike", axis->GetNbins(), axis->GetXbins()->fArray);
	histo_TTlike_ = new TH1F (TString ("histo_TTlike_") + VarLabel, "histo_TTlike", axis->GetNbins(), axis->GetXbins()->fArray);
	h_SysError_ = new TH1F (TString ("SystErrors") + VarLabel, "SystErrors", axis->GetNbins(), axis->GetXbins()->fArray);
	histo_Wlike_cvg_  = new TH1F(TString("histo_Wlike_cvg")+VarLabel,"histo_Wlike_cvg",NofIterations,1,NofIterations+1); 
	histo_TTlike_cvg_ = new TH1F(TString("histo_TTlike_cvg")+VarLabel,"histo_TTlike_cvg",NofIterations,1,NofIterations+1); 
	histo_Wlike_EstEff_ = 0;
	histo_TTlike_EstEff_ = 0;
	//Sumw2
	histo_inclusif_->Sumw2();
	histo_0bj_->Sumw2();
	histo_1bj_->Sumw2();
	histo_Wjets_->Sumw2();
	histo_ttjets_->Sumw2();
	histo_Wlike_->Sumw2();
	histo_TTlike_->Sumw2();

	NofIterations_ = NofIterations; 
	RatioNbEvts_Wlike_1bjet_ = 0.5;
	RatioNbEvts_TTlike_0bjet_ = 0.5;
	Nevts0bj_ = 0;
	Nevts1bj_ = 0;

	c_Wjets_ = new TCanvas("cWjets");
	c_ttjets_ = new TCanvas("cttjets");
	c_Wjets_Split = new TCanvas("cWjetsRatio");
	c_ttjets_Split = new TCanvas("cttjetsRatio");
	c_Wjets_Nbjet        = new TCanvas("cWjetsNbjet");
	c_ttjets_Nbjet       = new TCanvas("cttjetsNbjet");
	c_Wjets_Split_Nbjet  = new TCanvas("cWjetsRatio_Nbjet");
	c_ttjets_Split_Nbjet = new TCanvas("cttjetsRatio_Nbjet");
	l_Wjets_ = new TLegend(0.6,0.6,0.9,0.9);
	l_ttjets_ = new TLegend(0.6,0.6,0.9,0.9);
	l_Wjets_Nbjet  = new TLegend(0.6,0.6,0.9,0.9);
	l_ttjets_Nbjet = new TLegend(0.6,0.6,0.9,0.9);
}

WJetShapeEstimation::WJetShapeEstimation(int nbins, float xmin, float xmax, int NofIterations, string VarLabel){
	histo_inclusif_ = new TH1F (TString ("histo_inclusif_") + VarLabel, "histo_inclusif", nbins, xmin, xmax);
	histo_0bj_ = new TH1F (TString ("histo_0bj_") + VarLabel, "histo_0bj", nbins, xmin, xmax);
	histo_1bj_ = new TH1F (TString ("histo_1bj_") + VarLabel, "histo_1bj", nbins, xmin, xmax);
	histo_Wjets_ = new TH1F (TString ("histo_Wjets_") + VarLabel, "histo_Wjets", nbins, xmin, xmax);
	histo_ttjets_ = new TH1F (TString ("histo_ttjets_") + VarLabel, "histo_ttjets", nbins, xmin, xmax);
	histo_Wlike_ = new TH1F (TString ("histo_Wlike_") + VarLabel, "histo_Wlike", nbins, xmin, xmax);
	histo_TTlike_ = new TH1F (TString ("histo_TTlike_") + VarLabel, "histo_TTlike", nbins, xmin, xmax);
	h_SysError_ = new TH1F (TString ("SystErrors") + VarLabel, "SystErrors", nbins, xmin, xmax);
	histo_Wlike_cvg_  = new TH1F(TString("histo_Wlike_cvg")+VarLabel,"histo_Wlike_cvg",NofIterations,1,NofIterations+1); 
	histo_TTlike_cvg_ = new TH1F(TString("histo_TTlike_cvg")+VarLabel,"histo_TTlike_cvg",NofIterations,1,NofIterations+1); 
	histo_Wlike_EstEff_ = 0;
	histo_TTlike_EstEff_ = 0;
	//Sumw2
	histo_inclusif_->Sumw2();
	histo_0bj_->Sumw2();
	histo_1bj_->Sumw2();
	histo_Wjets_->Sumw2();
	histo_ttjets_->Sumw2();
	histo_Wlike_->Sumw2();
	histo_TTlike_->Sumw2();
	
	NofIterations_ = NofIterations; 
	RatioNbEvts_Wlike_1bjet_ = 0.5;
	RatioNbEvts_TTlike_0bjet_ = 0.5;
	Nevts0bj_ = 0;
	Nevts1bj_ = 0;

	c_Wjets_ = new TCanvas("cWjets");
	c_ttjets_ = new TCanvas("cttjets");
	c_Wjets_Split = new TCanvas("cWjetsRatio");
	c_ttjets_Split = new TCanvas("cttjetsRatio");
	c_Wjets_Nbjet        = new TCanvas("cWjetsNbjet");
	c_ttjets_Nbjet       = new TCanvas("cttjetsNbjet");
	c_Wjets_Split_Nbjet  = new TCanvas("cWjetsRatio_Nbjet");
	c_ttjets_Split_Nbjet = new TCanvas("cttjetsRatio_Nbjet");
	l_Wjets_ = new TLegend(0.6,0.6,0.9,0.9);
	l_ttjets_ = new TLegend(0.6,0.6,0.9,0.9);
	l_Wjets_Nbjet  = new TLegend(0.6,0.6,0.9,0.9);
	l_ttjets_Nbjet = new TLegend(0.6,0.6,0.9,0.9);
}

WJetShapeEstimation::~WJetShapeEstimation(){
	if(histo_Wlike_) delete histo_Wlike_;
	if(histo_TTlike_) delete histo_TTlike_;
	if(histo_Wlike_EstEff_) delete histo_Wlike_EstEff_;
	if(histo_TTlike_EstEff_) delete histo_TTlike_EstEff_;
	if(histo_Wjets_) delete histo_Wjets_;
	if(histo_ttjets_) delete histo_ttjets_;
	if(histo_Wlike_cvg_) delete histo_Wlike_cvg_;
	if(histo_TTlike_cvg_) delete histo_TTlike_cvg_;
	if(histo_0bj_) delete histo_0bj_;
	if(histo_1bj_) delete histo_1bj_;
	if(histo_inclusif_) delete histo_inclusif_;
	if(h_SysError_) delete h_SysError_;
	if(c_Wjets_) delete c_Wjets_;
	if(c_ttjets_) delete c_ttjets_;
	if(c_Wjets_Split) delete c_Wjets_Split;
	if(c_ttjets_Split) delete c_ttjets_Split;
	if(c_Wjets_Nbjet) delete c_Wjets_Nbjet;
	if(c_ttjets_Nbjet) delete c_ttjets_Nbjet;
	if(c_Wjets_Split_Nbjet) delete c_Wjets_Split_Nbjet;
	if(c_ttjets_Split_Nbjet) delete c_ttjets_Split_Nbjet;
	if(l_Wjets_) delete l_Wjets_;
	if(l_ttjets_) delete l_ttjets_;
	if(l_Wjets_Nbjet) delete l_Wjets_Nbjet;
	if(l_ttjets_Nbjet) delete l_ttjets_Nbjet;
}

void WJetShapeEstimation::Fill(int nofbjets, float variable, float weight){
	histo_inclusif_->Fill(variable, weight);
	if(nofbjets == 0){
		histo_0bj_->Fill(variable, weight);
		histo_Wlike_->Fill(variable, weight);
	}
	if(nofbjets == 1){
		histo_1bj_->Fill(variable, weight);
		histo_TTlike_->Fill(variable, weight);
	}
}

void WJetShapeEstimation::FillMC(bool isWJets, bool isTtJets, float variable, float weight){
	if(isWJets) histo_Wjets_->Fill(variable, weight);
	if(isTtJets) histo_ttjets_->Fill(variable, weight);
	
}

void WJetShapeEstimation::AddMCHitos(TH1F* hwjets, TH1F* httjets){
	histo_Wjets_ = (TH1F*) hwjets->Clone(histo_Wjets_->GetName());
	histo_ttjets_ = (TH1F*) httjets->Clone(histo_ttjets_->GetName());
}

void WJetShapeEstimation::Write(TFile* fout, string label){
	if(fout==0) return;
	fout->cd();
	string dirname = "VJetShape"+label;
	fout->mkdir(dirname.c_str());
	fout->cd(dirname.c_str());
	if(h_SysError_ !=0 ) h_SysError_->Write();
	if(histo_Wlike_ != 0) histo_Wlike_->Write();
	if(histo_TTlike_ != 0) histo_TTlike_->Write();
	if(histo_Wlike_EstEff_ != 0) histo_Wlike_EstEff_->Write();
	if(histo_TTlike_EstEff_ != 0) histo_TTlike_EstEff_->Write();
	if(histo_0bj_ != 0) histo_0bj_->Write();
	if(histo_1bj_ != 0) histo_1bj_->Write();
	if(histo_inclusif_ != 0) histo_inclusif_->Write();
	if(histo_Wlike_cvg_ != 0) histo_Wlike_cvg_->Write();
	if(histo_TTlike_cvg_ != 0) histo_TTlike_cvg_->Write();
	if(c_Wjets_ !=0) c_Wjets_->Write();
	if(c_ttjets_ !=0) c_ttjets_->Write();
	if(c_Wjets_Split !=0) c_Wjets_Split->Write();
	if(c_ttjets_Split !=0) c_ttjets_Split->Write();
	if(c_Wjets_Nbjet !=0) c_Wjets_Nbjet->Write();
	if(c_ttjets_Nbjet !=0) c_ttjets_Nbjet->Write();
	if(c_Wjets_Split_Nbjet !=0) c_Wjets_Split_Nbjet->Write();
	if(c_ttjets_Split_Nbjet !=0) c_ttjets_Split_Nbjet->Write();
}


void WJetShapeEstimation::SubstractHisto(TH1F*& histo, TH1F* histoToSubstract)
{
	if(histo->GetNbinsX() != histoToSubstract->GetNbinsX()) {cout<<"Nb of bins do not match!"<<endl;return;}
	histo->Add((TH1F*) histoToSubstract, -1);
}


void WJetShapeEstimation::IterativeHistoSubstraction(float RatioNbEvts_Wlike_1bjet, float RatioNbEvts_TTlike_0bjet, float Nevts0bj, float Nevts1bj)
{

        TH1F* temp_0bjet;
	TH1F* temp_1bjet;

	/******* Start looping ***********/
	for(int i = 0; i<NofIterations_; i++)
	{
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// 1st step : scale "W-like 0bjet" histogram to the estimated number of "W-like" events in the 1 b-jet bin
		histo_Wlike_->Scale( (RatioNbEvts_Wlike_1bjet*Nevts1bj)/histo_Wlike_->Integral());
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// 2nd step : substract scaled "W-like 0bjet" histogram from the "1 bjet" one (meant to become the histo_TTlike_1bjet)
		temp_1bjet = (TH1F*)histo_1bj_->Clone();
		
		SubstractHisto(temp_1bjet, histo_Wlike_);
		histo_TTlike_cvg_->SetBinContent(i+1, DiffEstimation(temp_1bjet,histo_TTlike_));
		histo_TTlike_ = (TH1F*)temp_1bjet->Clone();
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// 3rd step :  scale "TT-like 1 bjet" histogram (i.e. "1 bjet" histogram with its "W-like" component already substracted)
		//             to the estimated number of "TT-like" events in the 0 b-jet bin
		histo_TTlike_->Scale((RatioNbEvts_TTlike_0bjet*Nevts0bj)/histo_TTlike_->Integral());
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// 4th step :  substract scaled "TT-like 1bjet" histogram from the "0 bjet" one (meant to become the histo_Wlike)
		temp_0bjet = (TH1F*)histo_0bj_->Clone();
		SubstractHisto(temp_0bjet, histo_TTlike_);
		histo_Wlike_cvg_->SetBinContent(i+1, DiffEstimation(temp_0bjet,histo_Wlike_));
		histo_Wlike_ = (TH1F*)temp_0bjet->Clone();
	}
	delete temp_0bjet;
	delete temp_1bjet;
}


void WJetShapeEstimation::Scale(float Ntt, float Nw){
	if(histo_Wlike_!=0 && histo_Wlike_->Integral()>0) histo_Wlike_->Scale(Nw/histo_Wlike_->Integral());
	if(histo_TTlike_!=0 && histo_TTlike_->Integral()>0) histo_TTlike_->Scale(Ntt/histo_TTlike_->Integral());
}


void WJetShapeEstimation::PrintResults(){
	cout<<endl;
	cout<<"***********************************************"<<endl;
	cout<<"*****      WJets Shape  Estimation Results  ***"<<endl;
	cout<<"***********************************************"<<endl;
	cout<<"Difference between histo_0bj & histo_1bj: "<<Chi2Normalized(histo_0bj_, histo_1bj_,true)<<endl;
	cout<<"Difference between histo_Wjets & histo_ttjets: "<<Chi2Normalized(histo_Wjets_, histo_ttjets_,true)<<endl;
	cout<<"Histogram for Wjets"<<endl;
	cout<<"Chi2/nof initial: histo_0bj vs histo_WjetsMC: "<<Chi2Normalized(histo_Wjets_, histo_0bj_,true)<<endl;
	cout<<"Chi2/nof initial: histo_WjetsEstimated vs histo_WjetsMC: "<<Chi2Normalized(histo_Wjets_,histo_Wlike_,true)<<endl;
	cout<<"Histogram for ttjets"<<endl;
	cout<<"Chi2/nof initial: histo_1bj vs histo_TTjetsMC: "<<Chi2Normalized(histo_ttjets_, histo_0bj_,true)<<endl;
	cout<<"Chi2/nof initial: histo_TTjetsEstimated vs histo_TTjetsMC: "<<Chi2Normalized(histo_ttjets_,histo_TTlike_,true)<<endl;
	cout<<"***********************************************"<<endl<<endl;;
}


void WJetShapeEstimation::AddNormalisationError(float factor){
        factor = factor/100.;
	for(int i=1;i<histo_Wlike_->GetNbinsX()+1;i++){
                float error = 0;
                error = sqrt((histo_Wlike_->GetBinError(i)*histo_Wlike_->GetBinError(i)) + (factor*factor*histo_Wlike_->GetBinContent(i)*histo_Wlike_->GetBinContent(i)));
                histo_Wlike_->SetBinError(i,error);
        }
}


void WJetShapeEstimation::ComputeSystematicError(bool add){
        for(int i=1;i<histo_Wlike_->GetNbinsX()+1;i++){
                float error = 0;
                float diff = 0;
                diff = histo_Wlike_->GetBinContent(i)-histo_Wjets_->GetBinContent(i);
                error = sqrt((histo_Wlike_->GetBinError(i)*histo_Wlike_->GetBinError(i)) + (diff*diff));
                if(add) histo_Wlike_->SetBinError(i,error);
        	h_SysError_->SetBinContent(i,diff);
	}
}

void WJetShapeEstimation::AddSystematicError(TH1F* SystError){
        for(int i=1;i<histo_Wlike_->GetNbinsX()+1;i++){
		float error = 0;
		error = sqrt((histo_Wlike_->GetBinError(i)*histo_Wlike_->GetBinError(i)) + (h_SysError_->GetBinContent(i)*h_SysError_->GetBinContent(i)));
		histo_Wlike_->SetBinError(i,error);
	}
}
								

void WJetShapeEstimation::Draw(string AxisLabel){
	int lineWidth = 2;
	int colorMC = kBlack;
	int colorIncl = kGreen;
	int colorbj = kBlue;
	int colorEst = kRed;
	int fstyle = 3001;

	//Settings
	histo_Wlike_->SetLineWidth(lineWidth);
	histo_Wlike_->SetLineColor(colorEst);
	histo_Wjets_->SetLineWidth(lineWidth);
	histo_Wjets_->SetLineColor(colorMC);
	histo_Wjets_->SetFillStyle(fstyle);
	histo_Wjets_->SetFillColor(colorMC);
	histo_TTlike_->SetLineWidth(lineWidth);
	histo_TTlike_->SetLineColor(colorEst);
	histo_ttjets_->SetLineWidth(lineWidth);
	histo_ttjets_->SetLineColor(colorMC);
	histo_ttjets_->SetFillStyle(fstyle);
	histo_ttjets_->SetFillColor(colorMC);
	histo_0bj_->SetLineWidth(lineWidth);
	histo_0bj_->SetLineColor(colorbj);
	histo_1bj_->SetLineWidth(lineWidth);
	histo_1bj_->SetLineColor(colorbj);
	histo_inclusif_->SetLineWidth(lineWidth);
	histo_inclusif_->SetLineColor(colorIncl);
	
	histo_Wlike_EstEff_  = (TH1F*) PlotCumulDiffEstimation(histo_Wjets_,histo_Wlike_,"hCumulDiffWlikeEst");
	if(histo_Wlike_EstEff_) histo_Wlike_EstEff_->GetXaxis()->SetTitle(AxisLabel.c_str());
	histo_TTlike_EstEff_ = (TH1F*) PlotCumulDiffEstimation(histo_ttjets_,histo_TTlike_,"hCumulDiffTTlikeEst");
	if(histo_TTlike_EstEff_) histo_TTlike_EstEff_->GetXaxis()->SetTitle(AxisLabel.c_str());
	
	//Clones
	TH1F* histo_TTlike_Clone = (TH1F*) histo_TTlike_->Clone();
	TH1F* histo_ttjets_Clone = (TH1F*) histo_ttjets_->Clone();
	TH1F* histo_Wlike_Clone  = (TH1F*) histo_Wlike_->Clone();
	TH1F* histo_Wjets_Clone  = (TH1F*) histo_Wjets_->Clone();
	TH1F* histo_0bj_Clone    = (TH1F*) histo_0bj_->Clone();
	TH1F* histo_1bj_Clone    = (TH1F*) histo_1bj_->Clone();
	TH1F* histo_inclusif_Clone = (TH1F*) histo_inclusif_->Clone();
	//Scale them

	if(histo_TTlike_Clone->Integral()>0)   histo_TTlike_Clone->Scale(1/histo_TTlike_Clone->Integral());
	if(histo_ttjets_Clone->Integral()>0)   histo_ttjets_Clone->Scale(1/histo_ttjets_Clone->Integral());
	if(histo_Wlike_Clone->Integral()>0)    histo_Wlike_Clone->Scale(1/histo_Wlike_Clone->Integral());
	if(histo_Wjets_Clone->Integral()>0)    histo_Wjets_Clone->Scale(1/histo_Wjets_Clone->Integral());
	if(histo_0bj_Clone->Integral()>0)      histo_0bj_Clone->Scale(1/histo_0bj_Clone->Integral());
	if(histo_1bj_Clone->Integral()>0)      histo_1bj_Clone->Scale(1/histo_1bj_Clone->Integral());
	if(histo_inclusif_Clone->Integral()>0) histo_inclusif_Clone->Scale(1/histo_inclusif_Clone->Integral());

	//////////////
	///   Wjets
	//////////////
	//TLegend
	l_Wjets_->AddEntry(histo_Wjets_,"MC: Wjets","l");
	l_Wjets_->AddEntry(histo_Wlike_,"Estimation: Wjets","l");
	//Plotting	
	c_Wjets_->cd();
	histo_Wjets_->GetXaxis()->SetTitle(AxisLabel.c_str());
	histo_Wjets_->Draw("HISTF");
	histo_Wlike_->Draw("esame");
	l_Wjets_->Draw("esame");
	//TLegend
	l_Wjets_Nbjet->AddEntry(histo_Wjets_Clone,"MC: Wjets","l");
	l_Wjets_Nbjet->AddEntry(histo_0bj_Clone,"Data - 0 bjets","l");
	l_Wjets_Nbjet->AddEntry(histo_Wlike_Clone,"Estimation: Wjets","l");
	//Plotting	
	c_Wjets_Nbjet->cd();
	histo_Wjets_Clone->Draw("HISTF");
	//histo_inclusif_Clone->Draw("esame");
	histo_0bj_Clone->Draw("esame");
	histo_Wlike_Clone->Draw("esame");
	l_Wjets_Nbjet->Draw("esame");
    	
	//Ratios
	TH1F* histo_0bj_Ratio       = PlotRatio(histo_0bj_,        histo_Wjets_);
	TH1F* histo_Wlike_Ratio     = PlotRatio(histo_Wlike_,      histo_Wjets_);
	TH1F* histo_0bj_NormRatio   = PlotRatio(histo_0bj_Clone,   histo_Wjets_Clone);
	TH1F* histo_Wlike_NormRatio = PlotRatio(histo_Wlike_Clone, histo_Wjets_Clone);
	//setting
	histo_0bj_Ratio->SetLineWidth(lineWidth);  histo_0bj_NormRatio->SetLineWidth(lineWidth);
	histo_0bj_Ratio->SetLineColor(colorbj);    histo_0bj_NormRatio->SetLineColor(colorbj);
	histo_Wlike_Ratio->SetLineWidth(lineWidth);histo_Wlike_NormRatio->SetLineWidth(lineWidth);
	histo_Wlike_Ratio->SetLineColor(colorEst); histo_Wlike_NormRatio->SetLineColor(colorEst);
	//TLine
	TLine* line = new TLine();
	line->SetY1(1);line->SetY2(1);
	line->SetLineColor(kBlack);
	line->SetLineStyle(2);
	line->SetLineWidth(1);
	//vector
	vector<TH1F*> hWjetsTop;
	vector<TH1F*> hWjetsBottom;
	//
	hWjetsTop.push_back(histo_Wjets_);
	hWjetsTop.push_back(histo_Wlike_);
	//
	hWjetsBottom.push_back(histo_Wlike_Ratio);
	//canvas
	c_Wjets_Split = TCanvasCreator(hWjetsTop,hWjetsBottom, line, l_Wjets_, AxisLabel.c_str());	
	c_Wjets_Split->SetName("c_WjetsRatio");
	//
	hWjetsTop.clear();
	hWjetsTop.push_back(histo_Wjets_Clone);
	hWjetsTop.push_back(histo_0bj_Clone);
	hWjetsTop.push_back(histo_Wlike_Clone);
	//
	hWjetsBottom.clear();
	//hWjetsBottom.push_back(histo_inclusif_Ratio);
	hWjetsBottom.push_back(histo_0bj_NormRatio);
	hWjetsBottom.push_back(histo_Wlike_NormRatio);
	//canvas
	c_Wjets_Split_Nbjet = TCanvasCreator(hWjetsTop,hWjetsBottom, line, l_Wjets_Nbjet, AxisLabel.c_str());	
	c_Wjets_Split_Nbjet->SetName("c_WjetsRatio_Nbjet");

	//////////////
	///   ttjets
	//////////////
        //TLegend
	l_ttjets_->AddEntry(histo_ttjets_,"MC: ttjets","l");
	l_ttjets_->AddEntry(histo_TTlike_,"Estimation: ttjets","l");
	//Plotting	
	c_ttjets_->cd();
	histo_ttjets_->GetXaxis()->SetTitle(AxisLabel.c_str());
	histo_ttjets_->Draw("HISTF");
	histo_TTlike_->Draw("esame");
	l_ttjets_->Draw("esame");
	//TLegend
	l_ttjets_Nbjet->AddEntry(histo_ttjets_Clone,"MC: ttjets","l");
	l_ttjets_Nbjet->AddEntry(histo_1bj_Clone,"Data - 1 bjet","l");
	l_ttjets_Nbjet->AddEntry(histo_TTlike_Clone,"Estimation: ttjets","l");
	//Plotting	
	c_ttjets_Nbjet->cd();
	histo_ttjets_Clone->Draw("HISTF");
	histo_1bj_Clone->Draw("esame");
	histo_TTlike_Clone->Draw("esame");
	l_ttjets_Nbjet->Draw("esame");
	
	//Ratios
        //TH1F* histo_inclusif_Ratio2 = PlotRatio(histo_inclusif_Clone, histo_ttjets_Clone);
        TH1F* histo_1bj_Ratio        = PlotRatio(histo_1bj_,         histo_ttjets_);
        TH1F* histo_TTlike_Ratio     = PlotRatio(histo_TTlike_,      histo_ttjets_);
        TH1F* histo_1bj_NormRatio    = PlotRatio(histo_1bj_Clone,    histo_ttjets_Clone);
        TH1F* histo_TTlike_NormRatio = PlotRatio(histo_TTlike_Clone, histo_ttjets_Clone);
	//setting
	histo_1bj_Ratio->SetLineWidth(lineWidth);   histo_1bj_NormRatio->SetLineWidth(lineWidth);
	histo_1bj_Ratio->SetLineColor(colorbj);     histo_1bj_NormRatio->SetLineColor(colorbj);
	histo_TTlike_Ratio->SetLineWidth(lineWidth);histo_TTlike_NormRatio->SetLineWidth(lineWidth);
	histo_TTlike_Ratio->SetLineColor(colorEst); histo_TTlike_NormRatio->SetLineColor(colorEst);
        //vector
	vector<TH1F*> httjetsTop;
	vector<TH1F*> httjetsBottom;
	//
	httjetsTop.push_back(histo_ttjets_);
	httjetsTop.push_back(histo_TTlike_);
	//
	httjetsBottom.push_back(histo_TTlike_Ratio);
	//canvas
	c_ttjets_Split = TCanvasCreator(httjetsTop,httjetsBottom, line, l_ttjets_, AxisLabel.c_str());	
	c_ttjets_Split->SetName("c_ttjetsRatio");
	//
	httjetsTop.clear();
	httjetsTop.push_back(histo_ttjets_Clone);
	httjetsTop.push_back(histo_1bj_Clone);
	httjetsTop.push_back(histo_TTlike_Clone);
	//
	httjetsBottom.clear();
	//httjetsBottom.push_back(histo_inclusif_Ratio);
	httjetsBottom.push_back(histo_1bj_NormRatio);
	httjetsBottom.push_back(histo_TTlike_NormRatio);
	//canvas
	c_ttjets_Split_Nbjet = TCanvasCreator(httjetsTop,httjetsBottom, line, l_ttjets_Nbjet, AxisLabel.c_str());	
	c_ttjets_Split_Nbjet->SetName("c_ttjetsRatio_Nbjet");
}
