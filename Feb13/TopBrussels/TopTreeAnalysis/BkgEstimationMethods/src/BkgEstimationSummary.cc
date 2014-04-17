#include "../interface/BkgEstimationSummary.h"


BkgEstimationSummary::BkgEstimationSummary(){
  DataColor_ = 1;
  QCDColor_ = 50;
  WJetColor_ = 9;
  TtJetColor_ = 8;
  
  QCDLColor_ = 2;
  WJetLColor_ = 4;
  TtJetLColor_ = 3;
  LWidth_ = 2;
  hData_ = 0;
  hQCDEstim_ = 0;
  hWJetEstim_ = 0;
  hTtJetEstim_ = 0;
  hCompleteEstim_ = 0;
  hQCDMC_ = 0;
  hWJetMC_ = 0;
  hTtJetMC_ = 0;
  hCompleteMC_ = 0;
  legend = 0;
  SumCanvas = 0;
  MCStackHisto = 0;
}


BkgEstimationSummary::BkgEstimationSummary(TString Name, TString XaxisTitle, TString YaxisTitle, TH1F* hData, TH1F* hQCDEstim, TH1F* hWJetEstim, TH1F* hTtJetEstim){
  DataColor_ = 1;
  QCDColor_ = 50;
  WJetColor_ = 9;
  TtJetColor_ = 8;
  
  QCDLColor_ = 2;
  WJetLColor_ = 4;
  TtJetLColor_ = 3;
  
  LWidth_ = 2;
  Name_ = Name;
  XaxisTitle_ = XaxisTitle;
  YaxisTitle_ = YaxisTitle;
  hQCDMC_ = 0;
  hWJetMC_ = 0;
  hTtJetMC_ = 0;
  hCompleteMC_ = 0;
  hData_ = hData;
  hQCDEstim_ = hQCDEstim;
  hWJetEstim_ = hWJetEstim;
  hTtJetEstim_ = hTtJetEstim;
  //
  hCompleteEstim_ = (TH1F*) hQCDEstim->Clone(TString(hQCDEstim->GetName())+TString("_Complete"));
  hCompleteEstim_->Add(hWJetEstim);
  hCompleteEstim_->Add(hTtJetEstim);
  //
  legend = 0;
  SumCanvas = 0;
  MCStackHisto = 0;
}

BkgEstimationSummary::BkgEstimationSummary(TString Name, TString XaxisTitle, TString YaxisTitle, TH1F* hData, TH1F* hQCDEstim, TH1F* hQCDMC, TH1F* hWJetEstim, TH1F* hWJetMC, TH1F* hTtJetEstim, TH1F* hTtJetMC){
  DataColor_ = 1;
  QCDColor_ = 50;
  WJetColor_ = 9;
  TtJetColor_ = 8;
  
  QCDLColor_ = 2;
  WJetLColor_ = 4;
  TtJetLColor_ = 3;
  
  LWidth_ = 2;
  Name_ = Name;
  XaxisTitle_ = XaxisTitle;
  YaxisTitle_ = YaxisTitle;
  hQCDMC_ = hQCDMC;
  hWJetMC_ = hWJetMC;
  hTtJetMC_ = hTtJetMC;
  //
  hCompleteMC_ = (TH1F*) hQCDMC->Clone(TString(hQCDMC->GetName())+TString("_Complete"));
  hCompleteMC_->Add(hWJetMC);
  hCompleteMC_->Add(hTtJetMC);
  //
  hData_ = hData;
  hQCDEstim_ = hQCDEstim;
  hWJetEstim_ = hWJetEstim;
  hTtJetEstim_ = hTtJetEstim;
  //
  hCompleteEstim_ = (TH1F*) hQCDEstim->Clone(TString(hQCDEstim->GetName())+TString("_Complete"));
  hCompleteEstim_->Add(hWJetEstim);
  hCompleteEstim_->Add(hTtJetEstim);
  //
  legend = 0;
  SumCanvas = 0;
  MCStackHisto = 0;
}

void BkgEstimationSummary::Draw(){
  SumCanvas  = new TCanvas(Name_);
  
  // MC prediction
  vector<TH1F*> vec;
  TH1F* hQCDMCClone = 0;
  TH1F* hWJetMCClone = 0;
  TH1F* hTtJetMCClone = 0;
  
  if(hQCDMC_){
  	hQCDMCClone = (TH1F*) hQCDMC_->Clone();
  	hQCDMCClone->SetFillColor(QCDColor_);
  	vec.push_back(hQCDMCClone);
  }
  if(hWJetMC_){
  	hWJetMCClone = (TH1F*) hWJetMC_->Clone();
  	hWJetMCClone->SetFillColor(WJetColor_);
  	vec.push_back(hWJetMCClone);
  }
  if(hTtJetMC_){
  	hTtJetMCClone = (TH1F*) hTtJetMC_->Clone();
  	hTtJetMCClone->SetFillColor(TtJetColor_);
  	vec.push_back(hTtJetMCClone);
  }
  if(vec.size()>0) MCStackHisto = THStackCreator(vec); 
  
  // estimation
  TH1F* hQCDEstimClone = 0;
  TH1F* hWJetEstimClone = 0;
  TH1F* hTtJetEstimClone = 0;
  if(hQCDEstim_){
  	hQCDEstimClone = (TH1F*) hQCDEstim_->Clone();
  	hQCDEstimClone->SetLineColor(QCDLColor_);
  	hQCDEstimClone->SetLineWidth(LWidth_);
  }
  if(hWJetEstim_){
  	hWJetEstimClone = (TH1F*) hWJetEstim_->Clone();
  	hWJetEstimClone->SetLineColor(WJetLColor_);
  	hWJetEstimClone->SetLineWidth(LWidth_);
  }
  if(hTtJetEstim_){
  	hTtJetEstimClone = (TH1F*) hTtJetEstim_->Clone();
  	hTtJetEstimClone->SetLineColor(TtJetLColor_);
  	hTtJetEstimClone->SetLineWidth(LWidth_);
  }
  //summing histo
  if(hWJetEstimClone && hQCDEstimClone) hWJetEstimClone->Add(hQCDEstimClone);
  if(hTtJetEstimClone && hWJetEstimClone) hTtJetEstimClone->Add(hWJetEstimClone);

  //Plotting
  SumCanvas->cd();
  hData_->Sumw2();
  hData_->GetXaxis()->SetTitle(XaxisTitle_);
  hData_->GetYaxis()->SetTitle(YaxisTitle_);
  hData_->GetYaxis()->SetRangeUser(0., hData_->GetMaximum()*1.2);
  hData_->Draw("e");
  if(MCStackHisto) MCStackHisto->Draw("same");
  if(hQCDEstimClone) hQCDEstimClone->Draw("same");
  if(hWJetEstimClone) hWJetEstimClone->Draw("same");
  if(hTtJetEstimClone) hTtJetEstimClone->Draw("same");
  hData_->Draw("esame");
  //Tlegend
  legend = new TLegend(0.65,0.55,0.88,0.88);
  if(hData_) legend->AddEntry(hData_,"Data","l");
  if(hQCDMCClone) legend->AddEntry(hQCDMCClone,"QCD - MC","F");
  if(hWJetMCClone) legend->AddEntry(hWJetMCClone,"WJets - MC","F");
  if(hTtJetMCClone) legend->AddEntry(hTtJetMCClone,"TtJets - MC","F");
  //pred
  if(hQCDEstimClone) legend->AddEntry(hQCDEstimClone,"QCD - Estim","l");
  if(hWJetEstimClone) legend->AddEntry(hWJetEstimClone,"WJets - Estim","l");
  if(hTtJetEstimClone) legend->AddEntry(hTtJetEstimClone,"TtJets - Estim","l");
  legend->Draw("same");
}

void BkgEstimationSummary::Write(TFile* fout, string label){
 	fout->cd();
	string dirname = "BkgSummary"+label;
	fout->mkdir(dirname.c_str());
	fout->cd(dirname.c_str());				
	if(SumCanvas) SumCanvas->Write();
}

BkgEstimationSummary::~BkgEstimationSummary(){
  if(hData_) delete hData_;
  if(hQCDEstim_) delete hQCDEstim_;
  if(hWJetEstim_) delete hWJetEstim_;
  if(hTtJetEstim_) delete hTtJetEstim_;
  if(hCompleteEstim_) delete hCompleteEstim_;
  if(hQCDMC_) delete hQCDMC_;
  if(hWJetMC_) delete hWJetMC_;
  if(hTtJetMC_) delete hTtJetMC_;		  
  if(hCompleteMC_) delete hCompleteMC_;
  if(legend) delete legend;
  if(SumCanvas) delete SumCanvas;
  if(MCStackHisto) delete MCStackHisto;
}
