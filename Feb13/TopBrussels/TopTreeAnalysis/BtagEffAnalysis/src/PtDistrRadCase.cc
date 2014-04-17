//REMARK:
//
//This class will contain 2 pt distributions, One in the case there is a radiation jet and on in the case there is not a radiation jet in the selected jets. The class is called inside the myNTupleAnalyzer and writes the plots in NTupleAnalyzed.cc

#include "../interface/PtDistrRadCase.h"


PtDistrRadCase::PtDistrRadCase(TString inStr){

  TString strT=inStr; strT += "_RadInside";
  TString strF=inStr; strF += "_NoRadInside";
  trueName_ = new TString(strT);
  falseName_ = new TString(strF);

  Pt_true_ = new TH1D(*trueName_,*trueName_,100,0,400);
  Pt_false_ = new TH1D(*falseName_,*falseName_,100,0,400);
}

PtDistrRadCase::~PtDistrRadCase(){}

void PtDistrRadCase::FillPlots(bool ContainsR, double pt){

  if(ContainsR) Pt_true_->Fill(pt);
  if(!ContainsR) Pt_false_->Fill(pt);

}

void PtDistrRadCase::Write(){

  Pt_true_->Write();
  Pt_false_->Write();

}

void PtDistrRadCase::WritePS(TString inStr){

  double Pt_true_Scale = Pt_true_->Integral(0,Pt_true_->GetXaxis()->GetNbins()+1);
  double Pt_false_Scale = Pt_false_->Integral(0,Pt_false_->GetXaxis()->GetNbins()+1);

  Pt_true_->Scale(1/Pt_true_Scale);
  Pt_false_->Scale(1/Pt_false_Scale);

  TString fname = inStr; fname+=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();
 
  Pt_true_->GetYaxis()->SetRangeUser(0,0.1);
  Pt_true_->Draw();
  Pt_false_->Draw("same");
  Pt_false_->SetLineColor(kRed);

  TLegend leg(0.55,0.8,0.9,0.9);
  leg.AddEntry(Pt_true_,"Radiation in sel.jets","l");
  leg.AddEntry(Pt_false_,"No radiation in sel.jets","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");  

  Pt_true_->Scale(Pt_true_Scale);
  Pt_false_->Scale(Pt_false_Scale);

}
