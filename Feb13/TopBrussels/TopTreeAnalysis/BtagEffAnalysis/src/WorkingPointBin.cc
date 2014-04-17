#include "../interface/WorkingPointBin.h"

WorkingPointBin::WorkingPointBin(int nWP){

  TString str="WorkingPoint_"; str+=nWP;
  genericName_ = new TString(str);

}

WorkingPointBin::~WorkingPointBin(){}

void WorkingPointBin::DefineHistos(int nBinsEff,double lowRangeEff,double upRangeEff,int nBinsEffDiff,double lowRangeEffDiff,double upRangeEffDiff,int nBinsEffPull,double lowRangeEffPull,double upRangeEffPull,int nBinsFRatio,double lowRangeFRatio,double upRangeFRatio,int nBinsFRatioDiff,double lowRangeFRatioDiff,double upRangeFRatioDiff,int nBinsFRatioPull,double lowRangeFRatioPull,double upRangeFRatioPull,int nBinsPar0,double lowRangePar0,double upRangePar0,int nBinsPar1,double lowRangePar1,double upRangePar1){

  //  cout << "WorkingPointBin::DefineHistos - start " << endl;
  title_TH1D_Eff_=new TString(); title_TH1D_Eff_->Append(*genericName_); title_TH1D_Eff_->Append("_Eff");
  TH1D_Eff_ = new TH1D(*title_TH1D_Eff_,*title_TH1D_Eff_,nBinsEff,lowRangeEff,upRangeEff);

  title_TH1D_EffMC_=new TString(); title_TH1D_EffMC_->Append(*genericName_); title_TH1D_EffMC_->Append("_EffMC");
  TH1D_EffMC_ = new TH1D(*title_TH1D_EffMC_,*title_TH1D_EffMC_,nBinsEff,lowRangeEff,upRangeEff);

  title_TH1D_EffErr_=new TString(); title_TH1D_EffErr_->Append(*genericName_); title_TH1D_EffErr_->Append("_EffErr");
  TH1D_EffErr_ = new TH1D(*title_TH1D_EffErr_,*title_TH1D_EffErr_,nBinsEff,lowRangeEff,upRangeEff);

  title_TH1D_EffDiff_=new TString(); title_TH1D_EffDiff_->Append(*genericName_); title_TH1D_EffDiff_->Append("_EffDiff");
  TH1D_EffDiff_ = new TH1D(*title_TH1D_EffDiff_,*title_TH1D_EffDiff_,nBinsEffDiff,lowRangeEffDiff,upRangeEffDiff);

  title_TH1D_EffPull_=new TString(); title_TH1D_EffPull_->Append(*genericName_); title_TH1D_EffPull_->Append("_EffPull");
  TH1D_EffPull_ = new TH1D(*title_TH1D_EffPull_,*title_TH1D_EffPull_,nBinsEffPull,lowRangeEffPull,upRangeEffPull);

  title_TH1D_FRatio_=new TString(); title_TH1D_FRatio_->Append(*genericName_); title_TH1D_FRatio_->Append("_FRatio");
  TH1D_FRatio_ = new TH1D(*title_TH1D_FRatio_,*title_TH1D_FRatio_,nBinsFRatio,lowRangeFRatio,upRangeFRatio);

  title_TH1D_FRatioMC_=new TString(); title_TH1D_FRatioMC_->Append(*genericName_); title_TH1D_FRatioMC_->Append("_FRatioMC");
  TH1D_FRatioMC_ = new TH1D(*title_TH1D_FRatioMC_,*title_TH1D_FRatioMC_,nBinsFRatio,lowRangeFRatio,upRangeFRatio);

  title_TH1D_FRatioErr_=new TString(); title_TH1D_FRatioErr_->Append(*genericName_); title_TH1D_FRatioErr_->Append("_FRatioErr");
  TH1D_FRatioErr_ = new TH1D(*title_TH1D_FRatioErr_,*title_TH1D_FRatioErr_,nBinsFRatio,lowRangeFRatio,upRangeFRatio);

  title_TH1D_FRatioMCErr_=new TString(); title_TH1D_FRatioMCErr_->Append(*genericName_); title_TH1D_FRatioMCErr_->Append("_FRatioMCErr");
  TH1D_FRatioMCErr_ = new TH1D(*title_TH1D_FRatioMCErr_,*title_TH1D_FRatioMCErr_,nBinsFRatio,lowRangeFRatio,upRangeFRatio);

  title_TH1D_FRatioDiff_=new TString(); title_TH1D_FRatioDiff_->Append(*genericName_); title_TH1D_FRatioDiff_->Append("_FRatioDiff");
  TH1D_FRatioDiff_ = new TH1D(*title_TH1D_FRatioDiff_,*title_TH1D_FRatioDiff_,nBinsFRatioDiff,lowRangeFRatioDiff,upRangeFRatioDiff);

  title_TH1D_FRatioPull_=new TString(); title_TH1D_FRatioPull_->Append(*genericName_); title_TH1D_FRatioPull_->Append("_FRatioPull");
  TH1D_FRatioPull_ = new TH1D(*title_TH1D_FRatioPull_,*title_TH1D_FRatioPull_,nBinsFRatioPull,lowRangeFRatioPull,upRangeFRatioPull);

  title_TH1D_FitParam_[0]=new TString(); title_TH1D_FitParam_[0]->Append(*genericName_); title_TH1D_FitParam_[0]->Append("_FitParam1");
  TH1D_FitParam_[0] = new TH1D(*title_TH1D_FitParam_[0],*title_TH1D_FitParam_[0],nBinsPar0,lowRangePar0,upRangePar0);

  title_TH1D_FitParam_[1]=new TString(); title_TH1D_FitParam_[1]->Append(*genericName_); title_TH1D_FitParam_[1]->Append("_FitParam2");
  TH1D_FitParam_[1] = new TH1D(*title_TH1D_FitParam_[1],*title_TH1D_FitParam_[1],nBinsPar1,lowRangePar1,upRangePar1);

}

void WorkingPointBin::FillHistos(double* WParray){

  //cout << "WorkingPointBin::FillHistos : " << WParray[10] << endl;

  TH1D_Eff_->Fill(WParray[8]);
  TH1D_EffErr_->Fill(WParray[9]);
  TH1D_EffMC_->Fill(WParray[6]);
  TH1D_EffDiff_->Fill(WParray[6]-WParray[8]);
  TH1D_EffPull_->Fill((WParray[6]-WParray[8])/WParray[9]);

  TH1D_FRatio_->Fill(WParray[14]);
  TH1D_FRatioErr_->Fill(WParray[15]);
  TH1D_FRatioMC_->Fill(WParray[10]);
  TH1D_FRatioMCErr_->Fill(WParray[11]);

  TH1D_FRatioDiff_->Fill(WParray[10]-WParray[14]);
  TH1D_FRatioPull_->Fill((WParray[10]-WParray[14])/WParray[15]);

  TH1D_FitParam_[0]->Fill(WParray[18]);
  TH1D_FitParam_[1]->Fill(WParray[20]);

}

void WorkingPointBin::WriteHistos(TFile *fout){

  TString dirNameWP;

  dirNameWP=*genericName_;
  TDirectory *dirWP = fout->mkdir(dirNameWP);
  dirWP->cd();

  TH1D_Eff_->Write();
  TH1D_EffErr_->Write();
  TH1D_EffMC_->Write();
  TH1D_EffDiff_->Write();
  TH1D_EffPull_->Write();
  TH1D_FRatio_->Write();
  TH1D_FRatioErr_->Write();
  TH1D_FRatioMC_->Write();
  TH1D_FRatioMCErr_->Write();
  TH1D_FRatioDiff_->Write();
  TH1D_FRatioPull_->Write();
  TH1D_FitParam_[0]->Write();
  TH1D_FitParam_[1]->Write();
}
