#include "../interface/TtJetEstimation.h"


TtJetEstimation::TtJetEstimation (bool isMC)
{
  isMC_  = isMC;
  NbinsX_ = 100;
  Xmin_ = 0;
  Xmax_ = 1000;
  ObsName_ = string ("Observable");
  AxisLabel_ = string ("Observable");
  h_ObservableCR_ = 0;
  h_ObservableSR_ = 0;
  h_ObservableEstimated_ = 0;
  h_ObservableExpected_ = 0;
  h_ObservableObserved_ = 0;
  h_SysError_ = 0;
  scaled_ = false;
  normalized_ = false;
  binningFixed_ = true;
  canvasRatioEstim_Exp = 0;
  canvasRatioCR_SR = 0;
  canvasCR = 0;
  canvasSR = 0;
  canvasSign = 0;
  canvasSignRatio = 0;
  canvasLMXContamination = 0;
  NRegion_.first = 0;
  NRegion_.second = 0;
  binsX_ = 0;
  line = 0;
  leg = 0;
  legSign = 0;
  legSignRatio = 0;
  leg_LMXContamination = 0;
  func = 0;
}

TtJetEstimation::~TtJetEstimation ()
{
  if (h_ObservableCR_)
    delete h_ObservableCR_;
  if (h_ObservableSR_)
    delete h_ObservableSR_;
  if (h_ObservableExpected_)
    delete h_ObservableExpected_;
  if (h_ObservableObserved_)
    delete h_ObservableObserved_;
  if (h_SysError_)
    delete h_SysError_;
  if (line)
    delete line;
  if (leg)
    delete leg;
  if (legSign)
    delete legSign;
  if (legSignRatio)
    delete legSignRatio;
  if (func)
    delete func;
  if (plotRatioCR_SR)
    delete plotRatioCR_SR;
  if (plotRatioEstim_Exp)
    delete plotRatioEstim_Exp;
  if (canvasSR)
    delete canvasSR;
  if (canvasCR)
    delete canvasCR;
  if (canvasRatioCR_SR)
    delete canvasRatioCR_SR;
  if (canvasRatioEstim_Exp)
    delete canvasRatioEstim_Exp;
  if (canvasSign)
    delete canvasSign;
  if (canvasSignRatio)
    delete canvasSignRatio;
  if (leg_LMXContamination)
    delete leg_LMXContamination;
  if (canvasLMXContamination)
    delete canvasLMXContamination;
  //if(h_ObservableEstimated_) delete h_ObservableEstimated_;
}

TtJetEstimation::TtJetEstimation (const TtJetEstimation & ttj)
{
  canvasRatioEstim_Exp = 0;
  canvasRatioCR_SR = 0;
  canvasCR = 0;
  canvasSR = 0;
  canvasSign = 0;
  canvasSignRatio = 0;
  canvasLMXContamination = 0;
  line = 0;
  leg = 0;
  legSign = 0;
  legSignRatio = 0;
  leg_LMXContamination = 0;
  func = 0;
  h_ObservableCR_ = 0;
  h_ObservableSR_ = 0;
  h_ObservableEstimated_ = 0;
  h_ObservableExpected_ = 0;
  h_ObservableObserved_ = 0;
  h_SysError_ = 0;
  scaled_ = ttj.scaled_;
  normalized_ = ttj.normalized_;
  ObsName_ = ttj.ObsName_;
  AxisLabel_ = ttj.AxisLabel_;
  NbinsX_ = ttj.NbinsX_;
  binningFixed_ = ttj.binningFixed_;
  binsX_ = new float[ttj.NbinsX_];
  for (int i = 0; i < ttj.NbinsX_; i++)
    binsX_[i] = ttj.binsX_[i];
  Xmin_ = ttj.Xmin_;
  Xmax_ = ttj.Xmax_;
  listOfDatasets_ = ttj.listOfDatasets_;
  listOfDatasetsTtJets_ = ttj.listOfDatasetsTtJets_;
  listOfDatasetsBkg_ = ttj.listOfDatasetsBkg_;
  listOfDatasetsNP_ = ttj.listOfDatasetsNP_;
  NRegion_ = ttj.NRegion_;
  if (ttj.h_ObservableCR_)
    h_ObservableCR_ = (TH1F *) ttj.h_ObservableCR_->Clone ();
  if (ttj.h_ObservableSR_)
    h_ObservableSR_ = (TH1F *) ttj.h_ObservableSR_->Clone ();
  if (ttj.h_ObservableEstimated_)
    h_ObservableEstimated_ = (TH1F *) ttj.h_ObservableEstimated_->Clone ();
  if (ttj.h_ObservableExpected_)
    h_ObservableExpected_ = (TH1F *) ttj.h_ObservableExpected_->Clone ();
  if (ttj.h_ObservableObserved_)
    h_ObservableObserved_ = (TH1F *) ttj.h_ObservableObserved_->Clone ();
  if (ttj.h_SysError_)
    h_SysError_ = (TH1F *) ttj.h_SysError_->Clone ();

  for (unsigned int i = 0; i < ttj.h_ObservableByDatasetsTtJets_.size (); i++) {
    h_ObservableByDatasetsTtJets_.push_back ((TH1F *) ttj.h_ObservableByDatasetsTtJets_[i]->Clone ());
  }
  for (unsigned int i = 0; i < ttj.h_ObservableByDatasetsBkg_.size (); i++) {
    h_ObservableByDatasetsBkg_.push_back ((TH1F *) ttj.h_ObservableByDatasetsBkg_[i]->Clone ());
  }
  for (unsigned int i = 0; i < ttj.h_ObservableByDatasetsNP_.size (); i++) {
    h_ObservableByDatasetsNP_.push_back ((TH1F *) ttj.h_ObservableByDatasetsNP_[i]->Clone ());
  }
  for (unsigned int i = 0; i < ttj.h_ObservableCRByDatasetsTtJets_.size (); i++) {
    h_ObservableCRByDatasetsTtJets_.push_back ((TH1F *) ttj.h_ObservableCRByDatasetsTtJets_[i]->Clone ());
  }
  for (unsigned int i = 0; i < ttj.h_ObservableCRByDatasetsBkg_.size (); i++) {
    h_ObservableCRByDatasetsBkg_.push_back ((TH1F *) ttj.h_ObservableCRByDatasetsBkg_[i]->Clone ());
  }
  for (unsigned int i = 0; i < ttj.h_ObservableCRByDatasetsNP_.size (); i++) {
    h_ObservableCRByDatasetsNP_.push_back ((TH1F *) ttj.h_ObservableCRByDatasetsNP_[i]->Clone ());
  }
  for (unsigned int i = 0; i < ttj.h_ObservableSRByDatasetsTtJets_.size (); i++) {
    h_ObservableSRByDatasetsTtJets_.push_back ((TH1F *) ttj.h_ObservableSRByDatasetsTtJets_[i]->Clone ());
  }
  for (unsigned int i = 0; i < ttj.h_ObservableSRByDatasetsBkg_.size (); i++) {
    h_ObservableSRByDatasetsBkg_.push_back ((TH1F *) ttj.h_ObservableSRByDatasetsBkg_[i]->Clone ());
  }
  for (unsigned int i = 0; i < ttj.h_ObservableSRByDatasetsNP_.size (); i++) {
    h_ObservableSRByDatasetsNP_.push_back ((TH1F *) ttj.h_ObservableSRByDatasetsNP_[i]->Clone ());
  }
  for (unsigned int i = 0; i < ttj.h_SignEstimated_.size (); i++) {
    h_SignEstimated_.push_back ((TH1F *) ttj.h_SignEstimated_[i]->Clone ());
  }
  for (unsigned int i = 0; i < ttj.h_SignExpected_.size (); i++) {
    h_SignExpected_.push_back ((TH1F *) ttj.h_SignExpected_[i]->Clone ());
  }

  if (ttj.line)
    line = (TLine *) ttj.line->Clone ();
  if (ttj.leg)
    leg = (TLegend *) ttj.leg->Clone ();
  if (ttj.legSign)
    legSign = (TLegend *) ttj.legSign->Clone ();
  if (ttj.legSignRatio)
    legSignRatio = (TLegend *) ttj.legSignRatio->Clone ();
  if (ttj.leg_LMXContamination)
    leg_LMXContamination = (TLegend *) ttj.leg_LMXContamination->Clone ();
  if (ttj.func)
    func = (TF1 *) ttj.func->Clone ();
  if (ttj.plotRatioCR_SR)
    plotRatioCR_SR = (TH1F *) ttj.plotRatioCR_SR->Clone ();
  if (ttj.plotRatioEstim_Exp)
    plotRatioEstim_Exp = (TH1F *) ttj.plotRatioEstim_Exp->Clone ();
  if (ttj.canvasSR)
    canvasSR = (TCanvas *) ttj.canvasSR->Clone ();
  if (ttj.canvasCR)
    canvasCR = (TCanvas *) ttj.canvasCR->Clone ();
  if (ttj.canvasRatioCR_SR)
    canvasRatioCR_SR = (TCanvas *) ttj.canvasRatioCR_SR->Clone ();
  if (ttj.canvasRatioEstim_Exp)
    canvasRatioEstim_Exp = (TCanvas *) ttj.canvasRatioEstim_Exp->Clone ();
  if (ttj.canvasSign)
    canvasSign = (TCanvas *) ttj.canvasSign->Clone ();
  if (ttj.canvasSignRatio)
    canvasSignRatio = (TCanvas *) ttj.canvasSignRatio->Clone ();
  if (ttj.canvasLMXContamination)
    canvasLMXContamination = (TCanvas *) ttj.canvasLMXContamination->Clone ();
}

void
TtJetEstimation::ConfigHistos (int NbinsX, float Xmin, float Xmax, string ObsName, string AxisLabel, bool fillMode)
{
  NbinsX_ = NbinsX;
  Xmin_ = Xmin;
  Xmax_ = Xmax;
  ObsName_ = ObsName;
  AxisLabel_ = AxisLabel;
  binningFixed_ = true;
  char hname[100];
  char title[100];
  sprintf (hname, "h_%s_CR", ObsName_.c_str ());
  sprintf (title, "%s - Control Region", ObsName_.c_str ());
  h_ObservableCR_ = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
  h_ObservableCR_->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  sprintf (hname, "h_%s_SR", ObsName_.c_str ());
  sprintf (title, "%s - Signal Region", ObsName_.c_str ());
  h_ObservableSR_ = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
  h_ObservableSR_->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  sprintf (hname, "h_%s_Est", ObsName_.c_str ());
  sprintf (title, "%s - Estimated", ObsName_.c_str ());
  h_ObservableEstimated_ = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
  h_ObservableEstimated_->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  sprintf (hname, "h_%s_Exp", ObsName_.c_str ());
  sprintf (title, "%s - Expected", ObsName_.c_str ());
  h_ObservableExpected_ = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
  h_ObservableExpected_->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  sprintf (hname, "h_%s_Obs", ObsName_.c_str ());
  sprintf (title, "%s - Observed", ObsName_.c_str ());
  h_ObservableObserved_ = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
  h_ObservableObserved_->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  //inclusive
  for (unsigned int i = 0; i < h_ObservableByDatasetsTtJets_.size (); i++) {
    sprintf (hname, "h_%s_TtJets_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - TtJets_%d", ObsName_.c_str (), i);
    h_ObservableByDatasetsTtJets_[i] = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
    h_ObservableByDatasetsTtJets_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableByDatasetsBkg_.size (); i++) {
    sprintf (hname, "h_%s_Bkg_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - Bkg_%d", ObsName_.c_str (), i);
    h_ObservableByDatasetsBkg_[i] = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
    h_ObservableByDatasetsBkg_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableByDatasetsNP_.size (); i++) {
    sprintf (hname, "h_%s_NP_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - NP_%d", ObsName_.c_str (), i);
    h_ObservableByDatasetsNP_[i] = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
    h_ObservableByDatasetsNP_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  //CR
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsTtJets_.size (); i++) {
    sprintf (hname, "h_%s_CR_TtJets_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - CR - TtJets_%d", ObsName_.c_str (), i);
    h_ObservableCRByDatasetsTtJets_[i] = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
    h_ObservableCRByDatasetsTtJets_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsBkg_.size (); i++) {
    sprintf (hname, "h_%s_CR_Bkg_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - CR - Bkg_%d", ObsName_.c_str (), i);
    h_ObservableCRByDatasetsBkg_[i] = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
    h_ObservableCRByDatasetsBkg_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsNP_.size (); i++) {
    sprintf (hname, "h_%s_CR_NP_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - CR -NP_%d", ObsName_.c_str (), i);
    h_ObservableCRByDatasetsNP_[i] = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
    h_ObservableCRByDatasetsNP_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  //SR
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsTtJets_.size (); i++) {
    sprintf (hname, "h_%s_SR_TtJets_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - SR - TtJets_%d", ObsName_.c_str (), i);
    h_ObservableSRByDatasetsTtJets_[i] = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
    h_ObservableSRByDatasetsTtJets_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsBkg_.size (); i++) {
    sprintf (hname, "h_%s_SR_Bkg_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - SR - Bkg_%d", ObsName_.c_str (), i);
    h_ObservableSRByDatasetsBkg_[i] = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
    h_ObservableSRByDatasetsBkg_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsNP_.size (); i++) {
    sprintf (hname, "h_%s_SR_NP_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - SR -NP_%d", ObsName_.c_str (), i);
    h_ObservableSRByDatasetsNP_[i] = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
    h_ObservableSRByDatasetsNP_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  //
  for (unsigned int i = 0; i < h_SignEstimated_.size (); i++) {
    sprintf (hname, "h_SignEstimated_%d", i);
    sprintf (title, "h_SignEstimated_%d", i);
    h_SignEstimated_[i] = new TH1F (hname, title, NbinsX_, Xmin, Xmax_);
    h_SignEstimated_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
    sprintf (hname, "h_SignExpected_%d", i);
    sprintf (title, "h_SignExpected_%d", i);
    h_SignExpected_[i] = new TH1F (hname, title, NbinsX_, Xmin_, Xmax_);
    h_SignExpected_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
}

void
TtJetEstimation::ConfigHistos (int NbinsX, float *binsX, string ObsName, string AxisLabel, bool fillMode)
{
  NbinsX_ = NbinsX;
  binsX_ = binsX;
  ObsName_ = ObsName;
  AxisLabel_ = AxisLabel;
  binningFixed_ = false;
  char hname[100];
  char title[100];
  sprintf (hname, "h_%s_CR", ObsName_.c_str ());
  sprintf (title, "%s - Control Region", ObsName_.c_str ());
  h_ObservableCR_ = new TH1F (hname, title, NbinsX_, binsX_);
  h_ObservableCR_->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  sprintf (hname, "h_%s_SR", ObsName_.c_str ());
  sprintf (title, "%s - Signal Region", ObsName_.c_str ());
  h_ObservableSR_ = new TH1F (hname, title, NbinsX_, binsX_);
  h_ObservableSR_->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  sprintf (hname, "h_%s_Est", ObsName_.c_str ());
  sprintf (title, "%s - Estimated", ObsName_.c_str ());
  h_ObservableEstimated_ = new TH1F (hname, title, NbinsX_, binsX_);
  h_ObservableEstimated_->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  sprintf (hname, "h_%s_Exp", ObsName_.c_str ());
  sprintf (title, "%s - Expected", ObsName_.c_str ());
  h_ObservableExpected_ = new TH1F (hname, title, NbinsX_, binsX_);
  h_ObservableExpected_->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  sprintf (hname, "h_%s_Obs", ObsName_.c_str ());
  sprintf (title, "%s - Observed", ObsName_.c_str ());
  h_ObservableObserved_ = new TH1F (hname, title, NbinsX_, binsX_);
  h_ObservableObserved_->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  //inclusive
  for (unsigned int i = 0; i < h_ObservableByDatasetsTtJets_.size (); i++) {
    sprintf (hname, "h_%s_TtJets_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - TtJets_%d", ObsName_.c_str (), i);
    h_ObservableByDatasetsTtJets_[i] = new TH1F (hname, title, NbinsX_, binsX_);
    h_ObservableByDatasetsTtJets_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableByDatasetsBkg_.size (); i++) {
    sprintf (hname, "h_%s_Bkg_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - Bkg_%d", ObsName_.c_str (), i);
    h_ObservableByDatasetsBkg_[i] = new TH1F (hname, title, NbinsX_, binsX_);
    h_ObservableByDatasetsBkg_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableByDatasetsNP_.size (); i++) {
    sprintf (hname, "h_%s_NP_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - NP_%d", ObsName_.c_str (), i);
    h_ObservableByDatasetsNP_[i] = new TH1F (hname, title, NbinsX_, binsX_);
    h_ObservableByDatasetsNP_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  //CR
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsTtJets_.size (); i++) {
    sprintf (hname, "h_%s_CR_TtJets_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - CR - TtJets_%d", ObsName_.c_str (), i);
    h_ObservableCRByDatasetsTtJets_[i] = new TH1F (hname, title, NbinsX_, binsX_);
    h_ObservableCRByDatasetsTtJets_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsBkg_.size (); i++) {
    sprintf (hname, "h_%s_CR_Bkg_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - CR - Bkg_%d", ObsName_.c_str (), i);
    h_ObservableCRByDatasetsBkg_[i] = new TH1F (hname, title, NbinsX_, binsX_);
    h_ObservableCRByDatasetsBkg_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsNP_.size (); i++) {
    sprintf (hname, "h_%s_CR_NP_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - CR -NP_%d", ObsName_.c_str (), i);
    h_ObservableCRByDatasetsNP_[i] = new TH1F (hname, title, NbinsX_, binsX_);
    h_ObservableCRByDatasetsNP_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  //SR
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsTtJets_.size (); i++) {
    sprintf (hname, "h_%s_SR_TtJets_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - SR - TtJets_%d", ObsName_.c_str (), i);
    h_ObservableSRByDatasetsTtJets_[i] = new TH1F (hname, title, NbinsX_, binsX_);
    h_ObservableSRByDatasetsTtJets_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsBkg_.size (); i++) {
    sprintf (hname, "h_%s_SR_Bkg_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - SR - Bkg_%d", ObsName_.c_str (), i);
    h_ObservableSRByDatasetsBkg_[i] = new TH1F (hname, title, NbinsX_, binsX_);
    h_ObservableSRByDatasetsBkg_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsNP_.size (); i++) {
    sprintf (hname, "h_%s_SR_NP_%d", ObsName_.c_str (), i);
    sprintf (title, "%s - SR -NP_%d", ObsName_.c_str (), i);
    h_ObservableSRByDatasetsNP_[i] = new TH1F (hname, title, NbinsX_, binsX_);
    h_ObservableSRByDatasetsNP_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
  }
  //
  for (unsigned int i = 0; i < h_SignEstimated_.size (); i++) {
    sprintf (hname, "h_SignEstimated_%d", i);
    sprintf (title, "h_SignEstimated_%d", i);
    h_SignEstimated_[i] = new TH1F (hname, title, NbinsX_, binsX_);
    h_SignEstimated_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
    h_SignEstimated_[i]->SetLineColor (listOfDatasetsNP_[i].Color ());
    h_SignEstimated_[i]->SetLineStyle (1);
    sprintf (hname, "h_SignExpected_%d", i);
    sprintf (title, "h_SignExpected_%d", i);
    h_SignExpected_[i] = new TH1F (hname, title, NbinsX_, binsX_);
    h_SignExpected_[i]->GetXaxis ()->SetTitle (AxisLabel_.c_str ());
    h_SignExpected_[i]->SetLineColor (listOfDatasetsNP_[i].Color ());
    h_SignExpected_[i]->SetLineStyle (2);
  }
}

void
TtJetEstimation::Fill (float variable, bool trigged, Selection& selection, AnalysisEnvironment& anaEnv, bool isInSR, string DatasetName)
{
  //determine the nof b-jets
  vector < TRootJet* > selectedJets = selection.GetSelectedJets (anaEnv.JetsPtCutCR, anaEnv.JetsEtaCutCR, anaEnv.applyJetID);
  //////////////////////////
  //btag jets collection
  //////////////////////////
  vector < TRootJet* > btagJets = selection.GetSelectedBJets(selectedJets,anaEnv.BtagAlgo_ttjEst, anaEnv.BtagDiscriCut_ttjEst);
  ////////////////////


  //define if this events belongs to CR
  bool isInCR = false;
  //Control Region for TtJetEstimation
  if (trigged && selection.GetSelectedMuons (anaEnv.MuonPtCutCR, anaEnv.MuonEtaCutCR, anaEnv.MuonRelIsoCutCR, selectedJets).size () == 1
      && selection.GetSelectedMuons (anaEnv.MuonPtCutVetoSR, anaEnv.MuonEtaCutSR, anaEnv.MuonRelIsoCutCR, selectedJets).size () == 1
      && selection.GetSelectedElectrons (anaEnv.ElectronPtCut, anaEnv.ElectronEtaCut, anaEnv.ElectronRelIsoCut).size () == 0 && selection.GetSelectedJets (anaEnv.JetsPtCutCR, anaEnv.JetsEtaCutCR).size () >= (unsigned int) anaEnv.NofJets) {
    if (btagJets.size () > 1) {
      TRootMuon* muon = selection.GetSelectedMuons (anaEnv.MuonPtCutCR, anaEnv.MuonEtaCutCR, anaEnv.MuonRelIsoCutCR, selectedJets)[0];
      //compute Mbl
      float Mbl = -9999;
      float DRblMin = 9999;
      for (unsigned int bj = 0; bj < btagJets.size (); bj++) {
	if (bj > 1)
	  break;
	float DeltaR = ROOT::Math::VectorUtil::DeltaR (*btagJets[bj], *muon);
	if (DeltaR < DRblMin) {
	  DRblMin = DeltaR;
	  TLorentzVector BL = *btagJets[bj] + *selection.GetSelectedMuons (anaEnv.MuonPtCutCR, anaEnv.MuonEtaCutCR, anaEnv.MuonRelIsoCutCR, selectedJets)[0];
	  Mbl = BL.M ();
	}
      }
      float HTBB = btagJets[0]->Pt () + btagJets[1]->Pt ();
      float DRBB = ROOT::Math::VectorUtil::DeltaR (*btagJets[0], *btagJets[1]);
      if (Mbl < anaEnv.MblCut && HTBB < anaEnv.HTBBCut && DRBB > anaEnv.DRBBCut) {
	isInCR = true;
      }
    }
  }
  // isInCR is filled

  if(isMC_){
  //Call other Fill methods depending on the dataset name 
  //Bgk
  if (DatasetName == "QCD" || DatasetName == "WJets" || DatasetName == "SingleTop") {
    int idataset = -1;
    for (unsigned int dd = 0; dd < listOfDatasetsBkg_.size (); dd++)
      if (DatasetName == listOfDatasetsBkg_[dd].Name ())
	idataset = dd;
    FillByDatasetsBkg (variable, idataset);
    if (isInSR)
      FillSRByDatasetsBkg (variable, idataset);
    if (isInCR)
      FillCRByDatasetsBkg (variable, idataset);
  }
  //TtJets
  if (DatasetName == "TTJets") {
    FillByDatasetsTtJets (variable, 0);
    if (isInSR)
      FillSRByDatasetsTtJets (variable, 0);
    if (isInCR)
      FillCRByDatasetsTtJets (variable, 0);
  }
  //SUSY
  if (DatasetName.find ("LM") < DatasetName.size ()) {
    int idataset = 0;
    for (unsigned int dd = 0; dd < listOfDatasetsNP_.size (); dd++)
      if (DatasetName == listOfDatasetsNP_[dd].Name ())
	idataset = dd;
    FillByDatasetsNP (variable, idataset);
    if (isInSR)
      FillSRByDatasetsNP (variable, idataset);
    if (isInCR)
      FillCRByDatasetsNP (variable, idataset);
  }
  }
  else{
    if(isInCR) h_ObservableCR_->Fill (variable);
    if(isInCR) h_ObservableEstimated_->Fill (variable);
    if(isInSR) h_ObservableSR_->Fill (variable);
    if(isInSR) h_ObservableObserved_->Fill (variable);
  }

}

void
TtJetEstimation::DrawHisto (Dataset dataset, TH1F * histo)
{
  histo->SetFillColor (dataset.Color ());
  histo->SetLineColor (dataset.Color ());
  histo->SetLineWidth (dataset.LineWidth ());
  histo->SetLineStyle (dataset.LineStyle ());
}

void
TtJetEstimation::Draw (string label)
{
  //TtJets
  for (unsigned int i = 0; i < listOfDatasetsTtJets_.size (); i++) {
    DrawHisto (listOfDatasetsTtJets_[i], h_ObservableByDatasetsTtJets_[i]);
    DrawHisto (listOfDatasetsTtJets_[i], h_ObservableCRByDatasetsTtJets_[i]);
    DrawHisto (listOfDatasetsTtJets_[i], h_ObservableSRByDatasetsTtJets_[i]);

  }
  //Bkg
  for (unsigned int i = 0; i < listOfDatasetsBkg_.size (); i++) {
    DrawHisto (listOfDatasetsBkg_[i], h_ObservableByDatasetsBkg_[i]);
    DrawHisto (listOfDatasetsBkg_[i], h_ObservableCRByDatasetsBkg_[i]);
    DrawHisto (listOfDatasetsBkg_[i], h_ObservableSRByDatasetsBkg_[i]);
  }
  //NP
  for (unsigned int i = 0; i < listOfDatasetsNP_.size (); i++) {
    DrawHisto (listOfDatasetsNP_[i], h_ObservableByDatasetsNP_[i]);
    DrawHisto (listOfDatasetsNP_[i], h_ObservableCRByDatasetsNP_[i]);
    DrawHisto (listOfDatasetsNP_[i], h_ObservableSRByDatasetsNP_[i]);
  }
  for (unsigned int i = 0; i < listOfDatasets_.size (); i++) {
    DrawHisto (listOfDatasets_[i], PlotByDatasets ()[i]);
    DrawHisto (listOfDatasets_[i], PlotByDatasetsSR ()[i]);
    DrawHisto (listOfDatasets_[i], PlotByDatasetsCR ()[i]);
  }

  // order: fist the plots, then the canvas !!

  bool bfit = false;
  ///////////////////////////
  // Ratio CR SR for ttjets
  ///////////////////////////
  func = new TF1 ("func", "[0]*x+[1]", PlotCR ()->GetXaxis ()->GetXmin (), PlotCR ()->GetXaxis ()->GetXmax ());
  vector < TH1F * >vec;
  if (PlotByDatasetsCRTtJets ().size () > 0)
    vec.push_back ((TH1F *) PlotByDatasetsCRTtJets ()[0]->Clone ());
  if (PlotByDatasetsSRTtJets ().size () > 0)
    vec.push_back ((TH1F *) PlotByDatasetsSRTtJets ()[0]->Clone ());
  if (vec.size () > 0 && vec[0]->Integral () > 0) {
    vec[0]->Sumw2 ();
    vec[0]->Scale (1 / vec[0]->Integral ());
    vec[0]->SetLineColor (2);
    vec[0]->SetLineWidth (2);
    vec[0]->SetFillStyle (0);
  }
  if (vec.size () > 1 && vec[1]->Integral () > 0) {
    vec[1]->Sumw2 ();
    vec[1]->Scale (1 / vec[1]->Integral ());
    vec[1]->SetLineColor (4);
    vec[1]->SetLineWidth (2);
    vec[1]->SetFillStyle (0);
  }
  TLegend *legCR_SR = new TLegend (0.6, 0.6, 0.88, 0.88);
  if (vec.size () > 0)
    legCR_SR->AddEntry (vec[0], "Control Region", "l");
  if (vec.size () > 1)
    legCR_SR->AddEntry (vec[1], "Signal Region", "l");
  if (bfit) {
    plotRatioCR_SR = PlotRatio (vec[0], vec[1], func, string ("CR_SR_Ratio"));
    plotRatioCR_SR->GetFunction ("func")->SetLineColor (2);
  }
  else
    plotRatioCR_SR = PlotRatio (vec[0], vec[1], string ("CR_SR_Ratio"));
  plotRatioCR_SR->SetFillStyle (0);
  plotRatioCR_SR->SetLineColor (1);

  ///////////////////////////
  // Ratio Estimation - Expectation
  ///////////////////////////
  if (bfit) {
    plotRatioEstim_Exp = PlotRatio (PlotEstimated (), PlotExpected (), func, string ("Est_Exp_Ratio"));
    plotRatioEstim_Exp->GetFunction ("func")->SetLineColor (2);
  }
  else
    plotRatioEstim_Exp = PlotRatio (PlotEstimated (), PlotExpected (), string ("Est_Exp_Ratio"));
  plotRatioEstim_Exp->SetFillStyle (0);
  plotRatioEstim_Exp->SetLineColor (1);
  vector < TH1F * >vec2;
  vec2.clear ();
  vec2.push_back ((TH1F *) PlotEstimated ()->Clone ());
  vec2.push_back ((TH1F *) PlotExpected ()->Clone ());
  vec2[0]->SetLineColor (2);
  vec2[0]->SetLineWidth (2);
  vec2[1]->SetLineColor (4);
  vec2[1]->SetLineWidth (2);
  TLegend *legEstim_Exp = new TLegend (0.6, 0.6, 0.88, 0.88);
  if (vec2.size () > 0)
    legEstim_Exp->AddEntry (vec2[0], "Estimated", "l");
  if (vec2.size () > 1)
    legEstim_Exp->AddEntry (vec2[1], "Expected", "l");


  ///////////////////////////
  //  Canvas
  ///////////////////////////
  canvasRatioCR_SR = TCanvasCreator (vec, plotRatioCR_SR, line, legCR_SR, string (""), string ("CRDivSR"));
  canvasRatioEstim_Exp = TCanvasCreator (vec2, plotRatioEstim_Exp, line, legEstim_Exp, string (""), string ("EstDivExp"));
  canvasSR = TCanvasCreator (PlotByDatasetsSR (), ListOfDatasets (), leg, string ("l"), string ("canvasSR"));
  canvasCR = TCanvasCreator (PlotByDatasetsCR (), ListOfDatasets (), leg, string ("l"), string ("canvasCR"));

  ///////////////////////////
  //  PLots significance
  //////////////////////////
  legSign = new TLegend (0.6, 0.6, 0.88, 0.88);
  legSignRatio = new TLegend (0.6, 0.6, 0.88, 0.88);
  canvasSign = new TCanvas (TString ("Significance") + label);
  canvasSign->cd ();
  TH1F **histoSignRatio = new TH1F *[listOfDatasetsNP_.size ()];
  for (unsigned int l = 0; l < listOfDatasetsNP_.size (); l++) {
    if (l == 0)
      h_SignEstimated_[l]->Draw ();
    else
      h_SignEstimated_[l]->Draw ("same");
    string name = "Estimated - " + listOfDatasetsNP_[l].Name ();
    legSign->AddEntry (h_SignEstimated_[l], name.c_str (), "l");
    legSignRatio->AddEntry (h_SignEstimated_[l], listOfDatasetsNP_[l].Name ().c_str (), "l");
    h_SignExpected_[l]->Draw ("same");
    name = "Expected - " + listOfDatasetsNP_[l].Name ();
    legSign->AddEntry (h_SignExpected_[l], name.c_str (), "l");
    name = string (h_SignEstimated_[l]->GetName ()) + "_ratio";
    histoSignRatio[l] = PlotRatio (h_SignExpected_[l], h_SignEstimated_[l], name);
  }
  legSign->Draw ("same");
  canvasSignRatio = new TCanvas (TString ("SignificanceRatio") + label);
  canvasSignRatio->cd ();
  for (unsigned int l = 0; l < listOfDatasetsNP_.size (); l++) {
    histoSignRatio[l]->SetLineStyle (1);
    if (l == 0)
      histoSignRatio[l]->Draw ();
    histoSignRatio[l]->Draw ("same");
  }
  legSignRatio->Draw ("same");

}

void
TtJetEstimation::PrintResults (float METCut)
{
  cout << endl;
  cout << "***********************************************" << endl;
  cout << "*****      TtJetEstimation Results          ***" << endl;
  cout << "***********************************************" << endl;
  if (PlotByDatasetsSRTtJets ().size () > 0)
    cout << "Nof ttbar events in SR:\t " << PlotByDatasetsSRTtJets ()[0]->Integral () << endl;
  if (PlotByDatasetsCRTtJets ().size () > 0)
    cout << "Nof ttbar events in CR:\t " << PlotByDatasetsCRTtJets ()[0]->Integral () << endl;
  cout << "Nof events in SR:\t " << PlotSR ()->Integral () << endl;
  cout << "Nof events in CR:\t " << PlotCR ()->Integral () << endl;
  cout << "Nof events expected:\t " << PlotExpected ()->Integral () << endl;
  cout << "Nof events estimated:\t " << PlotEstimated ()->Integral () << endl;
  cout << "Nof events observed:\t " << PlotObserved ()->Integral () << endl;
  cout << "***********************************************" << endl;
  cout << "Shape Differences: " << endl;
  if (PlotByDatasetsCRTtJets ().size () > 0 && PlotByDatasetsSRTtJets ().size () > 0)
    cout << "Difference between ttbar in CR & in SR: " << Chi2Normalized (PlotByDatasetsCRTtJets ()[0], PlotByDatasetsSRTtJets ()[0], true) << endl;
  cout << "Difference between ttbar Estimated & Expected: " << Chi2Normalized (PlotEstimated (), PlotEstimated (), true) << endl;
  cout << "***********************************************" << endl;
  pair < float, float >signExpected, signEstimated;
  if (METCut < 0)
    METCut = NRegion_.second;
  signExpected = SignificanceExpected (METCut);
  signEstimated = SignificanceEstimated (METCut);
  cout << "Significance expected:  " << signExpected.first << " +/- " << signExpected.second << endl;
  cout << "Significance estimated: " << signEstimated.first << " +/- " << signEstimated.second << endl;
  cout << "***********************************************" << endl << endl;
}


void
TtJetEstimation::Write (TFile * fout, string label)
{
  fout->cd ();
  string dirname = "TtJetEstimation" + label;
  fout->mkdir (dirname.c_str ());
  fout->cd (dirname.c_str ());
  h_ObservableCR_->Write ();
  h_ObservableSR_->Write ();
  h_ObservableEstimated_->Write ();
  h_ObservableExpected_->Write ();
  h_ObservableObserved_->Write ();
  if (h_SysError_)
    h_SysError_->Write ();
  if (canvasRatioCR_SR)
    canvasRatioCR_SR->Write ();
  if (canvasRatioEstim_Exp)
    canvasRatioEstim_Exp->Write ();
  if (canvasSign)
    canvasSign->Write ();
  try {
    if (canvasSignRatio)
      canvasSignRatio->Write ();
  }
  catch (...) {;
  };
  if (canvasLMXContamination)
    canvasLMXContamination->Write ();
  //inclusive
  ((TDirectory *) (fout->Get (dirname.c_str ())))->mkdir ("Inclusive");
  char subdirName[100];
  sprintf (subdirName, "%s/Inclusive", dirname.c_str ());
  fout->cd (subdirName);
  for (unsigned int i = 0; i < h_ObservableByDatasetsTtJets_.size (); i++)
    h_ObservableByDatasetsTtJets_[i]->Write ();
  for (unsigned int i = 0; i < h_ObservableByDatasetsBkg_.size (); i++)
    h_ObservableByDatasetsBkg_[i]->Write ();
  for (unsigned int i = 0; i < h_ObservableByDatasetsNP_.size (); i++)
    h_ObservableByDatasetsNP_[i]->Write ();
  //CR
  ((TDirectory *) (fout->Get (dirname.c_str ())))->mkdir ("ControlRegion");
  sprintf (subdirName, "%s/ControlRegion", dirname.c_str ());
  fout->cd (subdirName);
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsTtJets_.size (); i++)
    h_ObservableCRByDatasetsTtJets_[i]->Write ();
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsBkg_.size (); i++)
    h_ObservableCRByDatasetsBkg_[i]->Write ();
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsNP_.size (); i++)
    h_ObservableCRByDatasetsNP_[i]->Write ();
  if (canvasCR)
    canvasCR->Write ();
  //SR
  ((TDirectory *) (fout->Get (dirname.c_str ())))->mkdir ("SignalRegion");
  sprintf (subdirName, "%s/SignalRegion", dirname.c_str ());
  fout->cd (subdirName);
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsTtJets_.size (); i++)
    h_ObservableSRByDatasetsTtJets_[i]->Write ();
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsBkg_.size (); i++)
    h_ObservableSRByDatasetsBkg_[i]->Write ();
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsNP_.size (); i++)
    h_ObservableSRByDatasetsNP_[i]->Write ();
  if (canvasSR)
    canvasSR->Write ();
}

  ////////////////////////////
  //      Compute           //
  ////////////////////////////

bool TtJetEstimation::Scale (float Lumi)
{
  if (scaled_) {
    cout << "Already scaled" << endl;
    return false;
  }
  else {
    for (unsigned int i = 0; i < listOfDatasetsTtJets_.size (); i++) {
      h_ObservableByDatasetsTtJets_[i]->Scale (listOfDatasetsTtJets_[i].NormFactor () * Lumi);
      h_ObservableCRByDatasetsTtJets_[i]->Scale (listOfDatasetsTtJets_[i].NormFactor () * Lumi);
      h_ObservableSRByDatasetsTtJets_[i]->Scale (listOfDatasetsTtJets_[i].NormFactor () * Lumi);
    }
    for (unsigned int i = 0; i < listOfDatasetsBkg_.size (); i++) {
      h_ObservableByDatasetsBkg_[i]->Scale (listOfDatasetsBkg_[i].NormFactor () * Lumi);
      h_ObservableCRByDatasetsBkg_[i]->Scale (listOfDatasetsBkg_[i].NormFactor () * Lumi);
      h_ObservableSRByDatasetsBkg_[i]->Scale (listOfDatasetsBkg_[i].NormFactor () * Lumi);
    }
    for (unsigned int i = 0; i < listOfDatasetsNP_.size (); i++) {
      h_ObservableByDatasetsNP_[i]->Scale (listOfDatasetsNP_[i].NormFactor () * Lumi);
      h_ObservableCRByDatasetsNP_[i]->Scale (listOfDatasetsNP_[i].NormFactor () * Lumi);
      h_ObservableSRByDatasetsNP_[i]->Scale (listOfDatasetsNP_[i].NormFactor () * Lumi);
    }
    scaled_ = true;

    return true;
  }
}

void
TtJetEstimation::ReScale (float factor)
{
  for (unsigned int i = 0; i < listOfDatasetsTtJets_.size (); i++) {
    h_ObservableByDatasetsTtJets_[i]->Scale (factor);
    h_ObservableCRByDatasetsTtJets_[i]->Scale (factor);
    h_ObservableSRByDatasetsTtJets_[i]->Scale (factor);
  }
  for (unsigned int i = 0; i < listOfDatasetsBkg_.size (); i++) {
    h_ObservableByDatasetsBkg_[i]->Scale (factor);
    h_ObservableCRByDatasetsBkg_[i]->Scale (factor);
    h_ObservableSRByDatasetsBkg_[i]->Scale (factor);
  }
  for (unsigned int i = 0; i < listOfDatasetsNP_.size (); i++) {
    h_ObservableByDatasetsNP_[i]->Scale (factor);
    h_ObservableCRByDatasetsNP_[i]->Scale (factor);
    h_ObservableSRByDatasetsNP_[i]->Scale (factor);
  }
  h_ObservableCR_->Scale (factor);
  h_ObservableSR_->Scale (factor);
  h_ObservableEstimated_->Scale (factor);
  h_ObservableExpected_->Scale (factor);
  h_ObservableObserved_->Scale (factor);
  if (h_SysError_)
    h_SysError_->Scale (factor);
}


void
TtJetEstimation::Compute (bool doTt, bool doBkg, bool doNP, int iNP)
{
  if(isMC_){
	ComputeObservableSR (doTt, doBkg, doNP, iNP);
	ComputeObservableCR (doTt, doBkg, doNP, iNP);
        ComputeObservableExpected(doTt, doBkg);  
	ComputeObservableObserved (doTt, doBkg, doNP, iNP);
  }
}

void
TtJetEstimation::ComputeObservableCR (bool doTt, bool doBkg, bool doNP, int iNP)
{
  if(!isMC_) return;
  for (int i = 1; i < h_ObservableCR_->GetNbinsX () + 1; i++) {
    float value = 0;
    float error = 0;
    if (doTt) {
      for (unsigned j = 0; j < listOfDatasetsTtJets_.size (); j++)
	value += h_ObservableCRByDatasetsTtJets_[j]->GetBinContent (i);
      for (unsigned j = 0; j < listOfDatasetsTtJets_.size (); j++)
	error += h_ObservableCRByDatasetsTtJets_[j]->GetBinError (i);
    }
    if (doBkg) {
      for (unsigned j = 0; j < listOfDatasetsBkg_.size (); j++)
	value += h_ObservableCRByDatasetsBkg_[j]->GetBinContent (i);
      for (unsigned j = 0; j < listOfDatasetsBkg_.size (); j++)
	error += h_ObservableCRByDatasetsBkg_[j]->GetBinError (i);
    }
    if (doNP) {
      if (iNP < (int) listOfDatasetsNP_.size ())
	value += h_ObservableCRByDatasetsNP_[iNP]->GetBinContent (i);
      if (iNP < (int) listOfDatasetsNP_.size ())
	error += h_ObservableCRByDatasetsNP_[iNP]->GetBinError (i);
    }
    h_ObservableCR_->SetBinContent (i, value);
    h_ObservableCR_->SetBinError (i, error);
  }
}

void
TtJetEstimation::ComputeObservableSR (bool doTt, bool doBkg, bool doNP, int iNP)
{
  if(!isMC_) return;
  for (int i = 1; i < h_ObservableSR_->GetNbinsX () + 1; i++) {
    float value = 0;
    float error = 0;
    if (doTt) {
      for (unsigned j = 0; j < listOfDatasetsTtJets_.size (); j++)
	value += h_ObservableSRByDatasetsTtJets_[j]->GetBinContent (i);
      for (unsigned j = 0; j < listOfDatasetsTtJets_.size (); j++)
	error += h_ObservableSRByDatasetsTtJets_[j]->GetBinError (i);
    }
    if (doBkg) {
      for (unsigned j = 0; j < listOfDatasetsBkg_.size (); j++)
	value += h_ObservableSRByDatasetsBkg_[j]->GetBinContent (i);
      for (unsigned j = 0; j < listOfDatasetsBkg_.size (); j++)
	error += h_ObservableSRByDatasetsBkg_[j]->GetBinError (i);
    }
    if (doNP) {
      if (iNP < (int) listOfDatasetsNP_.size ())
	value += h_ObservableSRByDatasetsNP_[iNP]->GetBinContent (i);
      if (iNP < (int) listOfDatasetsNP_.size ())
	error += h_ObservableSRByDatasetsNP_[iNP]->GetBinError (i);
    }
    h_ObservableSR_->SetBinContent (i, value);
    h_ObservableSR_->SetBinError (i, error);
  }
}

void
TtJetEstimation::ComputeObservableExpected (bool doTt, bool doBkg)
{
  if(!isMC_) return;
  for (int i = 1; i < h_ObservableSR_->GetNbinsX () + 1; i++) {
    float value = 0;
    float error = 0;
    if (doTt) {
      for (unsigned j = 0; j < listOfDatasetsTtJets_.size (); j++)
	value += h_ObservableSRByDatasetsTtJets_[j]->GetBinContent (i);
      for (unsigned j = 0; j < listOfDatasetsTtJets_.size (); j++)
	error += h_ObservableSRByDatasetsTtJets_[j]->GetBinError (i);
    }
    if (doBkg) {
      for (unsigned j = 0; j < listOfDatasetsBkg_.size (); j++)
	value += h_ObservableSRByDatasetsBkg_[j]->GetBinContent (i);
      for (unsigned j = 0; j < listOfDatasetsBkg_.size (); j++)
	error += h_ObservableSRByDatasetsBkg_[j]->GetBinError (i);
    }
    h_ObservableExpected_->SetBinContent (i, value);
    h_ObservableExpected_->SetBinError (i, error);
  }
}

void
TtJetEstimation::ComputeObservableObserved (bool doTt, bool doBkg, bool doNP, int iNP)
{
  if(!isMC_) return;
  for (int i = 1; i < h_ObservableSR_->GetNbinsX () + 1; i++) {
    float value = 0;
    float error = 0;
    if (doTt) {
      for (unsigned j = 0; j < listOfDatasetsTtJets_.size (); j++)
	value += h_ObservableSRByDatasetsTtJets_[j]->GetBinContent (i);
      for (unsigned j = 0; j < listOfDatasetsTtJets_.size (); j++)
	error += h_ObservableSRByDatasetsTtJets_[j]->GetBinError (i);
    }
    if (doBkg) {
      for (unsigned j = 0; j < listOfDatasetsBkg_.size (); j++)
	value += h_ObservableSRByDatasetsBkg_[j]->GetBinContent (i);
      for (unsigned j = 0; j < listOfDatasetsBkg_.size (); j++)
	error += h_ObservableSRByDatasetsBkg_[j]->GetBinError (i);
    }
    if (doNP) {
      if (iNP < (int) listOfDatasetsNP_.size ())
	value += h_ObservableSRByDatasetsNP_[iNP]->GetBinContent (i);
      if (iNP < (int) listOfDatasetsNP_.size ())
	error += h_ObservableSRByDatasetsNP_[iNP]->GetBinError (i);
    }
    h_ObservableObserved_->SetBinContent (i, value);
    h_ObservableObserved_->SetBinError (i, error);
  }
}

bool TtJetEstimation::ComputeNormalisation (TH1F * hQCDEstim, TH1F * hWJetEstim)
{
  if (normalized_) {
    cout << "Already normalized" << endl;
    //return false;
  }
  //else{
  //h_ObservableEstimated_ = (TH1F*) h_ObservableCR_->Clone("h_ObservableEstimated");             
  //h_ObservableEstimated_ = (TH1F*) h_ObservableCR_->Clone();            
  CopyHisto (h_ObservableCR_, h_ObservableEstimated_);
  float
    n1 = NofEvents (h_ObservableEstimated_, NRegion_).first;
  float
    n2 = NofEvents (h_ObservableSR_, NRegion_).first;
  float
    nQCD = 0;
  float
    nWJet = 0;
  if (hQCDEstim)
    nQCD = NofEvents (hQCDEstim, NRegion_).first;
  if (hWJetEstim)
    nWJet = NofEvents (hWJetEstim, NRegion_).first;
  float
    factor = (n2 - nQCD - nWJet) / n1;
  h_ObservableEstimated_->Scale (factor);
  return true;
  //}
}


void
TtJetEstimation::SetPlotByDatasets (vector < TFile * >files, string histoname, vector < TH1F * >histos, vector < Dataset > listRef)
{
  for (unsigned int i = 0; i < files.size (); i++) {
    if (files[i]->Get (histoname.c_str ()) == NULL)
      cout << "Histo " << histoname << " not found in file " << files[i]->GetName () << endl;
    histos.push_back ((TH1F *) files[i]->Get (histoname.c_str ()));
  }
  if (files.size () != listRef.size ())
    cout << "Inconsistency in SetPlotByDatasets" << endl;
}

pair < float, float >
TtJetEstimation::NofEvents (TH1F * histo, pair < float, float >edges) const
{
  int imax = histo->GetNbinsX ();
  int icut1 = 0;
  int icut2 = 0;
  for (int i = 1; i < imax + 1; i++) {
    if (edges.first >= histo->GetBinLowEdge (i) && edges.first < histo->GetBinLowEdge (i + 1)) {
      icut1 = i;
      break;
    }
  }
  for (int i = icut1; i < imax + 1; i++) {
    if (edges.second >= histo->GetBinLowEdge (i) && edges.second < histo->GetBinLowEdge (i + 1)) {
      icut2 = i;
      break;
    }
  }
  float value = histo->Integral (icut1, icut2);
  float error = 0;
  for (int i = icut1; i <= icut2; i++) {
    error += histo->GetBinError (i);
  }
  pair < float, float >out (value, error);
  return out;
}

pair < float, float >
TtJetEstimation::NofEventsObserved (float VarCut)  const
{
  int imax = h_ObservableObserved_->GetNbinsX ();
  int icut = 0;
  for (int i = 1; i < imax + 1; i++) {
    if (VarCut >= h_ObservableObserved_->GetBinLowEdge (i) && VarCut < h_ObservableObserved_->GetBinLowEdge (i + 1)) {
      icut = i;
      break;
    }
  }
  float value = h_ObservableObserved_->Integral (icut, imax);
  float error = 0;
  for (int i = icut; i < imax; i++) {
    error += h_ObservableObserved_->GetBinError (i);
  }
  pair < float, float >out (value, error);
  return (out);
}

pair < float, float >
TtJetEstimation::NofEventsEstimated (float VarCut) const
{
  int imax = h_ObservableEstimated_->GetNbinsX ();
  int icut = 0;
  for (int i = 1; i < imax + 1; i++) {
    if (VarCut >= h_ObservableEstimated_->GetBinLowEdge (i) && VarCut < h_ObservableEstimated_->GetBinLowEdge (i + 1)) {
      icut = i;
      break;
    }
  }
  float nofEvtsCR = h_ObservableCR_->Integral (icut, imax);
  float nofEvtsNorm = NofEvents (h_ObservableSR_, NRegion_).first;
  float value = h_ObservableEstimated_->Integral (icut, imax);
  float error = 0;
  error = sqrt (1 / nofEvtsCR + 1 / nofEvtsNorm) * value;
  /*
     for(int i=icut;i<imax;i++){
     error+=h_ObservableEstimated_->GetBinError(i);
     }
   */
  pair < float, float >out (value, error);
  return (out);
}


pair < float, float >
TtJetEstimation::NofEventsExpected (float VarCut) const
{
  int imax = h_ObservableExpected_->GetNbinsX ();
  int icut = 0;
  for (int i = 1; i < imax + 1; i++) {
    if (VarCut >= h_ObservableExpected_->GetBinLowEdge (i) && VarCut < h_ObservableExpected_->GetBinLowEdge (i + 1)) {
      icut = i;
      break;
    }
  }
  float value = h_ObservableExpected_->Integral (icut, imax);
  float error = 0;
  for (int i = icut; i < imax; i++) {
    error += h_ObservableExpected_->GetBinError (i);
  }
  pair < float, float >out (value, error);
  return out;
}

pair < float, float >
TtJetEstimation::SignificanceExpected (float VarCut) const
{
  float data = NofEventsObserved (VarCut).first;
  float dataError = NofEventsObserved (VarCut).second;
  float expected = NofEventsExpected (VarCut).first;
  float expectedError = NofEventsExpected (VarCut).second;
  float significance = -1;
  if (data > 0)
    significance = (data - expected) / sqrt (data);
  float error = -1;
  if (data > 0 && expected > 0)
    error = significance * sqrt ((dataError / data) * (dataError / data) + (expectedError / expected) * (expectedError / expected));
  pair < float, float >out (significance, error);
  return out;
}

pair < float, float >
TtJetEstimation::SignificanceEstimated (float VarCut) const
{
  float data = NofEventsObserved (VarCut).first;
  float dataError = sqrt (data);	//sqrt(n) and not .second;
  float estimated = NofEventsEstimated (VarCut).first;
  float estimatedError = NofEventsEstimated (VarCut).second;
  float significance = -1;
  if (data > 0)
    significance = (data - estimated) / sqrt (data);
  float error = -1;
  if (data > 0 && estimated > 0)
    error = significance * sqrt ((dataError / data) * (dataError / data) + (estimatedError / estimated) * (estimatedError / estimated));
  pair < float, float >out (significance, error);
  return out;
}

void
TtJetEstimation::SetListOfDatasets (vector < Dataset > listOfDatasetsBkg, vector < Dataset > listOfDatasetsTtJets, vector < Dataset > listOfDatasetsNP)
{
  listOfDatasetsTtJets_ = listOfDatasetsTtJets;
  listOfDatasetsNP_ = listOfDatasetsNP;
  listOfDatasetsBkg_ = listOfDatasetsBkg;
  for (unsigned int i = 0; i < listOfDatasetsTtJets_.size (); i++) {
    TH1F *f = NULL;
    TH1F *g = NULL;
    TH1F *h = NULL;
    h_ObservableCRByDatasetsTtJets_.push_back (h);
    h_ObservableSRByDatasetsTtJets_.push_back (g);
    h_ObservableByDatasetsTtJets_.push_back (f);
  }
  for (unsigned int i = 0; i < listOfDatasetsBkg_.size (); i++) {
    TH1F *f = NULL;
    TH1F *g = NULL;
    TH1F *h = NULL;
    h_ObservableCRByDatasetsBkg_.push_back (h);
    h_ObservableSRByDatasetsBkg_.push_back (g);
    h_ObservableByDatasetsBkg_.push_back (f);
  }
  for (unsigned int i = 0; i < listOfDatasetsNP_.size (); i++) {
    TH1F *f = NULL;
    TH1F *g = NULL;
    TH1F *h = NULL;
    TH1F *k = NULL;
    TH1F *l = NULL;
    h_ObservableCRByDatasetsNP_.push_back (h);
    h_ObservableSRByDatasetsNP_.push_back (g);
    h_ObservableByDatasetsNP_.push_back (f);
    h_SignEstimated_.push_back(k);
    h_SignExpected_.push_back(l);
  }
  //order Bkg - TTJets - SUSY
  for (unsigned int i = 0; i < listOfDatasetsBkg.size (); i++)
    listOfDatasets_.push_back (listOfDatasetsBkg[i]);
  for (unsigned int i = 0; i < listOfDatasetsTtJets.size (); i++)
    listOfDatasets_.push_back (listOfDatasetsTtJets[i]);
  for (unsigned int i = 0; i < listOfDatasetsNP.size (); i++)
    listOfDatasets_.push_back (listOfDatasetsNP[i]);
}

vector < TH1F * >TtJetEstimation::PlotByDatasets ()
{
  //order Bkg - TTJets - SUSY
  vector < TH1F * >vec;
  for (unsigned int i = 0; i < h_ObservableByDatasetsBkg_.size (); i++)
    vec.push_back (h_ObservableByDatasetsBkg_[i]);
  for (unsigned int i = 0; i < h_ObservableByDatasetsTtJets_.size (); i++)
    vec.push_back (h_ObservableByDatasetsTtJets_[i]);
  for (unsigned int i = 0; i < h_ObservableByDatasetsNP_.size (); i++)
    vec.push_back (h_ObservableByDatasetsNP_[i]);
  return vec;
}

vector < TH1F * >TtJetEstimation::PlotByDatasetsCR ()
{
  //order Bkg - TTJets - SUSY
  vector < TH1F * >vec;
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsBkg_.size (); i++)
    vec.push_back (h_ObservableCRByDatasetsBkg_[i]);
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsTtJets_.size (); i++)
    vec.push_back (h_ObservableCRByDatasetsTtJets_[i]);
  for (unsigned int i = 0; i < h_ObservableCRByDatasetsNP_.size (); i++)
    vec.push_back (h_ObservableCRByDatasetsNP_[i]);
  return vec;
}

vector < TH1F * >TtJetEstimation::PlotByDatasetsSR ()
{
  //order Bkg - TTJets - SUSY
  vector < TH1F * >vec;
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsBkg_.size (); i++)
    vec.push_back (h_ObservableSRByDatasetsBkg_[i]);
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsTtJets_.size (); i++)
    vec.push_back (h_ObservableSRByDatasetsTtJets_[i]);
  for (unsigned int i = 0; i < h_ObservableSRByDatasetsNP_.size (); i++)
    vec.push_back (h_ObservableSRByDatasetsNP_[i]);
  return vec;
}

vector < Dataset > TtJetEstimation::ListOfDatasets () const
{
  //order Bkg - TTJets - SUSY
  vector < Dataset > vec;
  for (unsigned int i = 0; i < listOfDatasetsBkg_.size (); i++)
    vec.push_back (listOfDatasetsBkg_[i]);
  for (unsigned int i = 0; i < listOfDatasetsTtJets_.size (); i++)
    vec.push_back (listOfDatasetsTtJets_[i]);
  for (unsigned int i = 0; i < listOfDatasetsNP_.size (); i++)
    vec.push_back (listOfDatasetsNP_[i]);
  return vec;
}

vector < string > TtJetEstimation::ListOfDatasetsName () const
{
  vector < string > a;
  for (unsigned int i = 0; i < ListOfDatasets ().size (); i++)
    a.push_back (ListOfDatasets ()[i].Name ());
  return a;
}

vector < int >
TtJetEstimation::ListOfDatasetsColor () const
{
  vector < int >a;
  for (unsigned int i = 0; i < ListOfDatasets ().size (); i++)
    a.push_back (ListOfDatasets ()[i].Color ());
  return a;
}

vector < float >
TtJetEstimation::ListOfDatasetsNormFactor () const
{
  vector < float >a;
  for (unsigned int i = 0; i < ListOfDatasets ().size (); i++)
    a.push_back (ListOfDatasets ()[i].NormFactor ());
  return a;
}


void
TtJetEstimation::ComputePlotSignificance (TH1F * qcdShapeEstim, TH1F * wjetShapeEstim)
{
  pair < float, float >NR = NRegion_;
  for (int i = 1; i < NbinsX_; i++) {
    ComputeObservableExpected (true, true);
    DefineNR (pair < float, float >(0, h_ObservableSR_->GetBinLowEdge(i+1)));
    for (unsigned int l = 0; l < listOfDatasetsNP_.size (); l++) {
      //recompute
      ComputeObservableSR (true, true, true, l);	// add this LMX
      ComputeObservableCR (true, true, true, l);	// add this LMX
      ComputeObservableObserved (true, true, true, l);
      DefineNR (pair < float, float >(0, h_ObservableSR_->GetBinLowEdge(i+1)));
      ComputeNormalisation (qcdShapeEstim, wjetShapeEstim);
      //
      pair < float, float >sign;
      //Estimated
      sign = SignificanceEstimated (h_ObservableSR_->GetBinLowEdge(i+1));
      h_SignEstimated_[l]->SetBinContent (i + 1, sign.first);
      h_SignEstimated_[l]->SetBinError (i + 1, sign.second);
      //Expected
      sign = SignificanceExpected (h_ObservableSR_->GetBinLowEdge(i+1));
      h_SignExpected_[l]->SetBinContent (i + 1, sign.first);
      //h_SignExpected_[l]->SetBinError(i,sign.second);
      h_SignExpected_[l]->SetBinError (i + 1, 0.);
    }
  }
  NRegion_ = NR;

}

void
TtJetEstimation::ComputeEstimationLMXContamination (TH1F * qcdShapeEstim, TH1F * wjetShapeEstim)
{
  ComputeObservableExpected (true, true);
  canvasLMXContamination = new TCanvas ("canvasLMXContamination");
  canvasLMXContamination->cd ();
  //PlotExpected()->SetLineColor(1);
  TH1F *h0 = (TH1F *) PlotExpected ()->Clone ();
  h0->Scale (1 / h0->Integral ());
  h0->Draw ();
  leg_LMXContamination = new TLegend (0.7, 0.7, 0.88, 0.88);
  leg_LMXContamination->AddEntry (h0, "Expected", "l");
  for (unsigned int l = 0; l < listOfDatasetsNP_.size (); l++) {
    //recompute
    ComputeObservableSR (true, true, true, l);	// add this LMX
    ComputeObservableCR (true, true, true, l);	// add this LMX
    ComputeObservableObserved (true, true, true, l);
    cout << "Inte " << PlotEstimated ()->Integral (1, 4);
    ComputeNormalisation (qcdShapeEstim, wjetShapeEstim);
    cout << "Normalisation: " << NRegion_.first << " " << NRegion_.second << endl;
    cout << "Inte " << PlotEstimated ()->Integral (1, 4);
    //
    char name[100];
    sprintf (name, "Exp_%d", l);
    TH1F *h = (TH1F *) PlotEstimated ()->Clone (name);
    h->Scale (1 / h->Integral ());
    h->SetLineColor (listOfDatasetsNP_[l].Color ());
    h->Draw ("same");
    leg_LMXContamination->AddEntry (h, listOfDatasetsNP_[l].Name ().c_str (), "l");
  }
  leg_LMXContamination->Draw ("same");
}



void
TtJetEstimation::AddNormalisationError (float factor)
{
  for (int i = 1; i < h_ObservableEstimated_->GetNbinsX () + 1; i++) {
    float error = 0;
    error = sqrt ((h_ObservableEstimated_->GetBinError (i) * h_ObservableEstimated_->GetBinError (i)) + (h_ObservableEstimated_->GetBinContent (i) * h_ObservableEstimated_->GetBinContent (i)));
    h_ObservableEstimated_->SetBinError (i, error);
  }
}


void
TtJetEstimation::ComputeSystematicError (bool add)
{
  TH1F *h_ref = 0;
  if (h_ObservableSRByDatasetsTtJets_.size () > 0)
    h_ref = (TH1F *) h_ObservableSRByDatasetsTtJets_[0]->Clone ("h_ref");
  if (h_ObservableSRByDatasetsTtJets_.size () >= 1) {
    for (unsigned int i = 1; i < h_ObservableSRByDatasetsTtJets_.size (); i++) {
      h_ref->Add (h_ObservableSRByDatasetsTtJets_[i]);
    }
  }
  for (int i = 1; i < h_ObservableEstimated_->GetNbinsX () + 1; i++) {
    float error = 0;
    float diff = 0;
    diff = h_ObservableEstimated_->GetBinContent (i) - h_ref->GetBinContent (i);
    error = sqrt ((h_ObservableEstimated_->GetBinError (i) * h_ObservableEstimated_->GetBinError (i)) + (diff * diff));
    if (add)
      h_ObservableEstimated_->SetBinError (i, error);
    h_SysError_->SetBinContent (i, diff);
  }
  delete h_ref;
}

void
TtJetEstimation::AddSystematicError (TH1F * SystError)
{
  for (int i = 1; i < h_ObservableEstimated_->GetNbinsX () + 1; i++) {
    float error = 0;
    error = sqrt ((h_ObservableEstimated_->GetBinError (i) * h_ObservableEstimated_->GetBinError (i)) + (h_SysError_->GetBinContent (i) * h_SysError_->GetBinContent (i)));
    h_ObservableEstimated_->SetBinError (i, error);
  }
}
