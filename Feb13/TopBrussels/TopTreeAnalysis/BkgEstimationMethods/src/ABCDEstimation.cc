#include "../interface/ABCDEstimation.h"

ABCDEstimation::ABCDEstimation()
{
/**
//This is my constructor
*/
  Name_ = "";
  XaxisTitle_ = "";
  YaxisTitle_ = "";
  //histo
  h1DQCD_ = 0;
  h2DQCD_ = 0;
  h3DQCD_ = 0 ;
  h2D_ = 0;
  h3D_ = 0 ;
  hXProfQCD_ = 0;
  hYProfQCD_ = 0;
  /* cXProfQCD_ = new TCanvas("cXProfQCD"); */
  /* cYProfQCD_ = new TCanvas("cYProfQCD"); */
  hXProf_ = 0;
  hYProf_ = 0;
  /* cXProf_ = new TCanvas("cXProf"); */
  /* cYProf_ = new TCanvas("cYProf"); */
  hShapePred_ = 0;
  hShapePredQCD_ = 0;
  hEstim_vs_Xcut = 0;
  hEstim_vs_Ycut = 0;
  //for drawing
  func = 0;
  line = 0;
  plotRatioEstim_Exp = 0;
  plotRatioQCDEstim_Exp = 0;
  /* canvasRatioEstim_Exp = 0; */
  /* canvasRatioQCDEstim_Exp = 0; */
  hErrorVsLuminosity_ = 0;
  hErrorVsLuminosityQCD_ = 0;
  /* cErrorVsLuminosity_ = new TCanvas("cErrorVsLuminosity"); */
  /* cErrorVsLuminosityQCD_ = new TCanvas("cErrorVsLuminosityQCD"); */
  Lumi_ = 100;
  //tab

  region_ = -1;

  estimate_ = new float*[4];
  estimateQCD_ = new float*[4];
  for(int i=0;i<4;i++)
  {
  	estimate_[i] = new float[4];
  	estimateQCD_[i] = new float[4];
	for(int j=0;j<4;j++)
	{
		estimate_[i][j] = -1;
		estimateQCD_[i][j] = -1;
	}
  }
  estimateZ_ = 0;
  estimateQCDZ_ = 0;
}

ABCDEstimation::ABCDEstimation(TString hName, TString hTitle, TString XaxisTitle, TString YaxisTitle, int NXbins, float Xmin, float Xmax, int NYbins, float Ymin, float Ymax)
{
//
// A test class. A more elaborate class description.
//
  
  Name_ = hName;
  XaxisTitle_ = XaxisTitle;
  YaxisTitle_ = YaxisTitle;

  h1DQCD_ = 0;
  h2DQCD_ = new TH2F(hName+"_QCD",hTitle,NXbins,Xmin,Xmax,NYbins,Ymin,Ymax);
  h3DQCD_ = 0 ;
  hShapePredQCD_ = 0;
  hXProfQCD_ = new TGraphAsymmErrors(NXbins); hXProfQCD_->SetNameTitle(hName+"_XProfQCD","Profile along the X axis");
  hYProfQCD_ = new TGraphAsymmErrors(NYbins); hYProfQCD_->SetNameTitle(hName+"_YProfQCD","Profile along the Y axis");
  /* cXProfQCD_ = new TCanvas("cXProfQCD"); */
  /* cYProfQCD_ = new TCanvas("cYProfQCD"); */

  h2D_ = new TH2F(hName,hTitle,NXbins,Xmin,Xmax,NYbins,Ymin,Ymax);
  h3D_ = 0 ;
  hShapePred_ = 0;
  hXProf_ = new TGraphAsymmErrors(NXbins); hXProf_->SetNameTitle(hName+"_XProf","Profile along the X axis");
  hYProf_ = new TGraphAsymmErrors(NYbins); hYProf_->SetNameTitle(hName+"_YProf","Profile along the Y axis");
  /* cXProf_ = new TCanvas("cXProf"); */
  /* cYProf_ = new TCanvas("cYProf"); */

  hEstim_vs_Xcut = new TH1F(Name_+"Estim_vs_XCut","",NXbins,Xmin,Xmax);
  hEstim_vs_Ycut = new TH1F(Name_+"Estim_vs_YCut","",NYbins,Ymin,Ymax);

  region_ = -1;

  estimate_ = new float*[4];
  estimateQCD_ = new float*[4];
  for(int i=0;i<4;i++)
  {
  	estimate_[i] = new float[4];
  	estimateQCD_[i] = new float[4];
	for(int j=0;j<4;j++)
	{
		estimate_[i][j] = -1;
		estimateQCD_[i][j] = -1;
	}
  }
  estimateZ_ = 0;
  estimateQCDZ_ = 0;
  
  func = 0;
  line = 0;
  plotRatioEstim_Exp = 0;
  plotRatioQCDEstim_Exp = 0;
  /* canvasRatioEstim_Exp = 0; */
  /* canvasRatioQCDEstim_Exp = 0; */
  hErrorVsLuminosity_ = 0;
  hErrorVsLuminosityQCD_ = 0;
  /* cErrorVsLuminosity_ = new TCanvas("cErrorVsLuminosity"); */
  /* cErrorVsLuminosityQCD_ = new TCanvas("cErrorVsLuminosityQCD"); */
  Lumi_ = 100;
}

ABCDEstimation::ABCDEstimation(TString hName, TString hTitle, TString XaxisTitle, TString YaxisTitle, int NXbins, float Xmin, float Xmax, int NYbins, float Ymin, float Ymax, int NZbins, float Zmin, float Zmax)
{
  Name_ = hName;
  XaxisTitle_ = XaxisTitle;
  YaxisTitle_ = YaxisTitle;

  h1DQCD_ = new TH1F("ShapeExp_QCD", "Shape expectation", NZbins, Zmin, Zmax);
  h2DQCD_ = new TH2F(hName+"_2D_QCD",hTitle,NXbins,Xmin,Xmax,NYbins,Ymin,Ymax);
  h3DQCD_ = new TH3F(hName+"_3D_QCD",hTitle,NXbins,Xmin,Xmax,NYbins,Ymin,Ymax,NZbins,Zmin,Zmax) ;
  hShapePredQCD_ = new TH1F("ShapePrd_QCD", "Shape estimate", NZbins, Zmin, Zmax);
  hXProfQCD_ = new TGraphAsymmErrors(NXbins); hXProfQCD_->SetNameTitle(hName+"_XProfQCD","Profile along the X axis");
  hYProfQCD_ = new TGraphAsymmErrors(NYbins); hYProfQCD_->SetNameTitle(hName+"_YProfQCD","Profile along the Y axis");
  /* cXProfQCD_ = new TCanvas("cXProfQCD"); */
  /* cYProfQCD_ = new TCanvas("cYProfQCD"); */

  h2D_ = new TH2F(hName+"_2D",hTitle,NXbins,Xmin,Xmax,NYbins,Ymin,Ymax);
  h3D_ = new TH3F(hName+"_3D",hTitle,NXbins,Xmin,Xmax,NYbins,Ymin,Ymax,NZbins,Zmin,Zmax) ;
  hShapePred_ = new TH1F("ShapePrd", "Shape estimate", NZbins, Zmin, Zmax);
  hXProf_ = new TGraphAsymmErrors(NXbins); hXProf_->SetNameTitle(hName+"_XProf","Profile along the X axis");
  hYProf_ = new TGraphAsymmErrors(NYbins); hYProf_->SetNameTitle(hName+"_YProf","Profile along the Y axis");
  /* cXProf_ = new TCanvas("cXProf"); */
  /* cYProf_ = new TCanvas("cYProf"); */

  hEstim_vs_Xcut = new TH1F(hName+"Estim_vs_XCut","",NXbins,Xmin,Xmax);
  hEstim_vs_Ycut = new TH1F(hName+"Estim_vs_YCut","",NYbins,Ymin,Ymax);

  region_ = -1;

  estimate_ = new float*[4];
  estimateQCD_ = new float*[4];
  for(int i=0;i<4;i++)
  {
  	estimate_[i] = new float[4];
  	estimateQCD_[i] = new float[4];
	for(int j=0;j<4;j++)
	{
		estimate_[i][j] = -1;
		estimateQCD_[i][j] = -1;
	}
  }
  estimateZ_ = new float*[NZbins];
  estimateQCDZ_ = new float*[NZbins];
  for(int i=0;i<NZbins;i++)
  {
  	estimateZ_[i] = new float[4];
  	estimateQCDZ_[i] = new float[4];
	for(int j=0;j<4;j++)
	{
		estimateZ_[i][j] = -1;
		estimateQCDZ_[i][j] = -1;
	}
  }
  func = 0;
  line = 0;
  plotRatioEstim_Exp = 0;
  plotRatioQCDEstim_Exp = 0;
  /* canvasRatioEstim_Exp = 0; */
  /* canvasRatioQCDEstim_Exp = 0; */
  hErrorVsLuminosity_ = 0;
  hErrorVsLuminosityQCD_ = 0;
  /* cErrorVsLuminosity_ = new TCanvas("cErrorVsLuminosity"); */
  /* cErrorVsLuminosityQCD_ = new TCanvas("cErrorVsLuminosityQCD"); */
  Lumi_ = 100;
}

ABCDEstimation::ABCDEstimation(TString hName, TString hTitle, TString XaxisTitle, TString YaxisTitle, int NXbins, float *binsX, int NYbins, float *binsY, int NZbins, float *binsZ)
{
  Name_ = hName;
  XaxisTitle_ = XaxisTitle;
  YaxisTitle_ = YaxisTitle;
  h1DQCD_ = new TH1F("ShapeExp", "Shape expectation", NZbins, binsZ);
  h2DQCD_ = new TH2F(hName+"_2D_QCD",hTitle,NXbins,binsX,NYbins,binsY);
  h3DQCD_ = new TH3F(hName+"_3D_QCD",hTitle,NXbins,binsX,NYbins,binsY,NZbins,binsZ) ;
  hShapePredQCD_ = new TH1F("ShapePrd_QCD", "Shape estimate", NZbins, binsZ);
  hXProfQCD_ = new TGraphAsymmErrors(NXbins); hXProfQCD_->SetNameTitle(hName+"_XProfQCD","Profile along the X axis");
  hYProfQCD_ = new TGraphAsymmErrors(NYbins); hYProfQCD_->SetNameTitle(hName+"_YProfQCD","Profile along the Y axis");
  /* cXProfQCD_ = new TCanvas("cXProfQCD"); */
  /* cYProfQCD_ = new TCanvas("cYProfQCD"); */

  h2D_ = new TH2F(hName+"_2D",hTitle,NXbins,binsX,NYbins,binsY);
  h3D_ = new TH3F(hName+"_3D",hTitle,NXbins,binsX,NYbins,binsY,NZbins,binsZ) ;
  hShapePred_ = new TH1F("ShapePrd", "Shape estimate", NZbins, binsZ);
  hXProf_ = new TGraphAsymmErrors(NXbins); hXProf_->SetNameTitle(hName+"_XProf","Profile along the X axis");
  hYProf_ = new TGraphAsymmErrors(NYbins); hYProf_->SetNameTitle(hName+"_YProf","Profile along the Y axis");
  /* cXProf_ = new TCanvas("cXProf"); */
  /* cYProf_ = new TCanvas("cYProf"); */

  hEstim_vs_Xcut = new TH1F(hName+"Estim_vs_XCut","",NXbins,binsX);
  hEstim_vs_Ycut = new TH1F(hName+"Estim_vs_YCut","",NYbins,binsY);
  
  region_ = -1;

  estimate_ = new float*[4];
  estimateQCD_ = new float*[4];
  for(int i=0;i<4;i++)
  {
  	estimate_[i] = new float[4];
  	estimateQCD_[i] = new float[4];
	for(int j=0;j<4;j++)
	{
		estimate_[i][j] = -1;
		estimateQCD_[i][j] = -1;
	}
  }
  estimateZ_ = new float*[NZbins];
  estimateQCDZ_ = new float*[NZbins];
  for(int i=0;i<NZbins;i++)
  {
  	estimateZ_[i] = new float[4];
  	estimateQCDZ_[i] = new float[4];
	for(int j=0;j<4;j++)
	{
		estimateZ_[i][j] = -1;
		estimateQCDZ_[i][j] = -1;
	}
  }
  func = 0;
  line = 0;
  plotRatioEstim_Exp = 0;
  plotRatioQCDEstim_Exp = 0;
  /* canvasRatioEstim_Exp = 0; */
  /* canvasRatioQCDEstim_Exp = 0; */
  hErrorVsLuminosity_ = 0;
  hErrorVsLuminosityQCD_ = 0;
  /* cErrorVsLuminosity_ = new TCanvas("cErrorVsLuminosity"); */
  /* cErrorVsLuminosityQCD_ = new TCanvas("cErrorVsLuminosityQCD"); */
  Lumi_ = 100;
}

ABCDEstimation::ABCDEstimation( TH2F* histo, TH2F* histo_QCD){
  Name_ = histo->GetName();
  XaxisTitle_ = histo->GetXaxis()->GetTitle();
  YaxisTitle_ = histo->GetYaxis()->GetTitle();

  h1DQCD_ = 0;
  if(histo_QCD){
  	h2DQCD_ = histo_QCD;
  	h2DQCD_->SetName(Name_+"_2D_QCD");
  }
  h3DQCD_ = 0;
  h2D_ = histo;
  h2D_->SetName(Name_+"_2D");
  h3D_ = 0;
  hShapePred_ = 0;
  hShapePredQCD_ = 0;
  hXProfQCD_ = 0;
  hYProfQCD_ = 0;
  /* cXProfQCD_ = new TCanvas("cXProfQCD"); */
  /* cYProfQCD_ = new TCanvas("cYProfQCD"); */

  hXProf_ = 0;
  hYProf_ = 0;
  /* cXProf_ = new TCanvas("cXProf"); */
  /* cYProf_ = new TCanvas("cYProf"); */

  const TArrayD* Xbins = h2D_->GetXaxis()->GetXbins();
  const TArrayD* Ybins = h2D_->GetYaxis()->GetXbins();
  double* XbinsTab = new double[Xbins->GetSize()];
  double* YbinsTab = new double[Ybins->GetSize()];
  for(int i=0;i<Xbins->GetSize();i++){
  	XbinsTab[i] = Xbins->At(i);
  }
  for(int i=0;i<Ybins->GetSize();i++){
  	YbinsTab[i] = Ybins->At(i);
  }
  hEstim_vs_Xcut = new TH1F(Name_+"Estim_vs_XCut","",h2D_->GetNbinsX(),XbinsTab);
  hEstim_vs_Ycut = new TH1F(Name_+"Estim_vs_YCut","",h2D_->GetNbinsY(),YbinsTab);

  region_ = -1;

  estimate_ = new float*[4];
  estimateQCD_ = new float*[4];
  for(int i=0;i<4;i++)
  {
  	estimate_[i] = new float[4];
  	estimateQCD_[i] = new float[4];
	for(int j=0;j<4;j++)
	{
		estimate_[i][j] = -1;
		estimateQCD_[i][j] = -1;
	}
  }
  estimateZ_ = 0;
  estimateQCDZ_ = 0;
  
  func = 0;
  line = 0;
  plotRatioEstim_Exp = 0;
  plotRatioQCDEstim_Exp = 0;
  /* canvasRatioEstim_Exp = 0; */
  /* canvasRatioQCDEstim_Exp = 0; */
  hErrorVsLuminosity_ = 0;
  hErrorVsLuminosityQCD_ = 0;
  /* cErrorVsLuminosity_ = new TCanvas("cErrorVsLuminosity"); */
  /* cErrorVsLuminosityQCD_ = new TCanvas("cErrorVsLuminosityQCD"); */
  Lumi_ = 100;
}

ABCDEstimation::ABCDEstimation( TH3F* histo, TH3F* histo_QCD){
  Name_ = histo->GetName();
  XaxisTitle_ = histo->GetXaxis()->GetTitle();
  YaxisTitle_ = histo->GetYaxis()->GetTitle();

  h1DQCD_ = new TH1F("ShapeExp", "Shape expectation", histo->GetNbinsZ()+1, (float*)histo->GetZaxis()->GetXbins());
  if(histo_QCD){
  	h2DQCD_ = (TH2F*)histo_QCD->Project3D("xy");
  	h2DQCD_->SetName(Name_+"_2D_QCD");
  	h3DQCD_ = histo_QCD;
  	h3DQCD_->SetName(Name_+"_3D_QCD");
  	hShapePredQCD_ = new TH1F("ShapePrd_QCD", "Shape estimate", histo_QCD->GetNbinsZ()+1, (float*)histo_QCD->GetZaxis()->GetXbins());
  }
  h2D_ = (TH2F*)histo->Project3D("xy");
  h2D_->SetName(Name_+"_2D");
  h3D_ = histo;
  h3D_->SetName(Name_+"_3D");
  hShapePred_ = new TH1F("ShapePrd", "Shape estimate", histo->GetNbinsZ()+1, (float*)histo->GetZaxis()->GetXbins());
  hXProfQCD_ = 0;
  hYProfQCD_ = 0;
  hXProf_ = 0;
  hYProf_ = 0;

  region_ = -1;

  estimate_ = new float*[4];
  estimateQCD_ = new float*[4];
  for(int i=0;i<4;i++)
  {
  	estimate_[i] = new float[4];
  	estimateQCD_[i] = new float[4];
	for(int j=0;j<4;j++)
	{
		estimate_[i][j] = -1;
		estimateQCD_[i][j] = -1;
	}
  }
  estimateZ_ = new float*[histo->GetNbinsZ()];
  estimateQCDZ_ = new float*[histo->GetNbinsZ()];
  for(int i=0;i<histo->GetNbinsZ();i++)
  {
  	estimateZ_[i] = new float[4];
  	estimateQCDZ_[i] = new float[4];
	for(int j=0;j<4;j++)
	{
		estimateZ_[i][j] = -1;
		estimateQCDZ_[i][j] = -1;
	}
  }
  func = 0;
  line = 0;
  plotRatioEstim_Exp = 0;
  plotRatioQCDEstim_Exp = 0;
  /* canvasRatioEstim_Exp = 0; */
  /* canvasRatioQCDEstim_Exp = 0; */
  hErrorVsLuminosity_ = 0;
  hErrorVsLuminosityQCD_ = 0;
  /* cErrorVsLuminosity_ = new TCanvas("cErrorVsLuminosity"); */
  /* cErrorVsLuminosityQCD_ = new TCanvas("cErrorVsLuminosityQCD"); */
  Lumi_ = 100;
}

ABCDEstimation::~ABCDEstimation()
{
  if(h1DQCD_) delete h1DQCD_;
  if(h2DQCD_) delete h2DQCD_;
  if(hShapePredQCD_) delete hShapePredQCD_;
  if(plotRatioQCDEstim_Exp) delete plotRatioQCDEstim_Exp;
  /* if(canvasRatioQCDEstim_Exp) delete canvasRatioQCDEstim_Exp; */
  if(h2D_) delete h2D_;
  if(plotRatioEstim_Exp) delete plotRatioEstim_Exp;
  /* if(canvasRatioEstim_Exp) delete canvasRatioEstim_Exp; */
  if(hErrorVsLuminosity_)    delete hErrorVsLuminosity_;
  if(hErrorVsLuminosityQCD_) delete hErrorVsLuminosityQCD_;
  /* if(cErrorVsLuminosity_) delete cErrorVsLuminosity_; */
  /* if(cErrorVsLuminosityQCD_) delete cErrorVsLuminosityQCD_; */

  if(hXProfQCD_) delete hXProfQCD_;
  if(hYProfQCD_) delete hYProfQCD_;
  /* if(cXProfQCD_) delete cXProfQCD_; */
  /* if(cYProfQCD_) delete cYProfQCD_; */

  if(hXProf_) delete hXProf_;
  if(hYProf_) delete hYProf_;
  /* if(cXProf_) delete cXProf_; */
  /* if(cYProf_) delete cYProf_; */

  if(func) delete func;
  if(line) delete line;
  for(int i=0;i<4;i++) delete estimate_[i];
  for(int i=0;i<4;i++) delete estimateQCD_[i];
  if(h3D_) for(int i=0;i<h3D_->GetNbinsZ();i++) delete estimateZ_[i];
  if(h3DQCD_) for(int i=0;i<h3D_->GetNbinsZ();i++) delete estimateQCDZ_[i];
  if(h3D_) delete h3D_;
  if(h3DQCD_) delete h3DQCD_;
  if(hEstim_vs_Xcut) delete hEstim_vs_Xcut;
  if(hEstim_vs_Ycut) delete hEstim_vs_Ycut;
  try{
  	if(hShapePred_) delete hShapePred_;
  	delete [] estimate_;
  	if(h3D_) delete [] estimateZ_;
  	delete [] estimateQCD_;
  	if(h3DQCD_) delete [] estimateQCDZ_;
  }
  catch(...){}
}


ABCDEstimation::ABCDEstimation(const ABCDEstimation& a){
  Name_ = a.Name_;
  XaxisTitle_ = a.XaxisTitle_;
  YaxisTitle_ = a.YaxisTitle_;
  hShapePred_ = 0;
  hShapePredQCD_ = 0;
  hXProfQCD_ = 0;
  hYProfQCD_ = 0;
  /* cXProfQCD_ = 0; */
  /* cYProfQCD_ = 0; */
  hXProf_ = 0;
  hYProf_ = 0;
  /* cXProf_ = 0; */
  /* cYProf_ = 0; */

  h1DQCD_ = 0;
  h2DQCD_ = 0;
  h3DQCD_ = 0;
  h2D_ = 0;
  h3D_ = 0;
  func = 0;
  line = 0;
  plotRatioEstim_Exp = 0;
  plotRatioQCDEstim_Exp = 0;
  /* canvasRatioEstim_Exp = 0; */
  /* canvasRatioQCDEstim_Exp = 0; */
  hErrorVsLuminosity_ = 0;
  hErrorVsLuminosityQCD_ = 0;
  Lumi_ = 100;

  hEstim_vs_Xcut = 0;
  hEstim_vs_Ycut = 0;

  if(a.h1DQCD_) h1DQCD_ = (TH1F*) a.h1DQCD_->Clone(); 
  if(a.h2DQCD_) h2DQCD_ = (TH2F*) a.h2DQCD_->Clone(); 
  if(a.h3DQCD_) h3DQCD_ = (TH3F*) a.h3DQCD_->Clone(); 

  if(a.hShapePredQCD_)  hShapePredQCD_ = (TH1F*) a.hShapePredQCD_;
  if(a.hXProfQCD_)  hXProfQCD_ = (TGraphAsymmErrors*) a.hXProfQCD_;
  if(a.hYProfQCD_)  hYProfQCD_ = (TGraphAsymmErrors*) a.hYProfQCD_;
  /* if(a.cXProfQCD_)  cXProfQCD_ = (TCanvas*) a.cXProfQCD_; */
  /* if(a.cYProfQCD_)  cYProfQCD_ = (TCanvas*) a.cYProfQCD_; */
  
  if(a.hEstim_vs_Xcut) hEstim_vs_Xcut = (TH1F*) a.hEstim_vs_Xcut->Clone();
  if(a.hEstim_vs_Ycut) hEstim_vs_Ycut = (TH1F*) a.hEstim_vs_Ycut->Clone();

  if(a.h2D_) h2D_ = (TH2F*) a.h2D_->Clone(); 
  if(a.h3D_) h3D_ = (TH3F*) a.h3D_->Clone(); 

  if(a.hShapePred_) hShapePred_ = (TH1F*) a.hShapePred_->Clone();
  if(a.hXProf_)  hXProf_ = (TGraphAsymmErrors*) a.hXProf_;
  if(a.hYProf_)  hYProf_ = (TGraphAsymmErrors*) a.hYProf_;
  /* if(a.cXProf_)  cXProf_ = (TCanvas*) a.cXProf_; */
  /* if(a.cYProf_)  cYProf_ = (TCanvas*) a.cYProf_; */

  region_ = a.region_;
  estimate_ = new float*[4];
  estimateQCD_ = new float*[4];
  for(int i=0;i<4;i++){
  	estimate_[i] = new float[4];
  	estimateQCD_[i] = new float[4];
  	for(int j=0;j<4;j++){
  		estimate_[i][j] = a.estimate_[i][j];
  		estimateQCD_[i][j] = a.estimateQCD_[i][j];
  	}
  }
  if(a.h3D_){
  	estimateZ_ = new float*[a.h3D_->GetNbinsZ()];
  	for(int i=0;i<a.h3D_->GetNbinsZ();i++){
  		estimateZ_[i] = new float[4];
  		estimateQCDZ_[i] = new float[4];
  		for(int j=0;j<4;j++){
  			estimateZ_[i][j] = a.estimateZ_[i][j];
  			estimateQCDZ_[i][j] = a.estimateQCDZ_[i][j];
		}
	}
  }
  if(a.func) func = (TF1*) a.func->Clone();
  if(a.line) line = (TLine*) a.line->Clone();
  /* if(a.canvasRatioEstim_Exp) canvasRatioEstim_Exp = (TCanvas*) a.canvasRatioEstim_Exp->Clone(); */
  /* if(a.canvasRatioQCDEstim_Exp) canvasRatioQCDEstim_Exp = (TCanvas*) a.canvasRatioQCDEstim_Exp->Clone(); */
  if(a.plotRatioEstim_Exp) plotRatioEstim_Exp = (TH1F*) a.plotRatioEstim_Exp->Clone();
  if(a.plotRatioQCDEstim_Exp) plotRatioQCDEstim_Exp = (TH1F*) a.plotRatioQCDEstim_Exp->Clone();
  if(a.hErrorVsLuminosity_) hErrorVsLuminosity_ = (TGraph*) a.hErrorVsLuminosity_->Clone();
  if(a.hErrorVsLuminosityQCD_) hErrorVsLuminosityQCD_= (TGraph*) a.hErrorVsLuminosityQCD_->Clone();
  /* if(a.cErrorVsLuminosity_) cErrorVsLuminosity_      = (TCanvas*) a.cErrorVsLuminosity_->Clone(); */
  /* if(a.cErrorVsLuminosityQCD_) cErrorVsLuminosityQCD_= (TCanvas*) a.cErrorVsLuminosityQCD_->Clone(); */
  if(a.Lumi_) Lumi_ = a.Lumi_;
}


void ABCDEstimation::Fill(float X, float Y, float weight, bool isQCD) { 
	h2D_->Fill(X,Y,weight); 
	if(isQCD) h2DQCD_->Fill(X,Y,weight);
};
  
  
void ABCDEstimation::Fill3D(float X, float Y, float Z, float weight, bool isQCD) { 
	h2D_->Fill(X,Y,weight); 
	h3D_->Fill(X,Y,Z,weight); 
	if(isQCD){
		h1DQCD_->Fill(Y,weight); 
		h2DQCD_->Fill(X,Y,weight); 
		h3DQCD_->Fill(X,Y,Z,weight); 
	} 
};
  

ABCDEstimation& ABCDEstimation::operator +=(const ABCDEstimation& a)
{
  h1DQCD_->Add(a.h1DQCD_);
  h2DQCD_->Add(a.h2DQCD_);
  h3DQCD_->Add(a.h3DQCD_);
  h2D_->Add(a.h2D_);
  h3D_->Add(a.h3D_);
  return *this;
}

void ABCDEstimation::Draw(bool dofit){
        ///////////////////////////
        // Expection & Estimation 
        ///////////////////////////
        if(hShapePred_ )func = new TF1("func","[0]*x+[1]", hShapePred_->GetXaxis()->GetXmin(), hShapePred_->GetXaxis()->GetXmax());
        vector<TH1F*> vec;
        if(hShapePred_) vec.push_back((TH1F*) hShapePred_->Clone());
        if(h1DQCD_) vec.push_back((TH1F*) h1DQCD_->Clone());
        if(vec.size()>0 && vec[0]->Integral()>0) { vec[0]->Scale(1/vec[0]->Integral()); vec[0]->SetLineColor(2); vec[0]->SetLineWidth(2); vec[0]->SetFillStyle(0);}
        if(vec.size()>1 && vec[1]->Integral()>0) { vec[1]->Scale(1/vec[1]->Integral()); vec[1]->SetLineColor(4); vec[1]->SetLineWidth(2); vec[1]->SetFillStyle(0);}
        TLegend* legEstim_Exp = new TLegend(0.6,0.6,0.88,0.88);
        if(vec.size()>0) legEstim_Exp->AddEntry(vec[0],"Prediction", "l");
        if(vec.size()>1) legEstim_Exp->AddEntry(vec[1],"QCD Expcted", "l");
        if(dofit){
                plotRatioEstim_Exp = PlotRatio(vec[0],vec[1], func, string ("Estim_Exp_Ratio"));
                plotRatioEstim_Exp->GetFunction("func")->SetLineColor(2);
        }
        else if(vec.size()>1) plotRatioEstim_Exp = PlotRatio(vec[0],vec[1], string ("Estim_Exp_Ratio"));
        if(plotRatioEstim_Exp){
		plotRatioEstim_Exp->SetFillStyle(0);
       	 	plotRatioEstim_Exp->SetLineColor(1);
	}
        
	///////////////////////////
        // Expection & Estimation for QCD
        ///////////////////////////
        if(hShapePredQCD_ )func = new TF1("func","[0]*x+[1]", hShapePredQCD_->GetXaxis()->GetXmin(), hShapePredQCD_->GetXaxis()->GetXmax());
        vector<TH1F*> vec2;
        if(hShapePredQCD_) vec2.push_back((TH1F*) hShapePredQCD_->Clone());
        if(h1DQCD_) vec2.push_back((TH1F*) h1DQCD_->Clone());
        if(vec2.size()>0 && vec2[0]->Integral()>0) { vec2[0]->Scale(1/vec2[0]->Integral()); vec2[0]->SetLineColor(2); vec2[0]->SetLineWidth(2); vec2[0]->SetFillStyle(0);}
        if(vec2.size()>1 && vec2[1]->Integral()>0) { vec2[1]->Scale(1/vec2[1]->Integral()); vec2[1]->SetLineColor(4); vec2[1]->SetLineWidth(2); vec2[1]->SetFillStyle(0);}
        TLegend* legQCDEstim_Exp = new TLegend(0.6,0.6,0.88,0.88);
        if(vec2.size()>0) legQCDEstim_Exp->AddEntry(vec2[0],"Prediction", "l");
        if(vec2.size()>1) legQCDEstim_Exp->AddEntry(vec2[1],"QCD Expcted", "l");
        if(dofit){
                plotRatioQCDEstim_Exp = PlotRatio(vec2[0],vec2[1], func, string ("QCDEstim_Exp_Ratio"));
                plotRatioQCDEstim_Exp->GetFunction("func")->SetLineColor(2);
        }
        else if(vec2.size()>1) plotRatioQCDEstim_Exp = PlotRatio(vec2[0],vec2[1], string ("Estim_Exp_Ratio"));
        if(plotRatioQCDEstim_Exp){
		plotRatioQCDEstim_Exp->SetFillStyle(0);
       	 	plotRatioQCDEstim_Exp->SetLineColor(1);
	}

	//canvas
	/* canvasRatioEstim_Exp = TCanvasCreator(vec, plotRatioEstim_Exp, line, legEstim_Exp, string(""), string("EstDivExp")); */
	/* canvasRatioQCDEstim_Exp = TCanvasCreator(vec2, plotRatioQCDEstim_Exp, line, legQCDEstim_Exp, string(""), string("QCDEstDivExp")); */

}

void ABCDEstimation::Write(TFile* fout, string label){
        if(fout==0) return;
	fout->cd();
	string dirname = "ABCD"+label;
	fout->mkdir(dirname.c_str());
	fout->cd(dirname.c_str());
	if(h1DQCD_ != 0) h1DQCD_->Write();
	if(h2DQCD_ != 0) h2DQCD_->Write();
	if(h3DQCD_ != 0) h3DQCD_->Write();
	if(hXProfQCD_ !=0){
		hXProfQCD_->Write();
		/* cXProfQCD_->cd();
		hXProfQCD_->Draw("AP*"); */
		hXProfQCD_->GetXaxis()->SetTitle(XaxisTitle_);
	}
	/* if(cXProfQCD_ !=0) cXProfQCD_->Write(); */

	if(hYProfQCD_ !=0){
		hYProfQCD_->Write();
		/* cYProfQCD_->cd();
		hYProfQCD_->Draw("AP*"); */
		hYProfQCD_->GetXaxis()->SetTitle(YaxisTitle_);
	}
	/* if(cYProfQCD_ !=0) cYProfQCD_->Write(); */

	if(h2D_ != 0) h2D_->Write();
	if(h3D_ != 0) h3D_->Write();
	if(hXProf_ !=0){
		hXProf_->Write();
		/* cXProf_->cd();
		hXProf_->Draw("AP*"); */
		hXProf_->GetXaxis()->SetTitle(XaxisTitle_);
	}
	/* if(cXProf_ !=0) cXProf_->Write(); */

	if(hYProf_ !=0){
		hYProf_->Write();
		/* cYProf_->cd();
		hYProf_->Draw("AP*"); */
		hYProf_->GetXaxis()->SetTitle(YaxisTitle_);
	}
	/* if(cYProf_ !=0) cYProf_->Write(); */

	if(hShapePred_!=0) hShapePred_->Write();
	if(hShapePredQCD_!=0) hShapePredQCD_->Write();
	/* if(canvasRatioEstim_Exp!=0) canvasRatioEstim_Exp->Write(); */
	/* if(canvasRatioQCDEstim_Exp!=0) canvasRatioQCDEstim_Exp->Write(); */
	/* if(cErrorVsLuminosity_!=0) cErrorVsLuminosity_->Write(); */
	/* if(cErrorVsLuminosityQCD_!=0) cErrorVsLuminosityQCD_->Write(); */
	if(hEstim_vs_Xcut!=0) hEstim_vs_Xcut->Write();
	if(hEstim_vs_Ycut!=0) hEstim_vs_Ycut->Write();
}

void ABCDEstimation::ComputeEstimate(float cutXmin, float cutX0, float cutX1, float cutXmax, float cutYmin, float cutY0, float cutY1, float cutYmax, int region)
{
	Xmin_ = cutXmin; X0_ = cutX0; X1_ = cutX1; Xmax_ = cutXmax;
	Ymin_ = cutYmin; Y0_ = cutY0; Y1_ = cutY1; Ymax_ = cutYmax;
	region_ = region;

	int NbinsX=0, NbinsY=0, NbinsZ=0;
	float Xaxis_min=0, Xaxis_max=0, Yaxis_min=0, Yaxis_max=0;
	int Xmin=0, X0=0, X1=0, Xmax=0;
	int Ymin=0, Y0=0, Y1=0, Ymax=0;
	float A=0, B=0, C=0, D=0;
	float LumiBins[7] = {Lumi_/10,Lumi_/5,Lumi_/2,Lumi_,Lumi_*2,Lumi_*5,Lumi_*10};
	float ErrorBins[7], ErrorBinsQCD[7];
	TString ErrorCalcTitle = "Relative error calculated for";
	ErrorCalcTitle += (int)Lumi_;
	//
  	float** estimateLocal = new float*[4];
  	for(int i=0;i<4;i++){
		estimateLocal[i] = new float[4];
		for(int j=0;j<4;j++) estimateLocal[i][j] = 0;
	}
	//
	if(h2D_ == 0){
		cout<<"Can't perform the ABCD method : 2D histogram is empty..."<<endl;
	}
	else{
		NbinsX = h2D_->GetNbinsX()+1;
		NbinsY = h2D_->GetNbinsY()+1;

		float Xaxis_min = h2D_->GetXaxis()->GetXmin();
		float Xaxis_max = h2D_->GetXaxis()->GetXmax();
		float Yaxis_min = h2D_->GetYaxis()->GetXmin();
		float Yaxis_max = h2D_->GetYaxis()->GetXmax();

		Xmin = (int)(NbinsX*(cutXmin-Xaxis_min)/(Xaxis_max-Xaxis_min));
		X0   = (int)(NbinsX*(cutX0-Xaxis_min)/(Xaxis_max-Xaxis_min));
		X1   = (int)(NbinsX*(cutX1-Xaxis_min)/(Xaxis_max-Xaxis_min));
		Xmax = (int)(NbinsX*(cutXmax-Xaxis_min)/(Xaxis_max-Xaxis_min));

		Ymin = (int)(NbinsY*(cutYmin-Yaxis_min)/(Yaxis_max-Yaxis_min));
		Y0   = (int)(NbinsY*(cutY0-Yaxis_min)/(Yaxis_max-Yaxis_min));
		Y1   = (int)(NbinsY*(cutY1-Yaxis_min)/(Yaxis_max-Yaxis_min));
		Ymax = (int)(NbinsY*(cutYmax-Yaxis_min)/(Yaxis_max-Yaxis_min));

		A = h2D_->Integral(Xmin,X0,Ymin,Y0);
		B = h2D_->Integral(X1,Xmax,Ymin,Y0);
		C = h2D_->Integral(Xmin,X0,Y1,Ymax);
		D = h2D_->Integral(X1,Xmax,Y1,Ymax);

		ABCDEstimator(A,B,C,D,estimate_);

		float Num=0, Den=0, Est=0;
		for(int i=1;i<NbinsX;i++){
                if(region==1 || region==2) Num = h2D_->Integral(i,i+1,Ymin,Y0);
                else                       Num = h2D_->Integral(i,i+1,Y1,Y0);
                Den = h2D_->Integral(i,i+1,Ymin,Y0)+h2D_->Integral(i,i+1,Y1,Y0);
		Est = (Den>0 ? Num/Den : 0);
		hXProf_->SetPoint(i-1,(Xaxis_max-Xaxis_min)*(i-1)/(NbinsX-1),Est);
		hXProf_->SetPointError(i-1,0,0,Est-WilsonScoreIntervalLow(Num,Den),WilsonScoreIntervalHigh(Num,Den)-Est);
		}
		for(int i=1;i<NbinsY;i++){
                if(region==1 || region==3) Num = h2D_->Integral(Xmin,X0,i,i+1);
                else                       Num = h2D_->Integral(X1,Xmax,i,i+1);
                Den = h2D_->Integral(Xmin,X0,i,i+1)+h2D_->Integral(X1,Xmax,i,i+1);
		Est = (Den>0 ? Num/Den : 0);
		hYProf_->SetPoint(i-1,(Yaxis_max-Yaxis_min)*(i-1)/(NbinsY-1),Est);
		hYProf_->SetPointError(i-1,0,0,Est-WilsonScoreIntervalLow(Num,Den),WilsonScoreIntervalHigh(Num,Den)-Est);
		}
		ErrorBins[0] = ABCDRelErrorCalculator(A/10,B/10,C/10,D/10,region);
		ErrorBins[1] = ABCDRelErrorCalculator(A/5, B/5, C/5, D/5, region);
		ErrorBins[2] = ABCDRelErrorCalculator(A/2, B/2, C/2, D/2, region);
		ErrorBins[3] = ABCDRelErrorCalculator(A,B,C,D,region);
		ErrorBins[4] = ABCDRelErrorCalculator(A*2, B*2, C*2, D*2, region);
		ErrorBins[5] = ABCDRelErrorCalculator(A*5, B*5, C*5, D*5, region);
		ErrorBins[6] = ABCDRelErrorCalculator(A*10,B*10,C*10,D*10,region);
		hErrorVsLuminosity_ = new TGraph(7,LumiBins,ErrorBins);
		hErrorVsLuminosity_ ->SetTitle(ErrorCalcTitle);
		hErrorVsLuminosity_ ->GetXaxis()->SetTitle("Integrated luminosity [/pb]");
		/* cErrorVsLuminosity_->cd();hErrorVsLuminosity_ ->Draw("AC*"); */

		//compute Estim vs cut
		//X axis
		for(int i=X0;i<Xmax;i++){
			//replace Y1 by Y0 ? both are possible
			A = h2D_->Integral(Xmin,X0,Ymin,Y0);
			B = h2D_->Integral(i,Xmax,Ymin,Y0);
			C = h2D_->Integral(Xmin,X0,Y1,Ymax);
			D = h2D_->Integral(i,Xmax,Y1,Ymax);
			ABCDEstimator(A,B,C,D,estimateLocal);
			if(region>0 && region<5){
				hEstim_vs_Xcut->SetBinContent(i,estimateLocal[region-1][2]);
				hEstim_vs_Xcut->SetBinError(i,estimateLocal[region-1][3]);
			}
		}
		//Y axis
		for(int i=Y0;i<Ymax;i++){
			//replace X1 by X0 ? both are possible
			A = h2D_->Integral(Xmin,X0,Ymin,Y0);
			B = h2D_->Integral(X1,Xmax,Ymin,Y0);
			C = h2D_->Integral(Xmin,X0,i,Ymax);
			D = h2D_->Integral(X1,Xmax,i,Ymax);
			ABCDEstimator(A,B,C,D,estimateLocal);
			if(region>0 && region<5){
				hEstim_vs_Ycut->SetBinContent(i,estimateLocal[region-1][2]);
				hEstim_vs_Ycut->SetBinError(i,estimateLocal[region-1][3]);
			}
		}
	}
	if(h3D_ == 0) {
		cout<<"Can't perform a 3D estimate : 3D histogram is empty..."<<endl;
	}
	else{
  		NbinsZ = h3D_->GetNbinsZ()+1;
		for(int i=1;i<NbinsZ;i++)
		{
    			A = h3D_->Integral(Xmin,X0,Ymin,Y0,i,i+1);
    			B = h3D_->Integral(X1,Xmax,Ymin,Y0,i,i+1);
    			C = h3D_->Integral(Xmin,X0,Y1,Ymax,i,i+1);
    			D = h3D_->Integral(X1,Xmax,Y1,Ymax,i,i+1);
			ABCDEstimator(A,B,C,D,region,i-1,estimateZ_);
			hShapePred_->SetBinContent(i,estimateZ_[i-1][2]);	
			hShapePred_->SetBinError(i,estimateZ_[i-1][3]);	
		}	
	}
	if(h2DQCD_ == 0){
		cout<<"Can't perform the ABCD method : 2D histogram is empty..."<<endl;
	}
	else{
		NbinsX = h2DQCD_->GetNbinsX()+1;
		NbinsY = h2DQCD_->GetNbinsY()+1;

		Xaxis_min = h2DQCD_->GetXaxis()->GetXmin();
		Xaxis_max = h2DQCD_->GetXaxis()->GetXmax();
		Yaxis_min = h2DQCD_->GetYaxis()->GetXmin();
		Yaxis_max = h2DQCD_->GetYaxis()->GetXmax();

		Xmin = (int)(NbinsX*(cutXmin-Xaxis_min)/(Xaxis_max-Xaxis_min));
		X0   = (int)(NbinsX*(cutX0-Xaxis_min)/(Xaxis_max-Xaxis_min));
		X1   = (int)(NbinsX*(cutX1-Xaxis_min)/(Xaxis_max-Xaxis_min));
		Xmax = (int)(NbinsX*(cutXmax-Xaxis_min)/(Xaxis_max-Xaxis_min));

		Ymin = (int)(NbinsY*(cutYmin-Yaxis_min)/(Yaxis_max-Yaxis_min));
		Y0   = (int)(NbinsY*(cutY0-Yaxis_min)/(Yaxis_max-Yaxis_min));
		Y1   = (int)(NbinsY*(cutY1-Yaxis_min)/(Yaxis_max-Yaxis_min));
		Ymax = (int)(NbinsY*(cutYmax-Yaxis_min)/(Yaxis_max-Yaxis_min));

		A = h2DQCD_->Integral(Xmin,X0,Ymin,Y0);
		B = h2DQCD_->Integral(X1,Xmax,Ymin,Y0);
		C = h2DQCD_->Integral(Xmin,X0,Y1,Ymax);
		D = h2DQCD_->Integral(X1,Xmax,Y1,Ymax);

		ABCDEstimator(A,B,C,D,estimateQCD_);
		float Num=0, Den=0, Est=0;
		for(int i=1;i<NbinsX;i++){
                if(region==1 || region==2) Num = h2DQCD_->Integral(i,i+1,Ymin,Y0);
                else                       Num = h2DQCD_->Integral(i,i+1,Y1,Y0);
                Den = h2DQCD_->Integral(i,i+1,Ymin,Y0)+h2DQCD_->Integral(i,i+1,Y1,Y0);
		Est = (Den>0 ? Num/Den : 0);
		hXProfQCD_->SetPoint(i-1,(Xaxis_max-Xaxis_min)*(i-1)/(NbinsX-1),Est);
		hXProfQCD_->SetPointError(i-1,0,0,Est-WilsonScoreIntervalLow(Num,Den),WilsonScoreIntervalHigh(Num,Den)-Est);
		}
		for(int i=1;i<NbinsY;i++){
                if(region==1 || region==3) Num = h2DQCD_->Integral(Xmin,X0,i,i+1);
                else                       Num = h2DQCD_->Integral(X1,Xmax,i,i+1);
                Den = h2DQCD_->Integral(Xmin,X0,i,i+1)+h2DQCD_->Integral(X1,Xmax,i,i+1);
		Est = (Den>0 ? Num/Den : 0);
		hYProfQCD_->SetPoint(i-1,(Yaxis_max-Yaxis_min)*(i-1)/(NbinsY-1),Est);
		hYProfQCD_->SetPointError(i-1,0,0,Est-WilsonScoreIntervalLow(Num,Den),WilsonScoreIntervalHigh(Num,Den)-Est);
		}
		ErrorBinsQCD[0] = ABCDRelErrorCalculator(A/10,B/10,C/10,D/10,region);
		ErrorBinsQCD[1] = ABCDRelErrorCalculator(A/5, B/5, C/5, D/5, region);
		ErrorBinsQCD[2] = ABCDRelErrorCalculator(A/2, B/2, C/2, D/2, region);
		ErrorBinsQCD[3] = ABCDRelErrorCalculator(A,B,C,D,region);
		ErrorBinsQCD[4] = ABCDRelErrorCalculator(A*2, B*2, C*2, D*2, region);
		ErrorBinsQCD[5] = ABCDRelErrorCalculator(A*5, B*5, C*5, D*5, region);
		ErrorBinsQCD[6] = ABCDRelErrorCalculator(A*10,B*10,C*10,D*10,region);

		ErrorCalcTitle += " (multijet only)";
		hErrorVsLuminosityQCD_ = new TGraph(7,LumiBins,ErrorBinsQCD);
		hErrorVsLuminosityQCD_ ->SetTitle(ErrorCalcTitle);
		hErrorVsLuminosityQCD_ ->GetXaxis()->SetTitle("Integrated luminosity [/pb]");
		/* cErrorVsLuminosityQCD_->cd();hErrorVsLuminosityQCD_ ->Draw("AC*"); */
        }
	if(h3DQCD_ == 0) {
		cout<<"Can't perform a 3D estimate : 3D histogram is empty..."<<endl;
	}
	else{
  		NbinsZ = h3DQCD_->GetNbinsZ()+1;
		for(int i=1;i<NbinsZ;i++)
		{
    			A = h3DQCD_->Integral(Xmin,X0,Ymin,Y0,i,i+1);
    			B = h3DQCD_->Integral(X1,Xmax,Ymin,Y0,i,i+1);
    			C = h3DQCD_->Integral(Xmin,X0,Y1,Ymax,i,i+1);
    			D = h3DQCD_->Integral(X1,Xmax,Y1,Ymax,i,i+1);
			ABCDEstimator(A,B,C,D,region,i-1,estimateQCDZ_);
			hShapePredQCD_->SetBinContent(i,estimateQCDZ_[i-1][2]);	
			hShapePredQCD_->SetBinError(i,estimateQCDZ_[i-1][3]);	
		}	
	}
}

void ABCDEstimation::Print2DEstimate()
{
	TString Regions[4] = {"A","B","C","D"};
	cout<<endl;
	cout<<"***********************************************"<<endl;
	cout<<"*****      ABCD Estimation Results          ***"<<endl;
	cout<<"***********************************************"<<endl;			
	cout<<"------------ For QCD only   : ------------"<<endl;
	cout<<"---- "<<endl;
	cout<<"--- Correlation between the X and Y variables : "<<setprecision(3)<<GetCorrelationQCD()*100<<" %"<<endl;
	for(unsigned int i=0;i<4;i++)
	{
		cout<<"Region "<<Regions[i]<<" = "<<estimateQCD_[i][0]<<"+/-"<<estimateQCD_[i][1]<<"|| Estimate = "<<estimateQCD_[i][2]<<" +/- "<<estimateQCD_[i][3]<<endl;
	}
	cout<<"\\begin{table}"<<endl;
	cout<<"	\\centering"<<endl;
	cout<<"	 \\begin{tabular}{|c|c||c|}"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<" & $"<<Xmin_<<"<"<<XaxisTitle_<<"<"<<X0_<<"$ & $"<<X1_<<"<"<<XaxisTitle_<<"<"<<Xmax_<<"$ \\\\"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<"$"<<Ymin_<<"<"<<YaxisTitle_<<Y0_<< "$ & $"<<estimateQCD_[0][0]<<"\\pm"<<estimateQCD_[0][1]<<"$ & $"<<estimateQCD_[1][0]<<"\\pm"<<estimateQCD_[1][1]<<"$ \\\\"<<endl;
	if(region_ == 1) cout<<" & $\\left("<<estimateQCD_[0][2]<<"\\pm"<<estimateQCD_[0][3]<<"\\right)$ & \\\\"<<endl;
	if(region_ == 2) cout<<" &&$\\left("<<estimateQCD_[1][2]<<"\\pm"<<estimateQCD_[1][3]<<"\\right)$   \\\\"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<"$"<<Y1_<<"<"<<YaxisTitle_<<Ymax_<< "$ & $"<<estimateQCD_[2][0]<<"\\pm"<<estimateQCD_[2][1]<<"$ & $"<<estimateQCD_[3][0]<<"\\pm"<<estimateQCD_[3][1]<<"$ \\\\"<<endl;
	if(region_ == 3) cout<<" & $\\left("<<estimateQCD_[2][2]<<"\\pm"<<estimateQCD_[2][3]<<"\\right)$ & \\\\"<<endl;
	if(region_ == 4) cout<<" &&$\\left("<<estimateQCD_[3][2]<<"\\pm"<<estimateQCD_[3][3]<<"\\right)$   \\\\"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<"  \\end{tabular}"<<endl;
	cout<<"  \\caption{}"<<endl;
	cout<<"  \\label{tab:}"<<endl;
	cout<<"\\end{table}"<<endl;
	cout<<"------------ With all data  : ------------"<<endl;
	cout<<"---- "<<endl;
	cout<<"--- Correlation between the X and Y variables : "<<setprecision(3)<<GetCorrelation()*100<<" %"<<endl;
	for(unsigned int i=0;i<4;i++)
	{
		cout<<"Region "<<Regions[i]<<" = "<<estimate_[i][0]<<"+/-"<<estimate_[i][1]<<"|| Estimate = "<<estimate_[i][2]<<" +/- "<<estimate_[i][3]<<endl;
	}
	cout<<"\\begin{table}"<<endl;
	cout<<"	\\centering"<<endl;
	cout<<"	 \\begin{tabular}{|c|c||c|}"<<endl;
	cout<<"	\\hline"<<endl;
	cout.precision(1);
	cout<< fixed ;
	cout<<" & $"<<Xmin_<<"<"<<XaxisTitle_<<"<"<<X0_<<"$ & $"<<X1_<<"<"<<XaxisTitle_<<"<"<<Xmax_<<"$ \\\\"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<"$"<<Ymin_<<"<"<<YaxisTitle_<<"<"<<Y0_<< "$ & $"<<estimate_[0][0]<<"\\pm"<<estimate_[0][1]<<"$ & $"<<estimate_[1][0]<<"\\pm"<<estimate_[1][1]<<"$ \\\\"<<endl;
	if(region_ == 1) cout<<" & $\\left("<<estimate_[0][2]<<"\\pm"<<estimate_[0][3]<<"\\right)$ & \\\\"<<endl;
	if(region_ == 2) cout<<" &&$\\left("<<estimate_[1][2]<<"\\pm"<<estimate_[1][3]<<"\\right)$   \\\\"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<"$"<<Y1_<<"<"<<YaxisTitle_<<"<"<<Ymax_<< "$ & $"<<estimate_[2][0]<<"\\pm"<<estimate_[2][1]<<"$ & $"<<estimate_[3][0]<<"\\pm"<<estimate_[3][1]<<"$ \\\\"<<endl;
	if(region_ == 3) cout<<" & $\\left("<<estimate_[2][2]<<"\\pm"<<estimate_[2][3]<<"\\right)$ & \\\\"<<endl;
	if(region_ == 4) cout<<" &&$\\left("<<estimate_[3][2]<<"\\pm"<<estimate_[3][3]<<"\\right)$   \\\\"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<"  \\end{tabular}"<<endl;
	cout<<"  \\caption{}"<<endl;
	cout<<"  \\label{tab:}"<<endl;
	cout<<"\\end{table}"<<endl;
	cout<<"***********************************************"<<endl<<endl;			
}

void ABCDEstimation::Print3DEstimate()
{
	cout<<endl;
	cout<<"***********************************************"<<endl;
	cout<<"*****      ABCD Estimation Results          ***"<<endl;
	cout<<"***********************************************"<<endl;			
	if(h3DQCD_){
	cout<<"------------ For QCD only   : ------------"<<endl;
		for(int i=0;i<h3DQCD_->GetNbinsZ();i++)
		{
			cout<<"Data = "<<estimateQCDZ_[i][0]<<"+/-"<<estimateQCDZ_[i][1]<<"|| Estimate = "<<estimateQCDZ_[i][2]<<" +/- "<<estimateQCDZ_[i][3]<<endl;
		}
	}
	if(h3D_){
		cout<<"------------ With all data  : ------------"<<endl;
		for(int i=0;i<h3D_->GetNbinsZ();i++)
		{
			cout<<"Data = "<<estimateZ_[i][0]<<"+/-"<<estimateZ_[i][1]<<"|| Estimate = "<<estimateZ_[i][2]<<" +/- "<<estimateZ_[i][3]<<endl;
		}
	}
	cout<<"***********************************************"<<endl<<endl;			
	if(h3DQCD_==0 && h3D_==0) cout<<"Cannot be performed"<<endl;

}
/*
void ABCDEstimation::ABCDErrorCalc(float a, float b, float c, float d, float** tab)
{
	float Ntotal = a+b+c+d;
	tab[0][0] = a; tab[0][1] = Sigma_X(a,Ntotal); tab[0][2] = b*c/d; tab[0][3] = sqrt(pow(Sigma_XY(b,c,Ntotal)/d,2)+pow(b*c*Sigma_X(d,Ntotal)/pow(d,2),2)-2*b*c*CoVar_XY_Z(b,c,d,Ntotal)/pow(d,3));
	tab[1][0] = b; tab[1][1] = Sigma_X(b,Ntotal); tab[1][2] = a*d/c; tab[1][3] = sqrt(pow(Sigma_XY(a,d,Ntotal)/c,2)+pow(a*d*Sigma_X(c,Ntotal)/pow(c,2),2)-2*a*d*CoVar_XY_Z(a,d,c,Ntotal)/pow(c,3));
	tab[2][0] = c; tab[2][1] = Sigma_X(c,Ntotal); tab[2][2] = d*a/b; tab[2][3] = sqrt(pow(Sigma_XY(d,a,Ntotal)/b,2)+pow(d*a*Sigma_X(b,Ntotal)/pow(b,2),2)-2*d*a*CoVar_XY_Z(d,a,b,Ntotal)/pow(b,3));
	tab[3][0] = d; tab[3][1] = Sigma_X(d,Ntotal); tab[3][2] = c*b/a; tab[3][3] = sqrt(pow(Sigma_XY(c,b,Ntotal)/a,2)+pow(c*b*Sigma_X(a,Ntotal)/pow(a,2),2)-2*c*b*CoVar_XY_Z(c,b,a,Ntotal)/pow(a,3));
}

void ABCDEstimation::ABCDErrorCalc(float a, float b, float c, float d, int region, int i, float** tabZ)
{
	float Ntotal = a+b+c+d;
	if(region == 1) tabZ[i][0] = a; tabZ[i][1] = Sigma_X(a,Ntotal); tabZ[i][2] = b*c/d; tabZ[i][3] = sqrt(pow(Sigma_XY(b,c,Ntotal)/d,2)+pow(b*c*Sigma_X(d,Ntotal)/pow(d,2),2)-2*b*c*CoVar_XY_Z(b,c,d,Ntotal)/pow(d,3));
	if(region == 2) tabZ[i][0] = b; tabZ[i][1] = Sigma_X(b,Ntotal); tabZ[i][2] = a*d/c; tabZ[i][3] = sqrt(pow(Sigma_XY(a,d,Ntotal)/c,2)+pow(a*d*Sigma_X(c,Ntotal)/pow(c,2),2)-2*a*d*CoVar_XY_Z(a,d,c,Ntotal)/pow(c,3));
	if(region == 3) tabZ[i][0] = c; tabZ[i][1] = Sigma_X(c,Ntotal); tabZ[i][2] = d*a/b; tabZ[i][3] = sqrt(pow(Sigma_XY(d,a,Ntotal)/b,2)+pow(d*a*Sigma_X(b,Ntotal)/pow(b,2),2)-2*d*a*CoVar_XY_Z(d,a,b,Ntotal)/pow(b,3));
	if(region == 4) tabZ[i][0] = d; tabZ[i][1] = Sigma_X(d,Ntotal); tabZ[i][2] = c*b/a; tabZ[i][3] = sqrt(pow(Sigma_XY(c,b,Ntotal)/a,2)+pow(c*b*Sigma_X(a,Ntotal)/pow(a,2),2)-2*c*b*CoVar_XY_Z(c,b,a,Ntotal)/pow(a,3));
}
*/
void ABCDEstimation::ABCDEstimator(float a, float b, float c, float d, float** tab)
{
	float Ntotal = a+b+c+d;
	tab[0][0] = a; tab[0][1] = Sigma_X(a,Ntotal); tab[0][2] = b*c/d; tab[0][3] = sqrt(pow(Sigma_XY(b,c,Ntotal)/d,2)+pow(b*c*Sigma_X(d,Ntotal)/pow(d,2),2)-2*b*c*CoVar_XY_Z(b,c,d,Ntotal)/pow(d,3));
	tab[1][0] = b; tab[1][1] = Sigma_X(b,Ntotal); tab[1][2] = a*d/c; tab[1][3] = sqrt(pow(Sigma_XY(a,d,Ntotal)/c,2)+pow(a*d*Sigma_X(c,Ntotal)/pow(c,2),2)-2*a*d*CoVar_XY_Z(a,d,c,Ntotal)/pow(c,3));
	tab[2][0] = c; tab[2][1] = Sigma_X(c,Ntotal); tab[2][2] = d*a/b; tab[2][3] = sqrt(pow(Sigma_XY(d,a,Ntotal)/b,2)+pow(d*a*Sigma_X(b,Ntotal)/pow(b,2),2)-2*d*a*CoVar_XY_Z(d,a,b,Ntotal)/pow(b,3));
	tab[3][0] = d; tab[3][1] = Sigma_X(d,Ntotal); tab[3][2] = c*b/a; tab[3][3] = sqrt(pow(Sigma_XY(c,b,Ntotal)/a,2)+pow(c*b*Sigma_X(a,Ntotal)/pow(a,2),2)-2*c*b*CoVar_XY_Z(c,b,a,Ntotal)/pow(a,3));
}

void ABCDEstimation::ABCDEstimator(float a, float b, float c, float d, int region, int i, float** tabZ)
{
	float Ntotal = a+b+c+d;
	if(region == 1) tabZ[i][0] = a; tabZ[i][1] = Sigma_X(a,Ntotal); tabZ[i][2] = b*c/d; tabZ[i][3] = sqrt(pow(Sigma_XY(b,c,Ntotal)/d,2)+pow(b*c*Sigma_X(d,Ntotal)/pow(d,2),2)-2*b*c*CoVar_XY_Z(b,c,d,Ntotal)/pow(d,3));
	if(region == 2) tabZ[i][0] = b; tabZ[i][1] = Sigma_X(b,Ntotal); tabZ[i][2] = a*d/c; tabZ[i][3] = sqrt(pow(Sigma_XY(a,d,Ntotal)/c,2)+pow(a*d*Sigma_X(c,Ntotal)/pow(c,2),2)-2*a*d*CoVar_XY_Z(a,d,c,Ntotal)/pow(c,3));
	if(region == 3) tabZ[i][0] = c; tabZ[i][1] = Sigma_X(c,Ntotal); tabZ[i][2] = d*a/b; tabZ[i][3] = sqrt(pow(Sigma_XY(d,a,Ntotal)/b,2)+pow(d*a*Sigma_X(b,Ntotal)/pow(b,2),2)-2*d*a*CoVar_XY_Z(d,a,b,Ntotal)/pow(b,3));
	if(region == 4) tabZ[i][0] = d; tabZ[i][1] = Sigma_X(d,Ntotal); tabZ[i][2] = c*b/a; tabZ[i][3] = sqrt(pow(Sigma_XY(c,b,Ntotal)/a,2)+pow(c*b*Sigma_X(a,Ntotal)/pow(a,2),2)-2*c*b*CoVar_XY_Z(c,b,a,Ntotal)/pow(a,3));
}

float ABCDEstimation::ABCDErrorCalculator(float a, float b, float c, float d, int region)
{
    float Error = -1;
	float Ntotal = a+b+c+d;
	if(region == 1) Error = sqrt(pow(Sigma_XY(b,c,Ntotal)/d,2)+pow(b*c*Sigma_X(d,Ntotal)/pow(d,2),2)-2*b*c*CoVar_XY_Z(b,c,d,Ntotal)/pow(d,3));
	if(region == 2) Error = sqrt(pow(Sigma_XY(a,d,Ntotal)/c,2)+pow(a*d*Sigma_X(c,Ntotal)/pow(c,2),2)-2*a*d*CoVar_XY_Z(a,d,c,Ntotal)/pow(c,3));
	if(region == 3) Error = sqrt(pow(Sigma_XY(d,a,Ntotal)/b,2)+pow(d*a*Sigma_X(b,Ntotal)/pow(b,2),2)-2*d*a*CoVar_XY_Z(d,a,b,Ntotal)/pow(b,3));
	if(region == 4) Error = sqrt(pow(Sigma_XY(c,b,Ntotal)/a,2)+pow(c*b*Sigma_X(a,Ntotal)/pow(a,2),2)-2*c*b*CoVar_XY_Z(c,b,a,Ntotal)/pow(a,3));
    return Error;
}

float ABCDEstimation::ABCDRelErrorCalculator(float a, float b, float c, float d, int region)
{
    float RelError = -1;
	float Ntotal = a+b+c+d;
	if(region == 1) RelError = ((b*c/d)>0 ? sqrt(pow(Sigma_XY(b,c,Ntotal)/d,2)+pow(b*c*Sigma_X(d,Ntotal)/pow(d,2),2)-2*b*c*CoVar_XY_Z(b,c,d,Ntotal)/pow(d,3))/(b*c/d) : 0);
	if(region == 2) RelError = ((a*d/c)>0 ? sqrt(pow(Sigma_XY(a,d,Ntotal)/c,2)+pow(a*d*Sigma_X(c,Ntotal)/pow(c,2),2)-2*a*d*CoVar_XY_Z(a,d,c,Ntotal)/pow(c,3))/(a*d/c) : 0);
	if(region == 3) RelError = ((d*a/b)>0 ? sqrt(pow(Sigma_XY(d,a,Ntotal)/b,2)+pow(d*a*Sigma_X(b,Ntotal)/pow(b,2),2)-2*d*a*CoVar_XY_Z(d,a,b,Ntotal)/pow(b,3))/(d*a/b) : 0);
	if(region == 4) RelError = ((c*b/a)>0 ? sqrt(pow(Sigma_XY(c,b,Ntotal)/a,2)+pow(c*b*Sigma_X(a,Ntotal)/pow(a,2),2)-2*c*b*CoVar_XY_Z(c,b,a,Ntotal)/pow(a,3))/(c*b/a) : 0);
    return RelError;
}

void ABCDEstimation::ReScale(float factor){
	if(h2DQCD_ && h2DQCD_->GetEntries()>0) h2DQCD_->Scale(factor);
	if(h3DQCD_ && h3DQCD_->GetEntries()>0) h3DQCD_->Scale(factor);
	if(h2D_ && h2D_->GetEntries()>0) h2D_->Scale(factor);
	if(h3D_ && h3D_->GetEntries()>0) h3D_->Scale(factor);
}

float ABCDEstimation::CoVar_XY_Z(float X, float Y, float Z, float Ntot)
{
	return (-2*Ntot*(Ntot-1)*X*Y*Z/(pow(Ntot,3)));
}

float ABCDEstimation::CoVar_XY(float X, float Y, float Ntot)
{
	return (-Ntot*X*Y/(pow(Ntot,2)));
}

float ABCDEstimation::Sigma_XY(float X, float Y, float Ntot)
{
	return (sqrt(pow(X*Sigma_X(Y,Ntot),2)+pow(Y*Sigma_X(X,Ntot),2)+2*X*Y*CoVar_XY(X,Y,Ntot)));
}

float ABCDEstimation::Sigma_X(float X, float Ntot)
{
	return (sqrt(X*(1-(X/Ntot))));
}
float ABCDEstimation::Sigma_EffX(float X, float Ntot)
{
	return (sqrt(X*(1-(X/Ntot)))/Ntot);
}
float ABCDEstimation::WilsonScoreIntervalHigh(float Non, float Ntot)
{
	double T = (Ntot>0 ? 1/Ntot : 0);
	double p_hat = (Ntot>0 && Non>=0 && Ntot>=Non ? Non/Ntot : 0);
	double Int_High = ((p_hat+(T/2))/(1+T))+(sqrt(p_hat*(1-p_hat)*T+pow(T/2,2))/(1+T));
	return Int_High;
}
float ABCDEstimation::WilsonScoreIntervalLow(float Non, float Ntot)
{
	double T = (Ntot>0 ? 1/Ntot : 0);
	double p_hat = (Ntot>0 && Non>=0 && Ntot>=Non ? Non/Ntot : 0);
	double Int_Low = ((p_hat+(T/2))/(1+T))-(sqrt(p_hat*(1-p_hat)*T+pow(T/2,2))/(1+T));
	return Int_Low;
}
