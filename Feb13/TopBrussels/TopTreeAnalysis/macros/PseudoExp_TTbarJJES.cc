#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

#include "../JESMeasurement/ROOTObjects/interface/Monster.h"

// Root stuff
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TBranch.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TRandom3.h"

// RooFit librairies
#include "RooArgSet.h"
#include "RooAddition.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooFitResult.h"

#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"

#include "Style.C"

using namespace std;
using namespace RooFit;

/// TGraphAsymmErrors
map<string,TGraphAsymmErrors*> graphAsymmErr;
map<string,TGraphErrors*> graphErr;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

vector<double> CalculateParabola(TH1D* parabola, int size, bool fit)
{ 
  // AN ORIGINAL PETRA VAN MULDERS!!!
  vector<double> getParabolaValues;
  Int_t minBin = parabola->GetMinimumBin();

  //cout << " DChi2 value of minbin-1 " << parabola->GetBinContent(minBin-1) << endl;
  //cout << " DChi2 value of minbin " << parabola->GetBinContent(minBin) << endl;
  //cout << " DChi2 value of minbin+1 " << parabola->GetBinContent(minBin+1) << endl;
  
  if (!fit)
  {
    double x1 = parabola->GetBinCenter(minBin);
    double x2 = parabola->GetBinCenter(minBin+size);
    double x3 = parabola->GetBinCenter(minBin-size);
    double y1 = parabola->GetBinContent(minBin);
    double y2 = parabola->GetBinContent(minBin+size);
    double y3 = parabola->GetBinContent(minBin-size);
    
    //y = ax^2 + bx + c
    double a = (y3-y1)/((x3-x1)*(x3-x2))-(y2-y1)/((x2-x1)*(x3-x2));
    double b = (y2-y1)/(x2-x1)-(y3-y1)*(x2+x1)/((x3-x1)*(x3-x2))+(y2-y1)*(x2+x1)/((x2-x1)*(x3-x2));
    double c = y1-(y3-y1)*(x1*x1)/((x3-x1)*(x3-x2))+(y2-y1)*(x1*x1)/((x2-x1)*(x3-x2))-(y2-y1)*(x1)/(x2-x1)+(y3-y1)*(x2+x1)*(x1)/((x3-x1)*(x3-x2))-(y2-y1)*(x2+x1)*(x1)/((x2-x1)*(x3-x2));
    
    double minim = - b/(2*a); //minimum gives JES correction 
    double minimchi2 = a*minim*minim + b*minim +c;
    double minplus = (-b + sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //upper bound right
    double minmin  = (-b - sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //lower bound left
    
    double uncert1 = minplus-minim;
    double uncert2 = minim-minmin;
    double uncert=0.;
    if(uncert1 > uncert2) { uncert = uncert1; } else { uncert = uncert2; }
    
    getParabolaValues.push_back(a);
    getParabolaValues.push_back(b);
    getParabolaValues.push_back(c);
    getParabolaValues.push_back(minim);
    getParabolaValues.push_back(uncert);

    //cout << "Analytic min " << minim << " minmin " << minmin << " minplus " << minplus << " uncert " << uncert << endl;
  }
  else
  {
    double lowLim = parabola->GetBinCenter(minBin) - parabola->GetBinCenter(minBin - size); //(size*0.02);
    double highLim = parabola->GetBinCenter(minBin) + parabola->GetBinCenter(minBin + size); //(size*0.02);

    string func_title = "Fitted";
    TF1* func = new TF1(func_title.c_str(),"pol2");

    func->SetRange(lowLim,highLim);
    //func->FixParameter(0,0);
    
    parabola->Fit(func,"RQ0");

    double a = func->GetParameter(2);
    double b = func->GetParameter(1);
    double c = func->GetParameter(0);

    double minim = - b/(2*a); //minimum gives JES correction 
    double minimchi2 = a*minim*minim + b*minim +c;
    
    double minplus = (-b + sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //upper bound right
    double minmin  = (-b - sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //lower bound left
    
    double uncert1 = minplus-minim;
    double uncert2 = minim-minmin;
    double uncert=0.;
    if(uncert1 > uncert2) { uncert = uncert1; } else { uncert = uncert2; }

    getParabolaValues.push_back(a);
    getParabolaValues.push_back(b);
    getParabolaValues.push_back(c);
    getParabolaValues.push_back(minim);
    getParabolaValues.push_back(uncert);

//    cout << "a " << a << "  b " << b << "  c " << c << endl;
//    cout << "Fitted min " << minim << "  minmin " << minmin << "  minplus " << minplus << "  uncert " << uncert << endl;
  }

  return vector<double>(getParabolaValues);
}

vector< vector<float> > GetEstimation(string supermonster, string pathforPNG, TH2F* histosupermonster, bool measureTopMass = false, bool writePlots = false) //get the estimated corrections
{
	if( histosupermonster->GetEntries() < 1 )
	{
	  cout << "Not enough entries in the SuperMonster: " << supermonster << "  So, not calculating the final results for this..." << endl;
	  vector< vector<float> > empty;
	  return empty;
	}
	
  float corrL = 0;
  float uncCorrL = 0;
  float corrB = 0;
  float uncCorrB = 0;
  float chi2 = -1;
  int ndf = -1;
	
  TH2F* histo = histosupermonster;
  Int_t binx = 0; Int_t biny = 0; Int_t binz = 0;
  Int_t bin = histo->GetMinimumBin(); histo->GetBinXYZ(bin,binx,biny,binz);
  float min = histo->GetBinContent(binx,biny);
  
  // Set errors to zero
  for(int i=0; i<=histo->GetNbinsX(); i++)
  {
    for(int j=0; j<=histo->GetNbinsY(); j++)
    {
//      histo->SetBinError(i, j, 1000000000000000.);
      histo->SetBinContent(i, j, histo->GetBinContent(i, j) - min );
    }
  }
  
	if(measureTopMass)
	{
    float minX = histo->GetXaxis()->GetBinCenter(binx);
    float minY = histo->GetYaxis()->GetBinCenter(biny);
    int nBins = 9; // 2D-fit on a nBins x nBins grid
    bool goLeft = ( histo->GetBinContent(binx-1,biny) < histo->GetBinContent(binx+1,biny) );
    bool goDown = ( histo->GetBinContent(binx,biny-1) < histo->GetBinContent(binx,biny+1) );
    int nBinsLeft = -1, nBinsRight = -1, nBinsBottom = -1, nBinsTop = -1;
    
    if( (nBins % 2) != 0 )
    {
      nBinsLeft = nBins/2;
      nBinsRight = nBins/2;
      nBinsBottom = nBins/2;
      nBinsTop = nBins/2;
    }
    else
    {
      if(goLeft)
      {
        nBinsLeft = nBins/2;
        nBinsRight = nBins/2 - 1;
      }
      else
      {
        nBinsLeft = nBins/2 - 1;
        nBinsRight = nBins/2;
      }
      if(goDown)
      {
        nBinsBottom = nBins/2;
        nBinsTop = nBins/2 - 1;
      }
      else
      {
        nBinsBottom = nBins/2 - 1;
        nBinsTop = nBins/2;
      }
    }
    
    float stepXLeft = histo->GetXaxis()->GetBinCenter(binx-nBinsLeft) - histo->GetXaxis()->GetBinWidth(binx)/3;
    float stepXRight = histo->GetXaxis()->GetBinCenter(binx+nBinsRight) + histo->GetXaxis()->GetBinWidth(binx)/3;
    float stepYBottom = histo->GetYaxis()->GetBinCenter(biny-nBinsBottom) - histo->GetYaxis()->GetBinWidth(biny)/3;
    float stepYTop = histo->GetYaxis()->GetBinCenter(biny+nBinsTop) + histo->GetYaxis()->GetBinWidth(biny)/3;
	  
	  string func_title = supermonster + "_2DFitted";
//	  TF2* func = new TF2(func_title.c_str(),"[0]+[1]*((x-[2])*cos([5])+(y-[4])*sin([5]))^2+[3]*(-(x-[2])*sin([5])+(y-[4])*cos([5]))^2");
//	  TF2* func = new TF2(func_title.c_str(),"[0]+[1]*(x*cos([5])+y*sin([5])-[2])^2+[3]*(-x*sin([5])+y*cos([5])-[4])^2"); extraction of the actual values is wrong for this function...
	  TF2* func = new TF2(func_title.c_str(),"[0]+[1]*(x-[2])^2+[3]*(y-[4])^2+[5]*(x-[2])*(y-[4])");
	  func->SetRange(stepXLeft, stepYBottom, stepXRight, stepYTop);
//	  func->SetParLimits(1,0,999999);
//	  func->SetParameter(1,1000);
//    func->SetParLimits(3,0,999999);
	  func->SetParLimits(2,0.4,1.5);
	  func->SetParameter(2,0.95);
	  func->SetParLimits(4,140,195);
	  func->SetParameter(4,172.5);
//	  histo->Fit(func,"RQ0 WW");
	  histo->Fit(func,"RQ0");
	  
	  histo->GetXaxis()->SetRangeUser(stepXLeft, stepXRight);
	  histo->GetYaxis()->SetRangeUser(stepYBottom, stepYTop);
	  
	  TH2F* histoRatio = PlotRatio(histo, func, func_title+"_Ratio");
	  
	  // create canvasses
	  if(writePlots)
	  {
      TCanvas* tempCanvas = TCanvasCreator(histo, func, func_title, "lego2"); //TCanvasCreator(pLight, func, func_title);
      tempCanvas->SaveAs( (pathforPNG+func_title+".png").c_str() );
      tempCanvas->Write();
      
      TCanvas* tempCanvas2 = TCanvasCreator(histo, supermonster+"_SmallRange");
      tempCanvas2->SaveAs( (pathforPNG+supermonster+"_SmallRange.png").c_str() );
      tempCanvas2->Write();
      
      TCanvas* tempCanvasRatio = TCanvasCreator(histoRatio, func_title+"_Ratio");
      tempCanvasRatio->SaveAs( (pathforPNG+func_title+"_Ratio.png").c_str() );
      tempCanvasRatio->Write();
    }
    
    // calculate the errors
    float p0 = func->GetParameter(0), p1 = func->GetParameter(1), p2 = func->GetParameter(2), p3 = func->GetParameter(3), p4 = func->GetParameter(4);
    float minChi2 = func->Eval(p2,p4);
    
    float Dx = -4*p1*( p0 - ( minChi2 + 1 ) );
    float xErrPlus = fabs( ( 2*p1*p2 + sqrt(Dx) ) / (2*p1) - p2 );
    float xErrMin = fabs( ( 2*p1*p2 - sqrt(Dx) ) / (2*p1) - p2 );
    float xErr = xErrPlus;
    if(xErrMin > xErrPlus) xErr = xErrMin;
        
    float Dy = -4*p3*( p0 - ( minChi2 + 1 ) );
    float yErrPlus = fabs( ( 2*p3*p4 + sqrt(Dy) ) / (2*p3) - p4 );
    float yErrMin = fabs( ( 2*p3*p4 - sqrt(Dy) ) / (2*p3) - p4 );
    float yErr = yErrPlus;
    if(yErrMin > yErrPlus) yErr = yErrMin;
    
	  // store the analysis results
    corrL = -(1-p2)*100;
    uncCorrL = xErr*100;
    corrB = p4;
    uncCorrB = yErr;
    chi2 = func->GetChisquare();
    ndf = func->GetNDF();

    //clean up
    delete func;
	}
	else
	{
    Double_t min = histo->GetBinContent(binx,biny);
    // project to light flavour correction and shift to 0
    TH1D* pLight = histo->ProjectionX((supermonster+"_projx").c_str(),biny,biny);
    pLight->GetXaxis()->SetTitle("#Delta E_{l}");
    pLight->GetYaxis()->SetTitle("#Delta #chi^{2}");
        
    double shift = pLight->GetBinContent(pLight->GetMinimumBin());

    for (int i=1; i<pLight->GetNbinsX()+1; i++)
    	pLight->SetBinContent(i,pLight->GetBinContent(i)-shift);
	  pLight->Write();

    // fit the parabola with a "pol2"
    vector<double> resultsF_l = CalculateParabola(pLight,5,true);
    double aF_l = resultsF_l[0];
    double bF_l = resultsF_l[1];
    double cF_l = resultsF_l[2];
        
    string func_title = string(pLight->GetName())+"_Fitted";

    TF1* func = new TF1(func_title.c_str(),"pol2",pLight->GetBinCenter(0)-(pLight->GetBinWidth(0)/2),pLight->GetBinCenter(pLight->GetNbinsX())+(pLight->GetBinWidth(0)/2));
    func->SetParameter(2,aF_l);
    func->SetParameter(1,bF_l);
    func->SetParameter(0,cF_l);

    // project to b flavour correction and shift to 0
    TH1D* pB = histo->ProjectionY((supermonster+"_projy").c_str(),binx,binx);
    pB->GetXaxis()->SetTitle("#Delta E_{b}");
    pB->GetYaxis()->SetTitle("#Delta #chi^{2}");

    shift = pB->GetBinContent(pB->GetMinimumBin());
    for (int i=1; i<pB->GetNbinsX()+1; i++)
     	pB->SetBinContent(i,pB->GetBinContent(i)-shift);
	  pB->Write();

    // fit the parabola with a "pol2"
    vector<double> resultsF_b = CalculateParabola(pB,5,true);
    double aF_b = resultsF_b[0];
    double bF_b = resultsF_b[1];
    double cF_b = resultsF_b[2];

    string func_titleb = string(pB->GetName())+"_Fitted";

    TF1* funcb = new TF1(func_titleb.c_str(),"pol2",pB->GetBinCenter(0)-(pB->GetBinWidth(0)/2),pB->GetBinCenter(pB->GetNbinsX())+(pB->GetBinWidth(0)/2));
    funcb->SetParameter(2,aF_b);
    funcb->SetParameter(1,bF_b);
    funcb->SetParameter(0,cF_b);

    // create canvasses
    if(writePlots)
    {
      TCanvas* tempCanvasLight = TCanvasCreator(pLight, func, func_title);
      tempCanvasLight->SaveAs( (pathforPNG+func_title+".png").c_str() );
      tempCanvasLight->Write();
      TCanvas* tempCanvasB = TCanvasCreator(pB, funcb, func_titleb);
      tempCanvasB->SaveAs( (pathforPNG+func_titleb+".png").c_str() );
      tempCanvasB->Write();
    }
    
    //clean up
    delete pLight;
    delete pB;
    delete func;
    delete funcb;
      
	  // show the analysis results :-)
    corrL = -(1-resultsF_l[3])*100;
    uncCorrL = resultsF_l[4]*100;
    corrB = -(1-resultsF_b[3])*100;
    uncCorrB = resultsF_b[4]*100;
  }
  
/*  if(writePlots)
  {
    cout << "Estimated light Correction = " << corrL << "+-" << uncCorrL << " %" << endl;
    if(measureTopMass)
    {
      cout << "Estimated Top Mass = " << corrB << "+-" << uncCorrB << " GeV"<< endl;
      cout << "ChiSquare / NDF of the fit = " << chi2 << " / " << ndf << " = " << chi2/ndf << endl;
    }
    else
      cout << "Estimated b Correction = " << corrB << "+-" << uncCorrB << " %"<< endl;
  }*/
  
  vector< vector<float> > result;
  vector<float> lResult, bResult;
  lResult.push_back(corrL);
  lResult.push_back(uncCorrL);
  bResult.push_back(corrB);
  bResult.push_back(uncCorrB);
  result.push_back(lResult);
  result.push_back(bResult);
  return result;
}


int main (int argc, char *argv[])
{
  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for TTbar JES Pseudo-experiments ! " << endl;
  cout << "*************************************************************" << endl;
 
  setTDRStyle();
//  setMyStyle();

  TFile *fout = new TFile ("PseudoExp_TTbarJES.root", "RECREATE");
  TDirectory* dir = fout->mkdir("SuperMonsters");
  
  string pathPNG = "PlotsPseudoExpJES/";
  mkdir(pathPNG.c_str(),0777);
  
  string inputMonsters[7] = {
    "Monsters/KinFit_Monsters_TopMass_TTbarJets_SemiMuon_Analysis.root",
    "Monsters/KinFit_Monsters_TopMass_TTbarJets_Other.root",
    "Monsters/KinFit_Monsters_TopMass_ST_tWChannel.root",
    "Monsters/KinFit_Monsters_TopMass_ST_tChannel.root",
    "Monsters/KinFit_Monsters_TopMass_WJets.root",
    "Monsters/KinFit_Monsters_TopMass_ZJets.root",
    "Monsters/KinFit_Monsters_TopMass_QCD_Mu15.root"
  };

  float eqLumi[7] = {5215.88, 7823.79, 46694.43, 23127.57, 352.84, 632.24, 163.18};
  
  float nEventToSelect[7] = {199.805, 27.3711, 5.74005, 5.89609, 155.571, 15.3895, 19.0206};
  
  // initialize histograms
  histo1D["maxMVA_probcuts_weighted"] = new TH1F("maxMVA_probcuts_weighted","maxMVA_probcuts_weighted",50,0,1);
  
  histo1D["nTotalEventsUsedMCweighted"] = new TH1F("nTotalEventsUsedMCweighted","nTotalEventsUsedMCweighted",400,199.5,599.5);

  histo1D["mTopPE"] = new TH1F("mTopPE","mTopPE",40,160,180);
  histo1D["mTopUncPE"] = new TH1F("mTopUncPE","mTopUncPE",40,0.6,1.4);
  histo1D["mTopPull"] = new TH1F("mTopPull","mTopPull",50,-10,10);
  histo1D["DElPE"] = new TH1F("DElPE","DElPE",40,-15,5);
  histo1D["DElUncPE"] = new TH1F("DElUncPE","DElUncPE",40,0.6,1.4);
  histo1D["DElPull"] = new TH1F("DElPull","DElPull",50,-10,10);
  
  // configuration
  bool measureTopMass = true;
  
  int nPseudoExperiments = 500;
  
//  float 
  
  vector<float> topMass, topMassUnc, DEl, DElUnc;

  TRandom3 *rndm = new TRandom3();
  rndm->SetSeed(0);

  for(unsigned int iExp=0; iExp<nPseudoExperiments; iExp++)
  {
    cout << "Performing pseudo-experiment nr: " << iExp << endl;
    
    stringstream ss; ss << iExp;
    string superMonsterTitle = "SumChi2KinFit_probcuts_"+ss.str();
    cout << superMonsterTitle << endl;
    
    float nTotalEventsUsedMCweighted = 0;
    
    for(unsigned int iDataSet=0; iDataSet<7; iDataSet++)
    {
      TFile* inFile = new TFile(inputMonsters[iDataSet].c_str(),"READ");
      
      TTree* inTree = (TTree*) inFile->Get("Monsters");
      TBranch* m_br = (TBranch*) inTree->GetBranch("Monsters");
      
      TClonesArray* monsters = new TClonesArray("Monster",0);
      m_br->SetAddress(&monsters);
      
      inTree->GetEvent(0);
      
      float nTotalEventsMCweighted = 0;
      for(unsigned int iEvt=0; iEvt<monsters->GetEntries(); iEvt++)
      {
        Monster* m =  (Monster*) monsters->At(iEvt);
        float eventWeight = m->eventWeight()*eqLumi[iDataSet] / 36.1389;
        nTotalEventsMCweighted += eventWeight;
      }
      
      float nEventsToUse = nEventToSelect[iDataSet];
      double PEweight = nEventsToUse/nTotalEventsMCweighted;
      
//      cout << "nTotalEventsMCweighted = " << nTotalEventsMCweighted << "  nEventsToUse = " << nEventsToUse << "  PEweight = " << PEweight;
      
      float nEventsUsedMCweighted = 0;
      
      for(unsigned int iEvt=0; iEvt<monsters->GetEntries(); iEvt++)
      {
        Monster* m =  (Monster*) monsters->At(iEvt);
//        cout << m->maxMVA() << " " << m->eventWeight() << endl;
        
        float eventWeight = m->eventWeight()*eqLumi[iDataSet] / 36.1389;
//        cout << "weight: " << eventWeight << endl;
        
        double random = rndm->Uniform(1);
        
        if(random < PEweight)
        {
          nEventsUsedMCweighted += eventWeight;
          TH2F monsterHisto = m->fitResults();
          
          Int_t binx = 0; Int_t biny = 0; Int_t binz = 0;
        	Int_t bin = monsterHisto.GetMaximumBin(); monsterHisto.GetBinXYZ(bin,binx,biny,binz);
        	float max = monsterHisto.GetBinContent(binx,biny);
        	float minMonster = monsterHisto.GetBinContent(monsterHisto.GetMinimumBin());
          float ProbNoCorr = monsterHisto.GetBinContent( (monsterHisto.GetNbinsX()/2)+1, (monsterHisto.GetNbinsX()/2)+1 );

          if( m->maxMVA() >= 0. && max > 0.98 ) // maxMVA cut and maxProb cut
				  {
	    			bool fitnotconverged = false;
					  //estimated corrections after probcuts
          	// convert to a chi2 map chi2 = -2*log[P(chi2)]
          	TH2F* chi2Histo2 = (TH2F*) monsterHisto.Clone();
					  Int_t maxBinx = 0; Int_t maxBiny = 0; Int_t maxBinz = 0;
            Int_t maxBin = monsterHisto.GetMaximumBin(); monsterHisto.GetBinXYZ(maxBin,maxBinx,maxBiny,maxBinz);
            float mTopMax = monsterHisto.GetYaxis()->GetBinCenter(maxBiny);
            float minDistHoleMax = 9999.;
            Int_t nHoles = 0;
          	for (unsigned int nBinX = 1; nBinX < chi2Histo2->GetNbinsX()+1; nBinX++) {
           		for (unsigned int nBinY = 1; nBinY < chi2Histo2->GetNbinsY()+1; nBinY++) {
           			chi2Histo2->SetBinContent(nBinX,nBinY,-2*log(pow(10.,-16.)));
             		if (monsterHisto.GetBinContent((Int_t)nBinX,(Int_t)nBinY) > 0)
             		{
               		chi2Histo2->SetBinContent(nBinX,nBinY,-2*log(monsterHisto.GetBinContent((Int_t)nBinX,(Int_t)nBinY)));
             		}
             		else
             		{
								  fitnotconverged = true;
								  nHoles++;
								  float distHoleMax = sqrt( pow( (float) nBinX - maxBinx ,2) + pow( (float) nBinY - maxBiny ,2) );
								  if( minDistHoleMax > distHoleMax ) minDistHoleMax = distHoleMax;
								  chi2Histo2->SetBinContent(nBinX,nBinY,minMonster);
							  }
           		}
          	}
          	
          	if( !fitnotconverged )
          	{
  //        	  histo1D["maxMVA_probcuts_weighted"]->Fill(m->maxMVA(), eventWeight);

          	  chi2Histo2->Scale(eventWeight);
          	  
      	  		if(histo2D.find(superMonsterTitle) == histo2D.end())
	          	{
	            	histo2D[superMonsterTitle] = (TH2F*) chi2Histo2->Clone();
	            	histo2D[superMonsterTitle]->SetNameTitle(superMonsterTitle.c_str(),superMonsterTitle.c_str());
	            	histo2D[superMonsterTitle]->GetXaxis()->SetTitle("Light Jet correction factor");
	            	histo2D[superMonsterTitle]->GetYaxis()->SetTitle("Top mass");
	            	histo2D[superMonsterTitle]->SetDirectory(gROOT);
	          	}
            	else
	            	histo2D[superMonsterTitle]->Add(chi2Histo2);

          	} // if noHoles
            delete chi2Histo2;
				  } // if maxMVA >= 0 en MaxProb
        } // if survived random cut        
      } // end loop over monsters
      
//      cout << "  nEventsUsedMCweighted = " << nEventsUsedMCweighted << endl;
      nTotalEventsUsedMCweighted += nEventsUsedMCweighted;
      
      delete monsters;
      
      inFile->Close();
      delete inFile;
    } // end loop over datasets
    
    fout->cd();
    dir->cd();
    string pathPNGsupermonsters = pathPNG + "SuperMonsters/";
    mkdir(pathPNGsupermonsters.c_str(),0777);
 	  vector< vector<float> > result = GetEstimation(superMonsterTitle, pathPNGsupermonsters, histo2D[superMonsterTitle], measureTopMass, true);
    fout->cd();
    
  	DEl.push_back(result[0][0]);
  	DElUnc.push_back(result[0][1]);
  	topMass.push_back(result[1][0]);
  	topMassUnc.push_back(result[1][1]);
  	cout << result[0][0] << " " << result[0][1] << " " << result[1][0] << " " << result[1][1] <<  endl;
  	  	
 	  histo1D["mTopPE"]->Fill(result[1][0]);
    histo1D["mTopUncPE"]->Fill(result[1][1]);
    histo1D["DElPE"]->Fill(result[0][0]);
    histo1D["DElUncPE"]->Fill(result[0][1]);
    
  	histo1D["nTotalEventsUsedMCweighted"]->Fill(nTotalEventsUsedMCweighted);
    
  } // end loop over pseudo-experiments
  
  // Calculate and fill the pull
  for(unsigned int i=0; i<topMass.size(); i++)
  {
    float mTopPull = ( topMass[i] - histo1D["mTopPE"]->GetMean() ) / topMassUnc[i];
    float DElPull = ( DEl[i] - histo1D["DElPE"]->GetMean() ) / DElUnc[i];
    histo1D["mTopPull"]->Fill(mTopPull);
    histo1D["DElPull"]->Fill(DElPull);
  }
  
  // Fit the pull
  TF1* gausTopMass = new TF1("gausTopMass","gaus");
  gausTopMass->SetLineWidth(2);
  histo1D["mTopPull"]->Fit(gausTopMass,"Q");
  
  TF1* gausDEl = new TF1("gausDEl","gaus");
  gausDEl->SetLineWidth(2);
  histo1D["mTopPull"]->Fit(gausDEl,"Q");
  
  // Write 1D histo's
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
	{
		TH1F *temp = it->second;
//		int N = temp->GetNbinsX();
//  	temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
//  	temp->SetBinContent(N+1,0);
//		temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}

  
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}




