#include "TString.h"
#include <stdio.h>
#include <time.h>
#include <ctime>
#include <iomanip>

#include "../interface/PtEtaBin.h"

#include <TFile.h>
#include <TObject.h>
#include <TH1.h>
#include <TProfile.h>
#include <TH2.h>
#include "TCanvas.h"
#include <TBranch.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TLine.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TClassTable.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TColor.h"
#include "TPaveStats.h"

#include "TFractionFitter.h";


using namespace std;


float round(float depth, float n)
{
    float d=((int)((n*depth)+0.5))/depth;
    	
    return d;
}

vector<double> calcStatUnc (bool silent, vector<double> values) {
	
	vector<double> res;
	
	if (values.size() < 11) return res;
	
	double ntt = values[0]; float sntt = values[1];
	double L = values[2]; float sL = values[3];
	double ebtag = values[4]; float sebtag = values[5];
	double emistag = values[6]; float semistag = values[7];
	double FracBoverTotal = values[8];
	double echi2 = values[9];
	double esel = values[10];
	
	double ebtagcut = ( FracBoverTotal * ebtag ) + ( ( 1 - FracBoverTotal ) * emistag );
	
	double sebtagcut = FracBoverTotal*sebtag;
	
	if (ebtag == 1) {
		
		ebtagcut = 1; 
		sebtagcut = 0; 
        
	}
	
	//cout << ntt << " " << sntt << " " << L << " " << ebtagcut << " " << sebtagcut << " " << emistag << " " << semistag << " " << esel << " " << echi2 << endl;
	if (!silent) cout << "calcStatUnc::UncBtagCut - " << sebtagcut << endl;
	
	res.push_back(sebtagcut);
	
	double exs_term1 = (pow(ntt*sebtagcut,2))/(pow(ebtagcut,4)*pow(echi2,2)*pow(esel,2)*pow(L,2));
	double exs_term2 = (pow(ntt*sL,2))/(pow(ebtagcut,2)*pow(echi2,2)*pow(esel,2)*pow(L,4));
	double exs_term3 = (pow(sntt,2))/(pow(ebtagcut,2)*pow(echi2,2)*pow(esel,2)*pow(L,2));
	
	//cout << exs_term1 << " " << exs_term2 << " " << exs_term3 << endl;
	
	double exs = sqrt(exs_term1+exs_term2+exs_term3);
	
	if (!silent) cout << "calcStatUnc::UncXS - " << exs << endl;
	
	res.push_back(exs);
	
	return res;
}

std::pair<double,double> combinedFitter (int iWP, string nSystematic, bool silent,string dir="FitOutput/") {

    
    stringstream nDisc; nDisc << iWP/3;
    stringstream nWP; nWP << iWP;
    
    string PrefixPlot = "";
    string btag = "";
    if (iWP%3==0) {
        PrefixPlot="Fit_bTagLCut";
        btag="bTagL";
    } else if ( (iWP-1) %3==0) {
        PrefixPlot="Fit_bTagMCut";
        btag="bTagM";
    } else if ( (iWP-2) %3==0) {
        PrefixPlot="Fit_bTagTCut";
        btag="bTagT";
    }
    
    // load templates
    
    std::map<string,TH1D*> h_mu, h_el, h_tmpl;

    TString histo_mu = dir+"/nSystematic_"+nSystematic+"_nDisc_"+nDisc.str()+"_ptbinlow_0_etabinlow_-9990__actual_fithistos_channel_Mu_"+PrefixPlot+".root";
    TString histo_el = dir+"/nSystematic_"+nSystematic+"_nDisc_"+nDisc.str()+"_ptbinlow_0_etabinlow_-9990__actual_fithistos_channel_El_"+PrefixPlot+".root";
    
    TFile* f_mu = new TFile(histo_mu,"READ");
    TFile* f_el = new TFile(histo_el,"READ");
    
    TFile* eff_mu = new TFile((dir+"/Efficiencies_nSystematic_"+nSystematic+"_Mu.root").c_str(),"READ");
    TFile* eff_el = new TFile((dir+"/Efficiencies_nSystematic_"+nSystematic+"_El.root").c_str(),"READ");
    
    if (!silent) cout << "*> Reading mu+jets templates from " << histo_mu << endl;
    
    h_mu["Data"] = (TH1D*) f_mu->Get(("TH1Data_"+btag).c_str());
    h_mu["TTbar"] = (TH1D*) f_mu->Get(("TH1Data_TTbar_"+btag).c_str());
    h_mu["Background"] = (TH1D*) f_mu->Get(("TH1Data_VVMC_"+btag).c_str());
    h_mu["QCD"] = (TH1D*) f_mu->Get(("TH1Data_multijet_"+btag).c_str());
    h_mu["eff"] = (TH1D*) eff_mu->Get(("eff_WP_"+nWP.str()).c_str());
    
    if (!silent) cout << "*> Reading el+jets templates from " << histo_el << endl;
    
    h_el["Data"] = (TH1D*) f_el->Get(("TH1Data_"+btag).c_str());
    h_el["TTbar"] = (TH1D*) f_el->Get(("TH1Data_TTbar_"+btag).c_str());
    h_el["Background"] = (TH1D*) f_el->Get(("TH1Data_VVMC_"+btag).c_str());
    h_el["QCD"] = (TH1D*) f_el->Get(("TH1Data_multijet_"+btag).c_str());
    h_el["eff"] = (TH1D*) eff_el->Get(("eff_WP_"+nWP.str()).c_str());

    //if (!silent) cout << h_mu["Data"]->GetBinContent(23) << endl;
    //if (!silent) cout << h_el["Data"]->GetBinContent(23) << endl;
    
    // make combined templates
    
    h_tmpl["Data"] = new TH1D("Data","Data",h_mu["Data"]->GetNbinsX()*2,0,h_mu["Data"]->GetNbinsX()*2);
    h_tmpl["TTbar"] = new TH1D("TTbar","TTbar",h_mu["Data"]->GetNbinsX()*2,0,h_mu["Data"]->GetNbinsX()*2);
    h_tmpl["Background_mu"] = new TH1D("Background_mu","Background",h_mu["Data"]->GetNbinsX()*2,0,h_mu["Data"]->GetNbinsX()*2);
    h_tmpl["Background_el"] = new TH1D("Background_el","Background",h_mu["Data"]->GetNbinsX()*2,0,h_mu["Data"]->GetNbinsX()*2);
    h_tmpl["QCD_mu"] = new TH1D("QCD_mu","QCD_mu",h_mu["Data"]->GetNbinsX()*2,0,h_mu["Data"]->GetNbinsX()*2);
    h_tmpl["QCD_el"] = new TH1D("QCD_el","QCD_el",h_mu["Data"]->GetNbinsX()*2,0,h_mu["Data"]->GetNbinsX()*2);
    
    if (!silent) cout << h_mu["eff"]->GetBinContent(1) << endl;
    if (!silent) cout << h_mu["eff"]->GetBinContent(2) << endl;
    if (!silent) cout << h_mu["eff"]->GetBinContent(3) << endl;
    if (!silent) cout << h_mu["eff"]->GetBinContent(4) << endl;
    if (!silent) cout << h_mu["eff"]->GetBinContent(5) << endl;
    if (!silent) cout << h_mu["eff"]->GetBinContent(6) << endl;
    if (!silent) cout << h_mu["eff"]->GetBinContent(7) << endl;
    
    //return 0;
    
    double lum_mu = h_mu["eff"]->GetBinContent(1)*1000;
    double g_mu = h_mu["eff"]->GetBinContent(5);
    double eb_mu = h_mu["eff"]->GetBinContent(7);
    double ueb_mu = h_mu["eff"]->GetBinError(7);
    double eq_mu = h_mu["eff"]->GetBinContent(6);
    
    double effbtag_mu = (g_mu*eb_mu)+((1-g_mu)*eq_mu);
    double effsel_mu = h_mu["eff"]->GetBinContent(2)*h_mu["eff"]->GetBinContent(3)*h_mu["eff"]->GetBinContent(4);
    
    double lum_el = h_el["eff"]->GetBinContent(1)*1000;
    double g_el = h_el["eff"]->GetBinContent(5);
    double eb_el = h_el["eff"]->GetBinContent(7);
    double eq_el = h_el["eff"]->GetBinContent(6);
    
    double effbtag_el = (g_el*eb_el)+((1-g_el)*eq_el);
    double effsel_el = h_el["eff"]->GetBinContent(2)*h_el["eff"]->GetBinContent(3)*h_el["eff"]->GetBinContent(4);
    
    cout << eq_mu << " " << eq_el << endl;
    
    //exit(1);
    // shift the W+jets shape to look like PDF rew one
    
    /*for (int b=0; b<h_mu["Background"]->GetNbinsX();b++) {
        double mlb=h_mu["Background"]->GetBinCenter(b);
        //double w = (0.000189*mlb)+0.9627;
        
        double warr[50]={0,0,0,0.922458,1.12541,0.960938,0.942082,1.00285,0.94403,0.969829,1.02156,0.975999,0.985114,0.993313,0.974077,0.96081,1.00626,1.0038,1.09399,1.07583,1.00075,0.978554,0.95174,0.984451,0.994625,0.965449,1.00505,0.965751,1.00737,0.963902,1.01811,0.963445,1.0532,1.28836,1.02454,0.998502,1.10814,1.06594,1.18937,1.0661,1.01729,1.16974,2.13829,0.992107,1.02643,1.12031,1.24993,0.912381,1.04525,1.15072};
        
        //cout << mlb << " " << w << endl;
        //h["TH1Data_WJets_bTagM"]->SetBinContent(b,w*h["TH1Data_WJets_bTagM"]->GetBinContent(b));
        h_mu["Background"]->SetBinContent(b,warr[b]*h_mu["Background"]->GetBinContent(b));
        
    }
    for (int b=0; b<h_el["Background"]->GetNbinsX();b++) {
        double mlb=h_el["Background"]->GetBinCenter(b);
        //double w = (0.000189*mlb)+0.9627;
        
        double warr[50]={0,0,0,0.922458,1.12541,0.960938,0.942082,1.00285,0.94403,0.969829,1.02156,0.975999,0.985114,0.993313,0.974077,0.96081,1.00626,1.0038,1.09399,1.07583,1.00075,0.978554,0.95174,0.984451,0.994625,0.965449,1.00505,0.965751,1.00737,0.963902,1.01811,0.963445,1.0532,1.28836,1.02454,0.998502,1.10814,1.06594,1.18937,1.0661,1.01729,1.16974,2.13829,0.992107,1.02643,1.12031,1.24993,0.912381,1.04525,1.15072};
        
        //cout << mlb << " " << w << endl;
        //h["TH1Data_WJets_bTagM"]->SetBinContent(b,w*h["TH1Data_WJets_bTagM"]->GetBinContent(b));
        h_el["Background"]->SetBinContent(b,warr[b]*h_el["Background"]->GetBinContent(b));
        
    }*/
       
    
    for (int i=1; i<h_tmpl["Data"]->GetNbinsX()+1; i++) {
    
        if (i <= h_mu["Data"]->GetNbinsX()) {

            //if (!silent) cout << "mu " << i << endl;
            
            //if (!silent) cout << "effbtag = " << effbtag << endl;
            //if (!silent) cout << "effsel = " << effsel << endl;
            
            double ttbarcont = h_mu["TTbar"]->GetBinContent(i);
            double wcont = h_mu["Background"]->GetBinContent(i);
            double datacont = h_mu["Data"]->GetBinContent(i);
            double qcdcont = h_mu["QCD"]->GetBinContent(i);
            
            h_tmpl["Data"]->SetBinContent(i,datacont);
            h_tmpl["TTbar"]->SetBinContent(i,ttbarcont);
            h_tmpl["Background_mu"]->SetBinContent(i,wcont);
            h_tmpl["QCD_mu"]->SetBinContent(i,qcdcont);
            
        } else {
            
            /*double attbarcont = h_mu["TTbar"]->GetBinContent(i-h_mu["Data"]->GetNbinsX());
            
            h_tmpl["Data"]->SetBinContent(i,h_mu["Data"]->GetBinContent(i-h_mu["Data"]->GetNbinsX()));
            h_tmpl["TTbar"]->SetBinContent(i,attbarcont);
            h_tmpl["Background_el"]->SetBinContent(i,h_mu["Background"]->GetBinContent(i-h_mu["Data"]->GetNbinsX()));
            h_tmpl["QCD_el"]->SetBinContent(i,h_mu["QCD"]->GetBinContent(i-h_mu["Data"]->GetNbinsX()));
            
            continue;*/
                                    
            //if (!silent) cout << "el " << i << endl;
            
            //double ttbarcont = h_el["TTbar"]->GetBinContent(i-h_mu["Data"]->GetNbinsX())*(effsel_mu/effsel_el)*(effbtag_mu/effbtag_el)*(lum_mu/lum_el);
            //double wcont = h_el["Background"]->GetBinContent(i-h_mu["Data"]->GetNbinsX())*(effsel_mu/effsel_el)*(effbtag_mu/effbtag_el)*(lum_mu/lum_el);
            //double datacont = h_el["Data"]->GetBinContent(i-h_mu["Data"]->GetNbinsX())*(effsel_mu/effsel_el)*(effbtag_mu/effbtag_el)*(lum_mu/lum_el);
            //double qcdcont = h_el["QCD"]->GetBinContent(i-h_mu["Data"]->GetNbinsX())*(effsel_mu/effsel_el)*(effbtag_mu/effbtag_el)*(lum_mu/lum_el);
            
            double ttbarcont = h_el["TTbar"]->GetBinContent(i-h_mu["Data"]->GetNbinsX())*(effsel_mu/effsel_el)*(effbtag_mu/effbtag_el);
            double wcont = h_el["Background"]->GetBinContent(i-h_mu["Data"]->GetNbinsX())*(effsel_mu/effsel_el)*(effbtag_mu/effbtag_el);
            double datacont = h_el["Data"]->GetBinContent(i-h_mu["Data"]->GetNbinsX())*(effsel_mu/effsel_el)*(effbtag_mu/effbtag_el)*(lum_mu/lum_el);
            double qcdcont = h_el["QCD"]->GetBinContent(i-h_mu["Data"]->GetNbinsX())*(effsel_mu/effsel_el)*(effbtag_mu/effbtag_el)*(lum_mu/lum_el);
            

            h_tmpl["Data"]->SetBinContent(i,datacont);
            h_tmpl["TTbar"]->SetBinContent(i,ttbarcont);
            
            h_tmpl["Background_el"]->SetBinContent(i,wcont);
            h_tmpl["QCD_el"]->SetBinContent(i,qcdcont);
            //h_tmpl["Background_mu"]->SetBinContent(i,wcont);
            //h_tmpl["QCD_mu"]->SetBinContent(i,qcdcont);
            
            /*double ttbarcont = h_el["TTbar"]->GetBinContent(i-h_mu["Data"]->GetNbinsX());
            
            h_tmpl["Data"]->SetBinContent(i,h_el["Data"]->GetBinContent(i-h_mu["Data"]->GetNbinsX()));
            h_tmpl["TTbar"]->SetBinContent(i,ttbarcont);
            h_tmpl["Background_el"]->SetBinContent(i,h_el["Background"]->GetBinContent(i-h_mu["Data"]->GetNbinsX()));
            h_tmpl["QCD_el"]->SetBinContent(i,h_el["QCD"]->GetBinContent(i-h_mu["Data"]->GetNbinsX()));*/
        }
    
    }
        
    TFile* fouta = new TFile("CombinedFitHistos_beforefit.root","RECREATE");
    
    fouta->cd();
    if (!silent) {
        for (std::map<string,TH1D*>::const_iterator it=h_tmpl.begin();it!=h_tmpl.end();++it)
            it->second->Write();
        
    }
    fouta->Close();
    
    // now do the Fit with TFractionFitter
    
    double fttb=0; double efttb=0;
    double fbkg_mu=0; double efbkg_mu=0;
    double fbkg_el=0; double efbkg_el=0;
    double fqcd_mu=0; double efqcd_mu=0;
    double fqcd_el=0; double efqcd_el=0;
    
    Int_t status = -1;
    
    // TFRACTIONFITTER
    
    TObjArray *mc = new TObjArray(5);        // MC histograms are put in this array
    mc->Add(h_tmpl["TTbar"]);
    mc->Add(h_tmpl["Background_mu"]);
    mc->Add(h_tmpl["Background_el"]);
    //mc->Add(h_tmpl["QCD_mu"]);
    //mc->Add(h_tmpl["QCD_el"]);
    
    TFractionFitter* fit = new TFractionFitter(h_tmpl["Data"], mc); // initialise
    fit->Constrain(0,0.,1.0);               // constrain fraction 1 to be between 0 and 1
    fit->Constrain(1,0.,1.0);               // constrain fraction 2 to be between 0 and 1
    fit->Constrain(2,0.,1.0);               // constrain fraction 2 to be between 0 and 1
    fit->Constrain(3,0.,1.0);               // constrain fraction 2 to be between 0 and 1
    fit->Constrain(4,0.,1.0);               // constrain fraction 2 to be between 0 and 1

    //fit->SetRangeX(1,50);
    
    status = fit->Fit();               // perform the fit
    
    fit->GetResult(0,fttb,efttb);
    fit->GetResult(1,fbkg_mu,efbkg_mu);
    fit->GetResult(2,fbkg_el,efbkg_el);
    fit->GetResult(3,fqcd_mu,efqcd_mu);
    fit->GetResult(4,fqcd_el,efqcd_el);
    
    // correct stat unc for 5 parameter fit
    
    double ep1=sqrt(fit->GetFitter()->GetCovarianceMatrixElement(0,0));
    double ep2=sqrt(fit->GetFitter()->GetCovarianceMatrixElement(1,1));
    double ep3=sqrt(fit->GetFitter()->GetCovarianceMatrixElement(2,2));
    double ep4=sqrt(fit->GetFitter()->GetCovarianceMatrixElement(3,3));
    double ep5=sqrt(fit->GetFitter()->GetCovarianceMatrixElement(4,4));
    
    double cov[10];
    
    cov[0]=fit->GetFitter()->GetCovarianceMatrixElement(0,1);
    cov[1]=fit->GetFitter()->GetCovarianceMatrixElement(0,2);
    cov[2]=fit->GetFitter()->GetCovarianceMatrixElement(0,3);
    cov[3]=fit->GetFitter()->GetCovarianceMatrixElement(0,4);
    cov[4]=fit->GetFitter()->GetCovarianceMatrixElement(1,2);
    cov[5]=fit->GetFitter()->GetCovarianceMatrixElement(1,3);
    cov[6]=fit->GetFitter()->GetCovarianceMatrixElement(1,4);
    cov[7]=fit->GetFitter()->GetCovarianceMatrixElement(2,3);
    cov[8]=fit->GetFitter()->GetCovarianceMatrixElement(2,4);
    cov[9]=fit->GetFitter()->GetCovarianceMatrixElement(3,4);
    
    double sum=fttb+fbkg_mu+fbkg_el+fqcd_mu+fqcd_el;
    
    double d1=(1/(sum))-(fttb/pow(sum,2));
    double d2=fttb/pow(sum,2);
    double d3=d2;
    double d4=d4;
    double d5=d2;
    
    double t1 = pow(d1,2) * pow(ep1,2);
    double t2 = pow(d2,2) * pow(ep2,2);
    double t3 = pow(d3,2) * pow(ep3,2);
    double t4 = pow(d4,2) * pow(ep4,2);
    double t5 = pow(d4,2) * pow(ep5,2);
    
    double tcov[10];

    tcov[0] = 2*cov[0]*d1*(d2*d2);
    tcov[1] = 2*cov[1]*d1*(d3*d3);
    tcov[2] = 2*cov[2]*d1*(d4*d4);
    tcov[3] = 2*cov[3]*d1*(d5*d5);    
    tcov[4] = 2*cov[4]*(d2*d2)*(d3*d3);
    tcov[5] = 2*cov[5]*(d2*d2)*(d4*d4);
    tcov[6] = 2*cov[6]*(d2*d2)*(d5*d5);    
    tcov[7] = 2*cov[7]*(d3*d3)*(d4*d4);
    tcov[8] = 2*cov[8]*(d3*d3)*(d5*d5);    
    tcov[9] = 2*cov[9]*(d4*d4)*(d5*d5);
    
    double sumtcov = 0;
    
    for (int t=0;t<10;t++) sumtcov+=tcov[t];

    double efttb_fixed = sqrt(t1+t2+t3+t4+t5+sumtcov);

    if (!silent) cout << "****** Fit results ******" << endl;
    
    if (!silent) cout << "Status: " << status << endl;
    
    if (!silent) cout << "fttb = " << fttb << " efttb = " << efttb << endl;
    if (!silent) cout << "fbkg_mu = " << fbkg_mu << " efbkg_mu = " << efbkg_mu << endl;
    if (!silent) cout << "fbkg_el = " << fbkg_el << " efbkg_el = " << efbkg_el << endl;
    if (!silent) cout << "fqcd_mu = " << fqcd_mu << " fqcd_mu = " << fqcd_mu << endl;
    if (!silent) cout << "fqcd_el = " << fqcd_el << " fqcd_el = " << fqcd_el << endl;
    
    double nTTbar_fitted=fttb*h_tmpl["Data"]->Integral()*0.5;
    double enTTbar_fitted=efttb_fixed*h_tmpl["Data"]->Integral()*0.5;
    
    if (!silent) cout << "nTTbar fitted = " << nTTbar_fitted  << " +- " << enTTbar_fitted <<  endl;
    
    if (!silent) cout << "Fit Chi2/NDF = " << fit->GetChisquare()/fit->GetNDF() << endl;
    
    if (!silent) cout << endl << "****** XS calculation parameters ******"<<endl;
    if (!silent) cout << "mu eff(b-tag cut) = " << effbtag_mu << endl;
    if (!silent) cout << "mu eff(sel) = " << effsel_mu << endl;
    if (!silent) cout << "mu lumi = " << lum_mu << endl;
    if (!silent) cout << "el eff(b-tag cut) = " << effbtag_el << endl;
    if (!silent) cout << "el eff(sel) = " << effsel_el << endl;
    if (!silent) cout << "el lumi = " << lum_el << endl;

    if (!silent) cout << endl << "****** XS result ******"<<endl;

    vector<double> uncvals;
	
	uncvals.push_back(nTTbar_fitted);
	uncvals.push_back(enTTbar_fitted);
	uncvals.push_back(lum_mu);
	
	uncvals.push_back(0);
	
	uncvals.push_back(eb_mu);
	uncvals.push_back(ueb_mu);
	uncvals.push_back(eq_mu);
	uncvals.push_back(0.);
	uncvals.push_back(g_mu);
	uncvals.push_back(h_mu["eff"]->GetBinContent(3)*h_mu["eff"]->GetBinContent(4));
	uncvals.push_back(h_mu["eff"]->GetBinContent(2));
	
	vector<double> UncForXSResults = calcStatUnc(true,uncvals);
    
    float xs=nTTbar_fitted/(lum_mu*effsel_mu*effbtag_mu);
    float uxs=UncForXSResults[1];
    
    if (!silent) cout << "sigma_ttbar = "<< xs << " +- " << uxs << " pb"<< endl << endl;
    if (!silent) cout << "sigma_ttbar = "<< round(10,xs) << " +- " << round(10,uxs) << " pb"<< endl << endl;

    //if (!silent) cout << h_mu["Data"]->Integral() << " " << h_tmpl["Data"]->Integral() << endl;
        
    //if (!silent) cout << 2300*effsel_mu*effbtag_mu << endl;
    //if (!silent) cout << effsel_mu*effbtag_mu << endl;
    
    // DRAW FIT PLOT
    
    h_tmpl["TTbar"]->Scale((h_tmpl["Data"]->Integral()*fttb)/h_tmpl["TTbar"]->Integral());
    h_tmpl["Background_mu"]->Scale((h_tmpl["Data"]->Integral()*fbkg_mu)/h_tmpl["Background_mu"]->Integral());
    h_tmpl["Background_el"]->Scale((h_tmpl["Data"]->Integral()*fbkg_el)/h_tmpl["Background_el"]->Integral());
    h_tmpl["QCD_mu"]->Scale((h_tmpl["Data"]->Integral()*fqcd_mu)/h_tmpl["QCD_mu"]->Integral());
    h_tmpl["QCD_el"]->Scale((h_tmpl["Data"]->Integral()*fqcd_el)/h_tmpl["QCD_el"]->Integral());
    
    //ttbar->Scale(ttbar->Integral()-150/ttbar->Integral());
    //vvmc->Scale(150);
    
    TH1D* result = (TH1D*) h_tmpl["TTbar"]->Clone();
    result->Add(h_tmpl["Background_mu"],1);
    result->Add(h_tmpl["Background_el"],1);
    result->Add(h_tmpl["QCD_mu"],1);
    result->Add(h_tmpl["QCD_el"],1);
    
    TCanvas* pom = new TCanvas("fitres","fitres");
    
    pom->cd();
    
    //gPad->SetLogy();
    
    h_tmpl["Data"]->GetXaxis()->SetTitle("m_{#mu j} (GeV/c2)");
    h_tmpl["TTbar"]->GetXaxis()->SetTitle("m_{#mu j} (GeV/c2)");
    
    h_tmpl["Data"]->GetXaxis()->SetTitle("m_{lb} (GeV/c2)");
    h_tmpl["TTbar"]->GetXaxis()->SetTitle("m_{lb} (GeV/c2)");
    
    h_tmpl["Data"]->GetYaxis()->SetTitle("# Events");
    h_tmpl["TTbar"]->GetYaxis()->SetTitle("# Events");
    
    h_tmpl["Data"]->SetLineColor(kBlack);
    h_tmpl["Data"]->SetMarkerColor(kBlack);    
    h_tmpl["Data"]->Draw("E");
    result->SetLineColor(kBlue);
    result->Draw("same hist");
    
    h_tmpl["TTbar"]->SetLineColor(kGreen);
    h_tmpl["TTbar"]->SetLineStyle(kDashed);
    h_tmpl["TTbar"]->Draw("same hist");
    
    h_tmpl["Background_mu"]->SetLineColor(kRed);
    h_tmpl["Background_mu"]->SetLineStyle(kDashed);
    h_tmpl["Background_mu"]->Draw("same hist");
    h_tmpl["Background_el"]->SetLineColor(kRed);
    h_tmpl["Background_el"]->SetLineStyle(kDashed);
    h_tmpl["Background_el"]->Draw("same hist");
    
    h_tmpl["QCD_mu"]->SetLineColor(kOrange);
    h_tmpl["QCD_mu"]->SetLineStyle(kDashed);
    h_tmpl["QCD_mu"]->Draw("same hist");
    h_tmpl["QCD_el"]->SetLineColor(kOrange);
    h_tmpl["QCD_el"]->SetLineStyle(kDashed);
    h_tmpl["QCD_el"]->Draw("same hist");
    
    TLegend *leg1 = new TLegend(0.65,0.73,0.998,0.99);
    
    leg1->AddEntry(h_tmpl["Data"],"Data", "P"); 
    leg1->AddEntry(result,"t#bar{t} + Background","L"); 
    leg1->AddEntry(h_tmpl["TTbar"],"t#bar{t}", "L"); 
    leg1->AddEntry(h_tmpl["Background_mu"],"Background", "L"); 
    leg1->AddEntry(h_tmpl["QCD_mu"],"Multi-jet", "L"); 
    
    stringstream f; f<<fit->GetChisquare()/fit->GetNDF();
    string fitProb = "Fit Chi2/ndf = "+f.str();
    TH1F* dummy = new TH1F("dummy","dummy",1,0,1); 
    dummy->SetLineColor(kWhite); dummy->SetMarkerColor(kWhite); 
    leg1->AddEntry(dummy,fitProb.c_str(), "P"); 
    
    leg1->Draw();
    
    TLatex* text = new TLatex(0.10,0.955,"CMS Preliminary 2.3 fb^{-1} at #sqrt{s}=8TeV");
        
    text->SetTextSize(0.05);
        
    text->SetNDC();
        
    text->Draw();   
    
    TFile* fout = new TFile("CombinedFitHistos.root","RECREATE");
    
    if (!silent) {
        for (std::map<string,TH1D*>::const_iterator it=h_tmpl.begin();it!=h_tmpl.end();++it)
            it->second->Write();
    
        pom->Write();

        fout->Close();
        
        pom->SaveAs("CombinedFitResult.png");
        
    }
    return std::pair<double,double>(xs,uxs);
    
}

int main (int argc, char* argv) {
    
/*
 int fitMode=0;
 string chi2cut="90";
 
 int iWP = 7;
 
 string nSystematic = "2";
*/
    int wp=7;
    
    bool doSyst=false;
    
    if (doSyst) {
        
        std::map<string,string> syst;
        
        syst["0"]="nominal";
        syst["1"]="JES-";
        syst["2"]="JES+";
        syst["3"]="JER-";
        syst["4"]="JER+";
        syst["5"]="PileUp-";
        syst["6"]="PileUp+";
        
        syst["17"]="WJetsXSDown-50perc";
        syst["18"]="WJetsXSUp-50per";
        syst["19"]="ZJetsXSDown-50perc";
        syst["20"]="ZJetsXSUp-50perc";
        
        /*
         
         syst["24"]="TopMass-163.5";
         syst["25"]="TopMass-181.5";
         
         syst["7"]="TTJets-ScaleDown";
         syst["8"]="TTJets-ScaleUp";
         syst["9"]="TTJets-MatchingDown";
         syst["10"]="TTJets-MatchingUp";
         
         syst["11"]="WJets-ScaleDown";
         syst["12"]="WJets-ScaleUp";
         syst["13"]="WJets-MatchingDown";
         syst["14"]="WJets-MatchingUp";
         
         */
        
        std::map<string,double> shift;
        
        double nomXS = 0;
        double unomXS = 0;
        
        stringstream wps; wps << wp;
        
        ofstream csv(("systematics/systematics_results_WP"+wps.str()+"_channel_Comb.csv").c_str(), ios::out | ios::trunc);
        
        for (std::map<string,string>::const_iterator it=syst.begin(); it != syst.end(); ++it) {
            
            std::pair<double,double> res = combinedFitter(wp,it->first,true);
            
            if (it->second == "nominal") {
                nomXS=res.first;
                unomXS=res.second;
            }
            
            double shift=(res.first-nomXS);
            
            //nominal ALL 1.325 0.035 0 0 -17.03 6.75 0 0 276.8 26 0 0 0.829618 0 
            
            string samp="ALL";
            
            if (it->second.find("TTJets") != string::npos) samp="TTJets";
            if (it->second.find("ZJets") != string::npos) samp="ZJets";
            if (it->second.find("WJets") != string::npos) samp="WJets";
            
            csv << it->second << " " << samp <<  " 0 0 0 0 0 0 0 0 " << res.first << " " << res.second << " " << shift << " 0 0 0 " << endl;
            
        }
        
        csv.close();
        
    } else {
        
        combinedFitter(wp,"-2",false);
        
    }
}
