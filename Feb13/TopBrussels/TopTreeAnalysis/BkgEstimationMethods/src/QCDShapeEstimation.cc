#include "../interface/QCDShapeEstimation.h"


QCDShapeEstimation::QCDShapeEstimation(){
	NofQCDEventsExpected_ = -9999.;
	h_QCDMCShape_SR_ = 0;
	h_QCDMCShape_CR_ = 0;
	h_Shape_CR_ = 0;
	h_QCDShapeEstimated_ = 0;
	h_SysError_ = 0;
	line = 0;
	func = 0;
	plotRatioCR_SR = 0;
	plotRatioEstim_Exp = 0;
	canvasRatioCR_SR = 0;
	canvasRatioEstim_Exp = 0;
	canvasControlShape = 0;
	h_ControlShapeProfile_ = 0;
	legControlShape_ = 0;
}

		
QCDShapeEstimation::QCDShapeEstimation (int Nbins, float Xmin, float Xmax, string xlabel, string varLabel ){
	NofQCDEventsExpected_ = -9999.;
	h_QCDMCShape_SR_ = new TH1F(TString("h_QCDMCShape_SR_")+varLabel,"QCD shape from SR",Nbins, Xmin, Xmax);
	h_QCDMCShape_SR_->GetXaxis()->SetTitle(xlabel.c_str());
	h_QCDMCShape_SR_->GetYaxis()->SetTitle("Nof events");
	h_QCDMCShape_CR_ = new TH1F(TString("h_QCDMCShape_CR_")+varLabel,"QCD shape from CR",Nbins, Xmin, Xmax);
	h_QCDMCShape_CR_->GetXaxis()->SetTitle(xlabel.c_str());
	h_QCDMCShape_CR_->GetYaxis()->SetTitle("Nof events");
	h_Shape_CR_ = new TH1F(TString("h_Shape_CR_")+varLabel,"QCD shape from CR",Nbins, Xmin, Xmax);
	h_Shape_CR_->GetXaxis()->SetTitle(xlabel.c_str());
	h_Shape_CR_->GetYaxis()->SetTitle("Nof events");
	h_QCDShapeEstimated_ = new TH1F(TString("h_QCDShapeEstimated_")+varLabel,"QCD Shape Estimated",Nbins, Xmin, Xmax);
	h_QCDShapeEstimated_->GetXaxis()->SetTitle(xlabel.c_str());
	h_QCDShapeEstimated_->GetYaxis()->SetTitle("Nof events");
	h_SysError_ = new TH1F(TString("h_SysError_")+varLabel,"Sys Error for QCD Shape Estimation",Nbins,Xmin,Xmax);
	h_SysError_->GetXaxis()->SetTitle(xlabel.c_str());
	h_SysError_->GetYaxis()->SetTitle("Nof events");
	h_Projection_ = 0;
	line = 0;
	func = 0;
	plotRatioCR_SR = 0;
	plotRatioEstim_Exp = 0;
	canvasRatioCR_SR = 0;
	canvasRatioEstim_Exp = 0;
	canvasControlShape = 0;
	h_ControlShapeProfile_ = 0;
	legControlShape_ = 0;
}
	
QCDShapeEstimation::QCDShapeEstimation (int Nbins, float* Xbins, string xlabel, string varLabel ){
	NofQCDEventsExpected_ = -9999.;
	h_QCDMCShape_SR_ = new TH1F(TString("h_QCDMCShape_SR_")+varLabel,"QCD shape from SR",Nbins, Xbins);
	h_QCDMCShape_SR_->GetXaxis()->SetTitle(xlabel.c_str());
	h_QCDMCShape_SR_->GetYaxis()->SetTitle("Nof events");
	h_QCDMCShape_CR_ = new TH1F(TString("h_QCDMCShape_CR_")+varLabel,"QCD shape from CR",Nbins, Xbins);
	h_QCDMCShape_CR_->GetXaxis()->SetTitle(xlabel.c_str());
	h_QCDMCShape_CR_->GetYaxis()->SetTitle("Nof events");
	h_Shape_CR_ = new TH1F(TString("h_Shape_CR_")+varLabel,"QCD shape from CR",Nbins, Xbins);
	h_Shape_CR_->GetXaxis()->SetTitle(xlabel.c_str());
	h_Shape_CR_->GetYaxis()->SetTitle("Nof events");
	h_QCDShapeEstimated_ = new TH1F(TString("h_QCDShapeEstimated_")+varLabel,"QCD Shape Estimated",Nbins, Xbins);
	h_QCDShapeEstimated_->GetXaxis()->SetTitle(xlabel.c_str());
	h_QCDShapeEstimated_->GetYaxis()->SetTitle("Nof events");
	h_SysError_ = new TH1F(TString("h_SysError_")+varLabel,"Sys Error for QCD Shape Estimation",Nbins,Xbins);
	h_SysError_->GetXaxis()->SetTitle(xlabel.c_str());
	h_SysError_->GetYaxis()->SetTitle("Nof events");
	h_Projection_ = 0;
	line = 0;
	func = 0;
	plotRatioCR_SR = 0;
	plotRatioEstim_Exp = 0;
	canvasRatioCR_SR = 0;
	canvasRatioEstim_Exp = 0;
	canvasControlShape = 0;
	h_ControlShapeProfile_ = 0;
}


QCDShapeEstimation::~QCDShapeEstimation(){
	//if(h_Projection_) {for(int i=0;i<h_QCDMCShape_SR_->GetNbinsX();i++) {delete h_Projection_[i];}; delete [] h_Projection_;}
	if(h_QCDMCShape_SR_)       delete h_QCDMCShape_SR_;
	if(h_QCDMCShape_CR_)       delete h_QCDMCShape_CR_;
	if(h_Shape_CR_)            delete h_Shape_CR_;
	if(h_QCDShapeEstimated_)   delete h_QCDShapeEstimated_;
	if(h_SysError_)            delete h_SysError_;
	if(line)                   delete line;
	if(func)                   delete func;
	if(plotRatioCR_SR)         delete plotRatioCR_SR;
	if(plotRatioEstim_Exp)     delete plotRatioEstim_Exp;
	//if(h_ControlShapeProfile_) delete h_ControlShapeProfile_;
	//if(legControlShape_)       delete legControlShape_;
	if(canvasRatioCR_SR)       delete canvasRatioCR_SR;
	if(canvasRatioEstim_Exp)   delete canvasRatioEstim_Exp;
	//if(canvasControlShape)     delete canvasControlShape;
}

void QCDShapeEstimation::SetControlShape(vector<float> Xbins, string hname){
	double* tab = new double[Xbins.size()];
	for(unsigned int i=0;i<Xbins.size();i++){
		tab[i] = Xbins[i];
		char name[100];
		sprintf(name,"h_ControlShape_%f",Xbins[i]);
		TH1F* h = (TH1F*) h_QCDMCShape_SR_->Clone((name+hname).c_str());
		h->Reset();
		pair<float, TH1F*> p(Xbins[i],h);
		h_ControlShape_.push_back(p);
	}
	//
	legControlShape_ = new TLegend(0.7,0.5,0.9,0.9);
	h_ControlShapeProfile_ = new TProfile(("h_ControlShapeProfile"+hname).c_str(),("Profile: obs vs "+hname).c_str(),Xbins.size()-1,tab);
	h_ControlShapeProfile_->GetXaxis()->SetTitle(hname.c_str());
	if(h_QCDMCShape_SR_) h_ControlShapeProfile_->GetYaxis()->SetTitle(h_QCDMCShape_SR_->GetXaxis()->GetTitle());
	else h_ControlShapeProfile_->GetYaxis()->SetTitle("Observable");
	//
	h_Projection_ = new TH1F*[h_QCDMCShape_SR_->GetNbinsX()];
	for(int i=0;i<h_QCDMCShape_SR_->GetNbinsX();i++){
		char name[100];
		sprintf(name,"h_Projection_%d",i);
		h_Projection_[i] = new TH1F((name+hname).c_str(),"",Xbins.size()-1,tab);
	}
}

void QCDShapeEstimation::SetControlShape(vector<float> Xbins){
	double* tab = new double[Xbins.size()];
	for(unsigned int i=0;i<Xbins.size();i++){
		tab[i] = Xbins[i];
		char name[100];
		sprintf(name,"h_ControlShape_%f",Xbins[i]);
		TH1F* h = (TH1F*) h_QCDMCShape_SR_->Clone(name);
		h->Reset();
		pair<float, TH1F*> p(Xbins[i],h);
		h_ControlShape_.push_back(p);
	}
	//
	legControlShape_ = new TLegend(0.7,0.5,0.9,0.9);
	h_ControlShapeProfile_ = new TProfile("h_ControlShapeProfile","Profile: obs vs RelIso",Xbins.size()-1,tab);
	h_ControlShapeProfile_->GetXaxis()->SetTitle("RelIso");
	if(h_QCDMCShape_SR_) h_ControlShapeProfile_->GetYaxis()->SetTitle(h_QCDMCShape_SR_->GetXaxis()->GetTitle());
	else h_ControlShapeProfile_->GetYaxis()->SetTitle("Observable");
	//
	h_Projection_ = new TH1F*[h_QCDMCShape_SR_->GetNbinsX()];
	for(int i=0;i<h_QCDMCShape_SR_->GetNbinsX();i++){
		char name[100];
		sprintf(name,"h_Projection_%d",i);
		h_Projection_[i] = new TH1F(name,"",Xbins.size()-1,tab);
	}
}


void QCDShapeEstimation::Fill(float value, float weight, bool isQCD, bool isInSR, bool isInCR ){
	if(isQCD){
		if(isInSR) h_QCDMCShape_SR_->Fill(value, weight);
		if(isInCR) h_QCDMCShape_CR_->Fill(value, weight);
	}
	if(isInCR){
		h_Shape_CR_->Fill(value, weight);
		h_QCDShapeEstimated_->Fill(value, weight);
	}
}

void QCDShapeEstimation::Fill(float value, float weight, float controlValue, bool isQCD, bool isInSR, bool isInCR ){
	if(isQCD){
		if(isInSR) h_QCDMCShape_SR_->Fill(value, weight);
		if(isInCR) h_QCDMCShape_CR_->Fill(value, weight);
	}
	if(isInCR){
		h_Shape_CR_->Fill(value, weight);
		h_QCDShapeEstimated_->Fill(value, weight);
	}
	//control of the shapes
	if(isQCD){
		if(h_ControlShapeProfile_) h_ControlShapeProfile_->Fill(controlValue,value);
		if(h_ControlShape_.size()>0){
	 		for(unsigned int i=0;i<h_ControlShape_.size()-1;i++){
				if(controlValue>h_ControlShape_[i].first && controlValue<h_ControlShape_[i+1].first)
				h_ControlShape_[i].second->Fill(value, weight);
	 		}
			if(controlValue>h_ControlShape_[h_ControlShape_.size()-1].first) 
				h_ControlShape_[h_ControlShape_.size()-1].second->Fill(value, weight);
		}
	}
}


void QCDShapeEstimation::Normalize(float nofQCDExp) { 
	NofQCDEventsExpected_ = nofQCDExp; 
	if(h_QCDShapeEstimated_ && h_QCDShapeEstimated_->Integral()>0) h_QCDShapeEstimated_->Scale(nofQCDExp/h_QCDShapeEstimated_->Integral());
}


void QCDShapeEstimation::AddNormalisationError(float factor){
	factor = factor/100.;
	for(int i=1;i<h_QCDShapeEstimated_->GetNbinsX()+1;i++){
		float error = 0;
		error = sqrt((h_QCDShapeEstimated_->GetBinError(i)*h_QCDShapeEstimated_->GetBinError(i)) + (factor*factor*h_QCDShapeEstimated_->GetBinContent(i)*h_QCDShapeEstimated_->GetBinContent(i)));
		h_QCDShapeEstimated_->SetBinError(i,error);
	}
}
                

void QCDShapeEstimation::ComputeSystematicError(bool add){
	for(int i=1;i<h_QCDShapeEstimated_->GetNbinsX()+1;i++){
		float error = 0;
		float diff = 0;
		diff = h_QCDShapeEstimated_->GetBinContent(i)-h_QCDMCShape_SR_->GetBinContent(i);
		error = sqrt((h_QCDShapeEstimated_->GetBinError(i)*h_QCDShapeEstimated_->GetBinError(i)) + (diff*diff));
		if(add) h_QCDShapeEstimated_->SetBinError(i,error);
		h_SysError_->SetBinContent(i,diff);
	}
}
		
void QCDShapeEstimation::AddSystematicError(TH1F* SystError){
	for(int i=1;i<h_QCDShapeEstimated_->GetNbinsX()+1;i++){
		float error = 0;
		error = sqrt((h_QCDShapeEstimated_->GetBinError(i)*h_QCDShapeEstimated_->GetBinError(i)) + (h_SysError_->GetBinContent(i)*h_SysError_->GetBinContent(i)));
		h_QCDShapeEstimated_->SetBinError(i,error);
	}
}


void QCDShapeEstimation::Draw(bool dofit, string label){
	///////////////////////////
	// Ratio CR SR for QCD
	///////////////////////////
	func = new TF1("func","[0]*x+[1]", h_QCDMCShape_CR_->GetXaxis()->GetXmin(), h_QCDMCShape_CR_->GetXaxis()->GetXmax());
	vector<TH1F*> vec;
	if(h_QCDMCShape_CR_) vec.push_back((TH1F*) h_QCDMCShape_CR_->Clone());
	if(h_QCDMCShape_SR_) vec.push_back((TH1F*) h_QCDMCShape_SR_->Clone());
	if(vec.size()>0 && vec[0]->Integral()>0) { vec[0]->Sumw2(); vec[0]->Scale(1/vec[0]->Integral()); vec[0]->SetLineColor(2); vec[0]->SetLineWidth(2); vec[0]->SetFillStyle(0);}
	if(vec.size()>1 && vec[1]->Integral()>0) { vec[1]->Sumw2(); vec[1]->Scale(1/vec[1]->Integral()); vec[1]->SetLineColor(4); vec[1]->SetLineWidth(2); vec[1]->SetFillStyle(0);}
	TLegend* legCR_SR = new TLegend(0.6,0.6,0.88,0.88);
	if(vec.size()>0) legCR_SR->AddEntry(vec[0],"Control Region", "l");
	if(vec.size()>1) legCR_SR->AddEntry(vec[1],"Signal Region", "l");
	if(dofit){
		plotRatioCR_SR = PlotRatio(vec[0],vec[1], func, string ("CR_SR_Ratio"));
		plotRatioCR_SR->GetFunction("func")->SetLineColor(2);
	}
	else plotRatioCR_SR = PlotRatio(vec[0],vec[1], string ("CR_SR_Ratio"));
	plotRatioCR_SR->SetFillStyle(0);
	plotRatioCR_SR->SetLineColor(1);
					
  	///////////////////////////
    	// Ratio Estimation - Expectation
     	///////////////////////////
        
	if(dofit){
	        plotRatioEstim_Exp = PlotRatio( h_QCDShapeEstimated_ , h_QCDMCShape_SR_, func, string ("Est_Exp_Ratio"));
		plotRatioEstim_Exp->GetFunction("func")->SetLineColor(2);
        }
	else plotRatioEstim_Exp = PlotRatio(h_QCDShapeEstimated_,h_QCDMCShape_SR_, string ("Est_Exp_Ratio"));
	plotRatioEstim_Exp->SetFillStyle(0);
	plotRatioEstim_Exp->SetLineColor(1);
	vector<TH1F*> vec2;
	vec2.clear();
	vec2.push_back((TH1F*) h_QCDShapeEstimated_->Clone());
	vec2.push_back((TH1F*) h_QCDMCShape_SR_->Clone());
	//vec2[0]->SetLineColor(2); vec2[0]->SetLineWidth(2);
	//vec2[1]->SetLineColor(4); vec2[1]->SetLineWidth(2);
	//if(vec2[0]->Integral()>0) { vec2[0]->Sumw2(); vec2[0]->Scale(1/vec2[0]->Integral()); vec2[0]->SetLineColor(2); vec2[0]->SetLineWidth(2); vec2[0]->SetFillStyle(0);}
	//if(vec2[1]->Integral()>0) { vec2[1]->Sumw2(); vec2[1]->Scale(1/vec2[1]->Integral()); vec2[1]->SetLineColor(4); vec2[1]->SetLineWidth(2); vec2[1]->SetFillStyle(0);}
	TLegend* legEstim_Exp = new TLegend(0.6,0.6,0.88,0.88);
	if(vec2.size()>0) legEstim_Exp->AddEntry(vec2[0],"Estimated", "l");
	if(vec2.size()>1) legEstim_Exp->AddEntry(vec2[1],"Expected", "l");

        //Create the canvas
	canvasRatioCR_SR = TCanvasCreator(vec, plotRatioCR_SR, line, legCR_SR, string(""), string("CRDivSR"));
	canvasRatioEstim_Exp = TCanvasCreator(vec2, plotRatioEstim_Exp, line, legEstim_Exp, string(""), string("EstDivExp"));
	canvasControlShape = new TCanvas(TString("canvasControlShape")+label);
	//
	canvasControlShape->cd();
	for(unsigned int i=0;i<h_ControlShape_.size();i++){
		h_ControlShape_[i].second->SetLineColor(1+i);
		h_ControlShape_[i].second->Sumw2();
		if(h_ControlShape_[i].second->Integral()>0) h_ControlShape_[i].second->Scale(1/h_ControlShape_[i].second->Integral());
		if(i==0) h_ControlShape_[i].second->Draw("e");
		else h_ControlShape_[i].second->Draw("esame");
		char legChar [100];
		if(i<h_ControlShape_.size()-1) sprintf(legChar,"[%4.1f,%4.1f]",h_ControlShape_[i].first,h_ControlShape_[i+1].first);
		else sprintf(legChar,">%f",h_ControlShape_[i].first);
		legControlShape_->AddEntry(h_ControlShape_[i].second,legChar,"l");
	}
	legControlShape_->Draw("same");
	if(h_Projection_){
		for(int i=0;i<h_QCDMCShape_SR_->GetNbinsX();i++){
			for(unsigned int j=0;j<h_ControlShape_.size();j++){
				h_Projection_[i]->SetBinContent(j+1,h_ControlShape_[j].second->GetBinContent(i+1));
			}
		}
	}
}

void QCDShapeEstimation::PrintResults(){
	cout<<endl;
	cout<<"***********************************************"<<endl;
        cout<<"*****      QCDShapeEstimation Results       ***"<<endl;
	cout<<"***********************************************"<<endl;
	cout<<"Difference between QCD in CR & in SR: "<<Chi2Normalized(h_QCDMCShape_SR_, h_QCDMCShape_CR_,true)<<endl;
        cout<<"Difference between QCD Estimated & Expected: "<<Chi2Normalized(h_QCDMCShape_SR_, h_QCDShapeEstimated_,true)<<endl;
	cout<<"***********************************************"<<endl<<endl;
	
}

void QCDShapeEstimation::Write(TFile* fout, string label){
	if(fout==0) return;
	fout->cd();
	string dirname = "QCDShapeEstimation"+label;
	fout->mkdir(dirname.c_str());
	fout->cd(dirname.c_str());
	//plots
	if(h_QCDMCShape_SR_) h_QCDMCShape_SR_->Write();
	if(h_QCDMCShape_CR_) h_QCDMCShape_CR_->Write();
	if(h_Shape_CR_) h_Shape_CR_->Write();
	if(h_QCDShapeEstimated_) h_QCDShapeEstimated_->Write();
	if(h_SysError_) h_SysError_->Write();
	if(h_ControlShapeProfile_) h_ControlShapeProfile_->Write();
	//canvas
	if(canvasRatioCR_SR) canvasRatioCR_SR->Write();
	if(canvasRatioEstim_Exp) canvasRatioEstim_Exp->Write();
	//if(canvasControlShape) canvasControlShape->Write();// crash .. don't know why !!!!
	if(h_Projection_){
		for(int i=0;i<h_QCDMCShape_SR_->GetNbinsX();i++){
			h_Projection_[i]->Write();
		}
	}
}


