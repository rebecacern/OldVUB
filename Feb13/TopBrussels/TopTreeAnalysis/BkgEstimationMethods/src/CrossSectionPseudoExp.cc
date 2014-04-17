#include "../interface/CrossSectionPseudoExp.h"

CrossSectionPseudoExp::CrossSectionPseudoExp(){
	 h_ABCD_Estimation = new TH1F("h_ABCD_Estimation","ABCD Estimation",10,0,20);
	 h_ABCD_Diff = new TH1F("h_ABCD_Diff","ABCD: Estimation - Expectation",10,-20,20);
	 h_ABCD_Pull = new TH1F("h_ABCD_Pull","ABCD: Pull distribution",10,-5,5);
	 h_VJets_Estimation = new TH1F("h_VJets_Estimation","V+jets Estimation",10,0,50);
	 h_VJets_Diff = new TH1F("h_VJets_Diff","V+jets: Estimation - Expectation",10,-50,50);
	 h_VJets_Pull = new TH1F("h_VJets_Pull","V+Jets: Pull distribution",10,-5,5);
	 h_TTJets_Estimation = new TH1F("h_TTJets_Estimation","TT+jets Estimation",10,0,50);
	 h_TTJets_Diff = new TH1F("h_TTJets_Diff","TT+jets: Estimation - Expectation",10,-50,50);
	 h_TTJets_Pull = new TH1F("h_TTJets_Pull","TT+Jets: Pull distribution",10,-5,5);
	 h_Xsection_Estimation = new TH1F("h_Xsection_Estimation","Xsection Estimation",10,50,150);
	 h_Xsection_Diff = new TH1F("h_Xsection_Diff","Xsection: Estimation - Expectation",10,-50,50);
	 h_Xsection_Pull = new TH1F("h_Xsection_Pull","Xsection: Pull distribution",10,-5,5);
	 h_Xsection_RelError = new TH1F("h_Xsection_RelError","Xsection: Statistical uncertainty",20,0,1);
}

CrossSectionPseudoExp::CrossSectionPseudoExp(int NofBins, float QCDMean, float VJetMean, float TTJetMean){
	 h_ABCD_Estimation = new TH1F("h_ABCD_Estimation","ABCD Estimation",NofBins,QCDMean*0.5,QCDMean*1.5);
	 h_ABCD_Diff = new TH1F("h_ABCD_Diff","ABCD: Estimation - Expectation",NofBins,-QCDMean*0.5,QCDMean*0.5);
	 h_ABCD_Pull = new TH1F("h_ABCD_Pull","ABCD: Pull distribution",NofBins,-5,5);
	 h_VJets_Estimation = new TH1F("h_VJets_Estimation","V+jets Estimation",NofBins,VJetMean*0.5,VJetMean*1.5);
	 h_VJets_Diff = new TH1F("h_VJets_Diff","V+jets: Estimation - Expectation",NofBins,-VJetMean*0.5,VJetMean*0.5);
	 h_VJets_Pull = new TH1F("h_VJets_Pull","V+Jets: Pull distribution",NofBins,-5,5);
	 h_TTJets_Estimation = new TH1F("h_TTJets_Estimation","TT+jets Estimation",NofBins,TTJetMean*0.5,TTJetMean*1.5);
	 h_TTJets_Diff = new TH1F("h_TTJets_Diff","TT+jets: Estimation - Expectation",NofBins,TTJetMean*0.5,TTJetMean*1.5);
	 h_TTJets_Pull = new TH1F("h_TTJets_Pull","TT+Jets: Pull distribution",NofBins,-5,5);
	 h_Xsection_Estimation = new TH1F("h_Xsection_Estimation","Xsection Estimation",NofBins,0,200);
	 h_Xsection_Diff = new TH1F("h_Xsection_Diff","Xsection: Estimation - Expectation",NofBins,-50,+50);
	 h_Xsection_Pull = new TH1F("h_Xsection_Pull","Xsection: Pull distribution",NofBins,-5,5);
	 h_Xsection_RelError = new TH1F("h_Xsection_RelError","Xsection: Statistical uncertainty",50,0,1);
}

CrossSectionPseudoExp::~CrossSectionPseudoExp(){
	delete h_ABCD_Estimation;
	delete h_ABCD_Diff;
	delete h_ABCD_Pull;
	delete h_VJets_Estimation;
	delete h_VJets_Diff;
	delete h_VJets_Pull;
	delete h_TTJets_Estimation;
	delete h_TTJets_Diff;
	delete h_TTJets_Pull;
	delete h_Xsection_Estimation;
	delete h_Xsection_Diff;
	delete h_Xsection_Pull;
	delete h_Xsection_RelError;
}

void CrossSectionPseudoExp::Fill(MCExpectation& mcexp, pair<float,float> QCDEstim, pair<float,float> VJetEstim,pair<float,float> TTJetEstim,pair<float,float> XsectionEstim){
	MCExp.push_back(mcexp);
	QCDEstimation.push_back(QCDEstim);
	VJetEstimation.push_back(VJetEstim);
	//fill plots
	h_ABCD_Estimation->Fill(QCDEstim.first);
	h_ABCD_Diff->Fill(QCDEstim.first-mcexp.NQCD.first);
	h_ABCD_Pull->Fill((QCDEstim.first-mcexp.NQCD.first)/QCDEstim.second);
	//change NWlj to all V+jets
	h_VJets_Estimation->Fill(VJetEstim.first);
	h_VJets_Diff->Fill(VJetEstim.first-mcexp.NWlj.first);
	h_VJets_Pull->Fill((VJetEstim.first-mcexp.NWlj.first)/VJetEstim.second);
	//Ntt
	h_TTJets_Estimation->Fill(TTJetEstim.first);
	h_TTJets_Diff->Fill(TTJetEstim.first-mcexp.NTtJets.first);
	h_TTJets_Pull->Fill((TTJetEstim.first-mcexp.NTtJets.first)/TTJetEstim.second);
	//Xsection
	h_Xsection_Estimation->Fill(XsectionEstim.first);
	h_Xsection_Diff->Fill(XsectionEstim.first-mcexp.NTtJets.first/mcexp.Luminosity);//divided by luminosity
	h_Xsection_Pull->Fill((XsectionEstim.first-mcexp.NTtJets.first)/XsectionEstim.second);
	if(XsectionEstim.first>0) h_Xsection_RelError->Fill(XsectionEstim.second/XsectionEstim.first);
	else h_Xsection_RelError->Fill(-9999);
}


void CrossSectionPseudoExp::Compute(){
	//nothing for the moment
}

void CrossSectionPseudoExp::Draw(){
	//nothing for the moment
}

void CrossSectionPseudoExp::Write(TFile* fout, string label){
  fout->cd ();
  string dirname = "CrossSectionMCPseudoExp" + label;
  fout->mkdir (dirname.c_str ());
  fout->cd (dirname.c_str ());
  h_ABCD_Estimation->Write();
  h_ABCD_Diff->Write();
  h_ABCD_Pull->Write();
  h_VJets_Estimation->Write();
  h_VJets_Diff->Write();
  h_VJets_Pull->Write();
  h_TTJets_Estimation->Write();
  h_TTJets_Diff->Write();
  h_TTJets_Pull->Write();
  h_Xsection_Estimation->Write();
  h_Xsection_Diff->Write();
  h_Xsection_Pull->Write();
  h_Xsection_RelError->Write();
}

