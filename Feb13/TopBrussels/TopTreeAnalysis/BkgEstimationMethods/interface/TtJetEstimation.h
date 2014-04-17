#ifndef TtJetEstimation_h
#define TtJetEstimation_h

#include <iostream>
#include <iomanip>
#include <Math/VectorUtil.h>

#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TopTreeAnalysis/Content/interface/Dataset.h"
#include "TopTreeAnalysis/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysis/Selection/interface/Selection.h"

using namespace std;



        /**
	//	Aim: Estimation of ttbar contribution in a given observable (template normalized).
	//	This is done using a Control Region which provide a template and then normalize
	//	it to into a Normalisation Region (low part of the distribution).
	//	
	//	CR = Control Region\n
	//	SR = Signal Region\n
	//	NR = Normalization Region\n
	//
	//	Example of usage:
	//	- call default constructor
	//	- SetLisOfDatasets()  // used for MC comparison
	//	- ConfigHistos()
	//	- all the fill methods data/MC all/CR/SR:
	//		- FillByDatasetsBkg / FillBySRDatasetsBkg /FillByCRDatasetsBkg (Bkg refers to all the backgrounds of ttbar)
	//		- FillByDatasetsTtJets / FillBySRDatasetsTtJets /FillByCRDatasetsTtJets
	//		- FillByDatasetsNP / FillBySRDatasetsNP /FillByCRDatasetsNP (if New physics scenarii are used )
	//	- Scale(lumi)
	//	- ComputeObservableSR()
	//	- ComputeObservableCR()
	//	- ComputeObservableExpected()
	//	- ComputeObservableObserved()
	//	- DefineNR()
	//	- ComputeNormalization()
	//	- PrintResults()
	//	- ComputePlotSignificance()	
	//	- Draw()
	//	- Write()
	//
	//	Missing pieces:
	//	Work only on MC, this as to change. To be done:
	//	- Fill()  ( data for SR and CR )
	//	- No compute ObservableSR/CR and Observed ...
	*/
  
class TtJetEstimation {

  public:
  
  TtJetEstimation(bool isMC = true);
  ~TtJetEstimation();
 
  //copy constructor;
  TtJetEstimation(const TtJetEstimation&);

  //ConfigHistos should be called after SetListOfDatasets to run properly
  void ConfigHistos(int NbinsX, float Xmin, float Xmax, string ObsName, string AxisLabel, bool fillMode = true);
  void ConfigHistos(int NbinsX, float* binsX, string ObsName, string AxisLabel, bool fillMode = true);
  
  void Draw(string label = string(""));
  void PrintResults(float METCut = -999.);
  void Write(TFile* fout, string label = string(""));

  ////////////////////////////
  //      Fill              //
  ////////////////////////////
  void Fill(float variable, bool trigged, Selection& sel, AnalysisEnvironment& an, bool isInSR, string DatasetName = string("") );/** This method is the main Fill method and will call the others*/
  void FillSR(float obs) {h_ObservableSR_->Fill(obs);};
  void FillCR(float obs) {h_ObservableCR_->Fill(obs);};
  //Inclusive
  void FillByDatasetsBkg(float obs, unsigned int dataset) { if(dataset<h_ObservableByDatasetsBkg_.size()) h_ObservableByDatasetsBkg_[dataset]->Fill(obs);};
  void FillByDatasetsTtJets(float obs, unsigned int dataset) { if(dataset<h_ObservableByDatasetsTtJets_.size()) h_ObservableByDatasetsTtJets_[dataset]->Fill(obs);};
  void FillByDatasetsNP(float obs, unsigned int dataset) { if(dataset<h_ObservableByDatasetsNP_.size()) h_ObservableByDatasetsNP_[dataset]->Fill(obs);};
  //CR
  void FillCRByDatasetsBkg(float obs, unsigned int dataset) { if(dataset<h_ObservableCRByDatasetsBkg_.size()) h_ObservableCRByDatasetsBkg_[dataset]->Fill(obs);};
  void FillCRByDatasetsTtJets(float obs, unsigned int dataset) { if(dataset<h_ObservableCRByDatasetsTtJets_.size()) h_ObservableCRByDatasetsTtJets_[dataset]->Fill(obs);};
  void FillCRByDatasetsNP(float obs, unsigned int dataset) { if(dataset<h_ObservableCRByDatasetsNP_.size()) h_ObservableCRByDatasetsNP_[dataset]->Fill(obs);};
  //SR
  void FillSRByDatasetsBkg(float obs, unsigned int dataset) { if(dataset<h_ObservableSRByDatasetsBkg_.size()) h_ObservableSRByDatasetsBkg_[dataset]->Fill(obs);};
  void FillSRByDatasetsTtJets(float obs, unsigned int dataset) { if(dataset<h_ObservableSRByDatasetsTtJets_.size()) h_ObservableSRByDatasetsTtJets_[dataset]->Fill(obs);};
  void FillSRByDatasetsNP(float obs, unsigned int dataset) { if(dataset<h_ObservableSRByDatasetsNP_.size()) h_ObservableSRByDatasetsNP_[dataset]->Fill(obs);};

  ////////////////////////////
  //      Compute           //
  ////////////////////////////

  //sequence to be executed in that order
  bool Scale(float Lumi=100.); // rescale histograms with respect to Xsection
  void ReScale(float Factor=1.); // rescale histograms
  void Compute(bool doTt=true, bool doBkg=true, bool doNP=true, int iNP=0);/* Method calling other ComputeX() ... methods expect compute Normalisation**/
  void ComputeObservableSR(bool doTt=true, bool doBkg=true, bool doNP=true, int iNP=0);
  void ComputeObservableCR(bool doTt=true, bool doBkg=true, bool doNP=true, int iNP=0);
  void ComputeObservableExpected(bool doTt=true, bool doBkg=true); // no NP
  void ComputeObservableObserved(bool doTt=true, bool doBkg=true, bool doNP=true, int iNP=0); // w/o NP  (normaly observed = SR)
  bool ComputeNormalisation(TH1F* hQCDEstim = 0, TH1F* hWJetEstim = 0); // rescale histograms from CR with respect to nof events in NR   //->Fill Estimated (w/o NP contamination)

  ////////////////////////////
  //         Get plots      //
  ////////////////////////////
  TH1F* PlotCR() const {return h_ObservableCR_; };
  TH1F* PlotSR() const {return h_ObservableSR_; };
  //Inclusive
  vector<TH1F*> PlotByDatasets();
  vector<TH1F*> PlotByDatasetsBkg() const { return h_ObservableByDatasetsBkg_;};
  vector<TH1F*> PlotByDatasetsTtJets() const { return h_ObservableByDatasetsTtJets_;};
  vector<TH1F*> PlotByDatasetsNP() const { return h_ObservableByDatasetsNP_;};
  TH1F* PlotByDatasetsBkg(unsigned int dataset) const { if(dataset<h_ObservableByDatasetsBkg_.size()) return h_ObservableByDatasetsBkg_[dataset]; return NULL;};
  TH1F* PlotByDatasetsTtJets(unsigned int dataset) const { if(dataset<h_ObservableByDatasetsTtJets_.size()) return h_ObservableByDatasetsTtJets_[dataset]; return NULL;};
  TH1F* PlotByDatasetsNP(unsigned int dataset) const { if(dataset<h_ObservableByDatasetsNP_.size()) return h_ObservableByDatasetsNP_[dataset]; return NULL;};
  //CR
  vector<TH1F*> PlotByDatasetsCR();
  vector<TH1F*> PlotByDatasetsCRBkg() const { return h_ObservableCRByDatasetsBkg_;};
  vector<TH1F*> PlotByDatasetsCRTtJets() const { return h_ObservableCRByDatasetsTtJets_;};
  vector<TH1F*> PlotByDatasetsCRNP() const { return h_ObservableCRByDatasetsNP_;};
  TH1F* PlotByDatasetsCRBkg(unsigned int dataset) const { if(dataset<h_ObservableCRByDatasetsBkg_.size()) return h_ObservableCRByDatasetsBkg_[dataset]; return NULL;};
  TH1F* PlotByDatasetsCRTtJets(unsigned int dataset) const { if(dataset<h_ObservableCRByDatasetsTtJets_.size()) return h_ObservableCRByDatasetsTtJets_[dataset]; return NULL;};
  TH1F* PlotByDatasetsCRNP(unsigned int dataset) const { if(dataset<h_ObservableCRByDatasetsNP_.size()) return h_ObservableCRByDatasetsNP_[dataset]; return NULL;};
  //SR
  vector<TH1F*> PlotByDatasetsSR();
  vector<TH1F*> PlotByDatasetsSRBkg() const { return h_ObservableSRByDatasetsBkg_;};
  vector<TH1F*> PlotByDatasetsSRTtJets() const { return h_ObservableSRByDatasetsTtJets_;};
  vector<TH1F*> PlotByDatasetsSRNP() const { return h_ObservableSRByDatasetsNP_;};
  TH1F* PlotByDatasetsSRBkg(unsigned int dataset) const { if(dataset<h_ObservableSRByDatasetsBkg_.size()) return h_ObservableSRByDatasetsBkg_[dataset]; return NULL;};
  TH1F* PlotByDatasetsSRTtJets(unsigned int dataset) const { if(dataset<h_ObservableSRByDatasetsTtJets_.size()) return h_ObservableSRByDatasetsTtJets_[dataset]; return NULL;};
  TH1F* PlotByDatasetsSRNP(unsigned int dataset) const { if(dataset<h_ObservableSRByDatasetsNP_.size()) return h_ObservableSRByDatasetsNP_[dataset]; return NULL;};


  
  ////////////////////////////
  //  Estimation            //
  ////////////////////////////

  TH1F* PlotEstimated() const {return h_ObservableEstimated_;};
  TH1F* PlotExpected() const {return h_ObservableExpected_;};
  TH1F* PlotObserved() const {return h_ObservableObserved_;};
  pair<float,float> NofEvents( TH1F* histo, pair<float,float> edges) const;
  pair<float,float> NofEventsEstimated( float VarCut) const;
  pair<float,float> NofEventsExpected( float VarCut) const;
  pair<float,float> NofEventsObserved( float VarCut) const;
  pair<float,float> SignificanceEstimated( float VarCut) const;
  pair<float,float> SignificanceExpected( float VarCut) const;

  void ComputePlotSignificance(TH1F* qcdShapeEstim = 0, TH1F* wjetShapeEstim = 0);
  void ComputeEstimationLMXContamination(TH1F* qcdShapeEstim = 0, TH1F* wjetShapeEstim = 0);

  void AddNormalisationError(float factor);//this factor is a fraction %
  void ComputeSystematicError(bool add = true);
  void AddSystematicError(TH1F* SystError); //coming from MC ...
  TH1F* GetSystematicError() {return h_SysError_;};
  ////////////////////////////
  // Set plots (from files) //
  ////////////////////////////

  void SetPlotByDatasets(vector<TFile*> files, string histoname, vector<TH1F*> histos, vector<Dataset> listRef);
  //Inclusive
  void SetPlotByDatasetsBkg(vector<TFile*> files, string histoname){ SetPlotByDatasets(files, histoname, h_ObservableByDatasetsBkg_, listOfDatasetsBkg_);};
  void SetPlotByDatasetsTtJets(vector<TFile*> files, string histoname){ SetPlotByDatasets(files, histoname, h_ObservableByDatasetsTtJets_, listOfDatasetsTtJets_);};
  void SetPlotByDatasetsNP(vector<TFile*> files, string histoname){ SetPlotByDatasets(files, histoname, h_ObservableByDatasetsNP_, listOfDatasetsNP_);};
  //CR
  void SetPlotByDatasetsCRBkg(vector<TFile*> files, string histoname){ SetPlotByDatasets(files, histoname, h_ObservableCRByDatasetsBkg_, listOfDatasetsBkg_);};
  void SetPlotByDatasetsCRTtJets(vector<TFile*> files, string histoname){ SetPlotByDatasets(files, histoname, h_ObservableCRByDatasetsTtJets_, listOfDatasetsTtJets_);};
  void SetPlotByDatasetsCRNP(vector<TFile*> files, string histoname){ SetPlotByDatasets(files, histoname, h_ObservableCRByDatasetsNP_, listOfDatasetsNP_);};
  //SR
  void SetPlotByDatasetsSRBkg(vector<TFile*> files, string histoname){ SetPlotByDatasets(files, histoname, h_ObservableSRByDatasetsBkg_, listOfDatasetsBkg_);};
  void SetPlotByDatasetsSRTtJets(vector<TFile*> files, string histoname){ SetPlotByDatasets(files, histoname, h_ObservableSRByDatasetsTtJets_, listOfDatasetsTtJets_);};
  void SetPlotByDatasetsSRNP(vector<TFile*> files, string histoname){ SetPlotByDatasets(files, histoname, h_ObservableSRByDatasetsNP_, listOfDatasetsNP_);};



  ////////////////////////////
  // Datasets & factor & colors //
  ////////////////////////////
  vector<Dataset> ListOfDatasets() const;
  vector<string> ListOfDatasetsName() const;
  vector<int> ListOfDatasetsColor() const;
  vector<float> ListOfDatasetsNormFactor() const;
  
  vector<Dataset> ListOfDatasetsTtJets() const {return listOfDatasetsTtJets_;};
  vector<string> ListOfDatasetsNameTtJets() const { vector<string> a; for(unsigned int i=0;i<listOfDatasetsTtJets_.size();i++) a.push_back(listOfDatasetsTtJets_[i].Name()); return a;};
  vector<int> ListOfDatasetsColorTtJets() const { vector<int> a; for(unsigned int i=0;i<listOfDatasetsTtJets_.size();i++) a.push_back(listOfDatasetsTtJets_[i].Color()); return a;};
  vector<float> ListOfDatasetsNormFactorTtJets() const { vector<float> a; for(unsigned int i=0;i<listOfDatasetsTtJets_.size();i++) a.push_back(listOfDatasetsTtJets_[i].NormFactor()); return a;};
  vector<Dataset> ListOfDatasetsBkg() const {return listOfDatasetsBkg_;};
  vector<string> ListOfDatasetsNameBkg() const {vector<string> a; for(unsigned int i=0;i<listOfDatasetsBkg_.size();i++) a.push_back(listOfDatasetsBkg_[i].Name()); return a;};
  vector<int> ListOfDatasetsColorBkg() const {vector<int> a; for(unsigned int i=0;i<listOfDatasetsBkg_.size();i++) a.push_back(listOfDatasetsBkg_[i].Color()); return a;};
  vector<float> ListOfDatasetsNormFactorBkg() const {vector<float> a; for(unsigned int i=0;i<listOfDatasetsBkg_.size();i++) a.push_back(listOfDatasetsBkg_[i].NormFactor()); return a;};
  vector<Dataset> ListOfDatasetsNP() const {return listOfDatasetsNP_;};
  vector<string> ListOfDatasetsNameNP() const {vector<string> a; for(unsigned int i=0;i<listOfDatasetsNP_.size();i++) a.push_back(listOfDatasetsNP_[i].Name()); return a;};
  vector<int> ListOfDatasetsColorNP() const {vector<int> a; for(unsigned int i=0;i<listOfDatasetsNP_.size();i++) a.push_back(listOfDatasetsNP_[i].Color()); return a;};
  vector<float> ListOfDatasetsNormFactorNP() const {vector<float> a; for(unsigned int i=0;i<listOfDatasetsNP_.size();i++) a.push_back(listOfDatasetsNP_[i].NormFactor()); return a;};
  
  void SetListOfDatasets( vector<Dataset> listOfDatasetsBkg, vector<Dataset> listOfDatasetsTtJets, vector<Dataset> listOfDatasetsNP);

  ////////////////////////////
  //          NR            //
  ////////////////////////////
  
  std::pair<float,float> NRegion() const { return NRegion_;};
  void DefineNR(std::pair<float,float> NRegion) {NRegion_ = NRegion;};
	
		//private methods
 private:
  void DrawHisto(Dataset dataset, TH1F* histo);

		//private data
 private:

  bool isMC_;
  bool scaled_;
  bool normalized_;

  string ObsName_;
  string AxisLabel_;
  int NbinsX_;
  bool binningFixed_;
  float* binsX_;
  float Xmin_;
  float Xmax_;

  vector<Dataset> listOfDatasets_;
  vector<Dataset> listOfDatasetsTtJets_;
  vector<Dataset> listOfDatasetsBkg_;
  vector<Dataset> listOfDatasetsNP_;

  std::pair<float,float> NRegion_;
  
  TH1F* h_ObservableCR_;
  TH1F* h_ObservableSR_;
  TH1F* h_ObservableEstimated_;
  TH1F* h_ObservableExpected_; // no NP
  TH1F* h_ObservableObserved_;
  
  TH1F* h_SysError_;
 
  //inclusive
  vector<TH1F*> h_ObservableByDatasetsTtJets_;
  vector<TH1F*> h_ObservableByDatasetsBkg_;
  vector<TH1F*> h_ObservableByDatasetsNP_;
  //CR
  vector<TH1F*> h_ObservableCRByDatasetsTtJets_;
  vector<TH1F*> h_ObservableCRByDatasetsBkg_;
  vector<TH1F*> h_ObservableCRByDatasetsNP_;
  //SR
  vector<TH1F*> h_ObservableSRByDatasetsTtJets_;
  vector<TH1F*> h_ObservableSRByDatasetsBkg_;
  vector<TH1F*> h_ObservableSRByDatasetsNP_;

  vector<TH1F*> h_SignEstimated_;
  vector<TH1F*> h_SignExpected_;


  TLine* line;
  TLegend* leg;
  TLegend* legSign;
  TLegend* legSignRatio;
  TLegend* leg_LMXContamination;
  TF1 *func;
  TH1F* plotRatioCR_SR; 
  TH1F* plotRatioEstim_Exp; 
  TCanvas* canvasSR;
  TCanvas* canvasCR;
  TCanvas* canvasRatioCR_SR;
  TCanvas* canvasRatioEstim_Exp;
  TCanvas* canvasSign;
  TCanvas* canvasSignRatio;
  TCanvas* canvasLMXContamination;

};

#endif
