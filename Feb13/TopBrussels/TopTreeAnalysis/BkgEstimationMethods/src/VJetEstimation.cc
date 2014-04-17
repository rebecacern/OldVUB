#include "../interface/VJetEstimation.h"

ClassImp(VJetEstimation)

VJetEstimation::VJetEstimation(){
  Njets_ = NULL;
	BtagWorkingPoint_ = NULL;
  
  Ntt_          = NULL;
  Ntt_err_      = NULL;
  Ntt_err_up_   = NULL;
  Ntt_err_down_ = NULL;
	Nv_           = NULL;
	Nv_err_       = NULL;
  Nv_err_up_    = NULL;
  Nv_err_down_  = NULL;
	Nvb_          = NULL;
	Nvb_err_      = NULL;
  
	ebq_  = NULL;
	e0bq_ = NULL;
	e1bq_ = NULL;
	e2bq_ = NULL;
  
  minValue_ = NULL;
  
  MyLeg = NULL;
}

/**_________________________________________________________________________________________________________________________________________________*/
VJetEstimation::VJetEstimation(UInt_t NofBtagWorkingPoint, Float_t* BtagWorkingPoint, UInt_t NofJets, UInt_t NofJetBins, std::vector<Dataset> datasets, std::vector<std::string> ttLikeDatasetNames, std::vector<std::string> vLikeDatasetNames, std::vector<std::string> vbLikeDatasetNames):TObject(){
	cout<<"Object from the class VJetEstimation being instantiated"<<endl;
  vDatasets_ = datasets;
	NbOfDatasets_ = datasets.size();
	NbOfBtagWorkingPoint_ = NofBtagWorkingPoint;
	BtagWorkingPoint_ = new Float_t[NbOfBtagWorkingPoint_];
	for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++) BtagWorkingPoint_[i] = BtagWorkingPoint[i];
	MCdata_ = true;
	iDatasetsTTLike_.clear();
	iDatasetsVLike_.clear();
	iDatasetsVbLike_.clear();
  this->SetProcesses(std::vector<Bool_t>(NbOfDatasets_, kTRUE), ttLikeDatasetNames, vLikeDatasetNames, vbLikeDatasetNames);	
	NbOfJetsBins_ = NofJetBins;
	Njets_ = new UInt_t[NbOfJetsBins_];
	for(UInt_t i=0;i<NbOfJetsBins_;i++) Njets_[i] = NofJets+i;
	NbOfBJetsBins_ = 4;
  
	init_Nttlike_ = vector<Double_t>(NbOfJetsBins_,0.);
	init_Nvlike_  = vector<Double_t>(NbOfJetsBins_,0.);
	init_Eb_      = vector< vector<Double_t> >(NbOfJetsBins_,vector<Double_t>(NbOfBtagWorkingPoint_,0.) );
	init_Eudsc_   = vector< vector<Double_t> >(NbOfJetsBins_,vector<Double_t>(NbOfBtagWorkingPoint_,0.) );
	init_Euds_    = vector< vector<Double_t> >(NbOfJetsBins_,vector<Double_t>(NbOfBtagWorkingPoint_,0.) );
  
	minValue_ = new Double_t[NbOfJetsBins_];
	for(UInt_t i=0;i<NbOfJetsBins_;i++) minValue_[i] = 0;
	
  RescaledTTLikeEstimation = 0;
  RescaledVLikeEstimation = 0;
  RescaledVbLikeEstimation = 0;
  tCanva_RescaledTTLikeEstimation = 0;
  tCanva_RescaledVLikeEstimation = 0;
  tCanva_RescaledVbLikeEstimation = 0;
  
	
	hNbjets_mc_       = vector< vector< vector<TH1F> > >(NbOfBtagWorkingPoint_,vector< vector<TH1F> >(NbOfJetsBins_,vector< TH1F >(3)));
	hNbjets_pdf_mc_   = vector< vector< vector<TH1F> > >(NbOfBtagWorkingPoint_,vector< vector<TH1F> >(NbOfJetsBins_,vector< TH1F >(3)));
	hNbjets_pdf_est_  = vector< vector< vector<TH1F> > >(NbOfBtagWorkingPoint_,vector< vector<TH1F> >(NbOfJetsBins_,vector< TH1F >(3)));
	
	hNbjetsEstSummary = vector< vector< vector<TH1F> > >(NbOfBtagWorkingPoint_,vector< vector<TH1F> >(NbOfJetsBins_+1,vector< TH1F >(3)));
	hNbjetsMCSummary  = vector< vector< vector<TH1F> > >(NbOfBtagWorkingPoint_,vector< vector<TH1F> >(NbOfJetsBins_+1,vector< TH1F >(3)));
	hNjetsEstSummary  = vector< vector< vector<TH1F> > >(NbOfBtagWorkingPoint_,vector< vector<TH1F> >(NbOfBJetsBins_+1,vector< TH1F >(3)));
	hNjetsMCSummary   = vector< vector< vector<TH1F> > >(NbOfBtagWorkingPoint_,vector< vector<TH1F> >(NbOfBJetsBins_+1,vector< TH1F >(3)));
	
	char name[100];
	char title[100];
  
	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
    
      //hNbjets_mc_[i]       = new TH1F**[NbOfJetsBins_];
      //hNbjets_pdf_mc_[i]   = new TH1F**[NbOfJetsBins_];
      //hNbjets_pdf_est_[i]  = new TH1F**[NbOfJetsBins_];
		
      //hNbjetsEstSummary[i] = new TH1F**[NbOfJetsBins_+1];
      //hNbjetsMCSummary[i]  = new TH1F**[NbOfJetsBins_+1];
      //hNjetsEstSummary[i]  = new TH1F**[NbOfBJetsBins_+1];
      //hNjetsMCSummary[i]   = new TH1F**[NbOfBJetsBins_+1];
		
      //tCanva_Nbjets_Summary[i] = new TCanvas*[NbOfJetsBins_+1];
      //tCanva_Njets_Summary[i]  = new TCanvas*[NbOfBJetsBins_+1];
    
		for(UInt_t j=0;j<NbOfJetsBins_;j++){
        // Booking histograms for VLike,VbLike and TTlike predictions/estimations
        // for different jet multiplicity
      
        // estimation
      
        //hNbjets_mc_[i][j]          = new TH1F*[3];
			sprintf(name,"hNbjets_mc_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][0]       = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][0].Sumw2();
			sprintf(name,"hNbjets_mc_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][1]       = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][1].Sumw2();
			sprintf(name,"hNbjets_mc_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][2]       = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][2].Sumw2();
      
        //hNbjets_pdf_mc_[i][j]      = new TH1F*[3];
			sprintf(name,"hNbjets_pdf_mc_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][0]   = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][0].Sumw2();
			sprintf(name,"hNbjets_pdf_mc_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][1]   = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][1].Sumw2();
			sprintf(name,"hNbjets_pdf_mc_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][2]   = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][2].Sumw2();
      
        //hNbjets_pdf_est_[i][j]     = new TH1F*[3];
			sprintf(name,"hNbjets_pdf_est_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][0]  = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][0].Sumw2();
			sprintf(name,"hNbjets_pdf_est_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][1]  = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][1].Sumw2();
			sprintf(name,"hNbjets_pdf_est_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][2]  = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][2].Sumw2();
      
        //hNbjetsEstSummary[i][j]    = new TH1F*[3];
			sprintf(name,"hNbjetsEstSummary_VLike_%d_jets_wp_%d",Njets_[j],i);
        //sprintf(title,"V-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][0] = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsEstSummary_VbLike_%d_jets_wp_%d",Njets_[j],i);
        //sprintf(title,"Vb-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][1] = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsEstSummary_TTlike_%d_jets_wp_%d",Njets_[j],i);
        //sprintf(title,"TT-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][2] = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);
      
        // MC prediction
        //hNbjetsMCSummary[i][j]     = new TH1F*[3];
			sprintf(name,"hNbjetsMCSummary_VLike_%d_jets_wp_%d",Njets_[j],i);
        //sprintf(title,"V-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][0]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsMCSummary_VbLike_%d_jets_wp_%d",Njets_[j],i);
        //sprintf(title,"Vb-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][1]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsMCSummary_TTlike_%d_jets_wp_%d",Njets_[j],i);
        //sprintf(title,"TT-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][2]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
      
		}
    
      // Estimation summary (inclusive)
      //hNbjetsEstSummary[i][NbOfJetsBins_]    = new TH1F*[3];
		sprintf(name ,"hNbjetsEstSummary_VLike_Inclusive_wp_%d",i);
		sprintf(title,"V-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][0] = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsEstSummary_VbLike_Inclusive_wp_%d",i);
		sprintf(title,"Vb-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][1] = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsEstSummary_TTlike_Inclusive_wp_%d",i);
		sprintf(title,"TT-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][2] = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
    
      // MC prediction summary (inclusive)
      //hNbjetsMCSummary[i][NbOfJetsBins_]     = new TH1F*[3];
		sprintf(name ,"hNbjetsMCSummary_VLike_Inclusive_wp_%d",i);
		sprintf(title,"V-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][0]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsMCSummary_VbLike_Inclusive_wp_%d",i);
		sprintf(title,"Vb-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][1]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsMCSummary_TTlike_Inclusive_wp_%d",i);
		sprintf(title,"TT-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][2]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
    
		for(UInt_t j=0;j<NbOfBJetsBins_;j++){
        // Booking histograms for VLike,VbLike and TTlike predictions/estimations
        // for different b-jet multiplicity
      
        // estimation
        ///      			hsNjets_MC[i][j]          = NULL;
        ///      			hsNjets_Est[i][j]         = NULL;
        //hNjetsEstSummary[i][j]    = new TH1F*[3];
			sprintf(name,"hNjetsEstSummary_VLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"V-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][0] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsEstSummary_VbLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"Vb-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][1] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsEstSummary_TTlike_%d_bjets_wp_%d",j,i);
			sprintf(title,"TT-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][2] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
      
        // MC prediction
        //hNjetsMCSummary[i][j]     = new TH1F*[3];
			sprintf(name,"hNjetsMCSummary_VLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"V-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][0]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsMCSummary_VbLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"Vb-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][1]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsMCSummary_TTlike_%d_bjets_wp_%d",j,i);
			sprintf(title,"TT-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][2]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		}
    
      // Estimation summary (b-inclusive)
      ///    		hsNjets_MC[i][NbOfBJetsBins_]          = NULL;
      ///    		hsNjets_Est[i][NbOfBJetsBins_]         = NULL;
      //hNjetsEstSummary[i][NbOfBJetsBins_]    = new TH1F*[3];
		sprintf(name ,"hNjetsEstSummary_VLike_bInclusive_wp_%d",i);
		sprintf(title,"V-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][0] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsEstSummary_VbLike_bInclusive_wp_%d",i);
		sprintf(title,"Vb-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][1] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsEstSummary_TTlike_bInclusive_wp_%d",i);
		sprintf(title,"TT-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][2] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
    
      // MC prediction summary (b-inclusive)
      //hNjetsMCSummary[i][NbOfBJetsBins_]     = new TH1F*[3];
		sprintf(name ,"hNjetsMCSummary_VLike_bInclusive_wp_%d",i);
		sprintf(title,"V-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][0]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsMCSummary_VbLike_bInclusive_wp_%d",i);
		sprintf(title,"Vb-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][1]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsMCSummary_TTlike_bInclusive_wp_%d",i);
		sprintf(title,"TT-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][2]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
	}
	cout<<" -- Summary histograms correctly instantiated"<<endl;
  
    //init efficiency estimators
  eudsc_        = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eudsc_err_    = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eudsc_err_up_ = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eudsc_err_down_ = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eudsc_mc_     = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eudsc_err_mc_ = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_         = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_err_     = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_err_up_  = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_err_down_= vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_mc_      = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_err_mc_  = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_           = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_err_       = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_err_up_    = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_err_down_  = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_mc_        = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_err_mc_    = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
	cout<<" -- Estimators correctly instantiated"<<endl;
  
    //init containers
	N_ = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_,vector< vector< vector< Double_t > > >(NbOfJetsBins_,vector< vector< Double_t > >(NbOfBJetsBins_,vector< Double_t >(NbOfDatasets_,0))));//           = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_);
	N_err_ = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_,vector< vector< vector< Double_t > > >(NbOfJetsBins_,vector< vector< Double_t > >(NbOfBJetsBins_,vector< Double_t >(NbOfDatasets_,0))));//        = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_);
	N_err_hi_ = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_,vector< vector< vector< Double_t > > >(NbOfJetsBins_,vector< vector< Double_t > >(NbOfBJetsBins_,vector< Double_t >(NbOfDatasets_,0))));//     = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_);
	N_err_lo_ = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_,vector< vector< vector< Double_t > > >(NbOfJetsBins_,vector< vector< Double_t > >(NbOfBJetsBins_,vector< Double_t >(NbOfDatasets_,0))));//     = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_);
	Nbjets_ = vector< vector< vector< Double_t > > > (NbOfBtagWorkingPoint_,vector< vector< Double_t > >(NbOfJetsBins_,vector< Double_t >(NbOfBJetsBins_,0)));//       = vector< vector< vector< Double_t > > >(NbOfBtagWorkingPoint_);
	MultiJet_Est_ = vector< vector< vector< Double_t > > > (NbOfBtagWorkingPoint_,vector< vector< Double_t > >(NbOfJetsBins_,vector< Double_t >(NbOfBJetsBins_,0)));// = vector< vector< vector< Double_t > > >(NbOfBtagWorkingPoint_);
	cout<<" -- Containers correctly instantiated"<<endl;
  /*
   for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
   
   eudsc_[i]        = new Double_t[NbOfJetsBins_];
   eudsc_err_[i]    = new Double_t[NbOfJetsBins_];
   eudsc_mc_[i]     = new Double_t[NbOfJetsBins_];
   eudsc_err_mc_[i] = new Double_t[NbOfJetsBins_];
   euds_[i]         = new Double_t[NbOfJetsBins_];
   euds_err_[i]     = new Double_t[NbOfJetsBins_];
   euds_mc_[i]      = new Double_t[NbOfJetsBins_];
   euds_err_mc_[i]  = new Double_t[NbOfJetsBins_];
   eb_[i]           = new Double_t[NbOfJetsBins_];
   eb_err_[i]       = new Double_t[NbOfJetsBins_];
   eb_mc_[i]        = new Double_t[NbOfJetsBins_];
   eb_err_mc_[i]    = new Double_t[NbOfJetsBins_];
   
   N_[i].resize(NbOfJetsBins_);
   N_err_[i].resize(NbOfJetsBins_);
   N_err_hi_[i].resize(NbOfJetsBins_);
   N_err_lo_[i].resize(NbOfJetsBins_);
   Nbjets_[i].resize(NbOfJetsBins_);
   MultiJet_Est_[i].resize(NbOfJetsBins_);
   
   for(UInt_t j=0;j<NbOfJetsBins_;j++){
   
   eudsc_[i][j]       = 0;
   eudsc_err_[i][j]   = 0;
   eudsc_mc_[i][j]    = 0;
   eudsc_err_mc_[i][j]= 0;
   euds_[i][j]        = 0;
   euds_err_[i][j]    = 0;
   euds_mc_[i][j]     = 0;
   euds_err_mc_[i][j] = 0;
   eb_[i][j]          = 0;
   eb_err_[i][j]      = 0;
   eb_mc_[i][j]       = 0;
   eb_err_mc_[i][j]   = 0;
   
   N_[i][j].resize(NbOfBJetsBins_);
   N_err_[i][j].resize(NbOfBJetsBins_);
   N_err_hi_[i][j].resize(NbOfBJetsBins_);
   N_err_lo_[i][j].resize(NbOfBJetsBins_);
   Nbjets_[i][j].resize(NbOfBJetsBins_);
   MultiJet_Est_[i][j].resize(NbOfBJetsBins_);
   
   for(UInt_t k=0;k<NbOfBJetsBins_;k++){
   N_[i][j][k].resize(NbOfDatasets_);
   N_err_[i][j][k].resize(NbOfDatasets_);
   N_err_hi_[i][j][k].resize(NbOfDatasets_);
   N_err_lo_[i][j][k].resize(NbOfDatasets_);
   Nbjets_[i][j][k]       = 0;
   MultiJet_Est_[i][j][k] = 0;
   
   for(UInt_t l=0;l<NbOfDatasets_;l++){
   N_[i][j][k][l]     = 0;
   N_err_[i][j][k][l] = 0;
   N_err_hi_[i][j][k][l] = 0;
   N_err_lo_[i][j][k][l] = 0;
   }
   }
   
   }
   
   }
   */
    //init histograms used to calculate the b/mis-tagging efficiencies.
    //	hNbOfBGenJets_                                 = new TH1F**[NbOfDatasets_];
	hNbOfBGenJets_                                 = vector< vector < TH1F > >(NbOfDatasets_,vector < TH1F >(NbOfJetsBins_) );
	hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_     = vector< TH3F >(NbOfDatasets_);
	hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_     = vector< TH3F >(NbOfDatasets_);
	hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_   = vector< TH3F >(NbOfDatasets_);
	hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_  = vector< TH3F >(NbOfDatasets_);
  
	hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_    = vector< vector < TH3F > >(NbOfDatasets_,vector < TH3F >(NbOfBtagWorkingPoint_) );
	hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_    = vector< vector < TH3F > >(NbOfDatasets_,vector < TH3F >(NbOfBtagWorkingPoint_) );
	hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_  = vector< vector < TH3F > >(NbOfDatasets_,vector < TH3F >(NbOfBtagWorkingPoint_) );
	hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_ = vector< vector < TH3F > >(NbOfDatasets_,vector < TH3F >(NbOfBtagWorkingPoint_) );
  
  bTagEff_vs_Njets_   = vector< vector < TEfficiency* > >(NbOfDatasets_,vector < TEfficiency* >(NbOfBtagWorkingPoint_,0) );
  cTagEff_vs_Njets_   = vector< vector < TEfficiency* > >(NbOfDatasets_,vector < TEfficiency* >(NbOfBtagWorkingPoint_,0) );
  udsTagEff_vs_Njets_ = vector< vector < TEfficiency* > >(NbOfDatasets_,vector < TEfficiency* >(NbOfBtagWorkingPoint_,0) );
  misTagEff_vs_Njets_ = vector< vector < TEfficiency* > >(NbOfDatasets_,vector < TEfficiency* >(NbOfBtagWorkingPoint_,0) );
  
  bTagEff_vs_Njets_TTlike_   = vector< TGraphAsymmErrors* >(NbOfBtagWorkingPoint_,0);
  misTagEff_vs_Njets_TTlike_ = vector< TGraphAsymmErrors* >(NbOfBtagWorkingPoint_,0);
  misTagEff_vs_Njets_Vlike_  = vector< TGraphAsymmErrors* >(NbOfBtagWorkingPoint_,0);
  
	const Int_t JetPtNbins = 4;             Double_t JetPtBins[JetPtNbins+1]   = {30,50,80,120,250};
	const Int_t JetEtaNbins= 4;             Double_t JetEtaBins[JetEtaNbins+1] = {0,0.4,0.8,1.3,2.4};
	const Int_t JetNbins   = NbOfJetsBins_; Double_t JetNBins[JetNbins+1]; for(Int_t i=0;i<JetNbins+1;i++) JetNBins[i] = Njets_[0]+i;
  
	for(UInt_t i=0;i<NbOfDatasets_;i++){
    
      //		hNbOfBGenJets_[i] = new TH1F*[NbOfJetsBins_];
		for(UInt_t j=0; j<NbOfJetsBins_;j++){
			sprintf(name ,"hNbOfBGenJets_%d_%d",i,j);
			hNbOfBGenJets_[i][j] = TH1F(name,"",4,0,4);
			hNbOfBGenJets_[i][j].Sumw2();
			hNbOfBGenJets_[i][j].GetXaxis()->CenterLabels();
			hNbOfBGenJets_[i][j].GetXaxis()->SetTitle("Nb. of b-quark jets");
		}
		sprintf(name ,"hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->SetTitle("Nb. of jets");
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->CenterLabels();
    
		sprintf(name ,"hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->SetTitle("Nb. of jets");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->CenterLabels();
    
		sprintf(name ,"hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->SetTitle("Nb. of jets");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->CenterLabels();
    
		sprintf(name ,"hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->SetTitle("Nb. of jets");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->CenterLabels();
    
      //		bTagEff_vs_Njets_[i]   = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
      //		cTagEff_vs_Njets_[i]   = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
      //		udsTagEff_vs_Njets_[i] = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
      //		misTagEff_vs_Njets_[i] = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
    
      //		hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]    = new TH3F*[NbOfBtagWorkingPoint_];
      //		hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]    = new TH3F*[NbOfBtagWorkingPoint_];
      //		hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]  = new TH3F*[NbOfBtagWorkingPoint_];
      //		hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F*[NbOfBtagWorkingPoint_];
    
		for(UInt_t j=0; j<NbOfBtagWorkingPoint_;j++){
			sprintf(name ,"hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->SetTitle("Nb. of jets");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->CenterLabels();
      
			sprintf(name ,"hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->SetTitle("Nb. of jets");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->CenterLabels();
      
			sprintf(name ,"hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->SetTitle("Nb. of jets");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->CenterLabels();
      
			sprintf(name ,"hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->SetTitle("Nb. of jets");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->CenterLabels();
      
        //			bTagEff_vs_Njets_[i][j]   = 0;
        //			cTagEff_vs_Njets_[i][j]   = 0;
        //			udsTagEff_vs_Njets_[i][j] = 0;
        //			misTagEff_vs_Njets_[i][j] = 0;
		}
	}
  
    //init efficiency estimators
  ebq_   = new Double_t[NbOfJetsBins_];
  e0bq_  = new Double_t[NbOfJetsBins_];
  e1bq_  = new Double_t[NbOfJetsBins_];
  e2bq_  = new Double_t[NbOfJetsBins_];
    //init Ntt/Nvb/Nv estimators
  Ntt_          = new Double_t[NbOfJetsBins_];
  Ntt_err_      = new Double_t[NbOfJetsBins_];
  Ntt_err_up_   = new Double_t[NbOfJetsBins_];
  Ntt_err_down_ = new Double_t[NbOfJetsBins_];
  Nv_           = new Double_t[NbOfJetsBins_];
  Nv_err_       = new Double_t[NbOfJetsBins_];
  Nv_err_up_    = new Double_t[NbOfJetsBins_];
  Nv_err_down_  = new Double_t[NbOfJetsBins_];
  Nvb_          = new Double_t[NbOfJetsBins_];
  Nvb_err_      = new Double_t[NbOfJetsBins_];
	for(UInt_t i=0;i<NbOfJetsBins_;i++){
		ebq_[i]      = 0;
		e0bq_[i]     = 0;
		e1bq_[i]     = 0;
		e2bq_[i]     = 0;
		Ntt_[i]          = 0;
		Ntt_err_[i]      = 0;
		Ntt_err_up_[i]   = 0;
		Ntt_err_down_[i] = 0;
		Nv_[i]           = 0;
		Nv_err_[i]       = 0;
		Nv_err_up_[i]    = 0;
		Nv_err_down_[i]  = 0;
		Nvb_[i]          = 0;
		Nvb_err_[i]      = 0;
	}
	cout<<"Object from the class VJetEstimation correctly instantiated"<<endl;
}

/**________________________________________________________________________________________________________________*/
VJetEstimation::VJetEstimation(UInt_t NofBtagWorkingPoint, Float_t* BtagWorkingPoint, UInt_t NofJets, UInt_t NofJetBins, Double_t** EffXbq, UInt_t NofDatasets, vector<Int_t> iDTTLike, vector<Int_t> iDVLike, vector<Int_t> iDVbLike):TObject(){
	cout<<"Object from the class VJetEstimation being instantiated"<<endl;
	NbOfDatasets_ = NofDatasets;
	NbOfBtagWorkingPoint_ = NofBtagWorkingPoint;
	BtagWorkingPoint_ = new Float_t[NbOfBtagWorkingPoint_];
	for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++) BtagWorkingPoint_[i] = BtagWorkingPoint[i];
	MCdata_ = true;
	iDatasetsTTLike_.clear();
	iDatasetsVLike_.clear();
	iDatasetsVbLike_.clear();
	iDatasetsTTLike_ = iDTTLike;
	iDatasetsVLike_  = iDVLike;
	iDatasetsVbLike_ = iDVbLike;
	
	NbOfJetsBins_ = NofJetBins;
	Njets_ = new UInt_t[NbOfJetsBins_];
	for(UInt_t i=0;i<NbOfJetsBins_;i++) Njets_[i] = NofJets+i;
	NbOfBJetsBins_ = 4;
  
	init_Nttlike_ = vector<Double_t>(NbOfJetsBins_,0.);
	init_Nvlike_  = vector<Double_t>(NbOfJetsBins_,0.);
	init_Eb_      = vector< vector<Double_t> >(NbOfJetsBins_,vector<Double_t>(NbOfBtagWorkingPoint_,0.) );
	init_Eudsc_   = vector< vector<Double_t> >(NbOfJetsBins_,vector<Double_t>(NbOfBtagWorkingPoint_,0.) );
	init_Euds_    = vector< vector<Double_t> >(NbOfJetsBins_,vector<Double_t>(NbOfBtagWorkingPoint_,0.) );
  
	minValue_ = new Double_t[NbOfJetsBins_];
	for(UInt_t i=0;i<NbOfJetsBins_;i++) minValue_[i] = 0;
	
  RescaledTTLikeEstimation = 0;
  RescaledVLikeEstimation = 0;
  RescaledVbLikeEstimation = 0;
  tCanva_RescaledTTLikeEstimation = 0;
  tCanva_RescaledVLikeEstimation = 0;
  tCanva_RescaledVbLikeEstimation = 0;
  
	hNbjets_mc_       = vector< vector< vector< TH1F > > >(NbOfBtagWorkingPoint_,vector< vector< TH1F > >(NbOfJetsBins_,vector< TH1F >(3)));
	hNbjets_pdf_mc_   = vector< vector< vector< TH1F > > >(NbOfBtagWorkingPoint_,vector< vector< TH1F > >(NbOfJetsBins_,vector< TH1F >(3)));
	hNbjets_pdf_est_  = vector< vector< vector< TH1F > > >(NbOfBtagWorkingPoint_,vector< vector< TH1F > >(NbOfJetsBins_,vector< TH1F >(3)));
	
	hNbjetsEstSummary = vector< vector< vector< TH1F > > >(NbOfBtagWorkingPoint_,vector< vector< TH1F > >(NbOfJetsBins_+1,vector< TH1F >(3)));
	hNbjetsMCSummary  = vector< vector< vector< TH1F > > >(NbOfBtagWorkingPoint_,vector< vector< TH1F > >(NbOfJetsBins_+1,vector< TH1F >(3)));
	hNjetsEstSummary  = vector< vector< vector< TH1F > > >(NbOfBtagWorkingPoint_,vector< vector< TH1F > >(NbOfBJetsBins_+1,vector< TH1F >(3)));
	hNjetsMCSummary   = vector< vector< vector< TH1F > > >(NbOfBtagWorkingPoint_,vector< vector< TH1F > >(NbOfBJetsBins_+1,vector< TH1F >(3)));
	
	char name[100];
	char title[100];
  
	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
    
      //hNbjets_mc_[i]       = new TH1F**[NbOfJetsBins_];
      //hNbjets_pdf_mc_[i]   = new TH1F**[NbOfJetsBins_];
      //hNbjets_pdf_est_[i]  = new TH1F**[NbOfJetsBins_];
		
      //hNbjetsEstSummary[i] = new TH1F**[NbOfJetsBins_+1];
      //hNbjetsMCSummary[i]  = new TH1F**[NbOfJetsBins_+1];
      //hNjetsEstSummary[i]  = new TH1F**[NbOfBJetsBins_+1];
      //hNjetsMCSummary[i]   = new TH1F**[NbOfBJetsBins_+1];
		
      //tCanva_Nbjets_Summary[i] = new TCanvas*[NbOfJetsBins_+1];
      //tCanva_Njets_Summary[i]  = new TCanvas*[NbOfBJetsBins_+1];
    
		for(UInt_t j=0;j<NbOfJetsBins_;j++){
        // Booking histograms for VLike,VbLike and TTlike predictions/estimations
        // for different jet multiplicity
      
      
        //hNbjets_mc_[i][j]          = new TH1F*[3];
			sprintf(name,"hNbjets_mc_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][0]       = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][0].Sumw2();
			sprintf(name,"hNbjets_mc_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][1]       = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][1].Sumw2();
			sprintf(name,"hNbjets_mc_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][2]       = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][2].Sumw2();
      
        //hNbjets_pdf_mc_[i][j]      = new TH1F*[3];
			sprintf(name,"hNbjets_pdf_mc_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][0]   = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][0].Sumw2();
			sprintf(name,"hNbjets_pdf_mc_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][1]   = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][1].Sumw2();
			sprintf(name,"hNbjets_pdf_mc_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][2]   = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][2].Sumw2();
      
        //hNbjets_pdf_est_[i][j]     = new TH1F*[3];
			sprintf(name,"hNbjets_pdf_est_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][0]  = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][0].Sumw2();
			sprintf(name,"hNbjets_pdf_est_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][1]  = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][1].Sumw2();
			sprintf(name,"hNbjets_pdf_est_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][2]  = TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][2].Sumw2();
      
        //hNbjetsEstSummary[i][j]    = new TH1F*[3];
			sprintf(name,"hNbjetsEstSummary_VLike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"V-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][0] = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsEstSummary_VbLike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"Vb-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][1] = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsEstSummary_TTlike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"TT-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][2] = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
      
        // MC prediction
        //hNbjetsMCSummary[i][j]     = new TH1F*[3];
			sprintf(name,"hNbjetsMCSummary_VLike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"V-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][0]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsMCSummary_VbLike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"Vb-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][1]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsMCSummary_TTlike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"TT-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][2]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		}
    
      // Estimation summary (inclusive)
      //hNbjetsEstSummary[i][NbOfJetsBins_]    = new TH1F*[3];
		sprintf(name ,"hNbjetsEstSummary_VLike_Inclusive_wp_%d",i);
		sprintf(title,"V-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][0] = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsEstSummary_VbLike_Inclusive_wp_%d",i);
		sprintf(title,"Vb-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][1] = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsEstSummary_TTlike_Inclusive_wp_%d",i);
		sprintf(title,"TT-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][2] = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
    
      // MC prediction summary (inclusive)
      //hNbjetsMCSummary[i][NbOfJetsBins_]     = new TH1F*[3];
		sprintf(name ,"hNbjetsMCSummary_VLike_Inclusive_wp_%d",i);
		sprintf(title,"V-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][0]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsMCSummary_VbLike_Inclusive_wp_%d",i);
		sprintf(title,"Vb-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][1]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsMCSummary_TTlike_Inclusive_wp_%d",i);
		sprintf(title,"TT-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][2]  = TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
    
		for(UInt_t j=0;j<NbOfBJetsBins_;j++){
        // Booking histograms for VLike,VbLike and TTlike predictions/estimations
        // for different b-jet multiplicity
      
        // estimation
        ///      			hsNjets_MC[i][j]          = NULL;
        ///      			hsNjets_Est[i][j]         = NULL;
        //hNjetsEstSummary[i][j]    = new TH1F*[3];
			sprintf(name,"hNjetsEstSummary_VLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"V-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][0] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsEstSummary_VbLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"Vb-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][1] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsEstSummary_TTlike_%d_bjets_wp_%d",j,i);
			sprintf(title,"TT-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][2] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
      
        // MC prediction
        //hNjetsMCSummary[i][j]     = new TH1F*[3];
			sprintf(name,"hNjetsMCSummary_VLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"V-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][0]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsMCSummary_VbLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"Vb-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][1]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsMCSummary_TTlike_%d_bjets_wp_%d",j,i);
			sprintf(title,"TT-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][2]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		}
    
      // Estimation summary (b-inclusive)
      ///    		hsNjets_MC[i][NbOfBJetsBins_]          = NULL;
      ///    		hsNjets_Est[i][NbOfBJetsBins_]         = NULL;
      //hNjetsEstSummary[i][NbOfBJetsBins_]    = new TH1F*[3];
		sprintf(name ,"hNjetsEstSummary_VLike_bInclusive_wp_%d",i);
		sprintf(title,"V-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][0] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsEstSummary_VbLike_bInclusive_wp_%d",i);
		sprintf(title,"Vb-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][1] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsEstSummary_TTlike_bInclusive_wp_%d",i);
		sprintf(title,"TT-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][2] = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
    
      // MC prediction summary (b-inclusive)
      //hNjetsMCSummary[i][NbOfBJetsBins_]     = new TH1F*[3];
		sprintf(name ,"hNjetsMCSummary_VLike_bInclusive_wp_%d",i);
		sprintf(title,"V-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][0]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsMCSummary_VbLike_bInclusive_wp_%d",i);
		sprintf(title,"Vb-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][1]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsMCSummary_TTlike_bInclusive_wp_%d",i);
		sprintf(title,"TT-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][2]  = TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
	}
	cout<<" -- Summary histograms correctly instantiated"<<endl;
  
    //init efficiency estimators
  eudsc_        = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eudsc_err_    = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eudsc_err_up_ = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eudsc_err_down_ = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eudsc_mc_     = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eudsc_err_mc_ = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_         = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_err_     = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_err_up_  = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_err_down_= vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_mc_      = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  euds_err_mc_  = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_           = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_err_       = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_err_up_    = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_err_down_  = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_mc_        = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
  eb_err_mc_    = vector<vector<Double_t> >(NbOfBtagWorkingPoint_,vector<Double_t>(NbOfJetsBins_,0));
	cout<<" -- Estimators correctly instantiated"<<endl;
  
    //init containers
	N_ = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_,vector< vector< vector< Double_t > > >(NbOfJetsBins_,vector< vector< Double_t > >(NbOfBJetsBins_,vector< Double_t >(NbOfDatasets_,0))));//           = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_);
	N_err_ = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_,vector< vector< vector< Double_t > > >(NbOfJetsBins_,vector< vector< Double_t > >(NbOfBJetsBins_,vector< Double_t >(NbOfDatasets_,0))));//        = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_);
	N_err_hi_ = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_,vector< vector< vector< Double_t > > >(NbOfJetsBins_,vector< vector< Double_t > >(NbOfBJetsBins_,vector< Double_t >(NbOfDatasets_,0))));//     = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_);
	N_err_lo_ = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_,vector< vector< vector< Double_t > > >(NbOfJetsBins_,vector< vector< Double_t > >(NbOfBJetsBins_,vector< Double_t >(NbOfDatasets_,0))));//     = vector< vector< vector< vector< Double_t > > > >(NbOfBtagWorkingPoint_);
	Nbjets_ = vector< vector< vector< Double_t > > > (NbOfBtagWorkingPoint_,vector< vector< Double_t > >(NbOfJetsBins_,vector< Double_t >(NbOfBJetsBins_,0)));//       = vector< vector< vector< Double_t > > >(NbOfBtagWorkingPoint_);
	MultiJet_Est_ = vector< vector< vector< Double_t > > > (NbOfBtagWorkingPoint_,vector< vector< Double_t > >(NbOfJetsBins_,vector< Double_t >(NbOfBJetsBins_,0)));// = vector< vector< vector< Double_t > > >(NbOfBtagWorkingPoint_);
	cout<<" -- Containers correctly instantiated"<<endl;
  /*
   
   for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
   
   eudsc_[i]        = new Double_t[NbOfJetsBins_];
   eudsc_err_[i]    = new Double_t[NbOfJetsBins_];
   eudsc_mc_[i]     = new Double_t[NbOfJetsBins_];
   eudsc_err_mc_[i] = new Double_t[NbOfJetsBins_];
   euds_[i]         = new Double_t[NbOfJetsBins_];
   euds_err_[i]     = new Double_t[NbOfJetsBins_];
   euds_mc_[i]      = new Double_t[NbOfJetsBins_];
   euds_err_mc_[i]  = new Double_t[NbOfJetsBins_];
   eb_[i]           = new Double_t[NbOfJetsBins_];
   eb_err_[i]       = new Double_t[NbOfJetsBins_];
   eb_mc_[i]        = new Double_t[NbOfJetsBins_];
   eb_err_mc_[i]    = new Double_t[NbOfJetsBins_];
   
   N_[i].resize(NbOfJetsBins_);
   N_err_[i].resize(NbOfJetsBins_);
   N_err_hi_[i].resize(NbOfJetsBins_);
   N_err_lo_[i].resize(NbOfJetsBins_);
   Nbjets_[i].resize(NbOfJetsBins_);
   MultiJet_Est_[i].resize(NbOfJetsBins_);
   
   for(UInt_t j=0;j<NbOfJetsBins_;j++){
   
   eudsc_[i][j]       = 0;
   eudsc_err_[i][j]   = 0;
   eudsc_mc_[i][j]    = 0;
   eudsc_err_mc_[i][j]= 0;
   euds_[i][j]        = 0;
   euds_err_[i][j]    = 0;
   euds_mc_[i][j]     = 0;
   euds_err_mc_[i][j] = 0;
   eb_[i][j]          = 0;
   eb_err_[i][j]      = 0;
   eb_mc_[i][j]       = 0;
   eb_err_mc_[i][j]   = 0;
   
   N_[i][j].resize(NbOfBJetsBins_);
   N_err_[i][j].resize(NbOfBJetsBins_);
   N_err_hi_[i][j].resize(NbOfBJetsBins_);
   N_err_lo_[i][j].resize(NbOfBJetsBins_);
   Nbjets_[i][j].resize(NbOfBJetsBins_);
   MultiJet_Est_[i][j].resize(NbOfBJetsBins_);
   
   for(UInt_t k=0;k<NbOfBJetsBins_;k++){
   N_[i][j][k].resize(NbOfDatasets_);
   N_err_[i][j][k].resize(NbOfDatasets_);
   N_err_hi_[i][j][k].resize(NbOfDatasets_);
   N_err_lo_[i][j][k].resize(NbOfDatasets_);
   Nbjets_[i][j][k]       = 0;
   MultiJet_Est_[i][j][k] = 0;
   
   for(UInt_t l=0;l<NbOfDatasets_;l++){
   N_[i][j][k][l]     = 0;
   N_err_[i][j][k][l] = 0;
   N_err_hi_[i][j][k][l] = 0;
   N_err_lo_[i][j][k][l] = 0;
   }
   }
   }
   }
   */
    //init histograms used to calculate the b/mis-tagging efficiencies.
    //	hNbOfBGenJets_                                 = new TH1F**[NbOfDatasets_];
	hNbOfBGenJets_                                 = vector< vector < TH1F > >(NbOfDatasets_,vector < TH1F >(NbOfJetsBins_) );
    //	hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_     = new TH3F* [NbOfDatasets_];
	hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_     = vector< TH3F >(NbOfDatasets_);
	hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_     = vector< TH3F >(NbOfDatasets_);
	hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_   = vector< TH3F >(NbOfDatasets_);
	hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_  = vector< TH3F >(NbOfDatasets_);
  
	hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_    = vector< vector < TH3F > >(NbOfDatasets_,vector < TH3F >(NbOfBtagWorkingPoint_) );
	hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_    = vector< vector < TH3F > >(NbOfDatasets_,vector < TH3F >(NbOfBtagWorkingPoint_) );
	hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_  = vector< vector < TH3F > >(NbOfDatasets_,vector < TH3F >(NbOfBtagWorkingPoint_) );
	hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_ = vector< vector < TH3F > >(NbOfDatasets_,vector < TH3F >(NbOfBtagWorkingPoint_) );
  
  bTagEff_vs_Njets_   = vector< vector < TEfficiency* > >(NbOfDatasets_,vector < TEfficiency* >(NbOfBtagWorkingPoint_,0) );
  cTagEff_vs_Njets_   = vector< vector < TEfficiency* > >(NbOfDatasets_,vector < TEfficiency* >(NbOfBtagWorkingPoint_,0) );
  udsTagEff_vs_Njets_ = vector< vector < TEfficiency* > >(NbOfDatasets_,vector < TEfficiency* >(NbOfBtagWorkingPoint_,0) );
  misTagEff_vs_Njets_ = vector< vector < TEfficiency* > >(NbOfDatasets_,vector < TEfficiency* >(NbOfBtagWorkingPoint_,0) );
  
  bTagEff_vs_Njets_TTlike_   = vector< TGraphAsymmErrors* >(NbOfBtagWorkingPoint_,0);
  misTagEff_vs_Njets_TTlike_ = vector< TGraphAsymmErrors* >(NbOfBtagWorkingPoint_,0);
  misTagEff_vs_Njets_Vlike_  = vector< TGraphAsymmErrors* >(NbOfBtagWorkingPoint_,0);
  
	const Int_t JetPtNbins = 4;             Double_t JetPtBins[JetPtNbins+1]   = {30,50,80,120,250};
	const Int_t JetEtaNbins= 4;             Double_t JetEtaBins[JetEtaNbins+1] = {0,0.4,0.8,1.3,2.4};
	const Int_t JetNbins   = NbOfJetsBins_; Double_t JetNBins[JetNbins+1]; for(Int_t i=0;i<JetNbins+1;i++) JetNBins[i] = Njets_[0]+i;
  
	for(UInt_t i=0;i<NbOfDatasets_;i++){
    
      //		hNbOfBGenJets_[i] = new TH1F*[NbOfJetsBins_];
		for(UInt_t j=0; j<NbOfJetsBins_;j++){
			sprintf(name ,"hNbOfBGenJets_%d_%d",i,j);
			hNbOfBGenJets_[i][j] = TH1F(name,"",4,0,4);
			hNbOfBGenJets_[i][j].Sumw2();
			hNbOfBGenJets_[i][j].GetXaxis()->CenterLabels();
			hNbOfBGenJets_[i][j].GetXaxis()->SetTitle("Nb. of b-quark jets");
		}
		sprintf(name ,"hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
    hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
    hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetYaxis()->SetTitle("Jet |#eta|");
    hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->SetTitle("Nb. of jets");
    hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->CenterLabels();
    
		sprintf(name ,"hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->SetTitle("Nb. of jets");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->CenterLabels();
    
		sprintf(name ,"hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->SetTitle("Nb. of jets");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->CenterLabels();
    
		sprintf(name ,"hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->SetTitle("Nb. of jets");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i].GetZaxis()->CenterLabels();
    
      //		bTagEff_vs_Njets_[i]   = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
      //		cTagEff_vs_Njets_[i]   = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
      //		udsTagEff_vs_Njets_[i] = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
      //		misTagEff_vs_Njets_[i] = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
    
      //		hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]    = new TH3F*[NbOfBtagWorkingPoint_];
      //		hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]    = new TH3F*[NbOfBtagWorkingPoint_];
      //		hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]  = new TH3F*[NbOfBtagWorkingPoint_];
      //		hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F*[NbOfBtagWorkingPoint_];
    
		for(UInt_t j=0; j<NbOfBtagWorkingPoint_;j++){
			sprintf(name ,"hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->SetTitle("Nb. of jets");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->CenterLabels();
      
			sprintf(name ,"hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->SetTitle("Nb. of jets");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->CenterLabels();
      
			sprintf(name ,"hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->SetTitle("Nb. of jets");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->CenterLabels();
      
			sprintf(name ,"hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->SetTitle("Nb. of jets");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j].GetZaxis()->CenterLabels();
      
        //			bTagEff_vs_Njets_[i][j]   = 0;
        //			cTagEff_vs_Njets_[i][j]   = 0;
        //			udsTagEff_vs_Njets_[i][j] = 0;
        //			misTagEff_vs_Njets_[i][j] = 0;
		}
	}
  
    //init efficiency estimators
  ebq_   = new Double_t[NbOfJetsBins_];
  e0bq_  = new Double_t[NbOfJetsBins_];//EffXbq[0];
  e1bq_  = new Double_t[NbOfJetsBins_];//EffXbq[1];
  e2bq_  = new Double_t[NbOfJetsBins_];//EffXbq[2];
                                       //init Ntt/Nvb/Nv estimators
  Ntt_          = new Double_t[NbOfJetsBins_];
  Ntt_err_      = new Double_t[NbOfJetsBins_];
  Ntt_err_up_   = new Double_t[NbOfJetsBins_];
  Ntt_err_down_ = new Double_t[NbOfJetsBins_];
  Nv_           = new Double_t[NbOfJetsBins_];
  Nv_err_       = new Double_t[NbOfJetsBins_];
  Nv_err_up_    = new Double_t[NbOfJetsBins_];
  Nv_err_down_  = new Double_t[NbOfJetsBins_];
  Nvb_          = new Double_t[NbOfJetsBins_];
  Nvb_err_      = new Double_t[NbOfJetsBins_];
	for(UInt_t i=0;i<NbOfJetsBins_;i++){
		ebq_[i]      = 0;
    e0bq_[i]     = EffXbq[i][0];
    e1bq_[i]     = EffXbq[i][1];
    e2bq_[i]     = EffXbq[i][2];
		Ntt_[i]          = 0;
		Ntt_err_[i]      = 0;
		Ntt_err_up_[i]   = 0;
		Ntt_err_down_[i] = 0;
		Nv_[i]           = 0;
		Nv_err_[i]       = 0;
		Nv_err_up_[i]    = 0;
		Nv_err_down_[i]  = 0;
		Nvb_[i]          = 0;
		Nvb_err_[i]      = 0;
	}
	cout<<"Object from the class VJetEstimation correctly instantiated"<<endl;
}

/**________________________________________________________________________________________________________________*/
  //copy constructor
VJetEstimation::VJetEstimation(const VJetEstimation& vjet){
	cout<<"Using copy constructor..."<<endl;
}

/**________________________________________________________________________________________________________________*/
VJetEstimation::~VJetEstimation(){
  /*
   delete BtagWorkingPoint_;
   for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
   for(UInt_t j=0;j<NbOfJetsBins_;j++){
   for(UInt_t k=0;k<NbOfBJetsBins_;k++){
   delete [] N_[i][j][k];
   delete [] N_err_[i][j][k];
   }
   delete [] N_[i][j];
   delete [] N_err_[i][j];
   delete [] Nbjets_[i][j];
   }
   delete [] N_[i];
   delete [] N_err_[i];
   delete [] Nbjets_[i];
   }
   delete [] N_;
   delete [] N_err_;
   delete [] Nbjets_;
   
   for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
   delete [] eudsc_[i];
   delete [] euds_[i];
   delete [] eb_[i];
   
   for(UInt_t j=0;j<NbOfJetsBins_;j++){
   for(UInt_t k=0;k<3;k++){
   printf("6.4[%d].1[%d].5[%d]\n", i,j,k);
   delete hNbjetsEstSummary[i][j][k];
   delete hNbjetsMCSummary[i][j][k];
   delete hNbjets_mc_[i][j][k];
   delete hNbjets_pdf_mc_[i][j][k];
   delete hNbjets_pdf_est_[i][j][k];
   }
   printf("6.4[%d].1[%d]\n", i,j);
   delete [] hNbjetsEstSummary[i][j];
   delete [] hNbjetsMCSummary[i][j];
   delete [] hNbjets_mc_[i][j];
   delete [] hNbjets_pdf_mc_[i][j];
   delete [] hNbjets_pdf_est_[i][j];
   if (hsNbjets_MC[i][j] != NULL)
   delete hsNbjets_MC[i][j];
   if (hsNbjets_Est[i][j] != NULL)
   delete hsNbjets_Est[i][j];
   delete tCanva_Nbjets_Summary[i][j];
   
   printf("6.4[%d].2[%d]\n", i,j);
   //			delete hScanMin_Eb[i][j];
   //			delete hScanMin_Eudsc[i][j];
   //			delete hScanMin_Euds[i][j];
   }
   printf("6.5[%d]\n", i);
   delete [] hNbjetsEstSummary[i][NbOfJetsBins_];
   printf("6.5_1[%d]\n", i);
   delete [] hNbjetsMCSummary[i][NbOfJetsBins_];
   printf("6.5_2[%d]\n", i);
   printf(" -> %lu\n", (unsigned long) hsNbjets_MC[i][NbOfJetsBins_]);
   printf("6.5_22[%d]\n", i);
   if (hsNbjets_MC[i][NbOfJetsBins_] != NULL)
   delete hsNbjets_MC[i][NbOfJetsBins_];
   printf("6.5_3[%d]\n", i);
   if (hsNbjets_Est[i][NbOfJetsBins_] != NULL)
   delete hsNbjets_Est[i][NbOfJetsBins_];
   printf("6.5_4[%d]\n", i);
   delete tCanva_Nbjets_Summary[i][NbOfJetsBins_];
   
   printf("6.6[%d]\n", i);
   for(UInt_t j=0;j<NbOfBJetsBins_;j++){
   for(UInt_t k=0;k<3;k++){
   delete hNjetsEstSummary[i][j][k];
   delete hNjetsMCSummary[i][j][k];
   }
   delete [] hNjetsEstSummary[i][j];
   delete [] hNjetsMCSummary[i][j];
   if (hsNjets_MC[i][j] != NULL)
   delete hsNjets_MC[i][j];
   if (hsNjets_Est[i][j] != NULL)
   delete hsNjets_Est[i][j];
   delete tCanva_Njets_Summary[i][j];
   }
   delete [] hNjetsEstSummary[i][NbOfBJetsBins_];
   delete [] hNjetsMCSummary[i][NbOfBJetsBins_];
   if (hsNjets_MC[i][NbOfBJetsBins_] != NULL)
   delete hsNjets_MC[i][NbOfBJetsBins_];
   if (hsNjets_Est[i][NbOfBJetsBins_] != NULL)
   delete hsNjets_Est[i][NbOfBJetsBins_];
   delete tCanva_Njets_Summary[i][NbOfBJetsBins_];
   
   delete [] hNbjetsEstSummary[i];
   delete [] hNbjetsMCSummary[i];
   delete [] hNbjets_mc_[i];
   delete [] hNbjets_pdf_mc_[i];
   delete [] hNbjets_pdf_est_[i];
   delete [] hNjetsEstSummary[i];
   delete [] hNjetsMCSummary[i];
   delete [] hsNbjets_MC[i];
   delete [] hsNbjets_Est[i];
   delete [] hsNjets_MC[i];
   delete [] hsNjets_Est[i];
   delete [] tCanva_Nbjets_Summary[i];
   delete [] tCanva_Njets_Summary[i];
   }
   printf("7\n");
   
   //	for(int i=0;i<NbOfJetsBins_;i++){
   //		delete hScanMin_Ntt[i];
   //		delete hScanMin_Nv[i];
   //		delete hScanMin_Nvb[i];
   //	}
   printf("8\n");
   delete [] hNbjetsEstSummary;
   delete [] hNbjetsMCSummary;
   delete [] hNbjets_mc_;
   delete [] hNbjets_pdf_mc_;
   delete [] hNbjets_pdf_est_;
   delete [] hNjetsEstSummary;
   delete [] hNjetsMCSummary;
   delete [] hsNbjets_MC;
   delete [] hsNbjets_Est;
   delete [] hsNjets_MC;
   delete [] hsNjets_Est;
   delete [] tCanva_Nbjets_Summary;
   delete [] tCanva_Njets_Summary;
   //	delete [] hScanMin_Eb;
   //	delete [] hScanMin_Eudsc;
   //	delete [] hScanMin_Euds;
   printf("9\n");
   
   delete [] eudsc_;
   delete [] euds_;
   delete [] eb_;
   delete [] ebq_;
   delete [] Ntt_ ;
   delete [] Nv_;
   delete [] Nvb_;
   delete [] Njets_;
   delete [] MultiJet_Est_;
   printf("10\n");
   
   //	for(int i=0;i<NbOfBtagWorkingPoint_-1;i++){
   //		for(int j=0;j<NbOfJetsBins_;j++){
   //			for(int k=0;k<NbOfDatasets_;k++){
   //				for(int l=0;l<NbOfBJetsBins_;l++){
   //					for(int n=0;n<NbOfBJetsBins_;n++) cout<<condProb_[i][j][k][l][n]<<endl;
   //					delete [] condProb_[i][j][k][l];
   //					cout<<"Here I am : "<<i<<"/"<<j<<"/"<<k<<"/"<<l<<endl;
   //				}
   //				delete [] condProb_[i][j][k];
   //				cout<<"Here I loop"<<endl;
   //			}
   //			delete [] condProb_[i][j];
   //			cout<<"Here I loop again"<<endl;
   //		}
   //		delete [] condProb_[i];
   //		cout<<"Here I loop again and again"<<endl;
   //	}
   //	delete [] condProb_;
   //	cout<<"Here I should be"<<endl;
   
   printf("11\n");
   for(UInt_t i=0;i<NbOfDatasets_;i++){
   for(UInt_t j=0;j<NbOfJetsBins_;j++){
   delete hNbOfBGenJets_[i][j];
   }
   for(UInt_t j=0;j<NbOfBtagWorkingPoint_;j++){
   delete hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j];
   delete hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j];
   delete hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j];
   delete hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j];
   delete bTagEff_vs_Njets_[i][j];
   delete cTagEff_vs_Njets_[i][j];
   delete udsTagEff_vs_Njets_[i][j];
   delete misTagEff_vs_Njets_[i][j];
   }
   delete [] hNbOfBGenJets_[i];
   delete [] hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
   delete [] hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
   delete [] hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
   delete [] hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
   delete [] bTagEff_vs_Njets_[i];
   delete [] cTagEff_vs_Njets_[i];
   delete [] udsTagEff_vs_Njets_[i];
   delete [] misTagEff_vs_Njets_[i];
   delete hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
   delete hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
   delete hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
   delete hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
   }
   printf("12\n");
   delete [] hNbOfBGenJets_;
   delete [] hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
   delete [] hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
   delete [] hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
   delete [] hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
   delete [] bTagEff_vs_Njets_;
   delete [] cTagEff_vs_Njets_;
   delete [] udsTagEff_vs_Njets_;
   delete [] misTagEff_vs_Njets_;
   delete [] hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
   delete [] hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
   delete [] hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
   delete [] hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
   printf("13\n");
   
   //	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
   //		for(int j=0;j<NbOfJetsBins_+1;j++){
   //			delete [] NvVar_[i][j];
   //			delete [] NvBias_[i][j];
   //			delete [] NvbVar_[i][j];
   //			delete [] NvbBias_[i][j];
   //			delete [] NttVar_[i][j];
   //			delete [] NttBias_[i][j];
   //		}
   //		delete [] NvVar_[i];
   //		delete [] NvBias_[i];
   //		delete [] NvbVar_[i];
   //		delete [] NvbBias_[i];
   //		delete [] NttVar_[i];
   //		delete [] NttBias_[i];
   //	}
   printf("14\n");
   //	delete [] NvVar_;
   //	delete [] NvBias_;
   //	delete [] NvbVar_;
   //	delete [] NvbBias_;
   //	delete [] NttVar_;
   //	delete [] NttBias_;
   printf("15\n");
   
   if(RescaledTTLikeEstimation) delete RescaledTTLikeEstimation;
   if(RescaledVLikeEstimation)  delete RescaledVLikeEstimation;
   if(RescaledVbLikeEstimation) delete RescaledVbLikeEstimation;
   if(tCanva_RescaledTTLikeEstimation) delete tCanva_RescaledTTLikeEstimation;
   if(tCanva_RescaledVLikeEstimation)  delete tCanva_RescaledVLikeEstimation;
   if(tCanva_RescaledVbLikeEstimation) delete tCanva_RescaledVbLikeEstimation;
   printf("16\n");
   */
}

/**________________________________________________________________________________________________________________*/
/*
 void VJetEstimation::FillInputs(Double_t**** n, Double_t*** MultiJets_Estimated_Nj){
 for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
 for(UInt_t j=0;j<NbOfJetsBins_;j++){
 for(UInt_t k=0;k<NbOfBJetsBins_;k++){
 Nbjets_[i][j][k]  = 0;
 for(UInt_t l=0;l<NbOfDatasets_;l++){
 N_[i][j][k][l] = n[i][j][k][l];
 Nbjets_[i][j][k]  += N_[i][j][k][l];
 }
 Nbjets_[i][j][k]   += -MultiJets_Estimated_Nj[i][j][k];
 }
 }
 }
 }
 */
/**________________________________________________________________________________________________________________*/
/*
 void VJetEstimation::FillInputs(Double_t**** n){
 for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
 for(UInt_t j=0;j<NbOfJetsBins_;j++){
 for(UInt_t k=0;k<NbOfBJetsBins_;k++){
 Nbjets_[i][j][k]  = 0;
 for(UInt_t l=0;l<NbOfDatasets_;l++){
 N_[i][j][k][l] = n[i][j][k][l];
 Nbjets_[i][j][k]  += N_[i][j][k][l];
 }
 }
 }
 }
 }
 */
/**________________________________________________________________________________________________________________*/
void VJetEstimation::Fill(vector<TopTree::TRootJet*> &SelectedJets, UInt_t idx, Double_t (*btag_algo)(TopTree::TRootJet*), Double_t weight){
	Double_t btagDisc = 0;
	
  
    //	unsigned int nbofjetsIdx  = ((SelectedJets.size()-Njets_[0])<(unsigned int)(NbOfJetsBins_-1) ? (SelectedJets.size()-Njets_[0]) : (NbOfJetsBins_-1));
  if (SelectedJets.size()<Njets_[0])
    return;
  UInt_t nbofjetsIdx  = (UInt_t) (NbOfJetsBins_-1);
  if ((((UInt_t) SelectedJets.size())-(Njets_[0]))<((UInt_t)(NbOfJetsBins_-1)))
    nbofjetsIdx = (UInt_t) (((Int_t)SelectedJets.size())-Njets_[0]) ;
  
	
	UInt_t nbofbjetsIdx[NbOfBtagWorkingPoint_];
	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) nbofbjetsIdx[i]=0;
  
	std::vector<TopTree::TRootJet*> SelectedBGenJets;
  
	for(UInt_t i=0;i<SelectedJets.size();i++){
    btagDisc = btag_algo(SelectedJets[i]);
    
		for(UInt_t j=0;j<NbOfBtagWorkingPoint_;j++){
			if(btagDisc>BtagWorkingPoint_[j]) nbofbjetsIdx[j]++;
		}
    
		if(fabs(SelectedJets[i]->partonFlavour()) > 6) continue;
    
		if(fabs(SelectedJets[i]->partonFlavour()) == 5)	     hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx].Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx], weight);
		else if(fabs(SelectedJets[i]->partonFlavour()) == 4) hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx].Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx], weight);
		else                                                 hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx].Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx], weight);
    
		if(fabs(SelectedJets[i]->partonFlavour()) != 5)	     hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx].Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx], weight);
		else SelectedBGenJets.push_back(SelectedJets[i]);
    
		for(UInt_t j=0;j<NbOfBtagWorkingPoint_;j++){
			if(btagDisc>BtagWorkingPoint_[j]){
				if(fabs(SelectedJets[i]->partonFlavour()) == 5)      hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j].Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx], weight);
				else if(fabs(SelectedJets[i]->partonFlavour()) == 4) hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j].Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx], weight);
				else                                                 hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j].Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx], weight);
        
				if(fabs(SelectedJets[i]->partonFlavour()) != 5)      hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j].Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx], weight);
			}
    }
	}
  for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {
    N_[i][nbofjetsIdx][(nbofbjetsIdx[i]<NbOfBJetsBins_ ? nbofbjetsIdx[i] : NbOfBJetsBins_-1)][idx] += weight;
      //   Nbjets_[i][nbofjetsIdx][(nbofbjetsIdx[i]<NbOfBJetsBins_ ? nbofbjetsIdx[i] : NbOfBJetsBins_-1)] += weight ;
  }
    //	for(UInt_t i=0;i<NbOfBtagWorkingPoint_-1;i++) condProb_[i][nbofjetsIdx][idx][(nbofbjetsIdx[i+1]<NbOfBJetsBins_ ? nbofbjetsIdx[i+1] : NbOfBJetsBins_-1)][(nbofbjetsIdx[i]<NbOfBJetsBins_ ? nbofbjetsIdx[i] : NbOfBJetsBins_-1)] += weight;
  
	hNbOfBGenJets_[idx][nbofjetsIdx].Fill(SelectedBGenJets.size(), weight);
  
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::FillInputs(vector <vector< vector< vector< Double_t > > > > Inputs){
    // 1. Make sure the number of vectors does match the number of datasets specified at the instanciation of the VJetEstimation obj.
	if(Inputs.size() != NbOfDatasets_){
		cerr << "FillInputs Error 1 : the vector size does not match the number of datasets." << endl;
		exit(1);
	}
    // 2. Loop over the different datasets
	for(unsigned idx=0; idx<NbOfDatasets_;idx++){
      // 3. Make sure that each vector contains a number of histos equal to the number of b-tagging working points (NbOfBtagWorkingPoint_).
		if(Inputs[idx].size() != (UInt_t)NbOfBtagWorkingPoint_){
			cerr << "FillInputs Error 2 : the vector size does not match the number of b-tagging working points." << endl;
			exit(1);
		}
      // 4. Loop over the different b-tagging working points.
		for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
			if(Inputs[idx][i].size()<NbOfJetsBins_){
				cerr << "FillInputs Error 3 : the vector size does not match the maximal number of jets." << endl;
				exit(1);
			}
        // 5. Loop over the number of jets.
			for(UInt_t j=0; j<NbOfJetsBins_; j++){
				if(Inputs[idx][i][j].size()<NbOfBJetsBins_){
					cerr << "FillInputs Error 4 : the vector size does not match the maximal number of b-jets." << endl;
					exit(1);
				}
          // 6. Loop over the number of b-jets.
				for(UInt_t k=0; k<NbOfBJetsBins_; k++){
					N_[i][j][k][idx] = Inputs[idx][i][j][k];
				}
			}
		}
	}
	SumOverAllInputs();
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::FillInputs(vector <vector< vector< vector< Double_t > > > > Inputs, vector <vector< vector<Double_t> > > Multijet_est_nj){
    // 1. Make sure the number of vectors does match the number of datasets specified at the instanciation of the VJetEstimation obj.
	if(Inputs.size() != NbOfDatasets_){
		cerr << "FillInputs Error 1 : the vector size does not match the number of datasets." << endl;
		exit(1);
	}
    // 2. Loop over the different datasets
	for(unsigned idx=0; idx<NbOfDatasets_;idx++){
      // 3. Make sure that each vector contains a number of histos equal to the number of b-tagging working points (NbOfBtagWorkingPoint_).
		if(Inputs[idx].size() != NbOfBtagWorkingPoint_){
			cerr << "FillInputs Error 2 : the vector size does not match the number of b-tagging working points." << endl;
			exit(1);
		}
      // 4. Loop over the different b-tagging working points.
		for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
			if(Inputs[idx][i].size()<NbOfJetsBins_){
				cerr << "FillInputs Error 3 : the vector size does not match the maximal number of jets." << endl;
				exit(1);
			}
        // 5. Loop over the number of jets.
			for(UInt_t j=0; j<NbOfJetsBins_; j++){
				if(Inputs[idx][i][j].size()<NbOfBJetsBins_){
					cerr << "FillInputs Error 4 : the vector size does not match the maximal number of b-jets." << endl;
					exit(1);
				}
          // 6. Loop over the number of b-jets.
				for(UInt_t k=0; k<NbOfBJetsBins_; k++){
					N_[i][j][k][idx] = Inputs[idx][i][j][k];
				}
			}
		}
	}
    // 7. Sum the input weighted contributions
	SumOverAllInputs();
    // 8. Subtract the multi-jet background
	if(Multijet_est_nj.size() != NbOfBtagWorkingPoint_){
		cerr << "FillInputs Error 5 : the vector size does not match the number of b-tagging working points." << endl;
		exit(1);
	}
	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
		if(Multijet_est_nj[i].size() != NbOfJetsBins_){
			cerr << "FillInputs Error 6 : the vector size does not match the number of jets." << endl;
			exit(1);
		}
		for(UInt_t j=0;j<NbOfJetsBins_;j++){
			if(Multijet_est_nj[i][j].size() != NbOfBJetsBins_){
				cerr << "FillInputs Error 7 : the vector size does not match the number of b-jets." << endl;
				exit(1);
			}
			for(UInt_t k=0;k<NbOfBJetsBins_;k++){
				Nbjets_[i][j][k]   += -Multijet_est_nj[i][j][k];
			}
		}
	}
}

/**________________________________________________________________________________________________________________*/
Float_t VJetEstimation::WilsonScoreIntervalHigh(Float_t Non, Float_t Ntot)
{
	Double_t T = (Ntot>0 ? 1/Ntot : 0);
	Double_t p_hat = (Ntot>0 && Non>=0 && Ntot>=Non ? Non/Ntot : 0);
	Double_t Int_High = ((p_hat+(T/2))/(1+T))+(sqrt(p_hat*(1-p_hat)*T+pow(T/2,2))/(1+T));
	return Int_High;
}

/**________________________________________________________________________________________________________________*/
Float_t VJetEstimation::WilsonScoreIntervalLow(Float_t Non, Float_t Ntot)
{
	Double_t T = (Ntot>0 ? 1/Ntot : 0);
	Double_t p_hat = (Ntot>0 && Non>=0 && Ntot>=Non ? Non/Ntot : 0);
	Double_t Int_Low = ((p_hat+(T/2))/(1+T))-(sqrt(p_hat*(1-p_hat)*T+pow(T/2,2))/(1+T));
	return Int_Low;
}

/**________________________________________________________________________________________________________________*/
Float_t VJetEstimation::WilsonScoreIntervalMean(Float_t Non, Float_t Ntot)
{
	Double_t Err_High = WilsonScoreIntervalHigh(Non, Ntot)-Non;
	Double_t Err_Low  = Non-WilsonScoreIntervalLow(Non, Ntot);
	return (Err_High+Err_Low)/2;
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::ReScaleInputs(UInt_t idx, Float_t Ntot, Double_t factor){
	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
		for(UInt_t j=0;j<NbOfJetsBins_;j++){
			for(UInt_t k=0;k<NbOfBJetsBins_;k++){
				N_err_[i][j][k][idx]     = Ntot*WilsonScoreIntervalMean(N_[i][j][k][idx],Ntot);
				N_err_[i][j][k][idx]    *= factor;
        
				N_err_hi_[i][j][k][idx]  = Ntot*WilsonScoreIntervalHigh(N_[i][j][k][idx],Ntot);
				N_err_lo_[i][j][k][idx]  = Ntot*WilsonScoreIntervalLow(N_[i][j][k][idx],Ntot);
				N_err_hi_[i][j][k][idx] *= factor;
				N_err_lo_[i][j][k][idx] *= factor;
        
				N_[i][j][k][idx]        *= factor;
          //cout<<"N_["<<i<<"]["<<j<<"]["<<k<<"]["<<idx<<"] = "<<N_[i][j][k][idx]<<endl;
          //cout<<"N_err_["<<i<<"]["<<j<<"]["<<k<<"]["<<idx<<"] = "<<N_err_[i][j][k][idx]<<endl;
          //cout<<"Ntot = "<<Ntot<<endl;
			}
		}
	}	
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetProcesses(std::vector<Bool_t> processMask, std::vector<std::string> ttLikeDatasetNames, std::vector<std::string> vLikeDatasetNames, std::vector<std::string> vbLikeDatasetNames)
{
  if(processMask.size()==vDatasets_.size()) {
    processMask_ = processMask;
  } else {
    fprintf(stderr, "VJetEstimation::SetProcesses : processMask has wrong size (%lu) compared to datasets(%lu) !\n", processMask.size(), vDatasets_.size());
    exit(1);
  }
  iDatasetsTTLike_.clear();
  iDatasetsVLike_.clear();
  iDatasetsVbLike_.clear();
  Int_t indexDataset = 0;
  for (std::vector<Dataset>::const_iterator it=vDatasets_.begin(); it!=vDatasets_.end(); it++) {
    for (std::vector<std::string>::const_iterator like=ttLikeDatasetNames.begin(); like!=ttLikeDatasetNames.end(); like++) {
      if (it->Name() == *like) {
        iDatasetsTTLike_.push_back(indexDataset);
      }
    }
    for (std::vector<std::string>::const_iterator like=vLikeDatasetNames.begin(); like!=vLikeDatasetNames.end(); like++) {
      if (it->Name() == *like) {
        iDatasetsVLike_.push_back(indexDataset);
      }
    }
    for (std::vector<std::string>::const_iterator like=vbLikeDatasetNames.begin(); like!=vbLikeDatasetNames.end(); like++) {
      if (it->Name() == *like) {
        iDatasetsVbLike_.push_back(indexDataset);
      }
    }
    indexDataset++;
  }
}

    
/**________________________________________________________________________________________________________________*/
void VJetEstimation::SumOverAllInputs(){
	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
		for(UInt_t j=0;j<NbOfJetsBins_;j++){
			for(UInt_t k=0;k<NbOfBJetsBins_;k++){
				Nbjets_[i][j][k] = 0; // Reset value
				for(UInt_t l=0;l<NbOfDatasets_;l++)
          if (processMask_[l]) {
            Nbjets_[i][j][k] += N_[i][j][k][l]; //Loop over all the datasets and sum their weighted contributions.
          }
			}
		}
	}	
}
/**________________________________________________________________________________________________________________*/
void VJetEstimation::ComputeEffbqFromMC(){
  Double_t num0b = 0.;
  Double_t num1b = 0.;
  Double_t num2b = 0.;
  Double_t denom = 0.;
  for(UInt_t i=0;i<NbOfJetsBins_;i++){
    for(UInt_t idx=0;idx<iDatasetsTTLike_.size();idx++){
      num0b += hNbOfBGenJets_[iDatasetsTTLike_[idx]][i].GetBinContent(1);
      num1b += hNbOfBGenJets_[iDatasetsTTLike_[idx]][i].GetBinContent(2);
      num2b += hNbOfBGenJets_[iDatasetsTTLike_[idx]][i].GetBinContent(3);
      denom += hNbOfBGenJets_[iDatasetsTTLike_[idx]][i].Integral();
    }
    e0bq_[i] = num0b/denom;
    e1bq_[i] = num1b/denom;
    e2bq_[i] = num2b/denom;
  }
}
/**________________________________________________________________________________________________________________*/
void VJetEstimation::ComputeEffFromMC(){
	double x = 0;
	char name[100];
	TH1D* histo_total = 0;
	TH1D* histo_passed = 0;
  
  {
    UInt_t idx=0;
    for(std::vector<Bool_t>::const_iterator iter=processMask_.begin() ; iter!=processMask_.end() ;iter++) {
      if ( (*iter) == kFALSE ) {
        idx++;
        continue;
      }
      for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
        
        histo_total  = (TH1D*)hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx].Project3D("z");
        histo_passed = (TH1D*)hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i].Project3D("z");
        if(histo_total->GetEntries()>0){
          sprintf(name,"bTagEff_vs_Njets_%d_wp_%d",idx,i);
          bTagEff_vs_Njets_[idx][i] = new TEfficiency(*histo_passed,*histo_total);
          bTagEff_vs_Njets_[idx][i]->SetName(name);
          bTagEff_vs_Njets_[idx][i]->SetTitle(";Nb. of jets;eff.");
        }
        histo_total  = (TH1D*)hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx].Project3D("z");
        histo_passed = (TH1D*)hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i].Project3D("z");
        if(histo_total->GetEntries()>0){
          sprintf(name,"cTagEff_vs_Njets_%d_wp_%d",idx,i);
          cTagEff_vs_Njets_[idx][i] = new TEfficiency(*histo_passed,*histo_total);
          cTagEff_vs_Njets_[idx][i]->SetName(name);
          cTagEff_vs_Njets_[idx][i]->SetTitle(";Nb. of jets;eff.");
        }
        histo_total  = (TH1D*)hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx].Project3D("z");
        histo_passed = (TH1D*)hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i].Project3D("z");
        if(histo_total->GetEntries()>0){
          sprintf(name,"udsTagEff_vs_Njets_%d_wp_%d",idx,i);
          udsTagEff_vs_Njets_[idx][i] = new TEfficiency(*histo_passed,*histo_total);
          udsTagEff_vs_Njets_[idx][i]->SetName(name);
          udsTagEff_vs_Njets_[idx][i]->SetTitle(";Nb. of jets;eff.");
        }
        histo_total  = (TH1D*)hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx].Project3D("z");
        histo_passed = (TH1D*)hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i].Project3D("z");
        if(histo_total->GetEntries()>0){
          sprintf(name,"misTagEff_vs_Njets_%d_wp_%d",idx,i);
          misTagEff_vs_Njets_[idx][i] = new TEfficiency(*histo_passed,*histo_total);
          misTagEff_vs_Njets_[idx][i]->SetName(name);
          misTagEff_vs_Njets_[idx][i]->SetTitle(";Nb. of jets;eff.");
        }
        
      }
      idx++ ;
    }
  }
  
    //TGraphAsymmErrors* Combine(TCollection* pList, Option_t* opt = "", Int_t n = 0, const Double_t* w = 0)
  UInt_t nEmpty_vlike = 0;
  UInt_t nEmpty_ttlike = 0;
  Double_t w_ttlike[iDatasetsTTLike_.size()];
  Double_t w_vlike[iDatasetsVLike_.size()];
  {
    UInt_t index=0;
    for(UInt_t i=0;i<iDatasetsTTLike_.size();i++) {
      w_ttlike[index] = GetPredN(iDatasetsTTLike_[i]);
      if (w_ttlike[index]==0.) {
        nEmpty_ttlike++;
      } else {
        index++;
      }
    }
    index=0;
    for(UInt_t i=0;i<iDatasetsVLike_.size();i++) {
      w_vlike[index] = GetPredN(iDatasetsVLike_[i]);
      if (w_vlike[index]==0.) {
        nEmpty_vlike++; //if empty, the weight==0 or -1 causes a crash ; invalid custom weight found w = -1.00
                        //The "index" game is to same only the non empty weights in the range [0;iDatasetsVLike_.size()-nEmpty_vlike[ of indices
      } else {
        index++;
      }
    }
  }
	TList *pColl_bTagEff = new TList();
	TList *pColl_misTagEff_TTlike = new TList();
	TList *pColl_misTagEff_Vlike  = new TList();
  
	for(UInt_t j=0;j<NbOfBtagWorkingPoint_;j++){
		for(UInt_t i=0;i<iDatasetsTTLike_.size();i++){
      if(!bTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]) continue;
      pColl_bTagEff->Add((TEfficiency*)bTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]);
      
      if(!misTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]) continue;
      pColl_misTagEff_TTlike->Add((TEfficiency*)misTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]);
		}
		for(UInt_t i=0;i<iDatasetsVLike_.size();i++){
      if(!misTagEff_vs_Njets_[iDatasetsVLike_[i]][j]) continue;
      pColl_misTagEff_Vlike->Add((TEfficiency*)misTagEff_vs_Njets_[iDatasetsVLike_[i]][j]);
		}
    bTagEff_vs_Njets_TTlike_[j]   = bTagEff_vs_Njets_[iDatasetsTTLike_[0]][j]->Combine(pColl_bTagEff,"",iDatasetsTTLike_.size()-nEmpty_ttlike, w_ttlike);
    misTagEff_vs_Njets_TTlike_[j] = misTagEff_vs_Njets_[iDatasetsTTLike_[0]][j]->Combine(pColl_misTagEff_TTlike,"",iDatasetsTTLike_.size()-nEmpty_ttlike, w_ttlike);
    misTagEff_vs_Njets_Vlike_[j]  = misTagEff_vs_Njets_[iDatasetsVLike_[0]][j]->Combine(pColl_misTagEff_Vlike,"",iDatasetsVLike_.size()-nEmpty_vlike, w_vlike);
    
		for(UInt_t k = 0 ; k < NbOfJetsBins_ ; k++)
		{
      bTagEff_vs_Njets_TTlike_[j]->GetPoint(k,x,eb_mc_[j][k]);
      eb_err_mc_[j][k] = bTagEff_vs_Njets_TTlike_[j]->GetErrorY(k);
      
      misTagEff_vs_Njets_TTlike_[j]->GetPoint(k,x,eudsc_mc_[j][k]);
      eudsc_err_mc_[j][k] = misTagEff_vs_Njets_TTlike_[j]->GetErrorY(k);
      
      misTagEff_vs_Njets_Vlike_[j]->GetPoint(k,x,euds_mc_[j][k]);
      euds_err_mc_[j][k] = misTagEff_vs_Njets_Vlike_[j]->GetErrorY(k);
		}
    pColl_bTagEff->Clear();
    pColl_misTagEff_TTlike->Clear();
    pColl_misTagEff_Vlike->Clear();
	}
  pColl_bTagEff->Delete();
  pColl_misTagEff_TTlike->Delete();
  pColl_misTagEff_Vlike->Delete();
}

/**________________________________________________________________________________________________________________*/
/*
 void VJetEstimation::Write(TFile* file, string label, bool verbose){
 cout<<"--Start writing histograms for the VJetEstimation --"<<endl;
 file->cd();
 string dirname = "VJetEstimation"+label;
 file->mkdir(dirname.c_str());
 file->cd(dirname.c_str());
 
 char name[100];char title[100];
 
 for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
 for(UInt_t j=0;j<NbOfJetsBins_;j++){
 for(UInt_t k=0;k<3;k++){
 hNbjets_mc_[i][j][k]->Write();
 hNbjets_pdf_mc_[i][j][k]->Write();
 hNbjets_pdf_est_[i][j][k]->Write();
 }
 }
 if(  bTagEff_vs_Njets_TTlike_[i]!=0)   bTagEff_vs_Njets_TTlike_[i]->Write();
 if(misTagEff_vs_Njets_TTlike_[i]!=0) misTagEff_vs_Njets_TTlike_[i]->Write();
 if(misTagEff_vs_Njets_Vlike_[i] !=0) misTagEff_vs_Njets_Vlike_[i]->Write();
 
 for(UInt_t j=0;j<=NbOfJetsBins_;j++){
 // For Monte Carlo prediction 
 MyLeg->Clear();
 tCanva_Nbjets_Summary[i][j]->cd();
 
 if(j!=NbOfJetsBins_) sprintf(name,"hNbjetsMCStackSummary_%d_jets_wp_%d",Njets_[j],i);
 else                 sprintf(name,"hNbjetsMCStackSummary_Inclusive_wp_%d",i);
 if(j!=NbOfJetsBins_) sprintf(title,"V-like and tt-like MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
 else                 sprintf(title,"V-like and tt-like MC prediction summary (inclusive, WP nr%d)",i);
 
 hsNbjets_MC[i][j] = new THStack(name,title);
 
 hNbjetsMCSummary[i][j][0]->SetFillStyle(3004);
 hNbjetsMCSummary[i][j][0]->SetFillColor(kGreen+1);
 hNbjetsMCSummary[i][j][1]->SetFillStyle(3005);
 hNbjetsMCSummary[i][j][1]->SetFillColor(kBlue+1);
 hNbjetsMCSummary[i][j][2]->SetFillStyle(3006);
 hNbjetsMCSummary[i][j][2]->SetFillColor(kRed+1);
 
 hsNbjets_MC[i][j]->Add(hNbjetsMCSummary[i][j][0]);
 //hsNbjets_MC[i][j]->Add(hNbjetsMCSummary[i][j][1]);
 hsNbjets_MC[i][j]->Add(hNbjetsMCSummary[i][j][2]);
 
 MyLeg->AddEntry(hNbjetsMCSummary[i][j][0],"V-like events","f");
 //MyLeg->AddEntry(hNbjetsMCSummary[i][j][1],"Vb-like events","f");
 MyLeg->AddEntry(hNbjetsMCSummary[i][j][2],"#bar{t}t-like events","f");
 
 hsNbjets_MC[i][j]->Draw();
 hsNbjets_MC[i][j]->GetXaxis()->SetTitle("Nb of b-jets");
 
 if(i!=NbOfJetsBins_) sprintf(name,"hNbjetsEstStackSummary_%d_jets_wp_%d",Njets_[j],i);
 else                 sprintf(name,"hNbjetsEstStackSummary_Inclusive_wp_%d",i);
 if(i!=NbOfJetsBins_) sprintf(title,"V-like and tt-like estimation summary for %d jets (WP nr%d)",Njets_[j],i);
 else                 sprintf(title,"V-like and tt-like estimation summary (inclusive, WP nr%d)",i);
 
 hsNbjets_Est[i][j] = new THStack(name,title);
 
 hNbjetsEstSummary[i][j][0]->SetMarkerColor(kYellow+1);
 //hNbjetsEstSummary[i][j][1]->SetMarkerColor(kMagenta+1);
 hNbjetsEstSummary[i][j][2]->SetMarkerColor(kBlue+1);
 hNbjetsEstSummary[i][j][0]->SetMarkerStyle(22);
 //hNbjetsEstSummary[i][j][1]->SetMarkerStyle(23);
 hNbjetsEstSummary[i][j][2]->SetMarkerStyle(24);
 hNbjetsEstSummary[i][j][0]->SetMarkerSize(1.5);
 //hNbjetsEstSummary[i][j][1]->SetMarkerSize(1.5);
 hNbjetsEstSummary[i][j][2]->SetMarkerSize(1.5);
 
 hsNbjets_Est[i][j]->Add(hNbjetsEstSummary[i][j][0]);
 //hsNbjets_Est[i][j]->Add(hNbjetsEstSummary[i][j][1]);
 hsNbjets_Est[i][j]->Add(hNbjetsEstSummary[i][j][2]);
 
 MyLeg->AddEntry(hNbjetsEstSummary[i][j][0],"V-like estimation","p");
 //MyLeg->AddEntry(hNbjetsEstSummary[i][j][1],"Vb-like estimation","p");
 MyLeg->AddEntry(hNbjetsEstSummary[i][j][2],"TT-like estimation","p");
 
 hsNbjets_Est[i][j]->Draw("PEsame");
 hsNbjets_Est[i][j]->GetXaxis()->SetTitle("Nb of b-jets");
 MyLeg->Draw("same");
 
 tCanva_Nbjets_Summary[i][j]->Update();
 if(verbose)cout<<"Writing summary histograms for "<<j+4<<" jets, wp nr "<<i<<endl;
 tCanva_Nbjets_Summary[i][j]->Write();
 }
 for(UInt_t j=0;j<=NbOfBJetsBins_;j++){
 
 // For Monte Carlo prediction 
 MyLeg->Clear();
 tCanva_Njets_Summary[i][j]->cd();
 
 if(j!=4) sprintf(name,"hNjetsMCStackSummary_%d_bjets_wp_%d",j,i);
 else     sprintf(name,"hNjetsMCStackSummary_bInclusive_wp_%d",i);
 if(j!=4) sprintf(title,"V-like and tt-like MC prediction summary for %d b-jets (WP nr%d)",j,i);
 else     sprintf(title,"V-like and tt-like MC prediction summary (b-inclusive,WP nr%d)",i);
 hsNjets_MC[i][j] = new THStack(name,title);
 
 hNjetsMCSummary[i][j][0]->SetFillStyle(3004);
 hNjetsMCSummary[i][j][0]->SetFillColor(kGreen+1);
 //hNjetsMCSummary[i][j][1]->SetFillStyle(3005);
 //hNjetsMCSummary[i][j][1]->SetFillColor(kBlue+1);
 hNjetsMCSummary[i][j][2]->SetFillStyle(3006);
 hNjetsMCSummary[i][j][2]->SetFillColor(kRed+1);
 
 hsNjets_MC[i][j]->Add(hNjetsMCSummary[i][j][0]);
 //hsNjets_MC[i][j]->Add(hNjetsMCSummary[i][j][1]);
 hsNjets_MC[i][j]->Add(hNjetsMCSummary[i][j][2]);
 
 MyLeg->AddEntry(hNjetsMCSummary[i][j][0],"V-like events","f");
 //MyLeg->AddEntry(hNjetsMCSummary[i][j][1],"Vb-like events","f");
 MyLeg->AddEntry(hNjetsMCSummary[i][j][2],"#bar{t}t-like events","f");
 
 hsNjets_MC[i][j]->Draw();
 hsNjets_MC[i][j]->GetXaxis()->SetTitle("Nb of jets");
 
 //hNjetsEstSummary[i][0]->Write();
 //hNjetsEstSummary[i][1]->Write();
 
 if(j!=4) sprintf(name,"hNjetsEstStackSummary_%d_bjets_wp_%d",j,i);
 else     sprintf(name,"hNjetsEstStackSummary_bInclusive_wp_%d",i);
 if(j!=4) sprintf(title,"V-like and tt-like estimation summary for %d b-jets (WP nr%d)",j,i);
 else     sprintf(title,"V-like and tt-like estimation summary (b-inclusive,WP nr%d)",i);
 hsNjets_Est[i][j] = new THStack(name,title);
 
 //sprintf(name,"hNjetsEstStackLegend_%d_bjets",i);
 //MyLeg->SetName(name);
 
 hNjetsEstSummary[i][j][0]->SetMarkerColor(kYellow+1);
 //hNjetsEstSummary[i][j][1]->SetMarkerColor(kMagenta+1);
 hNjetsEstSummary[i][j][2]->SetMarkerColor(kBlue+1);
 hNjetsEstSummary[i][j][0]->SetMarkerStyle(22);
 //hNjetsEstSummary[i][j][1]->SetMarkerStyle(23);
 hNjetsEstSummary[i][j][2]->SetMarkerStyle(24);
 hNjetsEstSummary[i][j][0]->SetMarkerSize(1.5);
 //hNjetsEstSummary[i][j][1]->SetMarkerSize(1.5);
 hNjetsEstSummary[i][j][2]->SetMarkerSize(1.5);
 
 hsNjets_Est[i][j]->Add(hNjetsEstSummary[i][j][0]);
 //hsNjets_Est[i][j]->Add(hNjetsEstSummary[i][j][1]);
 hsNjets_Est[i][j]->Add(hNjetsEstSummary[i][j][2]);
 
 MyLeg->AddEntry(hNjetsEstSummary[i][j][0],"V-like estimation","p");
 //MyLeg->AddEntry(hNjetsEstSummary[i][j][1],"Vb-like estimation","p");
 MyLeg->AddEntry(hNjetsEstSummary[i][j][2],"TT-like estimation","p");
 
 hsNjets_Est[i][j]->Draw("PE SAME");
 hsNjets_Est[i][j]->GetXaxis()->SetTitle("Nb of jets");
 MyLeg->Draw("same");
 
 tCanva_Njets_Summary[i][j]->Update();
 tCanva_Njets_Summary[i][j]->Write();
 if(verbose)cout<<"Writing summary histograms for "<<j<<" b-jets, wp nr "<<i<<endl;
 }
 }
 for(UInt_t idx=0;idx<NbOfDatasets_;idx++){
 for(UInt_t i=0;i<NbOfJetsBins_;i++){
 if(hNbOfBGenJets_[idx][i]==0) continue;
 if(hNbOfBGenJets_[idx][i]->Integral() > 0) hNbOfBGenJets_[idx][i]->Scale(1/hNbOfBGenJets_[idx][i]->Integral());
 hNbOfBGenJets_[idx][i]->Write();
 }
 
 hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Write();
 hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Write();
 hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Write();
 hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Write();
 
 for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
 
 hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i]->Write();
 hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i]->Write();
 hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i]->Write();
 hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i]->Write();
 
 if(  bTagEff_vs_Njets_[idx][i]!=0)   bTagEff_vs_Njets_[idx][i]->Write();
 if(  cTagEff_vs_Njets_[idx][i]!=0)   cTagEff_vs_Njets_[idx][i]->Write();
 if(udsTagEff_vs_Njets_[idx][i]!=0) udsTagEff_vs_Njets_[idx][i]->Write();
 if(misTagEff_vs_Njets_[idx][i]!=0) misTagEff_vs_Njets_[idx][i]->Write();
 }
 }
 if(RescaledTTLikeEstimation) RescaledTTLikeEstimation->Write();
 if(RescaledVLikeEstimation)  RescaledVLikeEstimation->Write();
 if(RescaledVbLikeEstimation) RescaledVbLikeEstimation->Write();
 if(tCanva_RescaledTTLikeEstimation) tCanva_RescaledTTLikeEstimation->Write();
 if(tCanva_RescaledVLikeEstimation)  tCanva_RescaledVLikeEstimation->Write();
 if(tCanva_RescaledVbLikeEstimation) tCanva_RescaledVbLikeEstimation->Write();
 cout<<"--Done with writing histograms for the VJetEstimation --"<<endl;
 }
 */
/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Nttlike(vector<Double_t> init){
	if(init_Nttlike_.size()==init.size())
	{
		init_Nttlike_ = init;
	}
	else
	{
		cerr<<"Init Nttlike : Vector sizes do not match"<<endl;
		exit(1);
	}
};
/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Nttlike(UInt_t njets,Double_t init){
	if(init_Nttlike_.size()>njets)
	{
		init_Nttlike_[njets] = init;
	}
	else
	{
		cerr<<"Init Nttlike : njets exceeds vector size"<<endl;
		exit(1);
	}
};

/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Nvlike(vector<Double_t> init){
	if(init_Nvlike_.size()==init.size())
	{
		init_Nvlike_ = init;
	}
	else
	{
		cerr<<"Init Nvlike : Vector sizes do not match"<<endl;
		exit(1);
	}
};
/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Nvlike(UInt_t njets,Double_t init){
	if(init_Nvlike_.size()>njets)
	{
		init_Nvlike_[njets] = init;
	}
	else
	{
		cerr<<"Init Nvlike : njets exceeds vector size"<<endl;
		exit(1);
	}
};

/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Eb(vector< vector<Double_t> > init){
	if(init_Eb_.size()==init.size())
	{
		init_Eb_ = init;
	}
	else
	{
		cerr<<"Init Eb : Vector sizes do not match"<<endl;
		exit(1);
	}
};
/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Eb(UInt_t njets,vector<Double_t> init){
	if(init_Eb_.size()>njets)
	{
		init_Eb_[njets] = init;
	}
	else
	{
		cerr<<"Init Eb : njets exceeds vector size"<<endl;
		exit(1);
	}
};
/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Eb(UInt_t njets,UInt_t btagIdx, Double_t init){
	if(init_Eb_.size()>njets && init_Eb_[njets].size()>btagIdx)
	{
		init_Eb_[njets][btagIdx] = init;
	}
	else
	{
		cerr<<"Init Eb : njets/btagIdx exceeds vector size"<<endl;
		exit(1);
	}
};

/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Eudsc(vector< vector<Double_t> > init){
	if(init_Eudsc_.size()==init.size())
	{
		init_Eudsc_ = init;
	}
	else
	{
		cerr<<"Init Eudsc : Vector sizes do not match"<<endl;
		exit(1);
	}
};
/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Eudsc(UInt_t njets,vector<Double_t> init){
	if(init_Eudsc_.size()>njets)
	{
		init_Eudsc_[njets] = init;
	}
	else
	{
		cerr<<"Init Eudsc : njets exceeds vector size"<<endl;
		exit(1);
	}
};
/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Eudsc(UInt_t njets,UInt_t btagIdx, Double_t init){
	if(init_Eudsc_.size()>njets && init_Eudsc_[njets].size()>btagIdx)
	{
		init_Eudsc_[njets][btagIdx] = init;
	}
	else
	{
		cerr<<"Init Eudsc : njets/btagIdx exceeds vector size"<<endl;
		exit(1);
	}
};
/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Euds(vector< vector<Double_t> > init){
	if(init_Euds_.size()==init.size())
	{
		init_Euds_ = init;
	}
	else
	{
		cerr<<"Init Euds : Vector sizes do not match"<<endl;
		exit(1);
	}
};
/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Euds(UInt_t njets,vector<Double_t> init){
	if(init_Euds_.size()>njets)
	{
		init_Euds_[njets] = init;
	}
	else
	{
		cerr<<"Init Euds : njets exceeds vector size"<<endl;
		exit(1);
	}
};
/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetInitialValues_Euds(UInt_t njets,UInt_t btagIdx, Double_t init){
	if(init_Euds_.size()>njets && init_Euds_[njets].size()>btagIdx)
	{
		init_Euds_[njets][btagIdx] = init;
	}
	else
	{
		cerr<<"Init Euds : njets/btagIdx exceeds vector size"<<endl;
		exit(1);
	}
};

  /////////////////////////////////
  //Main Methods to be executed ...
  /////////////////////////////////
/**________________________________________________________________________________________________________________*/
void VJetEstimation::UnBinnedMaximumLikelihoodEst(UInt_t btag_wp_idx, UInt_t njets, vector<Int_t> &FixedVarIdx, bool doMinos, bool doWS, bool Verbose){
  
    // Declare observables
  char name[100];
  
	RooRealVar Ntt("Ntt","N_{t#bar{t}-like}",init_Nttlike_[njets],init_Nttlike_[njets]/5,init_Nttlike_[njets]*5);
	RooRealVar Nv( "Nv", "N_{V-like}",       init_Nvlike_[njets],init_Nvlike_[njets]/5,init_Nvlike_[njets]*5);
  
	vector<Int_t>::iterator Idx;
  
	RooRealVar  eb(    "eb",  "#epsilon_{b-tag}",      init_Eb_[njets][btag_wp_idx],0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 0 );
	if(Idx != FixedVarIdx.end()) eb.setConstant(kTRUE) ;
	RooRealVar  eudsc("eudsc","#epsilon_{mis-tag}",    init_Eudsc_[njets][btag_wp_idx],0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 1 );
	if(Idx != FixedVarIdx.end()) eudsc.setConstant(kTRUE) ;
	RooRealVar euds( "euds",  "#epsilon^{,}_{mis-tag}",init_Euds_[njets][btag_wp_idx],0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 2 );
	if(Idx != FixedVarIdx.end()) euds.setConstant(kTRUE) ;
  
    // Declare constants
  
	RooConstVar n("n","number of selected jets",Njets_[njets]) ;
	RooConstVar e0bq("e0bq","e0bq",e0bq_[njets]);
	RooConstVar e1bq("e1bq","e1bq",e1bq_[njets]);
	RooConstVar e2bq("e2bq","e2bq",e2bq_[njets]);
    //RooConstVar e3bq("e3bq","e3bq",e3bq_[njets]);
  
    // C o n s t r u c t   a   c a t e g o r y   w i t h   l a b e l s    a n d   i n d e c e s
    // -----------------------------------------------------------------------------------------
  
	RooCategory nbjets("nbjets","Number of b-jets");
	nbjets.defineType("N0bjet", 0);
	nbjets.defineType("N1bjet", 1);
	nbjets.defineType("N2bjets",2);
	nbjets.defineType("N3bjets",3);
  
    // C o n s t r u c t   f o r m u l a s
    // -------------------------------------------------------------------------------
  
  RooFormulaVar p0bjets_tt("p0bjets_tt","p0bjets_tt","(1-eb)*(1-eb)*pow((1-eudsc),n-2)*e2bq+(1-eb)*pow((1-eudsc),n-1)*e1bq+pow((1-eudsc),n)*e0bq",RooArgList(eb,eudsc,n,e0bq,e1bq,e2bq));
  RooFormulaVar p1bjets_tt("p1bjets_tt","p1bjets_tt","(2*eb*(1-eb)*pow(1-eudsc,n-2)+(1-eb)*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3))*e2bq+(eb*pow(1-eudsc,n-1)+(1-eb)*(n-1)*eudsc*pow(1-eudsc,n-2))*e1bq+(n*eudsc*pow(1-eudsc,n-1))*e0bq",RooArgList(eb,eudsc,n,e0bq,e1bq,e2bq));
  RooFormulaVar p2bjets_tt("p2bjets_tt","p2bjets_tt","(eb*eb*pow(1-eudsc,n-2)+2*eb*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3)+(1-eb)*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4))*e2bq+(eb*(n-1)*eudsc*pow(1-eudsc,n-2)+(1-eb)*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*e1bq+((n*(n-1)/2)*eudsc*eudsc*pow(1-eudsc,n-2))*e0bq",RooArgList(eb,eudsc,n,e0bq,e1bq,e2bq));
  RooFormulaVar p3bjets_tt("p3bjets_tt","p3bjets_tt","(eb*eb*(n-2)*eudsc*pow(1-eudsc,n-3)+2*eb*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4)+(n>4 ? pow((1-eb),2)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow((1-eudsc),n-5) : 0 ))*e2bq+(eb*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3)+(1-eb)*((n-1)*(n-2)*(n-3)/6)*pow(eudsc,3)*pow(1-eudsc,n-4))*e1bq+((n*(n-1)*(n-2)/6)*pow(eudsc,3)*pow(1-eudsc,n-3))*e0bq",RooArgList(eb,eudsc,n,e0bq,e1bq,e2bq));
  
  RooFormulaVar p0bjets_v ("p0bjets_v", "p0bjets_v" ,"pow(1-euds,n)",RooArgList(euds,n));
  RooFormulaVar p1bjets_v ("p1bjets_v", "p1bjets_v" ,"n*euds*pow(1-euds,n-1)",RooArgList(euds,n));
  RooFormulaVar p2bjets_v ("p2bjets_v", "p2bjets_v" ,"(n*(n-1)/2)*euds*euds*pow(1-euds,n-2)",RooArgList(euds,n));
  RooFormulaVar p3bjets_v ("p3bjets_v", "p3bjets_v" ,"((n)*(n-1)*(n-2)/6)*pow(euds,3)*pow(1-euds,n-3)",RooArgList(euds,n));
  
    // C o n s t r u c t   p . d . f 's
    // -------------------------------------------------------------------------------
  
	RooGenericPdf pbjets_tt("pbjets_tt","pbjets_tt","(nbjets==0)*p0bjets_tt+(nbjets==1)*p1bjets_tt+(nbjets==2)*p2bjets_tt+(nbjets==3)*p3bjets_tt",RooArgList(nbjets,p0bjets_tt,p1bjets_tt,p2bjets_tt,p3bjets_tt));
	RooExtendPdf  pbjets_tt_ext("pbjets_tt_ext","pbjets_tt_ext",pbjets_tt,Ntt);
  
	RooGenericPdf pbjets_v("pbjets_v","pbjets_v","(nbjets==0)*p0bjets_v+(nbjets==1)*p1bjets_v+(nbjets==2)*p2bjets_v+(nbjets==3)*p3bjets_v",RooArgList(nbjets,p0bjets_v,p1bjets_v,p2bjets_v,p3bjets_v));
	RooExtendPdf  pbjets_v_ext("pbjets_v_ext","pbjets_v_ext",pbjets_v,Nv);
  
  sprintf(name,"model_%d_wp_%d_jets",btag_wp_idx,Njets_[njets]);
  RooAddPdf model(name,"model",RooArgList(pbjets_tt,pbjets_v),RooArgList(Ntt,Nv));
	
    // C r e a t e   d a t a s e t 
    // -------------------------------------------------------------------------------
    // Sample a dataset for tt+jets events
  
	RooDataSet data("data","data",RooArgSet(nbjets)) ;
  for (Int_t i=0 ; i<GetPredNtotal(btag_wp_idx,njets,0) ; i++) { nbjets.setLabel("N0bjet")  ; data.add(RooArgSet(nbjets));}
  for (Int_t i=0 ; i<GetPredNtotal(btag_wp_idx,njets,1) ; i++) { nbjets.setLabel("N1bjet")  ; data.add(RooArgSet(nbjets));}
  for (Int_t i=0 ; i<GetPredNtotal(btag_wp_idx,njets,2) ; i++) { nbjets.setLabel("N2bjets") ; data.add(RooArgSet(nbjets));}
  for (Int_t i=0 ; i<GetPredNtotal(btag_wp_idx,njets,3) ; i++) { nbjets.setLabel("N3bjets") ; data.add(RooArgSet(nbjets));}
  
    // F i t   t h e   d a t a   a n d   c o n s t r u c t   t h e   l i k e l i h o o d   f u n c t i o n
    // ----------------------------------------------------------------------------------------------
  
  RooAbsReal* nll = model.createNLL(data);//,Optimize(0));
  
	RooMinimizer minimizer(*nll);
  
	//minimizer.optimizeConst(0) ; DO NOT SET IT TO TRUE, WILL NOT CONVERGE OTHERWISE
	minimizer.setPrintLevel(-1);
    //minimizer.setNoWarn();
  
    // Set algorithm
	minimizer.minimize("Minuit2", "Combined");
    //minimizer.minimize("GSLMultiMin", "ConjugateFR");
    //minimizer.minimize("GSLMultiMin", "BFGS2");
  
	if(doMinos) minimizer.minos();
	
	RooFitResult* fit_result = minimizer.save();
	
	if(Verbose)fit_result->Print("v");
  
	Ntt_[njets]                = Ntt.getVal();
	Ntt_err_[njets]            = Ntt.getError();
	Ntt_err_up_[njets]         = Ntt.getErrorHi();
	Ntt_err_down_[njets]       = Ntt.getErrorLo();
	Nv_[njets]                 = Nv.getVal();
	Nv_err_[njets]             = Nv.getError();
	Nv_err_up_[njets]          = Nv.getErrorHi();
	Nv_err_down_[njets]        = Nv.getErrorLo();
  
	minValue_[njets]           = fit_result->minNll();
  
	eb_[btag_wp_idx][njets]             = eb.getVal();
	eb_err_[btag_wp_idx][njets]         = eb.getError();
	eb_err_up_[btag_wp_idx][njets]      = eb.getErrorHi();
	eb_err_down_[btag_wp_idx][njets]    = eb.getErrorLo();
	eudsc_[btag_wp_idx][njets]          = eudsc.getVal();
	eudsc_err_[btag_wp_idx][njets]      = eudsc.getError();
	eudsc_err_up_[btag_wp_idx][njets]   = eudsc.getErrorHi();
	eudsc_err_down_[btag_wp_idx][njets] = eudsc.getErrorLo();
	euds_[btag_wp_idx][njets]           = euds.getVal();
	euds_err_[btag_wp_idx][njets]       = euds.getError();
	euds_err_up_[btag_wp_idx][njets]    = euds.getErrorHi();
	euds_err_down_[btag_wp_idx][njets]  = euds.getErrorLo();
	
	delete fit_result;
	delete nll;
	
    // C r e a t e   w o r k s p a c e ,   i m p o r t   d a t a   a n d   m o d e l,   a n d   s a v e   w s   i n   f i l e
    // ----------------------------------------------------------------------------------------------------------------------
  if(doWS){
      // Create a new empty workspace
	  sprintf(name,"w_%d_wp_%d_jets",btag_wp_idx,Njets_[njets]);
	  RooWorkspace *w = new RooWorkspace(name,"workspace") ;
    
      // Import model and all its components into the workspace
	  w->import(model) ;
    
      // Import data into the workspace
	  w->import(data) ;
      // Print workspace contents
      //w->Print() ;
    
      // Save the workspace into a ROOT file
	  sprintf(name,"VJetEstimation_RooFit_WS_%d_wp_%d_jets.root",btag_wp_idx,Njets_[njets]);
	  w->writeToFile(name) ;
    
  	delete w;
  }
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::UnBinnedMaximumJointWPLikelihoodEst(UInt_t njets, vector<Int_t> &FixedVarIdx, bool doMinos, bool doWS, bool Verbose){
  
    // Declare observables
  char name[100];
	TString wp[3] = {"loose","medium","tight"};
	vector<Int_t>::iterator Idx;
  
	RooRealVar Ntt("Ntt","N_{t#bar{t}-like}",init_Nttlike_[njets],init_Nttlike_[njets]/5,init_Nttlike_[njets]*5);
	RooRealVar Nv( "Nv", "N_{V-like}",       init_Nvlike_[njets],init_Nvlike_[njets]/5,init_Nvlike_[njets]*5);
  
	RooRealVar* eb[NbOfBtagWorkingPoint_];
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 0 );
	for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
		eb[i] = new RooRealVar( "eb_"+wp[i],"#epsilon_{b-tag} ("+wp[i]+")", init_Eb_[njets][i],0.0,1.0);
		if(Idx != FixedVarIdx.end()) eb[i]->setConstant(kTRUE);
	}
	RooRealVar* eudsc[NbOfBtagWorkingPoint_];
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 1 );
	for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
		eudsc[i] = new RooRealVar( "eudsc_"+wp[i],"#epsilon_{mis-tag} ("+wp[i]+")", init_Eudsc_[njets][i],0.0,1.0);
		if(Idx != FixedVarIdx.end()) eudsc[i]->setConstant(kTRUE);
	}
	RooRealVar* euds[NbOfBtagWorkingPoint_];
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 2 );
	for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
		euds[i]  = new RooRealVar( "euds_"+wp[i],"#epsilon^{,}_{mis-tag} ("+wp[i]+")", init_Euds_[njets][i],0.0,1.0);//GetPredEuds( i,njets),0.0,1.0);
		if(Idx != FixedVarIdx.end()) euds[i]->setConstant(kTRUE);
	}
  
    // Declare constants
  
	RooConstVar n("n","number of selected jets",Njets_[njets]) ;
	RooConstVar e0bq("e0bq","e0bq",e0bq_[njets]);
	RooConstVar e1bq("e1bq","e1bq",e1bq_[njets]);
	RooConstVar e2bq("e2bq","e2bq",e2bq_[njets]);
    //RooConstVar e3bq("e3bq","e3bq",e3bq_[njets]);
  
    // C o n s t r u c t   a   c a t e g o r y   w i t h   l a b e l s    a n d   i n d e c e s
    // -----------------------------------------------------------------------------------------
  
	RooCategory nbjets("nbjets","Number of b-jets");
	nbjets.defineType("N0bjet", 0);
	nbjets.defineType("N1bjet", 1);
	nbjets.defineType("N2bjets",2);
	nbjets.defineType("N3bjets",3);
  
	RooCategory WPCat("WPCat","B-tagging working point") ;
	for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++) WPCat.defineType(wp[i]) ;
  
    // C o n s t r u c t   f o r m u l a s
    // -------------------------------------------------------------------------------
  
  RooFormulaVar *p0bjets_tt[NbOfBtagWorkingPoint_];
  RooFormulaVar *p1bjets_tt[NbOfBtagWorkingPoint_];
  RooFormulaVar *p2bjets_tt[NbOfBtagWorkingPoint_];
  RooFormulaVar *p3bjets_tt[NbOfBtagWorkingPoint_];
  
  for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
    p0bjets_tt[i] = new RooFormulaVar("p0bjets_tt_"+wp[i],"p0bjets_tt_"+wp[i],"(1-@0)*(1-@0)*pow((1-@1),@2-2)*e2bq+(1-@0)*pow((1-@1),@2-1)*e1bq+pow((1-@1),@2)*e0bq",RooArgList((*eb[i]),(*eudsc[i]),n,e0bq,e1bq,e2bq));
    p1bjets_tt[i] = new RooFormulaVar("p1bjets_tt_"+wp[i],"p1bjets_tt_"+wp[i],"(2*@0*(1-@0)*pow(1-@1,@2-2)+(1-@0)*(1-@0)*(@2-2)*@1*pow(1-@1,@2-3))*e2bq+(@0*pow(1-@1,@2-1)+(1-@0)*(@2-1)*@1*pow(1-@1,@2-2))*e1bq+(@2*@1*pow(1-@1,@2-1))*e0bq",RooArgList((*eb[i]),(*eudsc[i]),n,e0bq,e1bq,e2bq));
    p2bjets_tt[i] = new RooFormulaVar("p2bjets_tt_"+wp[i],"p2bjets_tt_"+wp[i],"(@0*@0*pow(1-@1,@2-2)+2*@0*(1-@0)*(@2-2)*@1*pow(1-@1,@2-3)+(1-@0)*(1-@0)*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4))*e2bq+(@0*(@2-1)*@1*pow(1-@1,@2-2)+(1-@0)*((@2-1)*(@2-2)/2)*@1*@1*pow(1-@1,@2-3))*e1bq+((@2*(@2-1)/2)*@1*@1*pow(1-@1,@2-2))*e0bq",RooArgList((*eb[i]),(*eudsc[i]),n,e0bq,e1bq,e2bq));
    p3bjets_tt[i] = new RooFormulaVar("p3bjets_tt_"+wp[i],"p3bjets_tt_"+wp[i],"(@0*@0*(@2-2)*@1*pow(1-@1,@2-3)+2*@0*(1-@0)*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4)+(@2>4 ? pow((1-@0),2)*((@2-2)*(@2-3)*(@2-4)/6)*pow(@1,3)*pow((1-@1),@2-5) : 0 ))*e2bq+(@0*((@2-1)*(@2-2)/2)*@1*@1*pow(1-@1,@2-3)+(1-@0)*((@2-1)*(@2-2)*(@2-3)/6)*pow(@1,3)*pow(1-@1,@2-4))*e1bq+((@2*(@2-1)*(@2-2)/6)*pow(@1,3)*pow(1-@1,@2-3))*e0bq",RooArgList((*eb[i]),(*eudsc[i]),n,e0bq,e1bq,e2bq));
  }
  
  RooFormulaVar *p0bjets_v[NbOfBtagWorkingPoint_]; 
  RooFormulaVar *p1bjets_v[NbOfBtagWorkingPoint_]; 
  RooFormulaVar *p2bjets_v[NbOfBtagWorkingPoint_]; 
  RooFormulaVar *p3bjets_v[NbOfBtagWorkingPoint_]; 
  
  for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
    p0bjets_v[i] = new RooFormulaVar("p0bjets_v_"+wp[i], "p0bjets_v_"+wp[i] ,"pow(1-@0,@1)",RooArgList((*euds[i]),n));
    p1bjets_v[i] = new RooFormulaVar("p1bjets_v_"+wp[i], "p1bjets_v_"+wp[i] ,"@1*@0*pow(1-@0,@1-1)",RooArgList((*euds[i]),n));
    p2bjets_v[i] = new RooFormulaVar("p2bjets_v_"+wp[i], "p2bjets_v_"+wp[i] ,"(@1*(@1-1)/2)*@0*@0*pow(1-@0,@1-2)",RooArgList((*euds[i]),n));
    p3bjets_v[i] = new RooFormulaVar("p3bjets_v_"+wp[i], "p3bjets_v_"+wp[i] ,"((@1)*(@1-1)*(@1-2)/6)*pow(@0,3)*pow(1-@0,@1-3)",RooArgList((*euds[i]),n));
  }
  
    // C o n s t r u c t   p . d . f 's
    // -------------------------------------------------------------------------------
  
  RooGenericPdf *pbjets_tt[NbOfBtagWorkingPoint_];
  RooGenericPdf *pbjets_v[NbOfBtagWorkingPoint_];
  RooExtendPdf  *pbjets_v_ext[NbOfBtagWorkingPoint_];
  RooExtendPdf  *pbjets_tt_ext[NbOfBtagWorkingPoint_];
  
  for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
    pbjets_tt[i] = new RooGenericPdf("pbjets_tt_"+wp[i],"pbjets_tt_"+wp[i],"(@0==0)*@1+(@0==1)*@2+(@0==2)*@3+(@0==3)*@4",RooArgList(nbjets,*p0bjets_tt[i],*p1bjets_tt[i],*p2bjets_tt[i],*p3bjets_tt[i]));
    pbjets_v[i]  = new RooGenericPdf("pbjets_v_"+wp[i],"pbjets_v_"+wp[i],"((@0)==0)*@1+(@0==1)*@2+(@0==2)*@3+(@0==3)*@4",RooArgList(nbjets,*p0bjets_v[i],*p1bjets_v[i],*p2bjets_v[i],*p3bjets_v[i]));
    
    pbjets_tt_ext[i] = new RooExtendPdf("pbjets_tt_ext_"+wp[i],"pbjets_tt_ext_"+wp[i],*pbjets_tt[i],Ntt);
    pbjets_v_ext[i]  = new RooExtendPdf("pbjets_v_ext_"+wp[i] ,"pbjets_v_ext_"+wp[i] ,*pbjets_v[i],Nv);
  }
  
  RooAddPdf* model[NbOfBtagWorkingPoint_];
  for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
	  sprintf(name,"model_%d_wp_%d_jets",i,Njets_[njets]);
    model[i] = new RooAddPdf(name,"model_"+wp[i],RooArgList(*pbjets_tt_ext[i],*pbjets_v_ext[i]));//,RooArgList(Ntt,Nv));
  }
  
	sprintf(name,"model_joint_wp_%d_jets",Njets_[njets]);
  RooSimultaneous simPdf(name,"simultaneous pdf",WPCat) ;
  for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++) simPdf.addPdf(*model[i],wp[i]) ;
	
    // C r e a t e   d a t a s e t 
    // -------------------------------------------------------------------------------
  
	RooDataSet data("data","data",RooArgSet(nbjets,WPCat)) ;
	for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
		WPCat.setLabel(wp[i]);
		for(Int_t j=0; j<GetPredNtotal(i,njets,0); j++){nbjets.setLabel("N0bjet") ; data.add(RooArgSet(nbjets,WPCat));}
		for(Int_t j=0; j<GetPredNtotal(i,njets,1); j++){nbjets.setLabel("N1bjet") ; data.add(RooArgSet(nbjets,WPCat));}
		for(Int_t j=0; j<GetPredNtotal(i,njets,2); j++){nbjets.setLabel("N2bjets"); data.add(RooArgSet(nbjets,WPCat));}
		for(Int_t j=0; j<GetPredNtotal(i,njets,3); j++){nbjets.setLabel("N3bjets"); data.add(RooArgSet(nbjets,WPCat));}
	}
  
    // F i t   t h e   d a t a   a n d   c o n s t r u c t   t h e   l i k e l i h o o d   f u n c t i o n
    // ----------------------------------------------------------------------------------------------
  
	RooAbsReal* nll = simPdf.createNLL(data);//,Optimize(0));
  
	RooMinimizer minimizer(*nll);
    //RooMinuit minimizer(*nll);
  
  //minimizer.optimizeConst(0) ; DO NOT SET IT TO TRUE, WILL NOT CONVERGE OTHERWISE
	minimizer.setPrintLevel(-1);
    //minimizer.setNoWarn();
  
    // Set algorithm
	minimizer.minimize("Minuit2", "Combined");
    //minimizer.minimize("GSLMultiMin", "ConjugateFR");
    //minimizer.minimize("GSLMultiMin", "BFGS2");
  
	if(doMinos) minimizer.minos();
	
	RooFitResult* fit_result = minimizer.save();
	
	if(Verbose)fit_result->Print("v");
  
	Ntt_[njets]                = Ntt.getVal();
	Ntt_err_[njets]            = Ntt.getError();
	Ntt_err_up_[njets]         = Ntt.getErrorHi();
	Ntt_err_down_[njets]       = Ntt.getErrorLo();
	Nv_[njets]                 = Nv.getVal();
	Nv_err_[njets]             = Nv.getError();
	Nv_err_up_[njets]          = Nv.getErrorHi();
	Nv_err_down_[njets]        = Nv.getErrorLo();
	minValue_[njets]           = fit_result->minNll();
	for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
		eb_[i][njets]        = eb[i]->getVal();
		eb_err_[i][njets]    = eb[i]->getError();
		eb_err_up_[i][njets] = eb[i]->getErrorHi();
		eb_err_down_[i][njets]    = eb[i]->getErrorLo();
		eudsc_[i][njets]     = eudsc[i]->getVal();
		eudsc_err_[i][njets] = eudsc[i]->getError();
		eudsc_err_up_[i][njets]   = eudsc[i]->getErrorHi();
		eudsc_err_down_[i][njets] = eudsc[i]->getErrorLo();
		euds_[i][njets]      = euds[i]->getVal();
		euds_err_[i][njets]  = euds[i]->getError();
		euds_err_up_[i][njets]    = euds[i]->getErrorHi();
		euds_err_down_[i][njets]  = euds[i]->getErrorLo();
	}
  
  
	for(UInt_t i=0; i<NbOfBtagWorkingPoint_; i++){
		delete eb[i];
		delete eudsc[i];
		delete euds[i];
		delete pbjets_tt[i];
		delete pbjets_tt_ext[i];
		delete pbjets_v[i];
		delete pbjets_v_ext[i];
	}
	delete fit_result;
  
    // C r e a t e   w o r k s p a c e ,   i m p o r t   d a t a   a n d   m o d e l,   a n d   s a v e   w s   i n   f i l e
    // ----------------------------------------------------------------------------------------------------------------------
  if(doWS){
      // Create a new empty workspace
  	sprintf(name,"w_joint_wp_%d_jets",Njets_[njets]);
  	RooWorkspace *w = new RooWorkspace(name,"workspace") ;
    
      // Import model and all its components into the workspace
  	w->import(simPdf) ;
    
      // Import data into the workspace
  	w->import(data) ;
      // Print workspace contents
      //w->Print() ;
    
      // Save the workspace into a ROOT file
  	sprintf(name,"VJetEstimation_RooFit_WS_JointWP_%d_jets.root",Njets_[njets]);
  	w->writeToFile(name) ;
    
  	delete w;
  }
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::PrintInputs(){
	for(UInt_t i=0;i<NbOfJetsBins_;i++) PrintInputs(i);
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::PrintResults(){
	cout<<endl;
	cout<<"***********************************************"<<endl;
	cout<<"*****      VJets Estimation Results         ***"<<endl;
	cout<<"***********************************************"<<endl;
  
	for(UInt_t i=0;i<NbOfJetsBins_;i++) PrintResults(i);
  
	cout<<"***********************************************************"<<endl;
	cout<<"************        ALL INCLUSIF             **************"<<endl; 
  cout<<std::fixed<<setprecision(4);
	cout<<"Btag working point : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<BtagWorkingPoint_[i]<<" || ";} cout<<endl;
	cout<<"Estimated/Predicted eb	 : ";	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstEb(i)   <<" +/- "<<GetEstEbErr(i)	<<" / "<<GetPredEb(i)   <<" || ";} cout<<endl;
	cout<<"Estimated/Predicted eudsc : ";	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstEudsc(i)<<" +/- "<<GetEstEudscErr(i)	<<" / "<<GetPredEudsc(i)<<" || ";} cout<<endl;
	cout<<"Estimated/Predicted euds  : ";	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstEuds(i) <<" +/- "<<GetEstEudsErr(i)	<<" / "<<GetPredEuds(i) <<" || ";} cout<<endl;
	cout<<std::fixed<<setprecision(5);
	cout<<"Estimated Ntt-like : "<<GetEstNtt(0)<<" / Predicted value from MC : "<< GetPredNtt(0)<<endl;
	cout<<"- 0 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNtt(i,-1,0)<<" +/- "<<GetEstNttErr(i,-1,0)<<" / "<< GetPredNtt(i,-1,0)<<" +/- "<<GetPredNttErr(i,-1,0) <<" || ";} cout<<endl;
	cout<<"- 1 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNtt(i,-1,1)<<" +/- "<<GetEstNttErr(i,-1,1)<<" / "<< GetPredNtt(i,-1,1)<<" +/- "<<GetPredNttErr(i,-1,1) <<" || ";} cout<<endl;
	cout<<"- 2 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNtt(i,-1,2)<<" +/- "<<GetEstNttErr(i,-1,2)<<" / "<< GetPredNtt(i,-1,2)<<" +/- "<<GetPredNttErr(i,-1,2) <<" || ";} cout<<endl;
	cout<<"- 3 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNtt(i,-1,3)<<" +/- "<<GetEstNttErr(i,-1,3)<<" / "<< GetPredNtt(i,-1,3)<<" +/- "<<GetPredNttErr(i,-1,3) <<" || ";} cout<<endl;
    // 	cout<<"Estimated Nvb-like : "<<GetEstNvb(0)<<" / Predicted value from MC : "<< GetPredNvb(0)<<endl;
    // 	cout<<"- 0 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNvb(i,-1,0)<<" / "<< GetPredNvb(i,-1,0) <<" || ";} cout<<endl;
    // 	cout<<"- 1 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNvb(i,-1,1)<<" / "<< GetPredNvb(i,-1,1) <<" || ";} cout<<endl;
    // 	cout<<"- 2 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNvb(i,-1,2)<<" / "<< GetPredNvb(i,-1,2) <<" || ";} cout<<endl;
    // 	cout<<"- 3 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNvb(i,-1,3)<<" / "<< GetPredNvb(i,-1,3) <<" || ";} cout<<endl;
	cout<<"Estimated Nv-like  : "<<GetEstNv(0)<<" / Predicted value from MC : "<< GetPredNv(0)<<endl;
	cout<<"- 0 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNv(i,-1,0)<<" +/- "<<GetEstNvErr(i,-1,0)<<" / "<<GetPredNv(i,-1,0)<<" +/- "<<GetPredNvErr(i,-1,0) <<" || ";} cout<<endl;
	cout<<"- 1 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNv(i,-1,1)<<" +/- "<<GetEstNvErr(i,-1,1)<<" / "<<GetPredNv(i,-1,1)<<" +/- "<<GetPredNvErr(i,-1,1) <<" || ";} cout<<endl;
	cout<<"- 2 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNv(i,-1,2)<<" +/- "<<GetEstNvErr(i,-1,2)<<" / "<<GetPredNv(i,-1,2)<<" +/- "<<GetPredNvErr(i,-1,2) <<" || ";} cout<<endl;
	cout<<"- 3 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNv(i,-1,3)<<" +/- "<<GetEstNvErr(i,-1,3)<<" / "<<GetPredNv(i,-1,3)<<" +/- "<<GetPredNvErr(i,-1,3) <<" || ";} cout<<endl;
	cout<<"***********************************************************"<<endl;
  cout << resetiosflags( std::ios_base::fixed ) << setprecision(-1);
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::PrintResults_LatexFormat(ostream &ofile, Int_t Eff_prec, Int_t N_prec){
  for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {
    ofile<<"Btag working point : "<<BtagWorkingPoint_[i]<<endl;
    ofile<<"\\begin{table}"<<endl;
    ofile<<"	\\centering"<<endl;
    ofile<<"		\\begin{supertabular}{l|";for(UInt_t j=0;j<NbOfJetsBins_+1;j++){ofile<<"c";};ofile<<"}"<<endl;
    ofile<<"			\\hline"<<endl;
    ofile<<" \\multicolumn{"<<NbOfJetsBins_+1<<"}{l}{Estimation / Monte Carlo prediction : }\\\\"<<endl;
    ofile<<"			      &";for(UInt_t j=0;j<NbOfJetsBins_;j++){ofile<<"$"<<Njets_[j]<<"$ jets & ";}ofile<<"Inclusive \\\\"<<endl;
    ofile<<"\\hline"<<endl;
    ofile<<" $\\epsilon_b$                &$"<<std::fixed<<setprecision(Eff_prec);for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstEb(   i,j)<<" \\pm "<<GetEstEbErr(   i,j)<<" $/$ "<<GetPredEb(   i,j)<<" \\pm "<<GetPredEbErr(   i,j)<<"$&$";}ofile<<GetEstEb(    i)<<" \\pm "<<GetEstEbErr(   i)<<"$/$"<<GetPredEb(   i)<<" \\pm "<<GetPredEbErr(   i)<<"$ \\\\"<<endl;
    ofile<<" $\\epsilon_{eudsc}$          &$"<<std::fixed<<setprecision(Eff_prec);for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstEudsc(i,j)<<"\\pm "<<GetEstEudscErr(i,j)<<" $/$ "<<GetPredEudsc(i,j)<<" \\pm "<<GetPredEudscErr(i,j)<<"$&$";}ofile<<GetEstEudsc( i)<<" \\pm "<<GetEstEudscErr(i)<<"$/$"<<GetPredEudsc(i)<<" \\pm "<<GetPredEudscErr(i)<<"$ \\\\"<<endl;
    ofile<<" $\\epsilon_{euds} $          &$"<<std::fixed<<setprecision(Eff_prec);for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstEuds( i,j)<<"\\pm "<<GetEstEudsErr( i,j)<<" $/$ "<<GetPredEuds( i,j)<<" \\pm "<<GetPredEudsErr( i,j)<<"$&$";}ofile<<GetEstEuds(  i)<<" \\pm "<<GetEstEudsErr( i)<<"$/$"<<GetPredEuds( i)<<" \\pm "<<GetPredEudsErr( i)<<"$ \\\\"<<endl;
    ofile<< std::fixed << setprecision(N_prec);        
    ofile<<"\\hline"<<endl;
    ofile<<" \\multicolumn{"<<NbOfJetsBins_+1<<"}{l}{$t\\bar{t}+jets$ and single top : }\\\\"<<endl;
    ofile<<"\\hline"<<endl;
    ofile<<" $0$ b-jet                    &$";for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNtt(i,j,0)<<" \\pm "<<GetEstNttErr(i,j,0)<<" $/$ "<<GetPredNtt(i,j,0)<<" \\pm "<<GetPredNttErr(i,j,0) <<"$&$";}ofile<<GetEstNtt(i,-1,0)<<" \\pm "<<GetEstNttErr(i,-1,0)<<" $/$ "<<GetPredNtt(i,-1,0)<<" \\pm "<<GetPredNttErr(i,-1,0)<<"$ \\\\"<<endl;
    ofile<<" $1$ b-jet                    &$";for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNtt(i,j,1)<<" \\pm "<<GetEstNttErr(i,j,1)<<" $/$ "<<GetPredNtt(i,j,1)<<" \\pm "<<GetPredNttErr(i,j,1) <<"$&$";}ofile<<GetEstNtt(i,-1,1)<<" \\pm "<<GetEstNttErr(i,-1,1)<<" $/$ "<<GetPredNtt(i,-1,1)<<" \\pm "<<GetPredNttErr(i,-1,1)<<"$ \\\\"<<endl;
    ofile<<" $2$ b-jets                   &$";for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNtt(i,j,2)<<" \\pm "<<GetEstNttErr(i,j,2)<<" $/$ "<<GetPredNtt(i,j,2)<<" \\pm "<<GetPredNttErr(i,j,2) <<"$&$";}ofile<<GetEstNtt(i,-1,2)<<" \\pm "<<GetEstNttErr(i,-1,2)<<" $/$ "<<GetPredNtt(i,-1,2)<<" \\pm "<<GetPredNttErr(i,-1,2)<<"$ \\\\"<<endl;
    ofile<<" $3$ b-jets                   &$";for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNtt(i,j,3)<<" \\pm "<<GetEstNttErr(i,j,3)<<" $/$ "<<GetPredNtt(i,j,3)<<" \\pm "<<GetPredNttErr(i,j,3) <<"$&$";}ofile<<GetEstNtt(i,-1,3)<<" \\pm "<<GetEstNttErr(i,-1,3)<<" $/$ "<<GetPredNtt(i,-1,3)<<" \\pm "<<GetPredNttErr(i,-1,3)<<"$ \\\\"<<endl;
    ofile<<"\\hline"<<endl;
    ofile<<" Inclusive                    &$";for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNtt(i,j)  <<" \\pm "<<GetEstNttErr(i,j)  <<" $/$ "<<GetPredNtt(i,j)  <<" \\pm "<<GetPredNttErr(i,j)   <<"$&$";}ofile<<GetEstNtt(i)     <<" \\pm "<<GetEstNttErr(i)     <<" $/$ "<<GetPredNtt(i)     <<" \\pm "<<GetPredNttErr(i)<<"$ \\\\"<<endl;
    ofile<<"\\hline"<<endl;
    ofile<<"\\hline"<<endl;
    ofile<<" \\multicolumn{"<<NbOfJetsBins_+1<<"}{l}{Multijet and $W/Z+jets$ : }\\\\"<<endl;
    ofile<<"\\hline"<<endl;
    ofile<<" $0$ b-jet                    &$";for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNv(i,j,0)<<" \\pm "<<GetEstNvErr(i,j,0)<<" $/$ "<<GetPredNv(i,j,0)<<" \\pm "<<GetPredNvErr(i,j,0) <<"$&$";}ofile<<GetEstNv(i,-1,0)<<" \\pm "<<GetEstNvErr(i,-1,0)<<" $/$ "<<GetPredNv(i,-1,0)<<" \\pm "<<GetPredNvErr(i,-1,0)<<"$ \\\\"<<endl;
    ofile<<" $1$ b-jet                    &$";for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNv(i,j,1)<<" \\pm "<<GetEstNvErr(i,j,1)<<" $/$ "<<GetPredNv(i,j,1)<<" \\pm "<<GetPredNvErr(i,j,1) <<"$&$";}ofile<<GetEstNv(i,-1,1)<<" \\pm "<<GetEstNvErr(i,-1,1)<<" $/$ "<<GetPredNv(i,-1,1)<<" \\pm "<<GetPredNvErr(i,-1,1)<<"$ \\\\"<<endl;
    ofile<<" $2$ b-jets                   &$";for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNv(i,j,2)<<" \\pm "<<GetEstNvErr(i,j,2)<<" $/$ "<<GetPredNv(i,j,2)<<" \\pm "<<GetPredNvErr(i,j,2) <<"$&$";}ofile<<GetEstNv(i,-1,2)<<" \\pm "<<GetEstNvErr(i,-1,2)<<" $/$ "<<GetPredNv(i,-1,2)<<" \\pm "<<GetPredNvErr(i,-1,2)<<"$ \\\\"<<endl;
    ofile<<" $3$ b-jets                   &$";for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNv(i,j,3)<<" \\pm "<<GetEstNvErr(i,j,3)<<" $/$ "<<GetPredNv(i,j,3)<<" \\pm "<<GetPredNvErr(i,j,3) <<"$&$";}ofile<<GetEstNv(i,-1,3)<<" \\pm "<<GetEstNvErr(i,-1,3)<<" $/$ "<<GetPredNv(i,-1,3)<<" \\pm "<<GetPredNvErr(i,-1,3)<<"$ \\\\"<<endl;
    ofile<<"\\hline"<<endl;
    ofile<<" Inclusive                    &$";for(UInt_t j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNv(i,j)  <<" \\pm "<<GetEstNvErr(i,j)  <<" $/$ "<<GetPredNv(i,j)  <<" \\pm "<<GetPredNvErr(i,j)   <<"$&$";}ofile<<GetEstNv(i)     <<" \\pm "<<GetEstNvErr(i)     <<" $/$ "<<GetPredNv(i)     <<" \\pm "<<GetPredNvErr(i)<<"$ \\\\"<<endl;
    ofile<<"\\hline"<<endl;
    ofile<<"\\hline"<<endl;
    ofile<<"		\\end{supertabular}"<<endl;
    ofile<<"	\\caption{}"<<endl;
    ofile<<"	\\label{tab:}"<<endl;
    ofile<<"\\end{table}"<<endl;
	}
  ofile << resetiosflags( std::ios_base::fixed ) << setprecision(-1);
}

void VJetEstimation::PrintResults_LatexFormat(Int_t Eff_prec, Int_t N_prec) {
  this->PrintResults_LatexFormat(std::cout, Eff_prec, N_prec);  
}

/**________________________________________________________________________________________________________________*/
bool VJetEstimation::CheckEstimation(Double_t threshold, UInt_t njets){ return (minValue_[njets]>threshold ? false : true);}

/**________________________________________________________________________________________________________________*/
bool VJetEstimation::CheckEstimation(Double_t threshold){
	for(UInt_t i=0;i<NbOfJetsBins_;i++)
		if(minValue_[i]>threshold) return false;
	return true;
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::CheckEstimationLinearity(vector<Int_t> &FixedVarIdx, bool doMinos, Int_t NbOfRescaleFact, Double_t **RescaleFact, Int_t Idx){
  
  RescaledTTLikeEstimation = new TGraphErrors(NbOfRescaleFact);
  RescaledVLikeEstimation  = new TGraphErrors(NbOfRescaleFact);
  tCanva_RescaledTTLikeEstimation = new TCanvas("tCanva_RescaledTTLikeEstimation","",-1);
  gROOT->GetListOfCanvases()->RemoveLast();
  tCanva_RescaledVLikeEstimation  = new TCanvas("tCanva_RescaledVLikeEstimation", "",-1);
  gROOT->GetListOfCanvases()->RemoveLast();
  
  /*
   vector<Int_t> Indeces;
   Double_t ****n_Init  = new Double_t***[NbOfBtagWorkingPoint_];
   Double_t ****n       = new Double_t***[NbOfBtagWorkingPoint_];
   for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
   n_Init[i] = new Double_t**[NbOfJetsBins_];
   n[i]      = new Double_t**[NbOfJetsBins_];
   for(UInt_t j=0;j<NbOfJetsBins_;j++){
   n_Init[i][j] = new Double_t*[NbOfDatasets_];
   n[i][j]      = new Double_t*[NbOfDatasets_];
   for(UInt_t k=0;k<NbOfBJetsBins_;k++){
   n_Init[i][j][k] = new Double_t[NbOfDatasets_];
   n[i][j][k]      = new Double_t[NbOfDatasets_];
   for(UInt_t l=0;l<NbOfDatasets_;l++) n_Init[i][j][k][l] = N_[i][j][k][l];
   }			
   }
   }
   cout<<"*****************Linearity check : results *****************"<<endl;
   for(Int_t nb=0;nb<NbOfRescaleFact;nb++){
   cout<<" --- ";
   for(UInt_t i=0;i<NbOfJetsBins_;i++){
   for(UInt_t j=0;j<NbOfBtagWorkingPoint_;j++){
   for(UInt_t k=0;k<NbOfBJetsBins_;k++){
   for(UInt_t l=0;l<NbOfDatasets_;l++){
   /// Multiply the original numbers of events (per datasets) 
   /// by a given factor "RescalFact"
   n[i][j][k][l] = n_Init[i][j][k][l]*RescaleFact[nb][l];
   }
   }
   }
   }
   //Estimate the number of tt-like and w-like events with the modified input numbers
   string method = "Minuit2";
   string option = "Combined";
   
   this->FillInputs(n);
   this->UnBinnedMaximumLikelihoodEst(FixedVarIdx, doMinos, false, false);
   
   RescaledTTLikeEstimation->SetPoint(nb,100*RescaleFact[nb][Idx],(Double_t)this->GetEstNtt(0));
   RescaledTTLikeEstimation->SetPointError(nb,0,RescaleFact[nb][Idx]*((Double_t)this->GetEstNttErr(0)));
   RescaledVLikeEstimation->SetPoint(nb,100*RescaleFact[nb][Idx],(Double_t)this->GetEstNv(0));
   RescaledVLikeEstimation->SetPointError(nb,0, RescaleFact[nb][Idx]*((Double_t)this->GetEstNvErr(0)));
   
   cout<<"-- Results :"<<endl;
   cout<<"--- V-like estimate  = "<<this->GetEstNv(0) <<"+/- "<<RescaleFact[nb][Idx]*((Double_t)this->GetEstNvErr(0))<<endl;
   cout<<"--- tt-like estimate = "<<this->GetEstNtt(0)<<"+/- "<<RescaleFact[nb][Idx]*((Double_t)this->GetEstNttErr(0))<<endl;
   cout<<"------------"<<endl;
   }
   
   RescaledTTLikeEstimation ->SetNameTitle("RescaledTTLikeEstimation","");
   RescaledVLikeEstimation  ->SetNameTitle("RescaledVLikeEstimation","");
   
   tCanva_RescaledTTLikeEstimation->cd();
   RescaledTTLikeEstimation->Draw("AC*");
   RescaledTTLikeEstimation->GetXaxis()->SetTitle("Rescaling factor (%)");
   RescaledTTLikeEstimation->GetYaxis()->SetTitle("Estimated nb of TT-like events");
   
   tCanva_RescaledVLikeEstimation->cd();
   RescaledVLikeEstimation->Draw("AC*");
   RescaledVLikeEstimation->GetXaxis()->SetTitle("Rescaling factor (%)");
   RescaledVLikeEstimation->GetYaxis()->SetTitle("Estimated nb of V-like events");
   
   // Put the number of events back to their original values
   this->FillInputs(n_Init);
   
   for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++){
   for(UInt_t j=0;j<NbOfJetsBins_;j++){
   for(UInt_t k=0;k<NbOfBJetsBins_;k++){
   delete [] n_Init[i][j][k];
   delete [] n[i][j][k];
   }
   delete [] n_Init[i][j];
   delete [] n[i][j];
   }
   delete[] n_Init[i];
   delete[] n[i];
   }
   */
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEb(Int_t wp) const{
	Double_t Effb = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) Effb += GetPredEb(wp,i)*GetPredNtotal(wp,i);
	if(GetPredNtotal(wp)>0) Effb /= GetPredNtotal(wp) ;
	return Effb;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEb(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eb_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEbErr(Int_t wp) const{
	Double_t EffbErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) EffbErr += pow(GetPredEbErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffbErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEbErr(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eb_err_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetPredEb(Int_t wp, Int_t njets, Double_t value) {
	eb_mc_[wp][njets] = value;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEudsc(Int_t wp) const{
	Double_t Effudsc = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) Effudsc += GetPredEudsc(wp,i)*GetPredNtotal(wp,i);
	Effudsc /= GetPredNtotal(wp);
	return Effudsc;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEudsc(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eudsc_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEudscErr(Int_t wp) const{
	Double_t EffudscErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) EffudscErr += pow(GetPredEudscErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffudscErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEudscErr(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eudsc_err_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetPredEudsc(Int_t wp, Int_t njets, Double_t value) {
	eudsc_mc_[wp][njets] = value;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEuds(Int_t wp) const{
	Double_t Effuds = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) Effuds += GetPredEuds(wp,i)*GetPredNtotal(wp,i);
	Effuds /= GetPredNtotal(wp);
	return Effuds;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEuds(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? euds_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEudsErr(Int_t wp) const{
	Double_t EffudsErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) EffudsErr += pow(GetPredEudsErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffudsErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredEudsErr(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? euds_err_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetPredEuds(Int_t wp, Int_t njets, Double_t value) {
	euds_mc_[wp][njets] = value;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNv(Int_t wp) const{
	Double_t Nvlike = 0;
	for(UInt_t i=0;i<iDatasetsVLike_.size();i++){
		for(UInt_t j=0; j<NbOfJetsBins_; j++){
			for(UInt_t k=0; k<NbOfBJetsBins_; k++) Nvlike  += N_[wp][j][k][iDatasetsVLike_[i]];
		}
	}
	return Nvlike;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNv(Int_t wp, Int_t njets) const{
	Double_t Nvlike = 0;
	for(UInt_t i=0;i<iDatasetsVLike_.size();i++){
		for(UInt_t j=0; j<NbOfBJetsBins_; j++) Nvlike  += N_[wp][njets][j][iDatasetsVLike_[i]];
	}
	return Nvlike;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNv(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t Nvlike = 0;
	if(njets == -1){
		for(UInt_t i=0;i<iDatasetsVLike_.size();i++){
			for(UInt_t j=0;j<NbOfJetsBins_;j++) Nvlike += N_[wp][j][nbjets][iDatasetsVLike_[i]];
		}
		return Nvlike;
	}
	else{
		for(UInt_t i=0;i<iDatasetsVLike_.size();i++){
			Nvlike  += N_[wp][njets][nbjets][iDatasetsVLike_[i]];
		}
		return Nvlike;
	}
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNvErr(Int_t wp) const{
	Double_t NvErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) NvErr += pow(GetPredNvErr(wp,i),2);
	return sqrt(NvErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNvErr(Int_t wp, Int_t njets) const{
	Double_t Nvlike = 0;
	for(UInt_t i=0;i<iDatasetsVLike_.size();i++){
		for(UInt_t j=0; j<NbOfBJetsBins_; j++) Nvlike  += pow(N_err_[wp][njets][j][iDatasetsVLike_[i]],2);
	}
	return sqrt(Nvlike);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNvErr(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t Nvlike = 0;
	if(njets == -1){
		for(UInt_t i=0;i<iDatasetsVLike_.size();i++){
			for(UInt_t j=0;j<NbOfJetsBins_;j++) Nvlike += pow(N_err_[wp][j][nbjets][iDatasetsVLike_[i]],2);
		}
		return sqrt(Nvlike);
	}
	else{
		for(UInt_t i=0;i<iDatasetsVLike_.size();i++){
			Nvlike  += pow(N_err_[wp][njets][nbjets][iDatasetsVLike_[i]],2);
		}
		return sqrt(Nvlike);
	}
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNvb(Int_t wp) const{
	Double_t Nvblike = 0;
	for(UInt_t i=0;i<iDatasetsVbLike_.size();i++){
		for(UInt_t j=0; j<NbOfJetsBins_; j++){
			for(UInt_t k=0; k<NbOfBJetsBins_; k++) Nvblike  += N_[wp][j][k][iDatasetsVbLike_[i]];
		}
	}
	return Nvblike;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNvb(Int_t wp, Int_t njets) const{
	Double_t Nvblike = 0;
	for(UInt_t i=0;i<iDatasetsVbLike_.size();i++){
		for(UInt_t j=0; j<NbOfBJetsBins_; j++) Nvblike  += N_[wp][njets][j][iDatasetsVbLike_[i]];
	}
	return Nvblike;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNvb(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t Nvblike = 0;
	if(njets == -1){
		for(UInt_t i=0;i<iDatasetsVbLike_.size();i++){
			for(UInt_t j=0;j<NbOfJetsBins_;j++) Nvblike += N_[wp][j][nbjets][iDatasetsVbLike_[i]];
		}
		return Nvblike;
	}
	else{
		for(UInt_t i=0;i<iDatasetsVbLike_.size();i++){
			Nvblike  += N_[wp][njets][nbjets][iDatasetsVbLike_[i]];
		}
		return Nvblike;
	}
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNvbErr(Int_t wp) const{
	Double_t NvbErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) NvbErr += pow(GetPredNvbErr(wp,i),2);
	return sqrt(NvbErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNvbErr(Int_t wp, Int_t njets) const{
	Double_t Nvblike = 0;
	for(UInt_t i=0;i<iDatasetsVbLike_.size();i++){
		for(UInt_t j=0; j<NbOfBJetsBins_; j++) Nvblike  += pow(N_err_[wp][njets][j][iDatasetsVbLike_[i]],2);
	}
	return sqrt(Nvblike);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNvbErr(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t Nvblike = 0;
	if(njets == -1){
		for(UInt_t i=0;i<iDatasetsVbLike_.size();i++){
			for(UInt_t j=0;j<NbOfJetsBins_;j++) Nvblike += pow(N_err_[wp][j][nbjets][iDatasetsVbLike_[i]],2);
		}
		return sqrt(Nvblike);
	}
	else{
		for(UInt_t i=0;i<iDatasetsVbLike_.size();i++){
			Nvblike  += pow(N_err_[wp][njets][nbjets][iDatasetsVbLike_[i]],2);
		}
		return sqrt(Nvblike);
	}
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNtt(Int_t wp) const{
	Double_t Nttlike = 0;
	for(UInt_t i=0;i<iDatasetsTTLike_.size();i++){
		for(UInt_t j=0; j<NbOfJetsBins_; j++){
			for(UInt_t k=0; k<NbOfBJetsBins_; k++) Nttlike  += N_[wp][j][k][iDatasetsTTLike_[i]];
		}
	}
	return Nttlike;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNtt(Int_t wp, Int_t njets) const{
	Double_t Nttlike = 0;
	for(UInt_t i=0;i<iDatasetsTTLike_.size();i++){
		for(UInt_t j=0; j<NbOfBJetsBins_; j++) Nttlike  += N_[wp][njets][j][iDatasetsTTLike_[i]];
	}
	return Nttlike;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNtt(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t Nttlike = 0;
	if(njets == -1){
		for(UInt_t i=0;i<iDatasetsTTLike_.size();i++){
			for(UInt_t j=0;j<NbOfJetsBins_;j++) Nttlike += N_[wp][j][nbjets][iDatasetsTTLike_[i]];
		}
		return Nttlike;
	}
	else{
		for(UInt_t i=0;i<iDatasetsTTLike_.size();i++){
			Nttlike  += N_[wp][njets][nbjets][iDatasetsTTLike_[i]];
		}
		return Nttlike;
	}
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNttErr(Int_t wp) const{
	Double_t NttErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) NttErr += pow(GetPredNttErr(wp,i),2);
	return sqrt(NttErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNttErr(Int_t wp, Int_t njets) const{
	Double_t Nttlike = 0;
	for(UInt_t i=0;i<iDatasetsTTLike_.size();i++){
		for(UInt_t j=0; j<NbOfBJetsBins_; j++) Nttlike  += pow(N_err_[wp][njets][j][iDatasetsTTLike_[i]],2);
	}
	return sqrt(Nttlike);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNttErr(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t Nttlike = 0;
	if(njets == -1){
		for(UInt_t i=0;i<iDatasetsTTLike_.size();i++){
			for(UInt_t j=0;j<NbOfJetsBins_;j++) Nttlike += pow(N_err_[wp][j][nbjets][iDatasetsTTLike_[i]],2);
		}
		return sqrt(Nttlike);
	}
	else{
		for(UInt_t i=0;i<iDatasetsTTLike_.size();i++){
			Nttlike  += pow(N_err_[wp][njets][nbjets][iDatasetsTTLike_[i]],2);
		}
		return sqrt(Nttlike);
	}
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNtotal(Int_t wp) const {
	Double_t NtotalMC = 0;
	if((UInt_t)wp<NbOfBtagWorkingPoint_)
	{
		for(UInt_t i=0;i<NbOfJetsBins_;i++){
      for(UInt_t j=0;j<NbOfBJetsBins_;j++) NtotalMC += Nbjets_[wp][i][j];
		}
	}
	return NtotalMC;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNtotal(Int_t wp, Int_t njets) const {
	Double_t NtotalMC = 0;
	if((UInt_t)wp<NbOfBtagWorkingPoint_)
		for(UInt_t i=0;i<NbOfBJetsBins_;i++) NtotalMC += Nbjets_[wp][njets][i];
	return NtotalMC;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredNtotal(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t NtotalMC = 0;
	if((UInt_t)wp<NbOfBtagWorkingPoint_)
		NtotalMC = Nbjets_[wp][njets][nbjets];
	return NtotalMC;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredN(Int_t idx) const {
	Double_t N_MC = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++){
    for(UInt_t j=0;j<NbOfBJetsBins_;j++) {
      N_MC += N_[0][i][j][idx];
    }
	}
	return N_MC;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetPredN(Int_t idx, Int_t njets) const {
	Double_t N_MC = 0;
	for(UInt_t i=0;i<NbOfBJetsBins_;i++) N_MC += N_[0][njets][i][idx];
	return N_MC;
}

  //////////////////////////////////////
  // Access to the parameter estimated values
  //////////////////////////////////////

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEb(Int_t wp) const{
	Double_t Effb = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) Effb += GetEstEb(wp,i)*GetPredNtotal(wp,i);
	Effb /= GetPredNtotal(wp);
	return Effb;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEb(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eb_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEbErr(Int_t wp) const{
	Double_t EffbErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) EffbErr += pow(GetEstEbErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffbErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEbErr(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eb_err_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEbErrHigh(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eb_err_up_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEbErrLow(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eb_err_down_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEudsc(Int_t wp) const{
	Double_t Effudsc = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) Effudsc += GetEstEudsc(wp,i)*GetPredNtotal(wp,i);
	Effudsc /= GetPredNtotal(wp);
	return Effudsc;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEudsc(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eudsc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEudscErr(Int_t wp) const{
	Double_t EffudscErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) EffudscErr += pow(GetEstEudscErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffudscErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEudscErr(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eudsc_err_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEudscErrHigh(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eudsc_err_up_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEudscErrLow(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? eudsc_err_down_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEuds(Int_t wp) const{
	Double_t Effuds = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) Effuds += GetEstEuds(wp,i)*GetPredNtotal(wp,i);
	Effuds /= GetPredNtotal(wp);
	return Effuds;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEuds(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? euds_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEudsErr(Int_t wp) const{
	Double_t EffudsErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) EffudsErr += pow(GetEstEudsErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffudsErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEudsErr(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? euds_err_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEudsErrHigh(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? euds_err_up_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstEudsErrLow(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ && (UInt_t)wp<NbOfBtagWorkingPoint_ ? euds_err_down_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNv(Int_t wp) const{
	Double_t Nv = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) Nv += GetEstNv(wp,i);
	return Nv;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNv(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ ? Nv_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNv(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t Nvlike = -9999;
	if((UInt_t)njets>=NbOfJetsBins_ || (UInt_t)wp>=NbOfBtagWorkingPoint_ ) return Nvlike;
	else if(njets == -1){
		Nvlike = 0;
		for(UInt_t i=0;i<NbOfJetsBins_;i++) Nvlike += Nv_bjets( Nv_[i], euds_[wp][i], (Int_t)Njets_[i], nbjets);
		return Nvlike;
	}
	else{
		Nvlike = Nv_bjets( Nv_[njets], euds_[wp][njets], (Int_t)Njets_[njets], nbjets);
		return Nvlike;
	}
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNvErr(Int_t wp) const{
	Double_t NvErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) NvErr += pow(GetEstNvErr(wp,i),2);
	return sqrt(NvErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNvErr(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ ? Nv_err_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNvErrUp(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ ? Nv_err_up_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNvErrDown(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ ? Nv_err_down_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNvErr(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t NvlikeErr = -9999;
	if(njets>=(Int_t)NbOfJetsBins_ || (UInt_t)wp>=NbOfBtagWorkingPoint_ ) return NvlikeErr;
	else if(njets == -1){
		NvlikeErr = 0;
		for(UInt_t i=0;i<NbOfJetsBins_;i++) NvlikeErr += pow(Nv_err_bjets( Nv_[i], Nv_err_[i], euds_[wp][i], euds_err_[wp][i], (Int_t)Njets_[i], nbjets),2);
		return sqrt(NvlikeErr);
	}
	else{
		NvlikeErr = Nv_err_bjets( Nv_[njets], Nv_err_[njets], euds_[wp][njets], euds_err_[wp][njets], (Int_t)Njets_[njets], nbjets);
		return NvlikeErr;
	}
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNvb(Int_t wp) const{
	Double_t Nvb = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) Nvb += GetEstNvb(wp,i);
	return Nvb;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNvb(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ ? Nvb_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNvb(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t Nvblike = -9999;
	if(njets>=(Int_t)NbOfJetsBins_ || (UInt_t)wp>=NbOfBtagWorkingPoint_ ) return Nvblike;
	else if(njets == -1){
		Nvblike = 0;
		for(UInt_t i=0;i<NbOfJetsBins_;i++) Nvblike += Nvb_bjets( Nvb_[i], eb_[wp][njets], euds_[wp][i], (Int_t)Njets_[i], nbjets);
		return Nvblike;
	}
	else{
		Nvblike = Nvb_bjets( Nvb_[njets], eb_[wp][njets], euds_[wp][njets], (Int_t)Njets_[njets], nbjets);
		return Nvblike;
	}
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNvbErr(Int_t wp) const{
	Double_t NvbErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) NvbErr += pow(GetEstNvbErr(wp,i),2);
	return sqrt(NvbErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNvbErr(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ ? Nvb_err_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNvbErr(Int_t wp, Int_t njets, Int_t nbjets) const{
	return -9999;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNtt(Int_t wp) const{
	Double_t Ntt = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) Ntt += GetEstNtt(wp,i);
	return Ntt;
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNtt(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ ? Ntt_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNtt(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t Nttlike = -9999;
	if(njets>=(Int_t)NbOfJetsBins_|| (UInt_t)wp>=NbOfBtagWorkingPoint_ ) return Nttlike;
	else if(njets == -1){
		Nttlike = 0;
		for(UInt_t i=0;i<NbOfJetsBins_;i++) Nttlike += Ntt_bjets( Ntt_[i], eb_[wp][i], eudsc_[wp][i], (Int_t)Njets_[i], nbjets);
		return Nttlike;
	}
	else{
		Nttlike = Ntt_bjets( Ntt_[njets], eb_[wp][njets], eudsc_[wp][njets], (Int_t)Njets_[njets], nbjets);
		return Nttlike;
	}
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNttErr(Int_t wp) const{
	Double_t NttErr = 0;
	for(UInt_t i=0;i<NbOfJetsBins_;i++) NttErr += pow(GetEstNttErr(wp,i),2);
	return sqrt(NttErr);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNttErr(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ ? Ntt_err_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNttErrUp(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ ? Ntt_err_up_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNttErrDown(Int_t wp, Int_t njets) const{
	return ((UInt_t)njets<NbOfJetsBins_ ? Ntt_err_down_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNttErr(Int_t wp, Int_t njets, Int_t nbjets) const{
	Double_t NttlikeErr = -9999;
	if(njets>=(Int_t)NbOfJetsBins_|| (UInt_t)wp>=NbOfBtagWorkingPoint_ ) return NttlikeErr;
	else if(njets == -1){
		NttlikeErr = 0;
		for(UInt_t i=0;i<NbOfJetsBins_;i++) NttlikeErr += pow(Ntt_err_bjets( Ntt_[i], Ntt_err_[i], eb_[wp][i], eb_err_[wp][i], eudsc_[wp][i], eudsc_err_[wp][i], (Int_t)Njets_[i], nbjets),2);
		return sqrt(NttlikeErr);
	}
	else{
		NttlikeErr = Ntt_err_bjets( Ntt_[njets], Ntt_err_[njets], eb_[wp][njets], eb_err_[wp][njets], eudsc_[wp][njets], eudsc_err_[wp][njets], (Int_t)Njets_[njets], nbjets);
		return NttlikeErr;
	}
}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNtotal(Int_t wp)                        const{ return GetEstNtt(wp)              + GetEstNv(wp)              + GetEstNvb(wp);}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNtotal(Int_t wp, Int_t njets)             const{ return GetEstNtt(wp,njets)        + GetEstNv(wp,njets)        + GetEstNvb(wp,njets);}

/**________________________________________________________________________________________________________________*/
Double_t VJetEstimation::GetEstNtotal(Int_t wp, Int_t njets, Int_t nbjets) const{ return GetEstNtt(wp,njets,nbjets) + GetEstNv(wp,njets,nbjets) + GetEstNvb(wp,njets,nbjets);}

void VJetEstimation::FillSummaryHistos(){
  
	TString NbjetsXLabel[4] = {"0","1","2","#geq 3"};
	Double_t N = 1;
	Double_t Nerr = 0;
  
  char name[100];char title[100];
  
  
    //init summary histos
  hsNbjets_MC  = vector< vector< THStack* > >(NbOfBtagWorkingPoint_,vector< THStack* >(NbOfJetsBins_+1, NULL));
  hsNbjets_Est = vector< vector< THStack* > >(NbOfBtagWorkingPoint_,vector< THStack* >(NbOfJetsBins_+1, NULL));
    //	hsNjets_MC   = vector< vector< THStack* > >(NbOfBtagWorkingPoint_,vector< THStack* >(NbOfBJetsBins_+1));
    //	hsNjets_Est  = vector< vector< THStack* > >(NbOfBtagWorkingPoint_,vector< THStack* >(NbOfBJetsBins_+1));
	MyLeg = new TLegend(0.6,0.6,0.9,0.9);
  
  
  
  tCanva_Nbjets_Summary = vector< vector< TCanvas* > >(NbOfBtagWorkingPoint_,vector< TCanvas* >(NbOfJetsBins_+1, NULL));
  tCanva_Njets_Summary  = vector< vector< TCanvas* > >(NbOfBtagWorkingPoint_,vector< TCanvas* >(NbOfBJetsBins_+1, NULL));  
  
	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {
      // MC prediction/EStimation summary (inclusive)
    sprintf(name ,"hNbjetsSummary_Inclusive_wp_%d",i);
    sprintf(title,"MC prediction/Estimation summary (inclusive, WP nr%d)",i);
    tCanva_Nbjets_Summary[i][NbOfJetsBins_] = new TCanvas(name,title,600,600);
    gROOT->GetListOfCanvases()->RemoveLast();  
    
      // MC prediction/Estimation summary (b-inclusive)
    sprintf(name ,"hNjetsSummary_bInclusive_wp_%d",i);
    sprintf(title,"MC prediction/Estimation summary (b-inclusive, WP nr%d)",i);
    tCanva_Njets_Summary[i][NbOfBJetsBins_] = new TCanvas(name,title,600,600);
    gROOT->GetListOfCanvases()->RemoveLast();
    
		for(UInt_t j=0;j<NbOfJetsBins_;j++){
			for(UInt_t k=0;k<NbOfBJetsBins_;k++){
          // For V-like and TT-like estimation
				hNbjetsEstSummary[i][j][0].SetBinContent(k+1,GetEstNv( i,j,k));
				hNbjetsEstSummary[i][j][0].SetBinError(k+1,GetEstNvErr( i,j,k));
          //hNbjetsEstSummary[i][j][1]->SetBinContent(k+1,GetEstNvb(i,j,k));
				hNbjetsEstSummary[i][j][2].SetBinContent(k+1,GetEstNtt(i,j,k));
				hNbjetsEstSummary[i][j][2].SetBinError(k+1,GetEstNttErr(i,j,k));
				
				hNbjetsEstSummary[i][j][0].GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);
          //hNbjetsEstSummary[i][j][1]->GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);
				hNbjetsEstSummary[i][j][2].GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);
        
          // For V-like and TT-like MC prediction
				hNbjetsMCSummary[i][j][0].SetBinContent(k+1,GetPredNv( i,j,k));
				hNbjetsMCSummary[i][j][0].SetBinError(k+1,GetPredNvErr( i,j,k));
          //hNbjetsMCSummary[i][j][1]->SetBinContent(k+1,GetPredNvb(i,j,k));
				hNbjetsMCSummary[i][j][2].SetBinContent(k+1,GetPredNtt(i,j,k));
				hNbjetsMCSummary[i][j][2].SetBinError(k+1,GetPredNttErr(i,j,k));
				
				hNbjetsMCSummary[i][j][0].GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);
          //hNbjetsMCSummary[i][j][1]->GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);
				hNbjetsMCSummary[i][j][2].GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);
				
				hNbjets_mc_[i][j][0].SetBinContent(k+1,GetPredNv( i,j,k));
				hNbjets_mc_[i][j][0].SetBinError(k+1,GetPredNvErr( i,j,k));
          //hNbjets_mc_[i][j][1]->SetBinContent(k+1,GetPredNvb(i,j,k));
          //hNbjets_mc_[i][j][1]->SetBinError(k+1,GetPredNvbErr( i,j,k));
				hNbjets_mc_[i][j][2].SetBinContent(k+1,GetPredNtt(i,j,k));
				hNbjets_mc_[i][j][2].SetBinError(k+1,GetPredNttErr( i,j,k));
				
				hNbjets_pdf_mc_[i][j][0].SetBinContent(k+1,Nv_bjets( N, euds_mc_[i][j], Njets_[j], k));
				hNbjets_pdf_mc_[i][j][0].SetBinError(k+1,Nv_err_bjets( N, Nerr, euds_mc_[i][j], euds_err_mc_[i][j], Njets_[j], k));
          //hNbjets_pdf_mc_[i][j][1]->SetBinContent(k+1,Nvb_bjets( N, eb_mc_[i][j], euds_mc_[i][j],  Njets_[j], k));
          //hNbjets_pdf_mc_[i][j][1]->SetBinError(k+1,Nvb_err_bjets(Nvb_[j], eb_mc_[i][j], euds_mc_[i][j],  Njets_[j], k));
				hNbjets_pdf_mc_[i][j][2].SetBinContent(k+1,Ntt_bjets(N, eb_mc_[i][j], eudsc_mc_[i][j], Njets_[j], k));
				hNbjets_pdf_mc_[i][j][2].SetBinError(k+1,Ntt_err_bjets( N, Nerr, eb_mc_[i][j], eb_err_mc_[i][j], eudsc_mc_[i][j], eudsc_err_mc_[i][j], Njets_[j], k));
				
				hNbjets_pdf_est_[i][j][0].SetBinContent(k+1,GetEstNv( i,j,k));
				hNbjets_pdf_est_[i][j][0].SetBinError(k+1,GetEstNvErr( i,j,k));
          //hNbjets_pdf_est_[i][j][1]->SetBinContent(k+1,GetEstNvb(i,j,k));
          //hNbjets_pdf_est_[i][j][1]->SetBinError(k+1,GetEstNvbErr(i,j,k));
				hNbjets_pdf_est_[i][j][2].SetBinContent(k+1,GetEstNtt(i,j,k));
				hNbjets_pdf_est_[i][j][2].SetBinError(k+1,GetEstNttErr(i,j,k));
			}
      
			hNbjets_mc_[i][j][0].Scale((hNbjets_mc_[i][j][0].Integral() != 0 ? 1/hNbjets_mc_[i][j][0].Integral() : 1));
			hNbjets_mc_[i][j][1].Scale((hNbjets_mc_[i][j][1].Integral() != 0 ? 1/hNbjets_mc_[i][j][1].Integral() : 1));
			hNbjets_mc_[i][j][2].Scale((hNbjets_mc_[i][j][2].Integral() != 0 ? 1/hNbjets_mc_[i][j][2].Integral() : 1));
			hNbjets_mc_[i][j][0].GetXaxis()->CenterLabels();hNbjets_mc_[i][j][0].GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_mc_[i][j][1].GetXaxis()->CenterLabels();hNbjets_mc_[i][j][1].GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_mc_[i][j][2].GetXaxis()->CenterLabels();hNbjets_mc_[i][j][2].GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_mc_[i][j][0].Scale((hNbjets_pdf_mc_[i][j][0].Integral() != 0 ? 1/hNbjets_pdf_mc_[i][j][0].Integral() : 1));
			hNbjets_pdf_mc_[i][j][1].Scale((hNbjets_pdf_mc_[i][j][1].Integral() != 0 ? 1/hNbjets_pdf_mc_[i][j][1].Integral() : 1));
			hNbjets_pdf_mc_[i][j][2].Scale((hNbjets_pdf_mc_[i][j][2].Integral() != 0 ? 1/hNbjets_pdf_mc_[i][j][2].Integral() : 1));
			hNbjets_pdf_mc_[i][j][0].GetXaxis()->CenterLabels();hNbjets_pdf_mc_[i][j][0].GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_mc_[i][j][1].GetXaxis()->CenterLabels();hNbjets_pdf_mc_[i][j][1].GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_mc_[i][j][2].GetXaxis()->CenterLabels();hNbjets_pdf_mc_[i][j][2].GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_est_[i][j][0].Scale((hNbjets_pdf_est_[i][j][0].Integral() != 0 ? 1/hNbjets_pdf_est_[i][j][0].Integral() : 1));
			hNbjets_pdf_est_[i][j][1].Scale((hNbjets_pdf_est_[i][j][1].Integral() != 0 ? 1/hNbjets_pdf_est_[i][j][1].Integral() : 1));
			hNbjets_pdf_est_[i][j][2].Scale((hNbjets_pdf_est_[i][j][2].Integral() != 0 ? 1/hNbjets_pdf_est_[i][j][2].Integral() : 1));
			hNbjets_pdf_est_[i][j][0].GetXaxis()->CenterLabels();hNbjets_pdf_est_[i][j][0].GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_est_[i][j][1].GetXaxis()->CenterLabels();hNbjets_pdf_est_[i][j][1].GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_est_[i][j][2].GetXaxis()->CenterLabels();hNbjets_pdf_est_[i][j][2].GetXaxis()->SetTitle("Nb. of b-tagged jets");
      
		}
		for(UInt_t j=0;j<NbOfBJetsBins_;j++){
        // For V-like and TT-like estimation
			hNbjetsEstSummary[i][NbOfJetsBins_][0].SetBinContent(j+1,GetEstNv( i,-1,j));
        //hNbjetsEstSummary[i][NbOfJetsBins_][1]->SetBinContent(j+1,GetEstNvb(i,-1,j));
			hNbjetsEstSummary[i][NbOfJetsBins_][2].SetBinContent(j+1,GetEstNtt(i,-1,j));
			hNbjetsEstSummary[i][NbOfJetsBins_][0].GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[j]);
			hNbjetsEstSummary[i][NbOfJetsBins_][1].GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[j]);
			hNbjetsEstSummary[i][NbOfJetsBins_][2].GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[j]);
        // For V-like and TT-like MC prediction
			hNbjetsMCSummary[i][NbOfJetsBins_][0].SetBinContent(j+1,GetPredNv( i,-1,j));
        //hNbjetsMCSummary[i][NbOfJetsBins_][1]->SetBinContent(j+1,GetPredNvb(i,-1,j));
			hNbjetsMCSummary[i][NbOfJetsBins_][2].SetBinContent(j+1,GetPredNtt(i,-1,j));
			hNbjetsMCSummary[i][NbOfJetsBins_][0].GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[i]);
        //hNbjetsMCSummary[i][NbOfJetsBins_][1]->GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[i]);
			hNbjetsMCSummary[i][NbOfJetsBins_][2].GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[i]);
		}
    /*******/
    /* Fill of THStack content */
    for (UInt_t j=0;j<=NbOfJetsBins_;j++) {
      MyLeg->Clear();
      
      if (j!=NbOfJetsBins_) {
          // MC prediction vs estimation
        sprintf(name,"hNbjetsSummary_%d_jets_wp_%d",Njets_[j],i);
        sprintf(title,"MC prediction/Estimation summary for %d jets (WP nr%d)",Njets_[j],i);
        tCanva_Nbjets_Summary[i][j] = new TCanvas(name,title,600,600);
        gROOT->GetListOfCanvases()->RemoveLast();
        
          // MC prediction vs estimation
        sprintf(name,"hNjetsSummary_%d_bjets_wp_%d",j,i);
        sprintf(title,"MC prediction/Estimation summary for %d b-jets (WP nr%d)",j,i);		
        tCanva_Njets_Summary[i][j] = new TCanvas(name,title,600,600);
        gROOT->GetListOfCanvases()->RemoveLast();
      }
      
      tCanva_Nbjets_Summary[i][j]->cd();
      
      if(j!=NbOfJetsBins_) sprintf(name,"hNbjetsMCStackSummary_%d_jets_wp_%d",Njets_[j],i);
      else                 sprintf(name,"hNbjetsMCStackSummary_Inclusive_wp_%d",i);
      if(j!=NbOfJetsBins_) sprintf(title,"V-like and tt-like MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
      else                 sprintf(title,"V-like and tt-like MC prediction summary (inclusive, WP nr%d)",i);
      
      hsNbjets_MC[i][j] = new THStack(name,title);
      
      hNbjetsMCSummary[i][j][0].SetFillStyle(3004);
      hNbjetsMCSummary[i][j][0].SetFillColor(kGreen+1);
      hNbjetsMCSummary[i][j][1].SetFillStyle(3005);
      hNbjetsMCSummary[i][j][1].SetFillColor(kBlue+1);
      hNbjetsMCSummary[i][j][2].SetFillStyle(3006);
      hNbjetsMCSummary[i][j][2].SetFillColor(kRed+1);
      
      hsNbjets_MC[i][j]->Add(& hNbjetsMCSummary[i][j][0]);
        //hsNbjets_MC[i][j]->Add(hNbjetsMCSummary[i][j][1]);
      hsNbjets_MC[i][j]->Add(& hNbjetsMCSummary[i][j][2]);
      
      MyLeg->AddEntry(& hNbjetsMCSummary[i][j][0],"V-like events","f");
        //MyLeg->AddEntry(hNbjetsMCSummary[i][j][1],"Vb-like events","f");
      MyLeg->AddEntry(& hNbjetsMCSummary[i][j][2],"#bar{t}t-like events","f");
      
      hsNbjets_MC[i][j]->Draw();
      hsNbjets_MC[i][j]->GetXaxis()->SetTitle("Nb of b-jets");
      
      if(i!=NbOfJetsBins_) sprintf(name,"hNbjetsEstStackSummary_%d_jets_wp_%d",Njets_[j],i);
      else                 sprintf(name,"hNbjetsEstStackSummary_Inclusive_wp_%d",i);
      if(i!=NbOfJetsBins_) sprintf(title,"V-like and tt-like estimation summary for %d jets (WP nr%d)",Njets_[j],i);
      else                 sprintf(title,"V-like and tt-like estimation summary (inclusive, WP nr%d)",i);
      
      hsNbjets_Est[i][j] = new THStack(name,title);
      
      hNbjetsEstSummary[i][j][0].SetMarkerColor(kYellow+1);
        //hNbjetsEstSummary[i][j][1]->SetMarkerColor(kMagenta+1);
      hNbjetsEstSummary[i][j][2].SetMarkerColor(kBlue+1);
      hNbjetsEstSummary[i][j][0].SetMarkerStyle(22);
        //hNbjetsEstSummary[i][j][1]->SetMarkerStyle(23);
      hNbjetsEstSummary[i][j][2].SetMarkerStyle(24);
      hNbjetsEstSummary[i][j][0].SetMarkerSize(1.5);
        //hNbjetsEstSummary[i][j][1]->SetMarkerSize(1.5);
      hNbjetsEstSummary[i][j][2].SetMarkerSize(1.5);
      
      hsNbjets_Est[i][j]->Add(& hNbjetsEstSummary[i][j][0]);
        //hsNbjets_Est[i][j]->Add(hNbjetsEstSummary[i][j][1]);
      hsNbjets_Est[i][j]->Add(& hNbjetsEstSummary[i][j][2]);
      
      MyLeg->AddEntry(& hNbjetsEstSummary[i][j][0],"V-like estimation","p");
        //MyLeg->AddEntry(hNbjetsEstSummary[i][j][1],"Vb-like estimation","p");
      MyLeg->AddEntry(& hNbjetsEstSummary[i][j][2],"TT-like estimation","p");
      
      hsNbjets_Est[i][j]->Draw("PEsame");
      hsNbjets_Est[i][j]->GetXaxis()->SetTitle("Nb of b-jets");
      MyLeg->Draw("same");
      
      tCanva_Nbjets_Summary[i][j]->Update();
        //    tCanva_Nbjets_Summary[i][j]->Write();
    }
    /*******/
  }
  for (std::vector<std::vector<TCanvas*> >::iterator iter_1=tCanva_Nbjets_Summary.begin(); iter_1!=tCanva_Nbjets_Summary.end(); iter_1++) {
    for (std::vector<TCanvas*>::iterator iter_2=iter_1->begin(); iter_2!=iter_1->end(); iter_2++) {
      if(true)cout<<"Writing summary histograms for X jets, wp nr Y"<<endl;
      sprintf(name, "macros/%s.pdf", (*iter_2)->GetName());
      (*iter_2)->Print(name);
    }
  }
}
void VJetEstimation::BckgdSubstraction(vector<MCObsExpectation*> &hists, vector<string> &name, Float_t lumi){
	if((UInt_t)NbOfBtagWorkingPoint_*NbOfJetsBins_ != hists.size()){
		cout<<"Mis-match between NbOfBtagWorkingPoint_*NbOfJetsBins_ and hists.size()"<<endl;
		return;
	}
	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {
		for(UInt_t j=0;j<NbOfJetsBins_;j++){
			for(UInt_t k=0;k<name.size();k++){
				if(!hists[i*NbOfJetsBins_+j]->GetHistogram(name[k])) continue;
				for(UInt_t l=0;l<NbOfBJetsBins_;l++){
					Nbjets_[i][j][l]  += - (hists[i*NbOfJetsBins_+j]->GetHistogram(name[k])->GetBinContent(l+1))*lumi;
				}
			}
		}
	}
}

void VJetEstimation::PrintInputs(UInt_t njets){
	cout<<"For N = "<<Njets_[njets]<<" jets _________________________________________________"<<endl;
	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {
		cout<<"B-tagging working point : "; cout<<this->BtagWorkingPoint_[i]<<endl;
		cout<<"Input values : Nbjets (0/1/2/3 b-jets) = "; cout<<GetPredNtotal(i,njets)<<" ("<<Nbjets_[i][njets][0]<<"/"<<Nbjets_[i][njets][1]<<"/"<<Nbjets_[i][njets][2]<<"/"<<Nbjets_[i][njets][3]<<")"<<endl;
		cout<<"Predicted eudsc : "<<eudsc_mc_[i][njets]<<" +/- "<<eudsc_err_mc_[i][njets]<<endl;
		cout<<"Predicted euds  : "<<euds_mc_[i][njets]<<" +/- "<<euds_err_mc_[i][njets]<<endl;
		cout<<"Predicted eb    : "<<eb_mc_[i][njets]<<" +/- "<<eb_err_mc_[i][njets]<<endl;
		cout<<"For tt-like datasets : "<<endl;
		for(UInt_t k=0;k<iDatasetsTTLike_.size();k++){
			if(N_[i][njets][0][iDatasetsTTLike_[k]]==0 && N_[i][njets][0][iDatasetsTTLike_[k]]==N_[i][njets][1][iDatasetsTTLike_[k]] && N_[i][njets][0][iDatasetsTTLike_[k]]==N_[i][njets][2][iDatasetsTTLike_[k]] && N_[i][njets][0][iDatasetsTTLike_[k]]==N_[i][njets][3][iDatasetsTTLike_[k]]) continue;
      /*			cout<<"Dataset "<<k+1<<" : N (0/1/2/3 jets) = "<<GetPredN(iDatasetsTTLike_[k],njets)<<" ("<<N_[i][njets][0][iDatasetsTTLike_[k]]<<"/"<<N_[i][njets][1][iDatasetsTTLike_[k]]<<"/"<<N_[i][njets][2][iDatasetsTTLike_[k]]<<"/"<<N_[i][njets][3][iDatasetsTTLike_[k]]<<")"<<endl; */
      printf("Dataset %s : N (0/1/2/3 jets) = %5lg (%5lg/%5lg/%5lg/%5lg)\n", vDatasets_[k].Name().c_str(), GetPredN(iDatasetsTTLike_[k],njets), N_[i][njets][0][iDatasetsTTLike_[k]], N_[i][njets][1][iDatasetsTTLike_[k]], N_[i][njets][2][iDatasetsTTLike_[k]], N_[i][njets][3][iDatasetsTTLike_[k]]);
		}
		cout<<"Total number : N (0/1/2/3 jets) = "<<GetPredNtt(i,njets)<<" ("<<GetPredNtt(i,njets,0)<<"/"<<GetPredNtt(i,njets,1)<<"/"<<GetPredNtt(i,njets,2)<<"/"<<GetPredNtt(i,njets,3)<<")"<<endl;
		cout<<"For V-like datasets : "<<endl;
		for(UInt_t k=0;k<iDatasetsVLike_.size();k++){
			if(N_[i][njets][0][iDatasetsVLike_[k]]==0 && N_[i][njets][0][iDatasetsVLike_[k]]==N_[i][njets][1][iDatasetsVLike_[k]] && N_[i][njets][0][iDatasetsVLike_[k]]==N_[i][njets][2][iDatasetsVLike_[k]] && N_[i][njets][0][iDatasetsVLike_[k]]==N_[i][njets][3][iDatasetsVLike_[k]]) continue;
			printf("Dataset %s : N (0/1/2/3 jets) = %lg (%lg/%lg/%lg/%lg)\n", vDatasets_[k].Name().c_str(), GetPredN(iDatasetsVLike_[k],njets), N_[i][njets][0][iDatasetsVLike_[k]], N_[i][njets][1][iDatasetsVLike_[k]], N_[i][njets][2][iDatasetsVLike_[k]], N_[i][njets][3][iDatasetsVLike_[k]]);
		}
		cout<<"Total number : N (0/1/2/3 jets) = "<<GetPredNv(i,njets)<<" ("<<GetPredNv(i,njets,0)<<"/"<<GetPredNv(i,njets,1)<<"/"<<GetPredNv(i,njets,2)<<"/"<<GetPredNv(i,njets,3)<<")"<<endl;
	}
}
void VJetEstimation::PrintResults(UInt_t njets){
	cout<<"***********************************************************"<<endl;
	cout<<"For N = "<<Njets_[njets]<<" jets"<<endl;
	cout<<std::fixed<<setprecision(4);
	cout<<"e0bq/e1bq/e2bq = : "<<e0bq_[njets]<<"/"<<e1bq_[njets]<<"/"<<e2bq_[njets]<<"/"<<endl;
	cout<<"Btag working point : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<BtagWorkingPoint_[i]<<" || ";} cout<<endl;
	cout<<"Estimated/Predicted eb	 : ";	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstEb(i,njets)	<<" +/- "<<GetEstEbErr(i,njets)		<<" // "<<GetPredEb(i,njets)   <<" +/- "<<GetPredEbErr(i,njets)<<" || ";}    cout<<endl;
	cout<<"Estimated/Predicted eudsc : ";	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstEudsc(i,njets)<<" +/- "<<GetEstEudscErr(i,njets)	<<" // "<<GetPredEudsc(i,njets)<<" +/- "<<GetPredEudscErr(i,njets)<<" || ";} cout<<endl;
	cout<<"Estimated/Predicted euds  : ";	for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstEuds(i,njets) <<" +/- "<<GetEstEudsErr(i,njets)	<<" // "<<GetPredEuds(i,njets) <<" +/- "<<GetPredEudsErr(i,njets)<<" || ";}  cout<<endl;
	cout<<std::fixed<<setprecision(5);
	cout<<"Estimated Ntt-like : ";cout<<GetEstNtt(0,njets)<<" +/- "<<GetEstNttErr(0,njets)<<" / Predicted value from MC : "<<GetPredNtt(0,njets)<<"	+/- "<<GetPredNttErr(0,njets)<<endl;
	cout<<"- 0 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNtt(i, njets, 0)<<" +/- "<<GetEstNttErr(i,njets,0)<<" // "<< GetPredNtt(i, njets, 0)<<" +/- "<<GetPredNttErr(i,njets,0)<<" || ";} cout<<endl;
	cout<<"- 1 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNtt(i, njets, 1)<<" +/- "<<GetEstNttErr(i,njets,1)<<" // "<< GetPredNtt(i, njets, 1)<<" +/- "<<GetPredNttErr(i,njets,1)<<" || ";} cout<<endl;
	cout<<"- 2 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNtt(i, njets, 2)<<" +/- "<<GetEstNttErr(i,njets,2)<<" // "<< GetPredNtt(i, njets, 2)<<" +/- "<<GetPredNttErr(i,njets,2)<<" || ";} cout<<endl;
	cout<<"- 3 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNtt(i, njets, 3)<<" +/- "<<GetEstNttErr(i,njets,3)<<" // "<< GetPredNtt(i, njets, 3)<<" +/- "<<GetPredNttErr(i,njets,3)<<" || ";} cout<<endl;
    // 	cout<<"Estimated Nvb-like : ";cout<<Nvb_[njets]<<" / Predicted value from MC : "<<GetPredNvb(0,njets)<<endl;
    // 	cout<<"- 0 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNvb(i, njets, 0)<<" // "<< GetPredNvb(i, njets, 0) <<" || ";} cout<<endl;
    // 	cout<<"- 1 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNvb(i, njets, 1)<<" // "<< GetPredNvb(i, njets, 1) <<" || ";} cout<<endl;
    // 	cout<<"- 2 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNvb(i, njets, 2)<<" // "<< GetPredNvb(i, njets, 2) <<" || ";} cout<<endl;
    // 	cout<<"- 3 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNvb(i, njets, 3)<<" // "<< GetPredNvb(i, njets, 3) <<" || ";} cout<<endl;
	cout<<"Estimated Nv-like  : ";cout<<GetEstNv(0,njets)<<" +/- "<<GetEstNvErr(0,njets)<<" / Predicted value from MC : "<<GetPredNv(0,njets)<<"	+/- "<<GetPredNvErr(0,njets)<<endl;
	cout<<"- 0 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNv(i, njets, 0)<<" +/- "<<GetEstNvErr(i,njets,0)<<" // "<< GetPredNv(i, njets, 0)<<" +/- "<<GetPredNvErr(i,njets,0)<<" || ";} cout<<endl;
	cout<<"- 1 b-jet          : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNv(i, njets, 1)<<" +/- "<<GetEstNvErr(i,njets,1)<<" // "<< GetPredNv(i, njets, 1)<<" +/- "<<GetPredNvErr(i,njets,1)<<" || ";} cout<<endl;
	cout<<"- 2 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNv(i, njets, 2)<<" +/- "<<GetEstNvErr(i,njets,2)<<" // "<< GetPredNv(i, njets, 2)<<" +/- "<<GetPredNvErr(i,njets,2)<<" || ";} cout<<endl;
	cout<<"- 3 b-jets         : ";for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) {cout<<GetEstNv(i, njets, 3)<<" +/- "<<GetEstNvErr(i,njets,3)<<" // "<< GetPredNv(i, njets, 3)<<" +/- "<<GetPredNvErr(i,njets,3)<<" || ";} cout<<endl;
	cout<<"***********************************************************"<<endl;
}

  ///////////////////////////////////////////////////////////////////	
  // Methods to calculates the estimators ...  (from equations)
  ///////////////////////////////////////////////////////////////////
Double_t VJetEstimation::Nbjets(Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t njets, Int_t nbjets) const{
	Double_t Ntt_bjets = this->Ntt_bjets(Ntt,eb,eudsc,njets, nbjets);
	Double_t Nv_bjets  = this->Nv_bjets( Nv, euds, njets, nbjets);
	Double_t Nvb_bjets = this->Nvb_bjets(Nvb,eb,euds, njets, nbjets);
	return (Ntt_bjets+Nv_bjets+Nvb_bjets);
}

Double_t VJetEstimation::Ntt_bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t njets, Int_t nbjets) const{
	if(nbjets == 0)      return Ntt_0bjet( Ntt, eb, eudsc, njets);
	else if(nbjets == 1) return Ntt_1bjet( Ntt, eb, eudsc, njets);
	else if(nbjets == 2) return Ntt_2bjets(Ntt, eb, eudsc, njets);
	else if(nbjets == 3) return Ntt_3bjets(Ntt, eb, eudsc, njets);
	else                 return -9999;
}

Double_t VJetEstimation::Ntt_err_bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb, Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t njets, Int_t nbjets) const{
	if(nbjets == 0)      return Ntt_err_0bjet( Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 1) return Ntt_err_1bjet( Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 2) return Ntt_err_2bjets(Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 3) return Ntt_err_3bjets(Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else                 return -9999;
}

Double_t VJetEstimation::Nvb_bjets(Double_t Nvb, Double_t eb, Double_t euds, Int_t njets, Int_t nbjets) const{
	if(nbjets == 0)      return Nvb_0bjet( Nvb, eb, euds, njets);
	else if(nbjets == 1) return Nvb_1bjet( Nvb, eb, euds, njets);
	else if(nbjets == 2) return Nvb_2bjets(Nvb, eb, euds, njets);
	else if(nbjets == 3) return Nvb_3bjets(Nvb, eb, euds, njets);
	else                 return -9999;
}

Double_t VJetEstimation::Nv_bjets(Double_t Nv, Double_t euds, Int_t njets, Int_t nbjets) const{
	if(nbjets == 0)      return Nv_0bjet( Nv, euds, njets);
	else if(nbjets == 1) return Nv_1bjet( Nv, euds, njets);
	else if(nbjets == 2) return Nv_2bjets(Nv, euds, njets);
	else if(nbjets == 3) return Nv_3bjets(Nv, euds, njets);
	else                 return -9999;
}

Double_t VJetEstimation::Nv_err_bjets(Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t njets, Int_t nbjets) const{
	if(nbjets == 0)      return Nv_err_0bjet( Nv, Nv_err, euds, euds_err, njets);
	else if(nbjets == 1) return Nv_err_1bjet( Nv, Nv_err, euds, euds_err, njets);
	else if(nbjets == 2) return Nv_err_2bjets(Nv, Nv_err, euds, euds_err, njets);
	else if(nbjets == 3) return Nv_err_3bjets(Nv, Nv_err, euds, euds_err, njets);
	else                 return -9999;
}

Double_t VJetEstimation::N0bjet(Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) const{
	Double_t Ntt_0bjet = this->Ntt_0bjet(Ntt,eb,eudsc,n);
	Double_t Nv_0bjet  = this->Nv_0bjet( Nv, euds, n);
	Double_t Nvb_0bjet = this->Nvb_0bjet(Nvb,eb,euds, n);
	return (Ntt_0bjet+Nv_0bjet+Nvb_0bjet);
}

Double_t VJetEstimation::Ntt_0bjet(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) const{
	Double_t  Ntt_0bjet = ((1-eb)*(1-eb)*pow((1-eudsc),n-2)*GetTTEff2bq(n-Njets_[0])
                         +         (1-eb)*pow((1-eudsc),n-1)*GetTTEff1bq(n-Njets_[0])
                         +                pow((1-eudsc),n)  *GetTTEff0bq(n-Njets_[0]))*Ntt;
	return (Ntt_0bjet < 0 ? 0 : Ntt_0bjet);
}
Double_t VJetEstimation::Ntt_err_0bjet(Double_t Ntt, Double_t Ntt_err, Double_t eb, Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) const{
	Double_t  Ntt_err_0bjet = pow((2*(eb-1)*pow((1-eudsc),n-2)*GetTTEff2bq(n-Njets_[0])
                                 +    (-1)*pow((1-eudsc),n-1)*GetTTEff1bq(n-Njets_[0]))*Ntt*eb_err,2)
  + pow((pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)*(1-eb)*(1-eb)*GetTTEff2bq(n-Njets_[0])
         +pow(-1.,n-1)*(n-1)*pow(eudsc-1,n-2)*       (1-eb)*GetTTEff1bq(n-Njets_[0])
         +pow(-1.,n  )*(n  )*pow(eudsc-1,n-1)*       (1   )*GetTTEff0bq(n-Njets_[0]))*Ntt*eudsc_err,2)
  + pow(Ntt_0bjet(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_0bjet < 0 ? 0 : sqrt(Ntt_err_0bjet));
}
Double_t VJetEstimation::Nvb_0bjet(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) const{
	Double_t  Nvb_0bjet = (1-eb)*pow((1-euds),n-1)*Nvb;
	return (Nvb_0bjet < 0 ? 0 : Nvb_0bjet);
}
Double_t VJetEstimation::Nv_0bjet (Double_t Nv, Double_t euds, Int_t n) const{
	Double_t  Nv_0bjet = pow((1-euds),n)*Nv;
	return (Nv_0bjet < 0 ? 0 : Nv_0bjet);
}
Double_t VJetEstimation::Nv_err_0bjet (Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) const{
	Double_t  Nv_err_0bjet = pow(pow(-1.,n)*n*pow(euds-1,n-1)*Nv*euds_err,2)
  + pow(pow((1-euds),n)*Nv_err,2);
	return (Nv_err_0bjet < 0 ? 0 : sqrt(Nv_err_0bjet));
}

Double_t VJetEstimation::N1bjet   (Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc, Double_t euds, Int_t n) const{
	Double_t Ntt_1bjet = this->Ntt_1bjet(Ntt,eb,eudsc,n);
	Double_t Nv_1bjet  = this->Nv_1bjet( Nv, euds, n);
	Double_t Nvb_1bjet = this->Nvb_1bjet(Nvb,eb,euds, n);
	return (Ntt_1bjet+Nv_1bjet+Nvb_1bjet);
}
Double_t VJetEstimation::Ntt_1bjet(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) const{
	Double_t  Ntt_1bjet =((2*eb*(1-eb)*pow(1-eudsc,n-2)+(1-eb)*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3))*GetTTEff2bq(n-Njets_[0])
                        +          (eb*pow(1-eudsc,n-1)+       (1-eb)*(n-1)*eudsc*pow(1-eudsc,n-2))*GetTTEff1bq(n-Njets_[0])
                        +                                                (n*eudsc*pow(1-eudsc,n-1))*GetTTEff0bq(n-Njets_[0]))*Ntt;
	return (Ntt_1bjet < 0 ? 0 : Ntt_1bjet);
}
Double_t VJetEstimation::Ntt_err_1bjet(Double_t Ntt, Double_t Ntt_err, Double_t eb,  Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) const{
	Double_t  Ntt_err_1bjet = pow(((2*(1-2*eb)*pow(1-eudsc,n-2)+2*(eb-1)*(n-2)*eudsc*pow(1-eudsc,n-3))*GetTTEff2bq(n-Njets_[0])
                                 +           (pow(1-eudsc,n-1)+    (-1)*(n-1)*eudsc*pow(1-eudsc,n-2))*GetTTEff1bq(n-Njets_[0]))*Ntt*eb_err,2)
  + pow(((2*eb*(1-eb)*(n-2)*pow(-1.,n-2)*pow(eudsc-1,n-3)+2*(eb-1)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(1-eudsc,n-4)))*GetTTEff2bq(n-Njets_[0])
         +(         eb*(n-1)*pow(-1.,n-1)*pow(eudsc-1,n-2)+  (1-eb)*(n-1)*(pow(1-eudsc,n-2)+eudsc*pow(-1.,n-2)*(n-2)*pow(1-eudsc,n-3)))*GetTTEff1bq(n-Njets_[0])
         +(                                                        (n  )*(pow(1-eudsc,n-1)+eudsc*pow(-1.,n-1)*(n-1)*pow(1-eudsc,n-2)))*GetTTEff0bq(n-Njets_[0]))*Ntt*eudsc_err,2)
  + pow(Ntt_1bjet(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_1bjet < 0 ? 0 : sqrt(Ntt_err_1bjet));
}
Double_t VJetEstimation::Nvb_1bjet(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) const{
	Double_t  Nvb_1bjet = (eb*pow(1-euds,n-1)+(1-eb)*(n-1)*euds*pow(1-euds,n-2))*Nvb;
	return (Nvb_1bjet < 0 ? 0 : Nvb_1bjet);
}
Double_t VJetEstimation::Nv_1bjet(Double_t Nv, Double_t euds, Int_t n) const{
	Double_t  Nv_1bjet = n*euds*pow(1-euds,n-1)*Nv;
	return (Nv_1bjet < 0 ? 0 : Nv_1bjet);
}
Double_t VJetEstimation::Nv_err_1bjet(Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) const{
	Double_t  Nv_err_1bjet = pow(n*(pow(1-euds,n-1)+euds*pow(-1.,n-1)*(n-1)*pow(euds-1,n-2))*Nv*euds_err,2)
  + pow(Nv_1bjet(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_1bjet < 0 ? 0 : sqrt(Nv_err_1bjet));
}

Double_t VJetEstimation::N2bjets(Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) const{
	Double_t  Ntt_2bjets = this->Ntt_2bjets(Ntt,eb,eudsc,n);
	Double_t  Nv_2bjets  = this->Nv_2bjets( Nv ,euds, n);
	Double_t  Nvb_2bjets = this->Nvb_2bjets(Nvb,eb,euds, n);
	return (Ntt_2bjets+Nv_2bjets+Nvb_2bjets);
}
Double_t VJetEstimation::Ntt_2bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) const{
	Double_t  Ntt_2bjets =((eb*eb*pow(1-eudsc,n-2)+2*eb*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3)+(1-eb)*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4))*GetTTEff2bq(n-Njets_[0])
                         +                                 (eb*(n-1)*eudsc*pow(1-eudsc,n-2)+       (1-eb)*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*GetTTEff1bq(n-Njets_[0])
                         +                                                                                   ((n*(n-1)/2)*eudsc*eudsc*pow(1-eudsc,n-2))*GetTTEff0bq(n-Njets_[0]))*Ntt;
	return (Ntt_2bjets<0 ? 0 : Ntt_2bjets);
}
Double_t VJetEstimation::Ntt_err_2bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb,  Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) const{
	Double_t  Ntt_err_2bjets  = pow(((2*eb*pow(1-eudsc,n-2)+2*(1-2*eb)*(n-2)*eudsc*pow(1-eudsc,n-3)+2*(eb-1)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4))*GetTTEff2bq(n-Njets_[0])
                                   +                                 ((n-1)*eudsc*pow(1-eudsc,n-2)+    (-1)*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*GetTTEff1bq(n-Njets_[0]))*Ntt*eb_err,2)
  + pow(((pow(eb,2)*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)+2*eb*(1-eb)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+pow(1-eb,2)*((n-2)*(n-3)/2)*(2*eudsc*pow(1-eudsc,n-4)+pow(eudsc,2)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5)))*GetTTEff2bq(n-Njets_[0])
         +                                                      (eb*(n-1)*(pow(1-eudsc,n-2)+eudsc*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3))+     (1-eb)*((n-1)*(n-2)/2)*(2*eudsc*pow(1-eudsc,n-3)+pow(eudsc,2)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))))*GetTTEff1bq(n-Njets_[0])
        +                                                                                                                 			(((n)*(n-1)/2)*(2*eudsc*pow(1-eudsc,n-2)+pow(eudsc,2)*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)))*GetTTEff0bq(n-Njets_[0])*Ntt*eudsc_err,2)
  + pow(Ntt_2bjets(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_2bjets < 0 ? 0 : sqrt(Ntt_err_2bjets));
}

Double_t VJetEstimation::Nvb_2bjets(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) const{
	Double_t  Nvb_2bjets = (eb*(n-1)*euds*pow(1-euds,n-2)+(1-eb)*((n-1)*(n-2)/2)*euds*euds*pow(1-euds,n-3))*Nvb;
	return (Nvb_2bjets<0 ? 0 : Nvb_2bjets);
}
Double_t VJetEstimation::Nv_2bjets (Double_t Nv, Double_t euds, Int_t n) const{
	Double_t  Nv_2bjets  = ((n*(n-1)/2)*euds*euds*pow((1-euds),n-2))*Nv;
	return (Nv_2bjets<0 ? 0 : Nv_2bjets);
}
Double_t VJetEstimation::Nv_err_2bjets(Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) const{
	Double_t  Nv_err_2bjets = pow((n*(n-1)/2)*(2*euds*pow(1-euds,n-2)+pow(euds,2)*pow(-1.,n-2)*(n-2)*pow(euds-1,n-3))*Nv*euds_err,2)
  + pow(Nv_2bjets(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_2bjets < 0 ? 0 : sqrt(Nv_err_2bjets));
}

Double_t VJetEstimation::N3bjets   (Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) const{
  Double_t Ntt_3bjets = this->Ntt_3bjets(Ntt, eb, eudsc,n);
	Double_t Nv_3bjets  = this->Nv_3bjets( Nv , euds, n);
	Double_t Nvb_3bjets = this->Nvb_3bjets(Nvb, eb, euds, n);
	return (Ntt_3bjets+Nv_3bjets+Nvb_3bjets);
}
Double_t VJetEstimation::Ntt_3bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) const{
	Double_t  Ntt_3bjets =((eb*eb*(n-2)*eudsc*pow(1-eudsc,n-3)+2*eb*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4)+(n>4 ? pow((1-eb),2)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow((1-eudsc),n-5) : 0 ))*GetTTEff2bq(n-Njets_[0])
                         +                                              (eb*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3)+ (1-eb)*((n-1)*(n-2)*(n-3)/6)*pow(eudsc,3)*pow(1-eudsc,n-4))*GetTTEff1bq(n-Njets_[0])
                         +                                                                                                          ((n*(n-1)*(n-2)/6)*pow(eudsc,3)*pow(1-eudsc,n-3))*GetTTEff0bq(n-Njets_[0]))*Ntt;
	return (Ntt_3bjets<0 ? 0 : Ntt_3bjets);
}
Double_t VJetEstimation::Ntt_err_3bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb,  Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) const{
	Double_t  Ntt_err_3bjets  = pow(((2*eb*(n-2)*eudsc*pow(1-eudsc,n-3)+2*(1-2*eb)*(n-2)*(n-3)/2*pow(eudsc,2)*pow(1-eudsc,n-4)+2*(eb-1)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow(1-eudsc,n-5))*GetTTEff2bq(n-Njets_[0])
                                   +                                             ((n-1)*(n-2)/2*pow(eudsc,2)*pow(1-eudsc,n-3)+    (-1)*((n-1)*(n-2)*(n-3)/6)*pow(eudsc,3)*pow(1-eudsc,n-4))*GetTTEff1bq(n-Njets_[0]))*Ntt*eb_err,2)
  + pow(((pow(eb,2)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+2*eb*(1-eb)*((n-2)*(n-3)/2)*(2*eudsc*pow(1-eudsc,n-4)+pow(eudsc,2)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5))+pow(1-eb,2)*((n-2)*(n-3)*(n-4)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-5)+pow(eudsc,3)*pow(-1.,n-5)*(n-5)*pow(eudsc-1,n-6)))*GetTTEff2bq(n-Njets_[0])
         +                                                     			            (eb*((n-1)*(n-2)/2)*(2*eudsc*pow(1-eudsc,n-3)+pow(eudsc,2)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+     (1-eb)*((n-1)*(n-2)*(n-3)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-4)+pow(eudsc,3)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5))))*GetTTEff1bq(n-Njets_[0])
        +                                                                                                                 				      			      		        (((n)*(n-1)*(n-2)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-3)+pow(eudsc,3)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4)))*GetTTEff0bq(n-Njets_[0])*Ntt*eudsc_err,2)
  + pow(Ntt_3bjets(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_3bjets < 0 ? 0 : sqrt(Ntt_err_3bjets));
}
Double_t VJetEstimation::Nvb_3bjets(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) const{
	Double_t  Nvb_3bjets = (eb*((n-1)*(n-2)/2)*euds*euds*pow(1-euds,n-3) + (1-eb)*((n-1)*(n-2)*(n-3)/6)*pow(euds,3)*pow(1-euds,n-4))*Nvb;
	return (Nvb_3bjets<0 ? 0 : Nvb_3bjets);
}
Double_t VJetEstimation::Nv_3bjets (Double_t Nv, Double_t euds, Int_t n) const{
	Double_t  Nv_3bjets  = ((n*(n-1)*(n-2)/6)*pow((euds),3)*pow((1-euds),n-3))*Nv;
	return (Nv_3bjets<0 ? 0 : Nv_3bjets);
}
Double_t VJetEstimation::Nv_err_3bjets(Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) const{
	Double_t  Nv_err_3bjets = pow((n*(n-1)*(n-2)/6)*(3*pow(euds,2)*pow(1-euds,n-3)+pow(euds,3)*pow(-1.,n-3)*(n-3)*pow(euds-1,n-4))*Nv*euds_err,2)
  + pow(Nv_3bjets(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_3bjets < 0 ? 0 : sqrt(Nv_err_3bjets));
}
