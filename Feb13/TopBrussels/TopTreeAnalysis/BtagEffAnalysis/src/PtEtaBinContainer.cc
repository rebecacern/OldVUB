//#include "/user/jmmaes/CMSSW/CMSSW_2_2_9_Analysis/CMSSW_2_2_9/src/TopBrussels/bTagAnalysis/interface/PtEtaBinContainer.h"
#include "../interface/PtEtaBinContainer.h"

 
//constructor
PtEtaBinContainer::PtEtaBinContainer(int debug, int nVar1, int nVar0, int nBdiscrAlgos ,bool varBinSize, bool doShift, bool doPtEtaBin, int inBin, int outBin, int totalbins){

  debug_=debug;
  nVar1_=nVar1;
  nVar0_=nVar0;
  //nBdiscrAlgos_=1;
  nBdiscrAlgos_=nBdiscrAlgos;

  ///  cout << nPtBins_ << " " << nEtaBins_ << endl;
  
  nPtBins_=-1; 
  nEtaBins_=-1; 
  //cout << nPtBins_ << " " << nEtaBins_ << endl;

  //double ptBins_[21];
  //cout << "length " << sizeof(ptBins_) << " " << sizeof(ptBins_[0]) << " " << (int) ( (double) sizeof(ptBins_)/ (double) sizeof(ptBins_[0])) << endl;
   

  //  for(int i=0; i<=nPtBins_; i++){
    //for(int i=0; i<=24; i++){
    //cout << "before: " << i << " " << nPtBins_ << " " << nEtaBins_ << endl;
    //cout << "length " << sizeof(ptBins_) << " " << sizeof(ptBins_[0]) << " " << (int) ( (double) sizeof(ptBins_)/ (double) sizeof(ptBins_[0])) << endl;
    //ptBins_[i]=-1.;
    //nPtBins_=20; 
    // nEtaBins_=1; 
    //cout << "after : " << i << " " << nPtBins_ << " " << nEtaBins_ << endl;
  // }

  //for(int i=0; i<nPtBins_; i++){
    //cout << ptBins_[i] << endl;
  //}  
  //cout << nPtBins_ << " " << nEtaBins_ << endl;
  
  //for(int i=0; i<nEtaBins_; i++){
  //etaBins_[i]=0;
  //}
  //cout << nPtBins_ << " " << nEtaBins_ << endl;

  if(!doPtEtaBin){
    nEtaBins_=0;
    nPtBins_=0; 
  }
  //cout << nPtBins_ << " " << nEtaBins_ << endl;

  
  double tmpptUpCut[51]={999999.  
			 ,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500};
  double tmpptLowCut[51]={0.      
			  ,0.,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490};
  
  double tmpptUpCutComplement[51]={0.   
				   ,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500};
  double tmpptLowCutComplement[51]={0.  
				    ,0.,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490};
  
  
  
  double tmpptUpCut2[101]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,345,350,355,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500};
  

  if(doPtEtaBin){
    cout << "!!!!!!!!doPtEtaBin!!!!!!!!! " << "-> start" << endl;

    double ptbinstmp[6]={30,45,60,80,110,999999};
    nPtBins_=5;    
    double etabinstmp[6]={0,0.48,0.96,1.44,1.92,2.4};
    nEtaBins_=5;

    //    double ptbinstmp[9]={30.,40.,50.,60.,75.,90.,110.,150.,999999.};
    /*double ptbinstmp[9]={30.,40.,50.,60.,75.,90.,110.,150.,300.};
    nPtBins_=8;    
    double etabinstmp[9]={0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4};
    nEtaBins_=8;*/

    /*double ptbinstmp[21]={30,34,38,42,46,49,53,58,62,66,71,76,81,86,92,100,108,118,132,155,999999}; 
    nPtBins_=20;
    double etabinstmp[2]={-997.,997.};
    nEtaBins_=1;*/

    /*double etabinstmp[21]={0,0.12,0.24,0.36,0.48,0.6,0.72,0.84,0.96,1.08,1.2,1.32,1.44,1.56,1.68,1.8,1.92,2.04,2.16,2.28,2.4};
    nEtaBins_=20;
    double ptbinstmp[2]={0,999998.};
    nPtBins_=1;*/
    
    
    //double ptbinstmp[4]={30,65,100,999999};
    //nPtBins_=3;    
    //double ptbinstmp[31]={30,32,35,38,40,43,46,48,51,53,56,59,62,65,68,71,74,77,81,84,88,92,97,102,108,114,122,132,145,168,999999};
    //nPtBins_=30;
    //double ptbinstmp[21]={30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,999999};
    //nPtBins_=20;
    
    //cout << nPtBins_ << " " << nEtaBins_ << endl;

    
    for(int i=0; i<=nPtBins_;i++){
      //cout << "ptdebug " << nPtBins_ << " " << i << " " << ptbinstmp[i] << " " << ptBins_[i]<< endl;
      ptBins_[i]=ptbinstmp[i];  
      //nPtBins_=20;
      //nEtaBins_=1;
      //cout << "ptdebug " << nPtBins_ << " " << i << " " << ptbinstmp[i] << " " << ptBins_[i]<< endl;
  
    }
    
    //cout << nPtBins_ << " " << nEtaBins_ << endl;
    for(int i=0; i<=nEtaBins_; i++){
      //      cout << " debug 5 " << i << endl;
      //cout << " debug 5 " << etabinstmp[i] << endl;
      etaBins_[i]=etabinstmp[i];
      //cout << " debug 6 " << endl;
      
    }
    
    
    cout << "!!!!!!!!doPtEtaBin!!!!!!!!! " << nPtBins_ << " " << nEtaBins_ << endl;
  }  
  cout << "nb of pt-eta bins " << nPtBins_ << " " << nEtaBins_ << endl;

  vector<double> TheCrasher;


  if(totalbins!=nEtaBins_+nPtBins_){
    cout << "totalbins " << totalbins << " != nEtaBins_+nPtBins_ " << nEtaBins_+nPtBins_ << " crashing in 3 2 1..." <<endl; //haha, I'm to lazy to implement this all in a better/safer way
    cout << TheCrasher[2] << endl;
  }
  
  ptUpCut_ = tmpptUpCut[inBin];//obsolete, but too lazy to remove everywhere in the code
  ptLowCut_ = tmpptLowCut[inBin];//obsolete
  
  ptUpCutComplement_ = tmpptUpCutComplement[outBin];//obsolete
  ptLowCutComplement_ = tmpptLowCutComplement[outBin];//obsolete
  
  cout << "Following pt-range will be included : " << ptLowCut_ << " " << ptUpCut_ << endl; //obsolete
  cout << "Following pt-range will be excluded : " << ptLowCutComplement_ << " " << ptUpCutComplement_ << endl;//obsolete

  binGlobal_ = new PtEtaBin*[nBdiscrAlgos_];
  for(int i=0; i<nPtBins_; i++){
    //cout << i << endl;
    binVector_[i] = new PtEtaBin*[nBdiscrAlgos_];
  }
  for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
    //cout << j << endl;
    binVector_[j] = new PtEtaBin*[nBdiscrAlgos_];
  }
  
  for (int k=0; k<nBdiscrAlgos_; k++){
    PtEtaBin *tempBin = new PtEtaBin(debug_, nVar1_, nVar0_, k, 0.,-999.,999999., 999999., varBinSize, doShift);
    binGlobal_[k]=tempBin;


    for(int i=0; i<nPtBins_; i++){
      //cout << "k: " << k << " - i: " << i << endl;
      PtEtaBin *tempBinV = new PtEtaBin(debug_, nVar1_, nVar0_, k, ptBins_[i] ,-998., ptBins_[i+1] ,998., varBinSize, doShift);
     binVector_[i][k]=tempBinV;
    }
    
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      int index=j-nPtBins_;
    //cout << "k: " << k << " - j: "  << j << endl;
      PtEtaBin *tempBinV = new PtEtaBin(debug_, nVar1_, nVar0_, k,0 ,etaBins_[index], 999998. ,etaBins_[index+1], varBinSize, doShift);
      binVector_[j][k]=tempBinV;
    }

  }

  cout << "done constructor PtEtaBinContainer" << endl;
  //cout << TheCrasher[2] << endl;

}

//destructor
PtEtaBinContainer::~PtEtaBinContainer(){

  //cout << "start deleting containers" << endl;
  
  for (int k=0; k<nBdiscrAlgos_; k++){ 
    //cout << "deleting container " << k  << endl;
    delete binGlobal_[k];
    //cout << "DONE deleting container " << k  << endl;
    for(int i=0; i<nPtBins_; i++){
      //cout << "deleting bin-container " << i << endl;
      delete binVector_[i][k];
      //cout << "DONE deleting bin-container " << i << endl;
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      //cout << "deleting bin-container " << j << endl;
      delete binVector_[j][k];
      //cout << "DONE deleting bin-container " << j << endl;
    }
  }

  //cout << "deleting container" << endl;
  delete [] binGlobal_;
  ///cout << "DONE deleting container global" << endl;
}

void PtEtaBinContainer::DefineSignalSamplePlots(int nBinsVar1, double lowRangVar1, double upRangeVar1, int nBinsVar2, double lowRangVar2, double upRangeVar2, int nBinsVar1SC, double lowRangVar1SC, double upRangeVar1SC, int nBinsVar2SC, double lowRangVar2SC, double upRangeVar2SC, int nBinsVar0, double lowRangeVar0, double upRangeVar0, int nBinsBtag[], double lowRangeBtag[], double upRangeBtag[]){

  vector<double> TheCrasher;
  //cout << TheCrasher[2] << endl;

  //if(debug_>1) cout << "PtEtaBinContainer::DefineSignalSamplePlots::START" << endl;
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      //if(debug_>1) cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << " i " << i << endl;
      //cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << " i " << i << endl;
      binVector_[i][k]->DefineSignalSamplePlots(k,nBinsVar1,lowRangVar1,upRangeVar1,nBinsVar2,lowRangVar2,upRangeVar2,nBinsVar1SC,lowRangVar1SC,upRangeVar1SC,nBinsVar2SC,lowRangVar2SC,upRangeVar2SC,nBinsVar0,lowRangeVar0,upRangeVar0,nBinsBtag,lowRangeBtag,upRangeBtag);
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      //if(debug_>1) cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << " j " << j << endl;
      //cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << " j " << j << endl;
      binVector_[j][k]->DefineSignalSamplePlots(k,nBinsVar1,lowRangVar1,upRangeVar1,nBinsVar2,lowRangVar2,upRangeVar2,nBinsVar1SC,lowRangVar1SC,upRangeVar1SC,nBinsVar2SC,lowRangVar2SC,upRangeVar2SC,nBinsVar0,lowRangeVar0,upRangeVar0,nBinsBtag,lowRangeBtag,upRangeBtag);
    }


    //if(debug_>1) cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << endl;
    //cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << endl;

    binGlobal_[k]->DefineSignalSamplePlots(k,nBinsVar1,lowRangVar1,upRangeVar1,nBinsVar2,lowRangVar2,upRangeVar2,nBinsVar1SC,lowRangVar1SC,upRangeVar1SC,nBinsVar2SC,lowRangVar2SC,upRangeVar2SC,nBinsVar0,lowRangeVar0,upRangeVar0,nBinsBtag,lowRangeBtag,upRangeBtag);
    //if(debug_>1) cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << endl;
    //cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << endl;

    //if(k>0)    cout << TheCrasher[2] << endl;

  }
  //if(debug_>1) cout << "PtEtaBinContainer::DefineSignalSamplePlots::END" << endl;

}

void PtEtaBinContainer::SetVarBins(std::map<int,vector<float> > rangesbTag) {

	for (int k=0; k<nBdiscrAlgos_; k++){
		for(int i=0; i<nPtBins_; i++){ 
			//if(debug_>1) cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << " i " << i << endl;
			//cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << " i " << i << endl;
			binVector_[i][k]->SetVarBins(rangesbTag[k]);
		}
		for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
			//if(debug_>1) cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << " j " << j << endl;
			//cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << " j " << j << endl;
			binVector_[j][k]->SetVarBins(rangesbTag[k]);
		}
		
		
		//if(debug_>1) cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << endl;
		//cout << "PtEtaBinContainer::DefineSignalSamplePlots " << " k "<< k << endl;
		
		binGlobal_[k]->SetVarBins(rangesbTag[k]);
	
	}
	
	
}

void PtEtaBinContainer::DefineControlSamplePlots(int nBinsControlVar1, double lowRangeControlVar1, double upRangeControlVar1, int nBinsControlVar2, double lowRangeControlVar2, double upRangeControlVar2, int nBinsControlVar, double lowRangeControlVar, double upRangeControlVar,int nBinsBtag[], double lowRangeBtag[], double upRangeBtag[]){
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->DefineControlSamplePlots(k,nBinsControlVar1,lowRangeControlVar1,upRangeControlVar1,nBinsControlVar2,lowRangeControlVar2,upRangeControlVar,nBinsControlVar,lowRangeControlVar,upRangeControlVar,nBinsBtag,lowRangeBtag,upRangeBtag);
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){	
      //cout << "PtEtaBinContainer::DefineControlSamplePlots " << " k "<< k << " i " << i << " j " << j << nBinsControlVar1 << " " << lowRangeControlVar1 << " " << upRangeControlVar1 << " " << nBinsControlVar2 << " " << lowRangeControlVar2 << " " << upRangeControlVar2 << " " << nBinsControlVar << " " << lowRangeControlVar << " " <<  upRangeControlVar << endl;
      binVector_[j][k]->DefineControlSamplePlots(k,nBinsControlVar1,lowRangeControlVar1,upRangeControlVar1,nBinsControlVar2,lowRangeControlVar2,upRangeControlVar,nBinsControlVar,lowRangeControlVar,upRangeControlVar,nBinsBtag,lowRangeBtag,upRangeBtag);
    }
    
    binGlobal_[k]->DefineControlSamplePlots(k,nBinsControlVar1,lowRangeControlVar1,upRangeControlVar1,nBinsControlVar2,lowRangeControlVar2,upRangeControlVar2,nBinsControlVar,lowRangeControlVar,upRangeControlVar,nBinsBtag,lowRangeBtag,upRangeBtag);
  }
}

void PtEtaBinContainer::SetErrorsSignalSamples(){
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->SetErrorsSignalSamples();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->SetErrorsSignalSamples();
    }
    binGlobal_[k]->SetErrorsSignalSamples();
  }
}

void PtEtaBinContainer::SetErrorsControlSamples(){
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->SetErrorsControlSamples();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->SetErrorsControlSamples();
    }
    binGlobal_[k]->SetErrorsControlSamples();
  }
}

//void PtEtaBinContainer::FillSignalSamplePlots(double weight, int partonFlavour, double *bTag, double var1, double var2, double var0, int partonFlavourControl, double bTagControl, double controlVar1, double controlVar2, double controlVar0, double lowCutVar0, double centralCutVar0, double upCutVar0){
void PtEtaBinContainer::FillSignalSamplePlots(double weight, double weight_nonrew, int partonFlavour, bool isW, bool isR, double chisq, double *bTag, std::map<int,vector<double> > WPMap, double var1, double var2, double var0, double var3, double lowCutVar0, double centralLowCutVar0, double centralUpCutVar0, double upCutVar0){
  for (int k=0; k<nBdiscrAlgos_; k++){
	  
	  double *bTagCuts = new double[3];
	  
	  if (WPMap.find(k) != WPMap.end()) {
		  if (WPMap[k].size() == 3) {
			  
			  bTagCuts[0]=WPMap[k][0];
			  bTagCuts[1]=WPMap[k][1];
			  bTagCuts[2]=WPMap[k][2];
			  
		  } else {
			  cout << "PtEtaBinContainer::FillXStemplates::ERROR - WMap["<<k<<"] contains != 3 cuts, code is build to accomodate 3 WP's (L,M,T)" << endl;
		  }
		  
	  } //else {
	  //cout << "PtEtaBinContainer::FillXStemplates::ERROR - WMap["<<k<<"] does not exist" << endl;
	  //}

	  
    for(int i=0; i<nPtBins_; i++){ 
      if(var1>ptBins_[i] && var1<ptBins_[i+1]){
	binVector_[i][k]->FillSignalSamplePlots(weight,weight_nonrew,partonFlavour,isW,isR,chisq,bTag[k],bTagCuts,var1,var2,var0,var3,lowCutVar0,centralLowCutVar0,centralUpCutVar0,upCutVar0,bTag[0]);
      }
    }

    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      int index=j-nPtBins_;
      if(var2>etaBins_[index] && var2<etaBins_[index+1]){// this is the absolute value of the pseudo-rapidity

	binVector_[j][k]->FillSignalSamplePlots(weight,weight_nonrew,partonFlavour,isW,isR,chisq,bTag[k],bTagCuts,var1,var2,var0,var3,lowCutVar0,centralLowCutVar0,centralUpCutVar0,upCutVar0,bTag[0]);
      } //else {cout << "WARNING:: this events falls outside the eta ranges " << var2 << endl;}
	
    }
    
    //if(var1>ptLowCut_ && var1<ptUpCut_){
    if(var1>ptLowCut_ && var1<ptUpCut_){
      if(!(var1>ptLowCutComplement_ && var1<ptUpCutComplement_)){
	binGlobal_[k]->FillSignalSamplePlots(weight,weight_nonrew,partonFlavour,isW,isR,chisq,bTag[k],bTagCuts,var1,var2,var0,var3,lowCutVar0,centralLowCutVar0,centralUpCutVar0,upCutVar0,bTag[0]);
      }
    }
  }
}

void PtEtaBinContainer::FillControlSamplePlots(double weight, int partonFlavour, bool isW, bool isR, double chisq, double *bTag, double controlVar1, double controlVar2, double controlVar, double lowCutVar0, double centralLowCutVar0, double centralUpCutVar0, double upCutVar0){
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      if(controlVar1>ptBins_[i] && controlVar1<ptBins_[i+1]){
	binVector_[i][k]->FillControlSamplePlots(weight,partonFlavour,isW,isR,chisq,bTag[k],controlVar1,controlVar2,controlVar,lowCutVar0,centralLowCutVar0,centralUpCutVar0,upCutVar0);
      }
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      int index=j-nPtBins_;
      if(controlVar2>etaBins_[index] && controlVar2<etaBins_[index+1]){
	binVector_[j][k]->FillControlSamplePlots(weight,partonFlavour,isW,isR,chisq,bTag[k],controlVar1,controlVar2,controlVar,lowCutVar0,centralLowCutVar0,centralUpCutVar0,upCutVar0);
      }
    }
    
    
    //if(controlVar1>ptLowCut_ && controlVar1<ptUpCut_){
    if(controlVar1>ptLowCut_ && controlVar1<ptUpCut_){
      if(!(controlVar1>ptLowCutComplement_ && controlVar1<ptUpCutComplement_)){
	binGlobal_[k]->FillControlSamplePlots(weight,partonFlavour,isW,isR,chisq,bTag[k],controlVar1,controlVar2,controlVar,lowCutVar0,centralLowCutVar0,centralUpCutVar0,upCutVar0);
      }
    }
 }
}

void PtEtaBinContainer::FillXStemplates(double weight, string datasetname, int partonFlavour, double *btag, std::map<int,vector<double> > WPMap, double controlVar0, double m3, double lowCutVar0, double centralLowCutVar0, double centralUpCutVar0, double upCutVar0) {
	for (int k=0; k<nBdiscrAlgos_; k++){
		
		double *btagCuts = new double[3];
		
		if (WPMap.find(k) != WPMap.end()) {
			if (WPMap[k].size() == 3) {
			
				btagCuts[0]=WPMap[k][0];
				btagCuts[1]=WPMap[k][1];
				btagCuts[2]=WPMap[k][2];
					
			} else {
				cout << "PtEtaBinContainer::FillXStemplates::ERROR - WMap["<<k<<"] contains != 3 cuts, code is build to accomodate 3 WP's (L,M,T)" << endl;
				exit(1);
			}

		} //else {
			//cout << "PtEtaBinContainer::FillXStemplates::ERROR - WMap["<<k<<"] does not exist" << endl;
		//}
		
		/*cout << "PtEtaBinContainer::FillXStemplates::bTag[" << k << "] " << btag[k] << endl;
		
		cout << "PtEtaBinContainer::FillXStemplates::bTagCuts[0] " << btagCuts[0] << endl;
		cout << "PtEtaBinContainer::FillXStemplates::bTagCuts[1] " << btagCuts[1] << endl;
		cout << "PtEtaBinContainer::FillXStemplates::bTagCuts[2] " << btagCuts[2] << endl;*/
		
		binGlobal_[k]->FillXStemplates(weight,datasetname,partonFlavour,btag[k],btagCuts,controlVar0,m3, lowCutVar0,centralLowCutVar0,centralUpCutVar0,upCutVar0);

	}
	
	//exit(1);
}

void PtEtaBinContainer::MakeSoverSBPlots(){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->MakeSoverSBPlots();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->MakeSoverSBPlots();
    }
    binGlobal_[k]->MakeSoverSBPlots();
  }
}

void PtEtaBinContainer::MakeMCEffPlots(){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->MakeMCEffPlots();
    } 
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->MakeMCEffPlots();
    }
    binGlobal_[k]->MakeMCEffPlots();
  }
}

void PtEtaBinContainer::MakeProfileXplots(){
  cout << "WARNING!! PtEtaBinContainer::MakeProfileXplots::WARNING  be careful with the profile plot. It is only made starting from 0 for the tagger, (this might not be the full range, nor the variable bins might not start at 0)" << endl;
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->MakeProfileXplots();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->MakeProfileXplots();
    }
    binGlobal_[k]->MakeProfileXplots();
  }
}

void PtEtaBinContainer::Make1DXplots(){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->Make1DXplots();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->Make1DXplots();
    }
    binGlobal_[k]->Make1DXplots();
  }
}

void PtEtaBinContainer::Make1DYplots(){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->Make1DYplots();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->Make1DYplots();
    }
    binGlobal_[k]->Make1DYplots();
  }
}

void PtEtaBinContainer::Make1DYVar2plots(){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->Make1DYVar2plots();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->Make1DYVar2plots();
    }
    binGlobal_[k]->Make1DYVar2plots();
  }
}

void PtEtaBinContainer::SetError1DXplots(){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->SetError1DXplots();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->SetError1DXplots();
    }
    binGlobal_[k]->SetError1DXplots();
  }
}

void PtEtaBinContainer::SetError1DYplots(){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->SetError1DYplots();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->SetError1DYplots();
    }
    binGlobal_[k]->SetError1DYplots();
  }
}

void PtEtaBinContainer::SetError1DYVar2plots(){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->SetError1DYVar2plots();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->SetError1DYVar2plots();
    }
    binGlobal_[k]->SetError1DYVar2plots();
  }
}

void PtEtaBinContainer::MakeXRatioPlot(bool useFit){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      //cout << "fit params for global bin, btagger " << k  << "  bin i: " << i << endl;
      binVector_[i][k]->MakeXRatioPlot(useFit);
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){ 
      //cout << "fit params for global bin, btagger " << k  << "  bin j: " << j << endl;
      binVector_[j][k]->MakeXRatioPlot(useFit);
    }
    //cout << "fit params for global bin, btagger " << k << endl;
    binGlobal_[k]->MakeXRatioPlot(useFit);
  }
}

void PtEtaBinContainer::ChangeLeftRightPars(){
  
  double *par0=new double(0);
  double *par1=new double(0);
  double *par2=new double(0);
  
  for (int k=0; k<nBdiscrAlgos_; k++){
    
    binGlobal_[k]->GetLeftRightPars(par0,par1,par2);
    
    /*cout << "PtEtaBinContainer::ChangeLeftRightPars "<< k << " par0 " << *par0 << endl; 
    cout << "PtEtaBinContainer::ChangeLeftRightPars "<< k << " par1 " << *par1 << endl; 
    cout << "PtEtaBinContainer::ChangeLeftRightPars "<< k << " par2 " << *par2 << endl; */
    
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->SetLeftRightPars(par0,par1,par2);
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->SetLeftRightPars(par0,par1,par2);
    }
    
  }
}

void PtEtaBinContainer::MakeSCVar12RatioPlot(){
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->MakeSCVar12RatioPlot();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->MakeSCVar12RatioPlot();
    }
    binGlobal_[k]->MakeSCVar12RatioPlot();
  }
}

void PtEtaBinContainer::MakeSCprojectionPlots(){
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->MakeSCprojectionPlots();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->MakeSCprojectionPlots();
    }
    binGlobal_[k]->MakeSCprojectionPlots();
  }
}

void PtEtaBinContainer::Make2DRatioPlot(){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->Make2DRatioPlot();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->Make2DRatioPlot();
    }
    binGlobal_[k]->Make2DRatioPlot();
  }
}

void PtEtaBinContainer::GetLRratio(double FMCBias,bool doNewF,double *newFpt, double *newFeta){

  cout << "PtEtaBinContainer::GetLRratio()::------------------------------------" << endl;
  
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      if(k==0) {
	cout << "bin: pt - eta " << i+1 << " - " << 0 << " new F_pt " << newFpt[i] << endl;
	binVector_[i][k]->GetLRratio(false,true,FMCBias,doNewF,newFpt[i]);
      }
      if(k!=0) binVector_[i][k]->GetLRratio(false,false,FMCBias,doNewF,newFpt[i]);
    }
    for(int j=nPtBins_; j<nEtaBins_+nPtBins_; j++){
      int index=j-nPtBins_;
      if(k==0) {
	cout << "bin: pt - eta " << 0 << " - " << j-nPtBins_+1 << " new F_eta " << newFeta[index] << endl;
	binVector_[j][k]->GetLRratio(false,true,FMCBias,doNewF,newFeta[index]);
      }
      if(k!=0) binVector_[j][k]->GetLRratio(false,false,FMCBias,doNewF,newFeta[index]);
    }


    if(k==0) {
      cout << "bin: global   " << endl;
      binGlobal_[k]->GetLRratio(true,true,FMCBias,false,-1);
    }
    if(k!=0) binGlobal_[k]->GetLRratio(false,false,FMCBias,false,-1);
  }

  cout << "PtEtaBinContainer::GetLRratio()::------------------------------------" << endl;

}

void PtEtaBinContainer::GetPercentiles(double left, double center, double right, int *returnval){

  cout << "PtEtaBinContainer::GetPercentiles()::------------------------------------" << endl;
  
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      //binVector_[i][k]->GetPercentiles(left,center,right);
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      //binVector_[j][k]->GetPercentiles(left,center,right);
    }
    binGlobal_[k]->GetPercentiles(left,center,right,returnval);
  }

  cout << "PtEtaBinContainer::GetPercentiles()::------------------------------------" << endl;

}

void PtEtaBinContainer::ReweighLeft(){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->ReweighLeft();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->ReweighLeft();
    }
    binGlobal_[k]->ReweighLeft();
  }
}

void PtEtaBinContainer::ReweighRight(){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->ReweighRight();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->ReweighRight();
    }
    
    binGlobal_[k]->ReweighRight();
  }
}

/*obsolete void PtEtaBinContainer::ReweighRightChangeFitParams(double par0, double par1){
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      for(int j=0; j<nEtaBins_; j++){
	binVector_[i][j][k]->ReweighRightChangeFitParams(par0, par1);
      }
    }
    binGlobal_[k]->ReweighRightChangeFitParams(par0, par1);
  }
  }*/

//void PtEtaBinContainer::FillReweighRight(bool do2D, bool useFit, double weight, int partonFlavour, double *bTag, double var1, double var2, double var0, int partonFlavourControl, double bTagControl, double controlVar1, double controlVar2, double controlVar0, double lowCutVar0, double centralCutVar0, double upCutVar0){
void PtEtaBinContainer::FillReweighRight(bool do2D, bool useFit, double weight, int partonFlavour, double *bTag, double var1, double var2, double var0, double lowCutVar0, double centralLowCutVar0, double centralUpCutVar0, double upCutVar0){
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      if(var1>ptBins_[i] && var1<ptBins_[i+1]){
	binVector_[i][k]->FillReweighRight(do2D, useFit, weight, partonFlavour, bTag[k], var1, var2, var0, lowCutVar0, centralLowCutVar0, centralUpCutVar0, upCutVar0, ptBins_[i], ptBins_[i+1],-1,-1);
      }
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      int index=j-nPtBins_;
      if(var2>etaBins_[index] && var2<etaBins_[index+1]){
	//cout << "index " << index << endl;
	binVector_[j][k]->FillReweighRight(do2D, useFit, weight, partonFlavour, bTag[k], var1, var2, var0, lowCutVar0, centralLowCutVar0, centralUpCutVar0, upCutVar0, -1 ,-1 , etaBins_[index], etaBins_[index+1]);
      }
    }
    


    //binGlobal_[k]->FillReweighRight(do2D, useFit, weight, partonFlavour, bTag[k], var1, var2, var0, partonFlavourControl,bTagControl,controlVar1,controlVar2,controlVar0, lowCutVar0, centralCutVar0, upCutVar0);
  

    //if(var1>ptLowCut_ && var1<ptUpCut_){
    if(var1>ptLowCut_ && var1<ptUpCut_){
      if(!(var1>ptLowCutComplement_ && var1<ptUpCutComplement_)){
	
	binGlobal_[k]->FillReweighRight(do2D, useFit, weight, partonFlavour, bTag[k], var1, var2, var0, lowCutVar0, centralLowCutVar0, centralUpCutVar0, upCutVar0, -1, -1, -1, -1);
      }
    }
  }
}

//void PtEtaBinContainer::FillReweighControl(bool do2D, bool useFit, double weight, int partonFlavour, double *bTag, double var1, double var2, double var0, int partonFlavourControl, double bTagControl, double controlVar1, double controlVar2, double controlVar0, double lowCutVar0, double centralCutVar0, double upCutVar0){
void PtEtaBinContainer::FillReweighControl(double* btag, double* btagCut, bool do2D, double weight, int partonFlavour, bool isW, bool isR, double controlVar1, double controlVar2, double controlVar0, double lowCutVar0, double centralLowCutVar0, double centralUpCutVar0, double upCutVar0, double chisq){
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      if(controlVar1>ptBins_[i] && controlVar1<ptBins_[i+1]){
	binVector_[i][k]->FillReweighControl(btag[k],btagCut,false, weight, partonFlavour, isW, isR, controlVar1,controlVar2,controlVar0, lowCutVar0, centralLowCutVar0, centralUpCutVar0, upCutVar0, chisq);//switch of the 2D reweigh to reduce the runtime 
      }
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      int index=j-nPtBins_;
      if(controlVar2>etaBins_[index] && controlVar2<etaBins_[index+1]){
	binVector_[j][k]->FillReweighControl(btag[k],btagCut,false, weight, partonFlavour, isW, isR, controlVar1,controlVar2,controlVar0, lowCutVar0, centralLowCutVar0, centralUpCutVar0, upCutVar0, chisq);//switch of the 2D reweigh to reduce the runtime 
      }
    }

    
    //binGlobal_[k]->FillReweighControl(do2D, useFit, weight, partonFlavour, bTag[k], var1, var2, var0, partonFlavourControl,bTagControl,controlVar1,controlVar2,controlVar0, lowCutVar0, centralCutVar0, upCutVar0);

    //if(controlVar1>ptLowCut_ && controlVar1<ptUpCut_){
    if(controlVar1>ptLowCut_ && controlVar1<ptUpCut_){
      if(!(controlVar1>ptLowCutComplement_ && controlVar1<ptUpCutComplement_)){
	
	binGlobal_[k]->FillReweighControl(btag[k],btagCut,do2D, weight, partonFlavour, isW, isR, controlVar1,controlVar2,controlVar0, lowCutVar0, centralLowCutVar0, centralUpCutVar0, upCutVar0, chisq);
      }
    }
  }
}

void PtEtaBinContainer::MakeReweighRatio(){
  for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->MakeReweighRatio();
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->MakeReweighRatio();
    }
    binGlobal_[k]->MakeReweighRatio();
  }
}

void PtEtaBinContainer::MeasureEff(bool doSCreweigh){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->MeasureEff(false);
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->MeasureEff(false);
    }
    binGlobal_[k]->MeasureEff(doSCreweigh);
  }
}

void PtEtaBinContainer::MeasureEffRR(bool doSCreweigh){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->MeasureEffRR(false);
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->MeasureEffRR(false);
    }
    binGlobal_[k]->MeasureEffRR(doSCreweigh);
  }
}

void PtEtaBinContainer::MeasureMistagEffRR(bool doSCreweigh){
    for (int k=0; k<nBdiscrAlgos_; k++){
        for(int i=0; i<nPtBins_; i++){ 
            binVector_[i][k]->MeasureMistagEffRR(false);
        }
        for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
            binVector_[j][k]->MeasureMistagEffRR(false);
        }
        binGlobal_[k]->MeasureMistagEffRR(doSCreweigh);
    }
}

void PtEtaBinContainer::GetWPEff(bool RRincl, bool doSCreweigh, double wp, int tagger, double *array, bool pseudo, bool more, int ptbin, int etabin, bool ptetabin, double runNb){
 
  /*for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      for(int j=0; j<nEtaBins_; j++){
	//binVector_[i][j][k]->GetWPEffRR(wp,eff,MCeff,FRatio,MCFRatio,FitPar0,FitPar1);
      }
    }
    //cout << "PtEtaBinContainer::GetWPEffRR " << k << endl;
  return binGlobal_[k]->GetWPEffRR(wp, array, pseudo, more);
  }*/
  
  if(ptetabin){ 
    if(ptbin==0 && etabin!=0){ 
      int index=ptbin+etabin;
      binVector_[index][tagger]->GetWPEff(RRincl, doSCreweigh, wp, array, pseudo, more, runNb);
    }
  }
  if(ptetabin){ 
    if(ptbin!=0 && etabin==0){ 
      int index=ptbin;
	  binVector_[index][tagger]->GetWPEff(RRincl, doSCreweigh, wp, array, pseudo, more, runNb);
    }
  }
  if(!ptetabin) binGlobal_[tagger]->GetWPEff(RRincl, doSCreweigh, wp, array, pseudo, more, runNb);
	
	//cout << "TAMAGOTSHIIIIIIIIIIIIIIIII" << endl;

	//cout << "First cout " << array[8]<< endl;
	

}

void PtEtaBinContainer::CoutWPEff(bool RRincl, bool doSCreweigh, double wp, int tagger, double *array, bool pseudo, bool more, int ptbin, int etabin, bool ptetabin,double runNb){
 
  for(int i=0; i<nPtBins_; i++){ 
    binVector_[i][tagger]->CoutWPEff(RRincl, doSCreweigh, wp, array, pseudo, more, i+1, 0, runNb);
  }
  for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
    int index=j-nPtBins_;
    binVector_[j][tagger]->CoutWPEff(RRincl, doSCreweigh, wp, array, pseudo, more, 0, index+1, runNb);
  }


  //cout << "PtEtaBinContainer::GetWPEffRR " << k << endl;
  binGlobal_[tagger]->CoutWPEff(RRincl, doSCreweigh, wp, array, pseudo, more, -1, -1, runNb);
    

  //if(ptetabin) return binVector_[ptbin][etabin][tagger]->CoutWPEff(RRincl, doSCreweigh, wp, array, pseudo, more);
  //if(!ptetabin) return binGlobal_[tagger]->CoutWPEff(RRincl, doSCreweigh, wp, array, pseudo, more);
}

std::map<std::string,float> PtEtaBinContainer::GetEffCalcDetails(float wp, int tagger, int ptbin, int etabin, bool ptetabin) {
	
	std::map<std::string,float> tmpres, res;
	
	if(ptetabin){ 
		if(ptbin==0 && etabin!=0){ 
			int index=ptbin+etabin;
			tmpres = binVector_[index][tagger]->GetEffCalcDetails();
		}
	}
	if(ptetabin){ 
		if(ptbin!=0 && etabin==0){ 
			int index=ptbin;
			tmpres = binVector_[index][tagger]->GetEffCalcDetails();
		}
	}
	if(!ptetabin) tmpres = binGlobal_[tagger]->GetEffCalcDetails();
	
	stringstream s; s<<wp;
	
	string search = "WP_"+s.str();
	
	if (tmpres.find(search+"_min") != tmpres.end()) {
		
		res["F"]=tmpres["F"];
		res["Fmc"]=tmpres["Fmc"];
		res["Fe"]=tmpres["Fe"];
		res["Fmce"]=tmpres["Fmce"];
		res["min"]=tmpres[search+"_min"];
		res["plus"]=tmpres[search+"_plus"];
		res["minMC"]=tmpres[search+"_minMC"];
		res["plusMC"]=tmpres[search+"_plusMC"];
		
		res["NLeft"]=tmpres["NLeft"];
		res["NRight"]=tmpres["NRight"];
		
	}

	return res;
}

std::vector<string> PtEtaBinContainer::GetFitPlotPaths(int tagger) {
	
	return binGlobal_[tagger]->GetFitPlotPaths();
	
	
}


void PtEtaBinContainer::MeasureEffLR(bool doSCreweigh){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      binVector_[i][k]->MeasureEffLR(doSCreweigh);
    }
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      binVector_[j][k]->MeasureEffLR(doSCreweigh);
    }
    binGlobal_[k]->MeasureEffLR(doSCreweigh);
  }
}

std::map<int,vector<double> > PtEtaBinContainer::doMLJTemplateFit(string chi2cut,int mode, string data_postfix,int nSystematic) {

	std::map<int,vector<double> > results;
	
	for (int k=0; k<nBdiscrAlgos_; k++){

		if (debug_ > 0) cout << endl << "*********************************************************************************************" << endl;
		if (debug_ > 0) cout << "PtEtaBinContainer::doMLJTemplateFit **************** TAGGER " << k << " *******************" << endl;
		if (debug_ > 0) cout << "*********************************************************************************************" << endl << endl;
		
		vector<double> tmp = binGlobal_[k]->doMLJTemplateFit(chi2cut,mode,data_postfix,nSystematic);
		
		for (unsigned int i=0; i<tmp.size(); i++) results[k].push_back(tmp[i]);
		
		if (debug_ > 0) cout << "PtEtaBinContainer::doMLJTemplateFit::results.size() = " << results[k].size() << endl;
	
	}
		
	return results;
}

void PtEtaBinContainer::WriteContainerToRootFile(TFile *fout, bool bins, bool w1, bool w2, bool w3, bool w4, bool w5, bool w6, bool w7, bool w8, bool w9, bool w10, bool w11){

  //fout->cd();
  TString dirNameVar;
  TString dirNameBDisc;
  TString dirNameBin;
    

  dirNameVar="Variable_"; dirNameVar+=nVar1_; dirNameVar+="_"; dirNameVar+=nVar0_;   //std::cout << dirNameVar << std::endl;
  TDirectory *dirVar = fout->mkdir(dirNameVar);
  dirVar->cd();

  // bool dobin1cout=false;
  bool dobin1cout=true;
 
  for (int k=0; k<nBdiscrAlgos_; k++){ 
   
    dirNameBDisc="Discriminator_"; dirNameBDisc+=k; //std::cout << dirNameBDisc << std::endl;
    TDirectory *dirBDisc = dirVar->mkdir(dirNameBDisc);
    dirBDisc->cd();
 
    //if(k==0) dobin1cout=true; else dobin1cout=false;

    //cout << "Writing the plots for bin/tagger nr " << k << endl;
    binGlobal_[k]->WriteHistosToRootFile(w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,dobin1cout);
        
    for(int i=0; i<nPtBins_; i++){ 
      dirNameBin="Bin_Pt_"; dirNameBin+=(int) ptBins_[i]; //std::cout << dirNameBin << std::endl;
      TDirectory *dirBin = dirBDisc->mkdir(dirNameBin);
      dirBin->cd();
      if(bins) binVector_[i][k]->WriteHistosToRootFile(w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,false);
    }
    
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){ 
      //cout << "Writing the plots for bin/tagger nr " << k << "   nPtBins_ "  << nPtBins_ << " nEtaBins_ " << nEtaBins_ <<  endl;
      int index=j-nPtBins_;
      dirNameBin= "Bin_Eta_"; dirNameBin+= (int) round(etaBins_[index]*10); //std::cout << dirNameBin << std::endl;
      TDirectory *dirBin = dirBDisc->mkdir(dirNameBin);
      dirBin->cd();
      if(bins) binVector_[j][k]->WriteHistosToRootFile(w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,false);
      //if(bins) cout << "writing bin i: " << i << " j: " << j << " for tagger " << k << endl;	
    }

  }
  
  //fout->Close();
  
}

void PtEtaBinContainer::WriteHistoToPSFile(TString* pathname){
 for (int k=0; k<nBdiscrAlgos_; k++){
    for(int i=0; i<nPtBins_; i++){ 
      if(k==0) binVector_[i][k]->WriteHistoToPSFile(pathname,true);
      if(k!=0) binVector_[i][k]->WriteHistoToPSFile(pathname,false);
    }
    
    for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
      if(k==0) binVector_[j][k]->WriteHistoToPSFile(pathname,true);
      if(k!=0) binVector_[j][k]->WriteHistoToPSFile(pathname,false);
    }

    if(k==0) binGlobal_[k]->WriteHistoToPSFile(pathname,true);
    if(k!=0) binGlobal_[k]->WriteHistoToPSFile(pathname,false);
  }
}

double PtEtaBinContainer::getmljMeanVal(){return binGlobal_[0]->getmljMeanVal();}
double PtEtaBinContainer::getmljMeannoRWVal(){return binGlobal_[0]->getmljMeannoRWVal();}
double PtEtaBinContainer::getmljMeanMCVal(){return binGlobal_[0]->getmljMeanMCVal();}
double PtEtaBinContainer::getmljSigmaVal(){return binGlobal_[0]->getmljSigmaVal();}
double PtEtaBinContainer::getmljSigmanoRWVal(){return binGlobal_[0]->getmljSigmanoRWVal();}
double PtEtaBinContainer::getmljSigmaMCVal(){return binGlobal_[0]->getmljSigmaMCVal();}

double PtEtaBinContainer::getmlj_W_MeanVal(){return binGlobal_[0]->getmlj_W_MeanVal();}
double PtEtaBinContainer::getmlj_W_MeannoRWVal(){return binGlobal_[0]->getmlj_W_MeannoRWVal();}
double PtEtaBinContainer::getmlj_W_MeanMCVal(){return binGlobal_[0]->getmlj_W_MeanMCVal();}
double PtEtaBinContainer::getmlj_W_SigmaVal(){return binGlobal_[0]->getmlj_W_SigmaVal();}
double PtEtaBinContainer::getmlj_W_SigmanoRWVal(){return binGlobal_[0]->getmlj_W_SigmanoRWVal();}
double PtEtaBinContainer::getmlj_W_SigmaMCVal(){return binGlobal_[0]->getmlj_W_SigmaMCVal();}

double PtEtaBinContainer::getmlj_R_MeanVal(){return binGlobal_[0]->getmlj_R_MeanVal();}
double PtEtaBinContainer::getmlj_R_MeannoRWVal(){return binGlobal_[0]->getmlj_R_MeannoRWVal();}
double PtEtaBinContainer::getmlj_R_MeanMCVal(){return binGlobal_[0]->getmlj_R_MeanMCVal();}
double PtEtaBinContainer::getmlj_R_SigmaVal(){return binGlobal_[0]->getmlj_R_SigmaVal();}
double PtEtaBinContainer::getmlj_R_SigmanoRWVal(){return binGlobal_[0]->getmlj_R_SigmanoRWVal();}
double PtEtaBinContainer::getmlj_R_SigmaMCVal(){return binGlobal_[0]->getmlj_R_SigmaMCVal();}

void PtEtaBinContainer::setLumi(double lumi) {
    
    for (int k=0; k<nBdiscrAlgos_; k++){

        binGlobal_[k]->setLumi(lumi);
        
        for(int i=0; i<nPtBins_; i++){ 
            binVector_[i][k]->setLumi(lumi);
        }
        
        for(int j=nPtBins_; j<nPtBins_+nEtaBins_; j++){
            binVector_[j][k]->setLumi(lumi);
        }

    
    }
}

