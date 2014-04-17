#include "../interface/MCPseudoExp.h"

MCPseudoExp::MCPseudoExp(){
      //initialize private members
	Luminosity_ = 100.0;
	Xsection_ = 100.0;
	PreSelEfficiency_ = 1.0;
	MeanNEvents_ = 100;
	NEvents_ = 100;
}

MCPseudoExp::~MCPseudoExp(){}

void MCPseudoExp::GenerateRandomNumbers(int nNumbers, int eventtree_nentries){
      //generate nNumbers random integers in interval ]0;eventtree_nentries];
	RandomNumbers_.clear();
	//old, risk of taking same events multiple time
	/*int ir;	
	TRandom3 randomgen(0);  //this random generator performs the best 
 	for(int i=0; i<nNumbers; i++){		
		ir = int(randomgen.Uniform(0,eventtree_nentries));		
		RandomNumbers_.push_back(ir);
	}*/
	int ir;	
	TRandom3 randomgen(0);  //this random generator performs the best 
 	for(int i=0; i<nNumbers; i++){		
		ir = int(randomgen.Uniform(0,eventtree_nentries));		
		bool uniquenumber = false;
		while(!uniquenumber){ //it has to be checked if the chosen number is taken already or not
		   if(RandomNumbers_.size() > 0){
		     for(int u=0; u<RandomNumbers_.size(); u++){
			if(RandomNumbers_[u]==ir){
			   ir = (ir + 1) % eventtree_nentries;
			   uniquenumber = false;
			   break;
			}
			else uniquenumber = true; //at least temporary unique!
		     }
		   }
		   else uniquenumber = true;
		}
		RandomNumbers_.push_back(ir); //only unique numbers should be stored in this vector
	}
}

void MCPseudoExp::CalculatePseudoExp(const Dataset* dataset, float& luminosity){
      ////////////////////////////////////////
      // Calculate amount of events to be picked out of a Monte Carlo dataset, and generate a list of integers to be used later on to do the pseudoexperiment 
      ////////////////////////////////////////
        int verbosity = 0;	
	TChain* eventTree = dataset->eventTree();
	
	cout << "dataset name: " << dataset->Name() << endl;
	cout << "dataset filenames: " << endl;
	for(unsigned int i=0; i<dataset->Filenames().size(); i++)  cout  << "	" << dataset->Filenames()[i]<< endl;
	cout << "Number of events: " << eventTree->GetEntries() << endl;

	Luminosity_ = luminosity;
	Xsection_ = dataset->Xsection();
	PreSelEfficiency_ = dataset->PreSelEfficiency();
	MeanNEvents_ = Xsection_ * Luminosity_ * PreSelEfficiency_;
	TRandom3 randomgenerator(0);
	NEvents_ = int(randomgenerator.Poisson(MeanNEvents_));  //the number of events in a pseudoexperiment follows Poisson statistics
	if(verbosity>0){
		cout << "Luminosity_: " << Luminosity_ << " 1/pb" << ", Xsection_: " << Xsection_ << " pb" << endl;
		cout << "MeanNEvents_: " << MeanNEvents_ << ", NEvents_: " << NEvents_ << endl;
	}	
	GenerateRandomNumbers(NEvents_,eventTree->GetEntries());
	RankVector(RandomNumbers_);
}

void MCPseudoExp::CalculatePseudoExp(const Dataset* dataset, float& luminosity, int& ievtmin, bool nPseudoSeries){
      ////////////////////////////////////////
      // Calculate amount of events to be picked out of a Monte Carlo dataset, and generate a list of integers to be used later on to do the pseudoexperiment 
      ////////////////////////////////////////
        int verbosity = 0;	
	TChain* eventTree = dataset->eventTree();
	
	cout << "dataset name: " << dataset->Name() << endl;
	cout << "dataset filenames: " << endl;
	for(unsigned int i=0; i<dataset->Filenames().size(); i++)  cout  << "	" << dataset->Filenames()[i]<< endl;
	cout << "Number of events: " << eventTree->GetEntries() << endl;

	Luminosity_ = luminosity;
	Xsection_ = dataset->Xsection();
	PreSelEfficiency_ = dataset->PreSelEfficiency();
	MeanNEvents_ = Xsection_ * Luminosity_ * PreSelEfficiency_;
	TRandom3 randomgenerator(0);
	NEvents_ = int(randomgenerator.Poisson(MeanNEvents_));  //the number of events in a pseudoexperiment follows Poisson statistics
	if(verbosity>0){
		cout << "Luminosity_: " << Luminosity_ << " 1/pb" << ", Xsection_: " << Xsection_ << " pb" << endl;
		cout << "MeanNEvents_: " << MeanNEvents_ << ", NEvents_: " << NEvents_ << endl;
	}
	cout<<ievtmin<<" "<<ievtmin+NEvents_<<" "<<dataset->eventTree()->GetEntries()<<endl;
	if(ievtmin>=0 && (ievtmin+NEvents_) < dataset->eventTree()->GetEntries() && !nPseudoSeries){//&& dataset->Name()!="TTJets"){
		RandomNumbers_.clear();
		for(int i=ievtmin;i<ievtmin+NEvents_;i++) RandomNumbers_.push_back(i); 
		ievtmin = ievtmin+NEvents_;
		cout<<"Take continuous events"<<endl;
	}
	else{
              
		cout<<"Take random events as  conditions are "<<nPseudoSeries<<endl;
		GenerateRandomNumbers(NEvents_,eventTree->GetEntries());
		RankVector(RandomNumbers_);
		if(ievtmin<=dataset->eventTree()->GetEntries()) ievtmin = dataset->eventTree()->GetEntries()+1;
	    
          }
}

void MCPseudoExp::CalculatePseudoExp(TTree* seleventtree, float& force_luminosity, const Dataset* dataset, AnalysisEnvironment anaEnv){
      ////////////////////////////////////////
      // Calculate amount of events to be picked out of a TTree of (selected) events, and generate a list of integers to be used later on to do the pseudoexperiment 
      ////////////////////////////////////////  
	int verbosity = 1;
	
	int neventsTTree = seleventtree->GetEntries();
	Luminosity_ = force_luminosity; //it is this lumi that will be used below, not necessarily the one in anaEnv !!

	if(anaEnv.Systematics == 0 || anaEnv.Systematics == 1){		
		MeanNEvents_ = neventsTTree * Luminosity_ / dataset->EquivalentLumi(); //this is actually the same as Xsection_ * Luminosity_ * PreSelEfficiency_ * SelEfficiency, which we need in this case
	}
	else if(anaEnv.Systematics == 2){
	        if(verbosity>0) cout<<"Smearing cross section..."<<endl;
		TRandom3 randomgeneratorGaus(0);	
		Xsection_ = randomgeneratorGaus.Gaus(dataset->Xsection(),dataset->XsectionError() * dataset->Xsection()); //smearing cross section
		MeanNEvents_ = neventsTTree * Luminosity_ * Xsection_ / (dataset->EquivalentLumi() * dataset->Xsection());
	}
	else cout << "WARNING!!! Check config file for 'Systematics' option, should be 0, 1 or 2!!!" << endl;


	TRandom3 randomgenerator(0);
	NEvents_ = int(randomgenerator.Poisson(MeanNEvents_));  //the number of events in a pseudoexperiment follows Poisson statistics
	if(verbosity>1){
		cout << "Luminosity_: " << Luminosity_ << " 1/pb" << endl;
		cout << "MeanNEvents_: " << MeanNEvents_ << ", NEvents_: " << NEvents_ << endl;
	}	
	GenerateRandomNumbers(NEvents_,neventsTTree);
	RankVector(RandomNumbers_);
}

void MCPseudoExp::RankVector(vector<int> &randomnumbers){
      // rank the vector from low (first element) to high
	rankstruct myComparison;
	sort(randomnumbers.begin(),randomnumbers.end(),myComparison);
}
