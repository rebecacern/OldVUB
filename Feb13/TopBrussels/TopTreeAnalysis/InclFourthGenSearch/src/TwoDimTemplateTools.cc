#include "../interface/TwoDimTemplateTools.h"

TwoDimTemplateTools::TwoDimTemplateTools()
{  
   nfixedbins_ = 80;
   xvariable_XaxisLabel_ = "Reconstructed top mass (GeV/c^{2})";
   xvariable_YaxisLabel_ = "#Events";
   xbins = NULL;   
   binningFileName_ = "";
}

TwoDimTemplateTools::TwoDimTemplateTools(string identifier, string xvariable, int nbinsxvariable, string yvariable, int nbinsyvariable)
{  
   identifier_ = identifier;
   nfixedbins_ = 80;
   xvariable_XaxisLabel_ = "Reconstructed top mass (GeV/c^{2})";
   xvariable_YaxisLabel_ = "#Events";
   xbins = NULL;
   binningFileName_ = ""; 
   xvariable_ = xvariable;
   yvariable_ = yvariable;
   nbinsxvariable_ = nbinsxvariable;
   nbinsyvariable_ = nbinsyvariable;     
   nbinsMap[xvariable_] = nbinsxvariable_;
   nbinsMap[yvariable_] = nbinsyvariable_;   
   historanges[xvariable_] = make_pair(0,3000);
   historanges[yvariable_] = make_pair(0,3000);
}

TwoDimTemplateTools::TwoDimTemplateTools(const TwoDimTemplateTools & t)
{
   xvariable_ = t.xvariable_;
   yvariable_ = t.yvariable_;
   nbinsxvariable_ = t.nbinsxvariable_;
   nbinsyvariable_ = t.nbinsyvariable_;
   xvariable_XaxisLabel_ = t.xvariable_XaxisLabel_;
   xvariable_YaxisLabel_ = t.xvariable_YaxisLabel_;
   nfixedbins_ = t.nfixedbins_;   
   xbins = t.xbins;   
   binningFileName_ = t.binningFileName_;   
   MSPlot_xvariableBins = t.MSPlot_xvariableBins;
   Histos_xvariableBins = t.Histos_xvariableBins;
   VariableValuesMap = t.VariableValuesMap;
   eventweightvector = t.eventweightvector;
   nbinsMap = t.nbinsMap;   
   historanges = t.historanges;      
   datasets_ = t.datasets_;
}

TwoDimTemplateTools::~TwoDimTemplateTools()
{  
}

void TwoDimTemplateTools::SetDatasets(vector < Dataset* > datasets)
{
   for(unsigned int d=0;d<datasets.size();d++)
   {
      datasets_.push_back(datasets[d]);
   }
}

void TwoDimTemplateTools::LoadTwoDimBinning(const string binningFileName)
{
   unsigned int xarraysize = nbinsxvariable_ + 2;
   xbins = NULL;
   
   //Binning for xvariable (originally "HT")
   TAxis *xaxis = NULL;
   char txAxisName[150];
   sprintf (txAxisName, "Binning_%s_SM", xvariable_.c_str ());
   
   TFile fBinning (binningFileName.c_str(), "READ"); //
   //cout<<" - Reading binning for variable "<<xvariable_<<endl;
   fBinning.GetObject (txAxisName, xaxis);
   fBinning.Close();
   xbins = xaxis->GetXbins ()->fArray; //this is the "HT" binning
   //cout<<""<<endl;
   for (unsigned int b = 0; b < xarraysize; b++)
   {        
	//cout<<" xbins["<<b<<"] = "<<xbins[b]<<endl;
   }
   
   //Binning for yvariable (originally "MTop")
   TAxis* axis = NULL;
   char tAxisName[150], histoName[150];
   for(int k = 1;k < nbinsxvariable_+1; k++)//k < nbinsxvariable_+1
   {    
     axis = NULL;
     TFile fBinning_bis (binningFileName.c_str(), "READ"); //
     
     sprintf (tAxisName, "Binning_%s_SM_%sbin%i", yvariable_.c_str (),xvariable_.c_str (),k);
     cout<<" - Getting "<<tAxisName<<"..."<<endl;
     fBinning_bis.GetObject (tAxisName, axis);     
     fBinning_bis.Close(); //Crucial! Otherwise histograms created in this file...    
     
     vector < pair < TH1F*,Dataset* > > Histos_datasets; // the histos are Mtop distributions per dataset
     cout << " - Creating histos..." << endl;
     if (axis == NULL)
     {
        cout<<"Warning: axis == NULL"<<endl;
	//cout<<"nbins = "<<nbins<<", range obs: "<<historanges[i].first<<" -> "<<historanges[i].second<<", variable "<<lstVar[i]<<endl;	
	for (unsigned int d = 0; d < datasets_.size (); d++)
	{  
	   sprintf (histoName, "h%s_%sbin%i_%s", yvariable_.c_str (),xvariable_.c_str (),k,(datasets_[d]->Name()).c_str());
	   Histos_datasets.push_back(make_pair (new TH1F (TString(histoName), yvariable_.c_str(), nfixedbins_, historanges[yvariable_].first, historanges[yvariable_].second),datasets_[d])); //to be fixed?
           (Histos_datasets.back()).first->GetXaxis()->SetTitle(xvariable_XaxisLabel_.c_str());
	   (Histos_datasets.back()).first->GetYaxis()->SetTitle(xvariable_YaxisLabel_.c_str());
	}
	Histos_xvariableBins[k] = Histos_datasets;
     }
     else
     {	  
	for (unsigned int d = 0; d < datasets_.size (); d++)
	{
	   sprintf (histoName, "h%s_%sbin%i_%s", yvariable_.c_str (),xvariable_.c_str (),k,(datasets_[d]->Name()).c_str());
	   //cout<<"histoName = "<<histoName<<endl;
	   Histos_datasets.push_back(make_pair (new TH1F (TString(histoName), yvariable_.c_str(), axis->GetNbins (),axis->GetXbins ()->fArray),datasets_[d]));
	   (Histos_datasets.back()).first->GetXaxis()->SetTitle(xvariable_XaxisLabel_.c_str());
	   (Histos_datasets.back()).first->GetYaxis()->SetTitle(xvariable_YaxisLabel_.c_str());	   
	   //cout<<" name = "<<(Histos_datasets[d].first)->GetName()<<endl;	
	}
	Histos_xvariableBins[k] = Histos_datasets;	
     } 
   } //end loop over "HT" bins

}

void TwoDimTemplateTools::Fill_for2DBinning(double xvariable_value,float yvariable_value,float fillweight)
{  
   VariableValuesMap[xvariable_].push_back(xvariable_value); //"HT"
   VariableValuesMap[yvariable_].push_back(yvariable_value); //"Mtop"	
   eventweightvector.push_back(fillweight);
}

void TwoDimTemplateTools::Fill(double xvariable_value,float yvariable_value,float fillweight,unsigned int datasetindex)
{ 
  for(int k = 1;k < nbinsxvariable_+1; k++)
  {    
    //cout<<"---> HT bin k = "<<k<<endl;
    double rightcriterion;
    if(k==nbinsxvariable_)
    {
       rightcriterion = xbins[k+1] + 1000000; //xbins is the "HT" binning read from a binning file
    }
    else
    {  
       rightcriterion = xbins[k+1];
    }       
    
    if(xvariable_value>=xbins[k] && xvariable_value<rightcriterion) //then the event is in the 1st (kth??) bin of the xvariable ("HT")
    {   
       (Histos_xvariableBins[k][datasetindex].first)->Fill(yvariable_value,fillweight);
       int N = (Histos_xvariableBins[k][datasetindex].first)->GetNbinsX();
       (Histos_xvariableBins[k][datasetindex].first)->SetBinContent(N,(Histos_xvariableBins[k][datasetindex].first)->GetBinContent(N)+(Histos_xvariableBins[k][datasetindex].first)->GetBinContent(N+1)); // same as in MultiSamplePlot
       (Histos_xvariableBins[k][datasetindex].first)->SetBinContent(N+1,0); // same as in MultiSamplePlot
       
    }
  }
}

void TwoDimTemplateTools::Write_for2DBinning(const string binningFileName)
{
   cout<<endl<<" ... Creating the binning"<<endl;	
   MakeBinning NewBins;
   NewBins.Binning_forTprimeAnalysis(xvariable_,yvariable_,VariableValuesMap,nbinsMap,eventweightvector,binningFileName);
}

void TwoDimTemplateTools::Convert2Dto1D(string postfix)
{ 
  cout<<" - Producing 2D to 1D converted templates"<<endl;
  vector<TH1F* > h2D_1Dconverted;//vector of dataset histograms, for each entry in the vector there is the histogram of a dataset which should be 1D converted from 2D
  //h2D_1Dconverted.clear();
  for(unsigned int d=0;d<datasets_.size();d++)
  {
     string datasetname = datasets_[d]->Name();
     string histoName = ""; //naming convention chosen the same as in macro for 1D histos (templates)...identifier
     //if(datasets_[d]->Name()=="Data" || datasets_[d]->Name()=="data" || datasets_[d]->Name()=="DATA")
       histoName = xvariable_ + "_" + yvariable_ + "_" + identifier_ + datasetname;
     //else 
     //  histoName = xvariable_ + "_" + yvariable_ + "_" + identifier_ + datasetname + postfix;
     TH1F* htemp = new TH1F(histoName.c_str(),histoName.c_str(),nbinsxvariable_*nbinsyvariable_,0,nbinsxvariable_*nbinsyvariable_);
     h2D_1Dconverted.push_back(htemp);
  }
 
  int b=1; //0 is underflowbin; k is the bin number of the 1D converted histogram
  int b_remember=1;
  int nbins_htemp = 0;
  
  for(int k = 1;k < nbinsxvariable_+1; k++)
  {
     for(unsigned int d=0;d<datasets_.size();d++)
     {
        //cout<<"d = "<<d<<", "<<datasets_[d]->Name()<<endl;
				TH1F* htemp =0;
				htemp = (TH1F*) Histos_xvariableBins[k][d].first;
				b=b_remember;
				//cout<<"b before loop over htemp bins = "<<b<<endl;
        nbins_htemp = htemp->GetNbinsX(); //should be nbinsyvariable_ + 1? (not sure... to be checked)
				//cout<<"nbins_htemp = "<<nbins_htemp<<endl;
				
				//cout<<"Integral (data) of top mass for HT bin "<<k<<": "<<htemp->Integral(0,htemp->GetNbinsX()+1)<<endl;		
				//for(int bini=0;bini<htemp->GetNbinsX()+2;bini++)
				//	cout<<" bini = "<<bini<<", content "<<htemp->GetBinContent(bini)<<endl;

			
				for(int j=2;j<nbins_htemp+1;j++)
				{
	    			float bincontentj = 0;
	    			bincontentj = htemp->GetBinContent(j);
	    			//cout<<" b = "<<b<<endl;
	    			//cout<<"      bincontent "<<j<<" = "<<bincontentj<<endl;
	    			h2D_1Dconverted[d]->SetBinContent(b,bincontentj);
	    			b++;
				}
				//cout<<"Updated integral of h2D_1Dconverted: "<<h2D_1Dconverted[d]->Integral(0,h2D_1Dconverted[d]->GetNbinsX()+1)<<endl;
     }
     
     //cout<<" b_remember before update = "<<b_remember<<endl;
     b_remember = b_remember + nbins_htemp - 1;
     //cout<<" b_remember after update = "<<b_remember<<endl;  
  }
  
  //fout->cd();
  
  for(unsigned int d=0;d<datasets_.size();d++)
  {
    string histoname;
    //cout<<" d = "<<d<<endl;
    h2D_1Dconverted[d]->GetYaxis()->SetTitle("#Events");
    templates1D.push_back(make_pair(h2D_1Dconverted[d],datasets_[d]));
    if(!(datasets_[d]->Name().find("NP") <= datasets_[d]->Name().size()) || (datasets_[d]->Name().find("NP_overlay") <= datasets_[d]->Name().size()))
    {
      histoname = h2D_1Dconverted[d]->GetName();
      //cout<<" histoname = "<<histoname<<endl;   
      templates1D_forMS.push_back(make_pair(h2D_1Dconverted[d],datasets_[d]));
    }

  }
  
}

void TwoDimTemplateTools::Write(TFile* fout,TDirectory* th1dir)
{
   //TDirectory* MTop_in_HTbins = fout->mkdir("MTop_in_HTbins");
   fout->cd();
   //MTop_in_HTbins->cd();
   for(int k = 1;k < nbinsxvariable_+1; k++)
   {  
	 			//cout<<"Integral (data) of top mass for HT bin "<<k<<": "<<((Histos_xvariableBins[k][0]).first)->Integral(0,((Histos_xvariableBins[k][0]).first)->GetNbinsX()+1)<<endl;
				MSPlot_xvariableBins[k] = new MultiSamplePlot(Histos_xvariableBins[k]);
        //for (unsigned int d = 0; d < datasets.size (); d++)
        //{
	//     (Histos_xvariableBins[k][d].first)->Write(); 
        //}	       
   }
   for(map<int,MultiSamplePlot*>::const_iterator it = MSPlot_xvariableBins.begin(); it != MSPlot_xvariableBins.end(); it++)
   {
        MultiSamplePlot *temp = it->second;
        int name_int = it->first;
				char HTbin[150];
				sprintf (HTbin, "%s_%s_%sbin%i",yvariable_.c_str (),identifier_.c_str (),xvariable_.c_str (),name_int);
        temp->Draw(false, HTbin, true, true, true, true, true, 5,false,true,true);//(bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal,bool addRatio, bool mergeVV, bool mergeTTV)
        temp->Write(fout, HTbin, false, "");
   }
   
   //2D to 1D converted
   fout->cd();
   th1dir->cd();
   cout<<" Writing the 2D to 1D converted template"<<endl; //not for MS, but the actual templates
   for(unsigned int d=0;d<datasets_.size();d++)
   {
        //cout<<"   -> Writing "<< templates1D[d]->GetName()<<"..."<<endl; 
        templates1D[d].first->Write();
   }
   fout->cd();
   cout<<" Writing the 2D to 1D converted template as MultiSampleplot"<<endl;
   MultiSamplePlot* templates1D_MS = new MultiSamplePlot(templates1D_forMS);
   char TemplateName[150];
   sprintf (TemplateName, "Template_%s_%s_%s",identifier_.c_str (),xvariable_.c_str (),yvariable_.c_str ());
   templates1D_MS->Draw(false, TemplateName, true, true, true, true, true,5,false,true,true);//(bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal,bool addRatio, bool mergeVV, bool mergeTTV)
	 templates1D_MS->Write(fout, TemplateName, false, "");
}


