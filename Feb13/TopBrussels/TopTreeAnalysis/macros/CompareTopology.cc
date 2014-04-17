// This is a macro to superimpose all equivalent histograms (i.e. with the same name) from multiple root files. The macro is based on hadd.C (a macro to add histogram files, From $ROOTSYS/tutorials/io/hadd.C)

#include <string.h>
#include <cmath>
#include <sstream>
#include "iomanip.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

using namespace std;

TList *FileList;
TFile *Target;

void SuperimposeRootfile( TDirectory *target, TList *sourcelist );


void ModifyHist(TH1F* &h, Color_t lcolor);


int main(int argc, char*argv[]){
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  
  Target = TFile::Open( "comparison.root", "RECREATE" );
  
  FileList = new TList();  // more than 2 source files: not tested
  FileList->Add( TFile::Open("SUSY_LM0_output.root") );
  FileList->Add( TFile::Open("Ttjets_output.root") );	
  
  SuperimposeRootfile( Target, FileList );

}   


void SuperimposeRootfile( TDirectory *target, TList *sourcelist ){ 
  cout << "	" << "========================================================" << endl; 
  cout << "	" << "This is a macro to superimpose plots of different root files." << endl;
  cout << "	" << "Only TH1F objects are superimposed." << endl;
   cout << "	" << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd( path );
  TDirectory *current_sourcedir = gDirectory;
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // loop over all keys in this directory
  TChain *globChain = 0;
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    // read object from first source file and create a canvas
    first_source->cd( path );
    TObject *obj = key->ReadObj();
    TCanvas *c1 = new TCanvas("c1",obj->GetName(),500,500);
    TLegend* legend_c1 = new TLegend(0.65,0.80,0.89,0.70);

    if ( obj->IsA()->InheritsFrom( "TH1F" ) ) {
      // descendant of TH1F -> prepare the histograms to be superimposed

      //      cout << "Modifying histogram " << obj->GetName() << endl;      
      TH1F *h1 = (TH1F*)obj;
      ModifyHist(h1,kRed);

      // loop over all source files and modify the correspondant
      // histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource ) {
        
        // make sure we are at the correct directory level by cd'ing to path
        nextsource->cd( path );
        TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
        if (key2) {
           TH1F *h2 = (TH1F*)key2->ReadObj();
	   ModifyHist(h2,kBlue);
           
	   double maxh1;
           double maxh2;
      	   maxh1 = h1->GetMaximum(10000.);
           maxh2= h2->GetMaximum(10000.);
           if(maxh1>maxh2){
   		h1->Draw();
   		h2->Draw("SAME");
           }
           else {
   		h2->Draw();
   		h1->Draw("SAME");
           }
	   legend_c1->SetTextFont(70);
   	   legend_c1->SetTextSize(0.03);
   	   legend_c1->AddEntry(h1,"susy","L");  // change this according to first source file
   	   legend_c1->AddEntry(h2,"ttbar semileptonic","L");    // change this according to second source file
   	   legend_c1->Draw("SAME");
        }

        nextsource = (TFile*)sourcelist->After( nextsource );
      }  // while ( nextsource )
    }
    else if ( obj->IsA()->InheritsFrom( "TTree" ) ) {	// not tested
      
      // loop over all source files create a chain of Trees "globChain"
      const char* obj_name= obj->GetName();

      globChain = new TChain(obj_name);
      globChain->Add(first_source->GetName());
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      //      const char* file_name = nextsource->GetName();
      // cout << "file name  " << file_name << endl;
     while ( nextsource ) {
     	  
       globChain->Add(nextsource->GetName());
       nextsource = (TFile*)sourcelist->After( nextsource );
     }

    } else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {  // not tested
      // it's a subdirectory

      cout << "Found subdirectory " << obj->GetName() << endl;

      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

      // newdir is now the starting point of another round of superimposing
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      SuperimposeRootfile( newdir, sourcelist );

    }     
    else {

      // object is of no type that we know or can handle
      cout << "Unknown object type, name: " 
           << obj->GetName() << ", object type: " << obj->ClassName() << endl;
    }

    // now draw and write the superimposed histograms to the target file
    // note that this will just store the canvas c1 in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if ( obj ) {
      target->cd();
           
      //!!if the object is a tree, it is stored in globChain...     
	if(obj->IsA()->InheritsFrom( "TTree" ))  // not tested
          globChain->Merge(target->GetFile(),0,"keep");
	else 
	  c1->Write( key->GetName() );
    }

  } // while ( ( TKey *key = (TKey*)nextkey() ) )

  // save modifications to target file
  target->SaveSelf(kTRUE);
  TH1F::AddDirectory(status);
  cout << "	" << "========================================================" << endl;
}


void ModifyHist (TH1F* &h, Color_t lcolor)
{
	double temp_integral;

	h->SetLineColor(lcolor);
	temp_integral = h->Integral();
	//cout << temp_integral << endl;
	if(temp_integral!=0) h->Scale(pow(temp_integral,-1));
}









/*
int main(int argc, char*argv[]){
   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);
   
   cout << "	" << "Comparison of plots" << endl;
   
   /*
   int numberofhistos;
   double maxh1;
   double maxh2;
   stringstream canvasname_stream;
   string canvasname("");
   const char* canvasname_char;	
  
   TFile s("SUSY_LM0_output.root");
   TList *shList = new TList();
   shList = s.GetList();
   numberofhistos = shList->GetEntries();
   //TH1F *susyhisto = new TH1F();
   double testing;
   testing = (shList->At(0))->GetBinContent(3);
   */
   
 /*  
   // get the number of histograms to loop over
   TFile s("SUSY_LM0_output.root");
   TObjArray *shArray = (TObjArray*) s.Get("histogramArray");
   numberofhistos = shArray->GetEntriesFast();
   s.Close();
   cout << numberofhistos << endl;
   
   // loop over the histograms to superimpose all equivalent histogram pairs
   for(unsigned int i=0;i<numberofhistos;i++){
   	TFile susy("SUSY_LM0_output.root");
	TObjArray *susyhArray = (TObjArray*) susy.Get("histogramArray");
   	TH1F *susyhisto = new TH1F();
	susyhisto = susyhArray->At(i);	// gives error... apparently susyhArray->At(i) is not a TH1F ???; tobefixed
   	ModifyHist(susyhisto,kRed);
   	maxh1 = susyhisto->GetMaximum(10000.);
   	TFile ttbar("Ttjets_output.root");
	TObjArray *ttbarhArray = (TObjArray*) ttbar.Get("histogramArray");
   	TH1F *ttbarhisto = ttbarhArray->At(i);
   	ModifyHist(ttbarhisto,kBlue);
   	maxh2= ttbarhisto->GetMaximum(10000.);
   
   	canvasname_char = susyhArray[i]->GetName();
   	canvasname = susyhArray[i]->GetName();	// gives error...; tobefixed
   	canvasname_stream << "c_" << canvasname;
   	canvasname = canvasname_stream.str();
   
   	TCanvas *c1 = new TCanvas("c1",canvasname.c_str(),500,500);
   	if(maxh1>maxh2){
   		susyhisto->Draw();
   		ttbarhisto->Draw("SAME");
   	}
   	else {
   		ttbarhisto->Draw();
   		susyhisto->Draw("SAME");
   	}
   
   	TLegend* legend_c1 = new TLegend(0.45,0.85,0.75,0.75);
   	legend_c1->SetTextFont(70);
   	legend_c1->SetTextSize(0.03);
   	legend_c1->AddEntry(susyhisto,"susy","L");
   	legend_c1->AddEntry(ttbarhisto,"ttbar semileptonic","L");
   	legend_c1->Draw("SAME");
   
   	TFile rootFile(outrootFileName.c_str());
   	c1->Write();
   	//rootFile.Close();	 
   }
   
   string outrootFileName("comparison.root");
   TFile* fout = new TFile(outrootFileName.c_str(),"RECREATE");
   cout << "	" << "Output written in "<< outrootFileName << endl;

} 
*/
