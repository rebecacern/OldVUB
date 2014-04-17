#include <iostream>
#include "TSystem.h"
//#include "freyaLooper.h"


void chainFreya(int nsel, int mode = 0, bool silent = false){  

  gSystem->CompileMacro("freyaLooper.C","k");

  double SFval = 0.95;
  bool SFplus = false;
  bool SFminus = false;
  // samples used
  double x_sec = 0.;
  char plotName[300];
  sprintf(plotName,"test");
  
  if (nsel == 0)                	{sprintf(plotName,"tt");}
  else if (nsel == 1)   		{sprintf(plotName,"twdr");}
  else if (nsel == -1)   		{sprintf(plotName,"twds");}
  else if (nsel == 2)			{sprintf(plotName, "st");}
  else if (nsel == 3)   		{sprintf(plotName,"wjets");}
  else if (nsel == 4)   		{sprintf(plotName,"zjets");}
  else if (nsel == 5)   		{sprintf(plotName,"dy");}
  else if (nsel == 6)   		{sprintf(plotName,"qcd_mu");}
  else if (nsel == 7)   		{sprintf(plotName,"di");}
  else if (nsel == 77)                	{sprintf(plotName,"non_tt");}
  else if (nsel == 66)                	{sprintf(plotName,"others");}
  else if (nsel == 666)                	{sprintf(plotName,"data");}
  else if (nsel == 555)                	{sprintf(plotName,"mc");}
  else if (nsel == 777)                 {sprintf(plotName,"all");}
  
  // if (mode != 0 &&  mode !=1 && mode !=2) mode = 0;
  if (!silent){
    if      (mode == 0) 	cout << " Electron-Muon Mixed channel " << endl;
    else if (mode == 1) 	cout << " Di-Muon channel " << endl;
    else if (mode == 2) 	cout << " Di-Electron channel " << endl;
    cout << "*************************************************" << endl;
  }
  
  bool nosf = true;
  char myRootFile[300];
  sprintf(myRootFile,"outputs/EleStudies_%d_%s.root", mode, plotName);
  
  TChain *myCh = new TChain("myTree","myTree");
  myCh->Add(myRootFile);

  freyaLooper an(myCh);
  std::cout << "*************************************************" << std::endl;
  std::cout << "starting loop! " << std::endl;
  std::cout << "chain: " << myCh->GetEntries() << " entries." << std::endl;
  std::cout << "*************************************************" << std::endl;
  an.myLoop(nsel,mode,silent);
  
}
