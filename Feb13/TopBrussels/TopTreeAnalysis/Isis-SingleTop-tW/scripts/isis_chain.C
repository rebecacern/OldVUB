/// Made faster Freya-style
#include <iostream>
#include "TSystem.h"


void isis_chain(int nsel = 0, int mode = 0, bool silent = false){  

  gSystem->CompileMacro("isis_looper.C","k");

 
  // samples used
  double x_sec = 0.;
  char plotName[300];
  sprintf(plotName,"test");
  
  if (nsel == 0)                	{sprintf(plotName,"tt");}
  else if (nsel == 1)   		{sprintf(plotName,"twdr");}
  else if (nsel == -1)   		{sprintf(plotName,"twds");}
  else if (nsel == 2)   		{sprintf(plotName,"zjetsall");}
  else if (nsel == 3)   		{sprintf(plotName,"di");}
  else if (nsel == 4)			{sprintf(plotName, "st");}
  else if (nsel == 5)   		{sprintf(plotName,"wjets");}
  else if (nsel == 6)   		{sprintf(plotName,"qcd_mu");}
  else if (nsel == 7)                	{sprintf(plotName,"others");}

  else if (nsel == 555)                	{sprintf(plotName,"mc");}
  
  else if (nsel == 666)                	{sprintf(plotName,"data");}
  
    
  else if (nsel == -10)                   {sprintf(plotName,"tt");}
  else if (nsel ==  10)                   {sprintf(plotName,"tt");}
  
  if (mode != 0 &&  mode !=1 && mode !=2) mode = 0;
  if (!silent){
    cout << "[Info:]" ;
    if      (mode == 0) 	cout << " emu channel, " ;
    else if (mode == 1) 	cout << " mumu channel, " ;
    else if (mode == 2) 	cout << " ee channel, " ;
  }
  
  char myRootFile[300];
  sprintf(myRootFile,"outputs/out_%d_%s.root", mode, plotName);
  
  if(nsel == -10){
  sprintf(myRootFile,"outputs/JERsysDown_%d_%s.root", mode, plotName);
  }else if(nsel == 10 ){
   sprintf(myRootFile,"outputs/JERsysUp_%d_%s.root", mode, plotName);
  }
  
  TChain *myCh = new TChain("myTree","myTree");
  myCh->Add(myRootFile);

  isis_looper an(myCh);
  std::cout << plotName << " (" << myRootFile << ",  " << myCh->GetEntries() << " entries )" << std::endl;
  an.myLoop(nsel,mode,silent);
  
}
