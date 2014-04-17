/// Made faster Freya-style
#include <iostream>
#include "TSystem.h"


void dis_chain(int nsel = 0, int mode = 0, bool silent = false){  

  gSystem->CompileMacro("Discriminator.C","k");

 
  // samples used
  double x_sec = 0.;
  char plotName[300];
  sprintf(plotName,"test");
  
  if (nsel == 0)                	{sprintf(plotName,"tt");}
  else if (nsel == 1)   		{sprintf(plotName,"twdr");}
 
  
  if (mode != 0 &&  mode !=1 && mode !=2) mode = 0;
  if (!silent){
    cout << "[Info:]" ;
    if      (mode == 0) 	cout << " emu channel, " ;
    else if (mode == 1) 	cout << " mumu channel, " ;
    else if (mode == 2) 	cout << " ee channel, " ;
  }
  
  char myRootFile[300];
  sprintf(myRootFile,"outputs/out_%d_%s.root", mode, plotName);
  
  TChain *myCh = new TChain("myTree","myTree");
  myCh->Add(myRootFile);

  Discriminator an(myCh);
  std::cout << plotName << " (" << myRootFile << ",  " << myCh->GetEntries() << " entries )" << std::endl;
  an.myLoop(nsel,mode,silent);
  
}
