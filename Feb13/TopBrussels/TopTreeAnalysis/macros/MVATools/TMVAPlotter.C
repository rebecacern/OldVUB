#include <iostream>
#include <vector>

#include "TList.h"
#include "TROOT.h"
#include "TKey.h"
#include "TString.h"
#include "TControlBar.h"
#include "TObjString.h"

#include "tmvaglob.C"


// main GUI
void TMVAPlotter( string fName = "TMVA.root" ) 
{   
   // Use this script in order to run the various individual macros
   // that plot the output of TMVA (e.g. running TMVAClassification.C),
   // stored in the file "TMVA.root"

   TString curMacroPath(gROOT->GetMacroPath());
   // uncomment next line for macros submitted to next root version
   gROOT->SetMacroPath(curMacroPath+":../TMVA/macros/:.:MVATools/");
   
   // for the sourceforge version, including $ROOTSYS/tmva/test in the
   // macro path is a mistake, especially if "./" was not part of path
   // add ../macros to the path (comment out next line for the ROOT version of TMVA)
   // gROOT->SetMacroPath(curMacroPath+":../macros:");
   
   cout << "--- Launch TMVA Plotter to save monitoring plots from: " << fName << endl;

   // s+(s+b)

   cout << "--- Plotting s/(s+b) plots" << endl;

   string cmd = ".x soversb.C(\""+fName+"\",\"InputVariables_Id\",\"Input variables\")";
   string cmd2 = ".x soversb.C(\""+fName+"\",\"InputVariables_Deco\",\"Input variables (Decorrelated)\")";

   gROOT->ProcessLine(cmd.c_str());
   gROOT->ProcessLine(cmd2.c_str());

   // input vars

   cout << "--- Plotting input variables" << endl;

   string cmd = ".x variables.C(\""+fName+"\",\"InputVariables_Id\",\"Input variables\")";
   string cmd2 = ".x variables.C(\""+fName+"\",\"InputVariables_Deco\",\"Input variables (Decorrelated)\")";

   gROOT->ProcessLine(cmd.c_str());
   gROOT->ProcessLine(cmd2.c_str());

   // correlations

   cout << "--- Plotting linear correlation plots" << endl;

   cmd = ".x correlations.C(\""+fName+"\")";

   gROOT->ProcessLine(cmd.c_str());

   // MVA output

   cout << "--- Plotting MVA output" << endl;

   cmd = ".x mvas.C(\""+fName+"\",0)";
   cmd2 = ".x mvas.C(\""+fName+"\",3)";

   gROOT->ProcessLine(cmd.c_str());
   gROOT->ProcessLine(cmd2.c_str());

   cout << "--- Done" << endl;

   gApplication->Terminate();
}
