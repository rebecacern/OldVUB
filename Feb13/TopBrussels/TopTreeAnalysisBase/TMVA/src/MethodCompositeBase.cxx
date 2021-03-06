// @(#)root/tmva $Id: MethodCompositeBase.cxx,v 1.1.2.1 2012/12/13 14:49:07 fblekman Exp $   
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss,Or Cohen 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodCompositeBase                                                   *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Virtual base class for all MVA method                                     *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *      Or Cohen        <orcohenor@gmail.com>    - Weizmann Inst., Israel         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      U. of Victoria, Canada                                                    * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

//_______________________________________________________________________
//
// This class is virtual class meant to combine more than one classifier//
// together. The training of the classifiers is done by classes that are//
// derived from this one, while the saving and loading of weights file  //
// and the evaluation is done here.                                     //
//_______________________________________________________________________

#include <algorithm>
#include <iomanip>
#include <vector>

#include "Riostream.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TObjString.h"

#include "TMVA/MethodCompositeBase.h"
#include "TMVA/MethodBoost.h"
#include "TMVA/MethodBase.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/Factory.h"
#include "TMVA/ClassifierFactory.h"

using std::vector;

ClassImp(TMVA::MethodCompositeBase)

//_______________________________________________________________________
TMVA::MethodCompositeBase::MethodCompositeBase( const TString& jobName, 
                                                Types::EMVA methodType,
                                                const TString& methodTitle,
                                                DataSetInfo& theData,
                                                const TString& theOption,
                                                TDirectory* theTargetDir )
   : TMVA::MethodBase( jobName, methodType, methodTitle, theData, theOption, theTargetDir )
{}

//_______________________________________________________________________
TMVA::MethodCompositeBase::MethodCompositeBase( Types::EMVA methodType,
                                                DataSetInfo& dsi,
                                                const TString& weightFile, 
                                                TDirectory* theTargetDir )
   : TMVA::MethodBase( methodType, dsi, weightFile, theTargetDir )
{}

//_______________________________________________________________________
TMVA::IMethod* TMVA::MethodCompositeBase::GetMethod( const TString &methodTitle ) const
{
   // returns pointer to MVA that corresponds to given method title
   vector<IMethod*>::const_iterator itrMethod    = fMethods.begin();
   vector<IMethod*>::const_iterator itrMethodEnd = fMethods.end();

   for (; itrMethod != itrMethodEnd; itrMethod++) {
      MethodBase* mva = dynamic_cast<MethodBase*>(*itrMethod);    
      if ( (mva->GetMethodName())==methodTitle ) return mva;
   }
   return 0;
}

//_______________________________________________________________________
TMVA::IMethod* TMVA::MethodCompositeBase::GetMethod( const Int_t index ) const
{
   // returns pointer to MVA that corresponds to given method index
   vector<IMethod*>::const_iterator itrMethod = fMethods.begin()+index;
   if (itrMethod<fMethods.end()) return *itrMethod;
   else                          return 0;
}


//_______________________________________________________________________
void TMVA::MethodCompositeBase::AddWeightsXMLTo( void* parent ) const 
{
   void* wght = gTools().xmlengine().NewChild(parent, 0, "Weights");
   gTools().AddAttr( wght, "NMethods",   fMethods.size()   );
   for (UInt_t i=0; i< fMethods.size(); i++) 
   {
      void* methxml = gTools().xmlengine().NewChild( wght, 0, "Method" );
      MethodBase* method = dynamic_cast<MethodBase*>(fMethods[i]);
      gTools().AddAttr(methxml,"Index",          i ); 
      gTools().AddAttr(methxml,"Weight",         fMethodWeight[i]); 
      gTools().AddAttr(methxml,"MethodSigCut",   method->GetSignalReferenceCut());
      gTools().AddAttr(methxml,"MethodTypeName", method->GetMethodTypeName());
      gTools().AddAttr(methxml,"MethodName",     method->GetMethodName()   ); 
      gTools().AddAttr(methxml,"JobName",        method->GetJobName());
      gTools().AddAttr(methxml,"Options",        method->GetOptions()); 
      method->AddWeightsXMLTo(methxml);
   }
}

//_______________________________________________________________________
TMVA::MethodCompositeBase::~MethodCompositeBase( void )
{
   // delete methods
   vector<IMethod*>::iterator itrMethod = fMethods.begin();
   for (; itrMethod != fMethods.end(); itrMethod++) {
      Log() << kVERBOSE << "Delete method: " << (*itrMethod)->GetName() << Endl;    
      delete (*itrMethod);
   }
   fMethods.clear();
}

//_______________________________________________________________________
void TMVA::MethodCompositeBase::ReadWeightsFromXML( void* wghtnode ) 
{
   // XML streamer
   UInt_t nMethods;
   TString methodName, methodTypeName, jobName, optionString;

   for (UInt_t i=0;i<fMethods.size();i++) delete fMethods[i];
   fMethods.clear();
   fMethodWeight.clear();
   gTools().ReadAttr( wghtnode, "NMethods",  nMethods );
   void* ch = gTools().xmlengine().GetChild(wghtnode);
   for (UInt_t i=0; i< nMethods; i++) {
      Double_t methodWeight, methodSigCut;
      gTools().ReadAttr( ch, "Weight",   methodWeight   );
      gTools().ReadAttr( ch, "MethodSigCut", methodSigCut);
      fMethodWeight.push_back(methodWeight);
      gTools().ReadAttr( ch, "MethodTypeName",  methodTypeName );
      gTools().ReadAttr( ch, "MethodName",  methodName );
      gTools().ReadAttr( ch, "JobName",  jobName );
      gTools().ReadAttr( ch, "Options",  optionString );

      if (i==0){
         // the cast on MethodBoost is ugly, but a similar line is also in ReadWeightsFromFile --> needs to be fixed later
         ((TMVA::MethodBoost*)this)->BookMethod( Types::Instance().GetMethodType( methodTypeName), methodName,  optionString );
      }
      fMethods.push_back(ClassifierFactory::Instance().Create(
         std::string(methodTypeName),jobName, methodName,DataInfo(),optionString));

      fMethodWeight.push_back(methodWeight);
      MethodBase* meth = dynamic_cast<MethodBase*>(fMethods.back());

      void* methXML = gTools().xmlengine().GetChild(ch);
      meth->SetupMethod();
      meth->ReadWeightsFromXML(methXML);
      meth->SetMsgType(kWARNING);
      meth->ParseOptions();
      meth->ProcessSetup();
      meth->CheckSetup();
      meth->SetSignalReferenceCut(methodSigCut);

      ch = gTools().xmlengine().GetNext(ch);
   }
   //Log() << kINFO << "Reading methods from XML done " << Endl;
}

//_______________________________________________________________________
void  TMVA::MethodCompositeBase::ReadWeightsFromStream( istream& istr )
{
   // text streamer
   TString var, dummy;
   TString methodName, methodTitle=GetMethodName(),
    jobName=GetJobName(),optionString=GetOptions();
   UInt_t methodNum; Double_t methodWeight;
   // and read the Weights (BDT coefficients)  
   istr >> dummy >> methodNum;
   Log() << kINFO << "Read " << methodNum << " Classifiers" << Endl;
   for (UInt_t i=0;i<fMethods.size();i++) delete fMethods[i];
   fMethods.clear();
   fMethodWeight.clear();
   for (UInt_t i=0; i<methodNum; i++) {
      istr >> dummy >> methodName >>  dummy >> fMethodIndex >> dummy >> methodWeight;
      if ((UInt_t)fMethodIndex != i) {
         Log() << kFATAL << "Error while reading weight file; mismatch MethodIndex=" 
               << fMethodIndex << " i=" << i 
               << " MethodName " << methodName
               << " dummy " << dummy
               << " MethodWeight= " << methodWeight 
               << Endl;
      }
      if (GetMethodType() != Types::kBoost || i==0) {
         istr >> dummy >> jobName;
         istr >> dummy >> methodTitle;
         istr >> dummy >> optionString;
         if (GetMethodType() == Types::kBoost)
            ((TMVA::MethodBoost*)this)->BookMethod( Types::Instance().GetMethodType( methodName), methodTitle,  optionString );
      }
      else methodTitle=Form("%s (%04i)",GetMethodName().Data(),fMethodIndex);
      fMethods.push_back(ClassifierFactory::Instance().Create( std::string(methodName), jobName, 
                                                               methodTitle,DataInfo(), optionString) );
      fMethodWeight.push_back( methodWeight );
      dynamic_cast<MethodBase*>(fMethods.back())->ReadWeightsFromStream(istr);
   }
}

//_______________________________________________________________________
Double_t TMVA::MethodCompositeBase::GetMvaValue( Double_t* err )
{
   // return composite MVA response
   Double_t mvaValue = 0;
   for (UInt_t i=0;i< fMethods.size(); i++) mvaValue+=fMethods[i]->GetMvaValue()*fMethodWeight[i];

   // cannot determine error
   if (err != 0) *err = -1;

   return mvaValue;
}
