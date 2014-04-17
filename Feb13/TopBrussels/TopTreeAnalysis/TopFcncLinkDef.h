#ifdef __CINT__

#include "TopFCNC/interface/TopFCNC_Evt.h"
//#include "TopFCNC/interface/TopFCNC_KinFit.h"
#include "Content/interface/AnalysisEnvironment.h"
#include "Content/interface/Dataset.h"
//#include "KinFitter/interface/TKinFitter.h"
//#include "Reconstruction/interface/Observables.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeProducer/interface/TRootGenEvent.h"
#include "../TopTreeProducer/interface/TRootNPGenEvent.h"
#include "../TopTreeProducer/interface/TRootHLTInfo.h"
#include "../TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootVertex.h"
#include "../TopTreeProducer/interface/TRootElectron.h"
#include "../TopTreeProducer/interface/TRootMuon.h"
#include "../TopTreeProducer/interface/TRootGenJet.h"
#include "../TopTreeProducer/interface/TRootJet.h"
#include "../TopTreeProducer/interface/TRootPFJet.h"
#include "../TopTreeProducer/interface/TRootMET.h"
#include "../TopTreeProducer/interface/TRootPFMET.h"
#include "../TopTreeProducer/interface/TRootTrackMET.h"
#include "../TopTreeProducer/interface/TRootParticle.h"
#include "../TopTreeProducer/interface/TRootMCParticle.h"
#include "../TopTreeProducer/interface/TRootGenTop.h"

#else

#include "TopFCNC/interface/TopFCNC_Evt.h"
//#include "TopFCNC/interface/TopFCNC_KinFit.h"
#include "Content/interface/AnalysisEnvironment.h"
#include "Content/interface/Dataset.h"
//#include "KinFitter/interface/TKinFitter.h"
//#include "Reconstruction/interface/Observables.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeProducer/interface/TRootGenEvent.h"
#include "../TopTreeProducer/interface/TRootNPGenEvent.h"
#include "../TopTreeProducer/interface/TRootHLTInfo.h"
#include "../TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootVertex.h"
#include "../TopTreeProducer/interface/TRootElectron.h"
#include "../TopTreeProducer/interface/TRootMuon.h"
#include "../TopTreeProducer/interface/TRootGenJet.h"
#include "../TopTreeProducer/interface/TRootJet.h"
#include "../TopTreeProducer/interface/TRootPFJet.h"
#include "../TopTreeProducer/interface/TRootMET.h"
#include "../TopTreeProducer/interface/TRootPFMET.h"
#include "../TopTreeProducer/interface/TRootTrackMET.h"
#include "../TopTreeProducer/interface/TRootParticle.h"
#include "../TopTreeProducer/interface/TRootMCParticle.h"
#include "../TopTreeProducer/interface/TRootGenTop.h"

#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class TopFCNC_Evt;
//#pragma link C++ class TopFCNC_KinFit;
#pragma link C++ class Dataset;
//#pragma link C++ class TKinFitter;
#pragma link C++ class TRootEvent;
#pragma link C++ class TRootGenEvent;
#pragma link C++ class TRootNPGenEvent;
#pragma link C++ class TRootHLTInfo;
#pragma link C++ class TRootRun;
#pragma link C++ class TRootVertex;
#pragma link C++ class TRootElectron;
#pragma link C++ class TRootMuon;
#pragma link C++ class TRootGenJet;
#pragma link C++ class TRootJet;
#pragma link C++ class TRootPFJet;
#pragma link C++ class TRootMET;
#pragma link C++ class TRootPFMET;
#pragma link C++ class TRootTrackMET;
#pragma link C++ class TRootParticle;
#pragma link C++ class TRootMCParticle;
#pragma link C++ class TRootGenTop;
#pragma link C++ class AnalysisEnvironment;
//#pragma link C++ class Observables;

#pragma link C++ struct TopTree::triggeredObject;

#endif
