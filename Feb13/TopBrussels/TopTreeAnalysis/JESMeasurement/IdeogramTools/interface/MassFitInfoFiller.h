#ifndef MassFitInfoFiller_h
#define MassFitInfoFiller_h

#include <iostream>

#include <TopQuarkAnalysis/TopMassIdeogram/interface/GenericTool.h>
#include <TopQuarkAnalysis/TopMassIdeogram/interface/MassFitInfoHeader.h>
#include <TopQuarkAnalysis/TopMassIdeogram/interface/MassFitInfoEvent.h>
#include <TopQuarkAnalysis/TopMassIdeogram/interface/MassFitInfo.h>
#include <TopQuarkAnalysis/TopMassIdeogram/interface/MassLikelihoodTool.h>

#include "TopTreeAnalysis/JESMeasurement/interface/LightMonster.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"

void fillMassFitInfoHeader(Dataset* dataSet, ideogram::MassFitInfoHeader& header, bool isSemiMu, bool containsPDFsyst);

void fillMassFitInfoEvent(LightMonster* monster, const ideogram::MassFitInfoHeader& header, ideogram::MassFitInfoEvent& event, float PUweight, float PUup, float PUdown, bool onlyHighestWeight, int bestCombi);

void fillMassFitResult(LightMonster* monster, ideogram::MassFitInfoEvent& event, bool onlyHighestWeight, int bestCombi);

int highestWeightCombi(LightMonster* monster, double chi2Cut, double bTagCut, double bTagEff, double misTagRate); // gives index of highest weight combi

#endif
