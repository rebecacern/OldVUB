#ifndef MonsterTools_h
#define MonsterTools_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"

#include "TopTreeAnalysis/JESMeasurement/interface/LightMonster.h"

using namespace std;

bool hadrJetsMVAMatched(LightMonster* monster, int jetCombi=0);

vector<double> CalculateParabola(TH1D* parabola, int size, bool fit);

#endif
