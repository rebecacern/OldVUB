#ifndef MakerJetOrigin_h
#define MakerJetOrigin_h

#include <iostream>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TH3.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TTree.h>
#include <sstream>
#include <vector>
#include "TKey.h"
#include "TLegend.h"
#include "TRandom3.h"

using namespace std;

class MakerJetOrigin{

 public:
  MakerJetOrigin();
  ~MakerJetOrigin();

  void add(bool, bool, bool, bool, double);
  void print();
  double get_b_Counter();
  double get_nonb_Counter();
  double get_hadqq_Counter();
  double get_radq_Counter();

 private:

  double is_b_Counter;
  double is_nonb_Counter;
  double is_hadqq_Counter;
  double is_radq_Counter;

};

#endif
