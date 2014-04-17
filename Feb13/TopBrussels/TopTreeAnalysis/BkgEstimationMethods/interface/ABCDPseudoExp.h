#ifndef _ABCDPseudoExp_h_
#define _ABCDPseudoExp_h_

#include <cmath>
#include <string>

#include "TH1D.h"
#include "TRandom.h"

/**
 *  AnalysisSW ABCDPseudoExp.h
 *
 *  Created by Thierry Caebergs () on 13/04/10.
 * 
 */

class ABCDPseudoExp
{
public:

  ABCDPseudoExp();
  
  ABCDPseudoExp(double a, double b, double c, double d, double errA, double errEstA);

  ABCDPseudoExp(double a, double b, double c, double d, double errA, unsigned int seed);

  ABCDPseudoExp(const ABCDPseudoExp& rhs);

  ~ABCDPseudoExp();

  const ABCDPseudoExp& operator=(const ABCDPseudoExp& rhs);

  void GenExp();

  void RetrieveExp(double& aPseudoExp, double& bPseudoExp, double& cPseudoExp, double& dPseudoExp);
  
  void Write(std::string prefixName=std::string(""));
  
  
private:
  double a,b,c,d;
  double errA, errEstA;
  double aExp,bExp,cExp,dExp;
  TH1D *pull, *pullErrExp, *distrA;
  TRandom *rndm;
};

#endif
