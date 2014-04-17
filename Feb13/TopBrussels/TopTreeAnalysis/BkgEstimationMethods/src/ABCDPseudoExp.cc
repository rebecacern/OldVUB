/*
 *  ABCDPseudoExp.cc
 *  AnalysisSW
 *
 *  Created by Thierry Caebergs on 13/04/10.
 *
 */

#include "../interface/ABCDPseudoExp.h"

#include "TRandom.h"

ABCDPseudoExp::ABCDPseudoExp()
{
  a = b = c = d = 0.;
  aExp = bExp = cExp = dExp = 0.;
  errA = errEstA = 0.;
  pull = pullErrExp = distrA = NULL;
}

ABCDPseudoExp::ABCDPseudoExp(double a, double b, double c, double d, double errA, double errEstA)
{
  this->a = a;
  this->b = b;
  this->c = c;
  this->d = d;
  this->errA = errA;
  this->errEstA = errEstA;
  this->rndm = new TRandom();
  this->aExp = this->bExp = this->cExp = this->dExp = 0.;
  pull = new TH1D("pull","pull of A region distibution", 20, -6., 6.);
  pull->SetDirectory(NULL);
  pullErrExp = new TH1D("pullErrExp","pull of A region distibution, using pseudo-experimental error", 20, -6., 6.);
  pullErrExp->SetDirectory(NULL);
  distrA = new TH1D("distrA", "Distribution of the population of the A region", 20, a/3., 4*a);
  distrA->SetDirectory(NULL);
}

ABCDPseudoExp::ABCDPseudoExp(double a, double b, double c, double d, double errA, unsigned int seed)
{
  ABCDPseudoExp();
  delete rndm;
  rndm = new TRandom(seed);
}

ABCDPseudoExp::ABCDPseudoExp(const ABCDPseudoExp& rhs)
{
  a = rhs.a;
  b = rhs.b;
  c = rhs.c;
  d = rhs.d;
  errA = rhs.errA;
  errEstA = rhs.errEstA;
  this->rndm = (TRandom *) rhs.rndm->Clone();
  this->aExp = rhs.aExp;
  this->bExp = rhs.bExp;
  this->cExp = rhs.cExp;
  this->dExp = rhs.dExp;
  pull = (TH1D*) rhs.pull->Clone();
  pullErrExp = (TH1D*) rhs.pullErrExp->Clone();
  distrA = (TH1D*) rhs.distrA->Clone();
}

ABCDPseudoExp::~ABCDPseudoExp(){ delete this->rndm; };

const ABCDPseudoExp& ABCDPseudoExp::operator=(const ABCDPseudoExp& rhs)
{
 if (&rhs != this) {
  this->~ABCDPseudoExp();
  ABCDPseudoExp(rhs);
 }
 return *this;
}

void ABCDPseudoExp::GenExp()
{
  this->aExp = rndm->PoissonD(this->a);
  this->bExp = rndm->PoissonD(this->b);
  this->cExp = rndm->PoissonD(this->c);
  this->dExp = rndm->PoissonD(this->d);
  double nTotExp = this->aExp + this->bExp + this->cExp + this->dExp;
  double aErrExp = sqrt(this->aExp*(1-this->aExp/nTotExp));
  if (dExp!=0)
  {
    double est = bExp*cExp/dExp;
    distrA->Fill(a);
    pull->Fill((est-a)/errA);
    pullErrExp->Fill((est-a)/aErrExp);
  }
}

void ABCDPseudoExp::RetrieveExp(double& aPseudoExp, double& bPseudoExp, double& cPseudoExp, double& dPseudoExp)
{
  aPseudoExp = this->aExp;
  bPseudoExp = this->bExp;
  cPseudoExp = this->cExp;
  dPseudoExp = this->dExp;
}
              
void ABCDPseudoExp::Write(std::string prefixName)
{
  pull->Write((prefixName+"ABCDPseudoExp_pull").c_str());
  pullErrExp->Write((prefixName+"ABCDPseudoExp_pullErrExp").c_str());
  distrA->Write((prefixName+"ABCDPseudoExp_distrA").c_str());
}
/*
if(anaEnv.nPseudoExp==0)
{
    //double a = PoissonD(abcdEstim.GetData(1));
  double b = PoissonD(abcdEstim.GetData(2));
  double c = PoissonD(abcdEstim.GetData(3));
  double d = PoissonD(abcdEstim.GetData(4));
  double est=0.;
  if (d!=0)
    est = b*c/d;
  pull->Fill((a-abcdEstim.GetData(1))/errA)
*/
