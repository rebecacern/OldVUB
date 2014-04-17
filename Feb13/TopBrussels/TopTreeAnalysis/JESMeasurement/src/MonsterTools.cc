#include "../interface/MonsterTools.h"

bool hadrJetsMVAMatched(LightMonster* monster, int jetCombi)
{
  bool goodCombi = false;
  unsigned int* MVAres = monster->mvaResult(jetCombi);
  int hadrBJetIndex = MVAres[2], lightJet1Index = MVAres[0], lightJet2Index = MVAres[1];
  if( hadrBJetIndex == monster->hadrBJet() && ( ( lightJet1Index == monster->hadrLJet1() && lightJet2Index == monster->hadrLJet2() )
    || ( lightJet2Index == monster->hadrLJet1() && lightJet1Index == monster->hadrLJet2() ) ) )
    goodCombi = true;
  return goodCombi;
}

vector<double> CalculateParabola(TH1D* parabola, int size, bool fit)
{ 
  vector<double> getParabolaValues;
  Int_t minBin = parabola->GetMinimumBin();

  //cout << " DChi2 value of minbin-1 " << parabola->GetBinContent(minBin-1) << endl;
  //cout << " DChi2 value of minbin " << parabola->GetBinContent(minBin) << endl;
  //cout << " DChi2 value of minbin+1 " << parabola->GetBinContent(minBin+1) << endl;
  
  if (!fit)
  {
    double x1 = parabola->GetBinCenter(minBin);
    double x2 = parabola->GetBinCenter(minBin+size);
    double x3 = parabola->GetBinCenter(minBin-size);
    double y1 = parabola->GetBinContent(minBin);
    double y2 = parabola->GetBinContent(minBin+size);
    double y3 = parabola->GetBinContent(minBin-size);
    
    //y = ax^2 + bx + c
    double a = (y3-y1)/((x3-x1)*(x3-x2))-(y2-y1)/((x2-x1)*(x3-x2));
    double b = (y2-y1)/(x2-x1)-(y3-y1)*(x2+x1)/((x3-x1)*(x3-x2))+(y2-y1)*(x2+x1)/((x2-x1)*(x3-x2));
    double c = y1-(y3-y1)*(x1*x1)/((x3-x1)*(x3-x2))+(y2-y1)*(x1*x1)/((x2-x1)*(x3-x2))-(y2-y1)*(x1)/(x2-x1)+(y3-y1)*(x2+x1)*(x1)/((x3-x1)*(x3-x2))-(y2-y1)*(x2+x1)*(x1)/((x2-x1)*(x3-x2));
    
    double minim = - b/(2*a); //minimum gives JES correction 
    double minimchi2 = a*minim*minim + b*minim +c;
    double minplus = (-b + sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //upper bound right
    double minmin  = (-b - sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //lower bound left
    
    double uncert1 = minplus-minim;
    double uncert2 = minim-minmin;
    double uncert=0.;
    if(uncert1 > uncert2) { uncert = uncert1; } else { uncert = uncert2; }
    
    getParabolaValues.push_back(a);
    getParabolaValues.push_back(b);
    getParabolaValues.push_back(c);
    getParabolaValues.push_back(minim);
    getParabolaValues.push_back(uncert);

    //cout << "Analytic min " << minim << " minmin " << minmin << " minplus " << minplus << " uncert " << uncert << endl;
  }
  else
  {
    double lowLim = parabola->GetBinCenter(minBin) - parabola->GetBinCenter(minBin - size); //(size*0.02);
    double highLim = parabola->GetBinCenter(minBin) + parabola->GetBinCenter(minBin + size); //(size*0.02);

    string func_title = "Fitted";
    TF1* func = new TF1(func_title.c_str(),"pol2");

    func->SetRange(lowLim,highLim);
    //func->FixParameter(0,0);
    
    parabola->Fit(func,"RQ0");

    double a = func->GetParameter(2);
    double b = func->GetParameter(1);
    double c = func->GetParameter(0);

    double minim = - b/(2*a); //minimum gives JES correction 
    double minimchi2 = a*minim*minim + b*minim +c;
    
    double minplus = (-b + sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //upper bound right
    double minmin  = (-b - sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //lower bound left
    
    double uncert1 = minplus-minim;
    double uncert2 = minim-minmin;
    double uncert=0.;
    if(uncert1 > uncert2) { uncert = uncert1; } else { uncert = uncert2; }

    getParabolaValues.push_back(a);
    getParabolaValues.push_back(b);
    getParabolaValues.push_back(c);
    getParabolaValues.push_back(minim);
    getParabolaValues.push_back(uncert);

//    cout << "a " << a << "  b " << b << "  c " << c << endl;
//    cout << "Fitted min " << minim << "  minmin " << minmin << "  minplus " << minplus << "  uncert " << uncert << endl;
  }

  return vector<double>(getParabolaValues);
}

