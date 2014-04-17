#include "../interface/MakerJetOrigin.h"

MakerJetOrigin::MakerJetOrigin(){
  is_b_Counter = 0;
  is_nonb_Counter = 0;
  is_hadqq_Counter = 0;
  is_radq_Counter = 0;
}

MakerJetOrigin::~MakerJetOrigin(){

}

void MakerJetOrigin::add(bool is_b, bool is_nonb, bool is_hadqq, bool is_radq, double weight){
  if(is_b) is_b_Counter += weight;
  if(is_nonb) is_nonb_Counter += weight;
  if(is_hadqq) is_hadqq_Counter += weight;
  if(is_radq) is_radq_Counter += weight;
}

void MakerJetOrigin::print(){
  double totalEntries = is_b_Counter + is_nonb_Counter;
  double totalnonbEntries = is_hadqq_Counter + is_radq_Counter;

  cout << "total b+nonb:              " << is_b_Counter + is_nonb_Counter <<endl;
  cout << "is_b_Counter:              " << is_b_Counter <<  "   " << (double) is_b_Counter/totalEntries << endl;
  cout << "is_nonb_Counter:           " << is_nonb_Counter <<  "   " << (double) is_nonb_Counter/totalEntries << "     MCmatching deficit (-" << is_nonb_Counter-(is_hadqq_Counter + is_radq_Counter) << "=" << (is_hadqq_Counter + is_radq_Counter) << " (" << (double) (is_nonb_Counter-(is_hadqq_Counter + is_radq_Counter))/is_nonb_Counter *100 << "%) )" << endl;

  if(totalnonbEntries!=0) {
    cout << endl << "is_hadqq_Counter:          " << is_hadqq_Counter  << "   " << (double) is_hadqq_Counter/totalnonbEntries << "     (" << (double) is_hadqq_Counter/totalEntries << ")" << endl;
    cout << "is_radq_Counter:           " << is_radq_Counter  << "   " << (double) is_radq_Counter/totalnonbEntries << "     (" << (double) is_radq_Counter/totalEntries << ")" <<  endl;
  }


} 

double MakerJetOrigin::get_b_Counter(){return is_b_Counter;}
double MakerJetOrigin::get_nonb_Counter(){return is_nonb_Counter;}
double MakerJetOrigin::get_hadqq_Counter(){return is_hadqq_Counter;}
double MakerJetOrigin::get_radq_Counter(){return is_radq_Counter;}

