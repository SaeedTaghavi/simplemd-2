#include "ExchangePT.hpp"


ExchangePT::ExchangePT( const RandomD01Handle& r ) : Exchange( r )
{
}



bool
ExchangePT::Trial( ExchangePT& target )
{
  const vector<double>& vars2   = target.GetVars();
  const vector<double>& consts2 = target.GetConsts();
  double energy2 = consts2[0];
  double volume2 = consts2[1];
  double beta2   = 1.0/vars2[0];
  double pext2   = vars2[1];

  double energy1 = consts[0];
  double volume1 = consts[1];
  double beta1   = 1.0/vars[0];
  double pext1   = vars[1];

  //taken from 19319382
  double delta   = (beta2-beta1)*(energy1-energy2) + (beta2*pext2-beta1*pext1)*(volume1-volume2);
  cerr << delta << "(" << (beta2-beta1)*(energy1-energy2) << ") delta\n";

  bool accept = false;
  if ( delta < 0 ){
    accept = true;
  }
  else{
    double ratio = exp( - delta );
    if ( random->get() < ratio ){
      accept = true;
    }
  }
  if ( accept ){
    cerr << beta1 << ":" << beta2 << " beta\n";
    cerr << energy1 << ":" << energy2 << " energy\n";
    cerr << pext1 << ":" << pext2 << " pressure\n";
    cerr << volume1 << ":" << volume2 << " volume\n";
    cerr << "accepted\n";
  }
  return accept;
}



