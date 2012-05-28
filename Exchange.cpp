#include "Exchange.hpp"


Exchange::Exchange( const RandomD01Handle& r ) : random( r )
{
}



void
Exchange::Set( vector<double>& c, vector<double>& v )
{
  vars = v;
  consts = c;
}



vector<double>&
Exchange::GetVars()
{
  return vars;
}



const vector<double>&
Exchange::GetConsts() const
{
  return consts;
}



bool
Exchange::Trial( Exchange& target )
{
  const vector<double>& consts2 = target.GetConsts();
  const vector<double>& vars2 = target.GetVars();
  double energy2 = consts2[0];
  double kT2     = vars2[0];
  double energy1 = consts[0];
  double kT1     = vars[0];
  double acceptance0 = -( energy1 / kT1 + energy2 / kT2 );
  double acceptance1 = -( energy1 / kT2 + energy2 / kT1 );
  bool accept = false;
  cerr << "(" << acceptance0 << "," << acceptance1 << ")\n";
  double beta1   = 1.0/vars[0];
  double beta2   = 1.0/vars2[0];
  double delta   = (beta2-beta1)*(energy1-energy2);
  cerr << beta1 << ":" << beta2 << " beta\n";
  cerr << energy1 << ":" << energy2 << " energy\n";
  cerr << delta << " delta\n";
  if ( acceptance0 < acceptance1 ){
    accept = true;
  }
  else{
    double ratio = exp( acceptance1 - acceptance0 );
    cerr << ratio << " ratio\n";
    ratio = exp( -delta );
    double r = random->get();
    cerr << r << " rand\n";
    if ( r < ratio ){
      accept = true;
    }
  }
  if ( accept ){
    cerr << "accepted\n";
  }
  return accept;
}



void
Exchange::exchange( Exchange& target )
{
  //体積やエネルギーは、交換されない。温度や圧力が交換される。
  vector<double>& targetVars = target.GetVars();
  vars.swap( targetVars );
}
