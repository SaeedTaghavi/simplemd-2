#include <cmath>
#include <cstdio>
#include "NosePoincare.hpp"

void
NosePoincare::ProgressL3( double tau )
{
  double coeff = ( 1 + tau * ps / (2 * Q ) );
  s  *= coeff * coeff;
  ps /= coeff;
  //cerr << ps << " " << coeff << " " << s << " L3\n";
}


void
NosePoincare::ProgressL2( double tau, double U )
{
  ps -= tau * U;
  //cerr << ps << " " << U << " L2\n";
}



void
NosePoincare::ProgressL1( double tau, double Ek, int dof )
{
  double gkT = kT * ( dof+0 );
  //ps -= tau * ( -Ek * 2 / (s*s) + gkT * ( 1 + log(s) ) - H0 );
  ps -= tau * ( -Ek * 1 / (s*s) + gkT * ( 1 + log(s) ) - H0 );
  //Sturgeonの式とてらして、ここは-ek/ssじゃないと合わないはず。
  //cerr << ps << " " << Ek * 2 / (s*s) << " " << gkT << " " << H0 << " L1\n";
  //(): -2ekss + gkT + gkTlog(s) + U - H0 + eks
}



double
NosePoincare::Ek()
{
  return ps*ps / ( 2*Q);
}


double
NosePoincare::Ep( int dof )
{
  return (dof + 0) * kT * log(s);
}



void NosePoincare::ReadNOPO( const Unit& unit, FILE* file )
{
  char buf[1000];
  fgets( buf, sizeof ( buf ), file );
  sscanf( buf, "%lf %lf %lf %lf %lf\n", &kT, &s, &ps, &Q, &H0);
  ps *= ( SI::pico * unit.sec ) * ( SI::kilo * unit.J / mol );
  Q  *= ( SI::pico * unit.sec ) * ( SI::pico * unit.sec ) * ( SI::kilo * unit.J / mol );
  H0 *= ( SI::kilo * unit.J / mol );
  kT *= ( unit.K() );
}



void NosePoincare::WriteNOPO( const Unit& unit, ostream& to )
{
  to << "@NOPO\n"
     << kT / ( unit.K() )
     << " " << s
     << " " << ps / ( SI::pico * unit.sec * SI::kilo * unit.J / mol )
     << " " << Q  / ( SI::pico * unit.sec * SI::pico * unit.sec * SI::kilo * unit.J / mol )
     << " " << H0 / ( SI::kilo * unit.J / mol ) << "\n";
}


  
void
NosePoincare::ResetH0( double ek, double ep, int dof )
{
  H0 = ek + ep + Ek() + Ep( dof );
}
