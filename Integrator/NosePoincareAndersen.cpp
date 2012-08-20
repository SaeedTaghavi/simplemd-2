#include <cmath>
#include <cstdio>
#include "NosePoincareAndersen.hpp"

void
NosePoincareAndersen::ProgressPv( double tau, double pressure )
{
  //Eq. (9b,h)
  pv += tau * s * ( pressure - Pext );
  //cerr << pressure << " Pv\n";
}



void
NosePoincareAndersen::ProgressPs( double tau, double ekss, double U, int dof, double volume )
{
  if ( ! disableTempControl ){
    //Eq. (9c,g)
    double dH = ekss + U + Ekv() + Eks() + Epv( volume ) + Eps( dof ) - Hzero();
    double gkT = ( dof+1 ) * kT;
    ps += tau * ( ekss*2 - pv - gkT - dH );
    //cerr << ekss << " " << pv << " " << gkT << " " << dH << " Ps\n";
    //(): ekss*2 - gkT - ekss - U - Eks - gkTlog(s) + Hzero;
  }
}



double
NosePoincareAndersen::ProgressSV( double tau )
{
  //Eq. (9d,e)
  double deltav;
  deltav  = tau * s * ( pv / Qv );
  if ( ! disableTempControl ){
    s      *= ( 1 + tau * ps / ( 2 * Qs ) ) / ( 1 - tau * ps / ( 2 * Qs ) );
  }
  deltav += tau * s * ( pv / Qv );
  //cerr  << pv <<" " << deltav << " " << " SV\n";
  return deltav;
}






double
NosePoincareAndersen::Eks()
{
  return ps*ps / ( 2*Qs);
}



double
NosePoincareAndersen::Ekv()
{
  return pv*pv / ( 2*Qv);
}


double
NosePoincareAndersen::Eps( int dof )
{
  return (dof + 1) * kT * log(s);
}



double
NosePoincareAndersen::Epv( double volume )
{
  return Pext * volume;
}



void NosePoincareAndersen::ReadNOPA( const Unit& unit, FILE* file )
{
  char buf[1000];
  fgets( buf, sizeof ( buf ), file );
  sscanf( buf, "%lf %lf %lf %lf %lf %lf %lf %lf\n", &kT, &Pext, &s, &ps, &Qs, &pv, &Qv, &H0);
  ps *= ( pico * unit.sec ) * ( kilo * unit.J / mol );
  Qs *= ( pico * unit.sec ) * ( pico * unit.sec ) * ( kilo * unit.J / mol );
  H0 *= ( kilo * unit.J / mol );
  kT *= ( unit.K() );
  Pext *= ( unit.atm() );
  pv *= ( pico * unit.sec ) * ( unit.Pa );
  Qv *= ( pico * unit.sec ) * ( pico * unit.sec ) * ( unit.Pa ) / pow( Angstro * unit.m, 3 );
}



void NosePoincareAndersen::WriteNOPA( const Unit& unit, ostream& to )
{
  to << "@NOPA\n"
     << kT / ( unit.K() )
     << " " << Pext / ( unit.atm() )
     << " " << s
     << " " << ps / ( pico * unit.sec * kilo * unit.J / mol )
     << " " << Qs / ( pico * unit.sec * pico * unit.sec * kilo * unit.J / mol )
     << " " << pv / ( pico * unit.sec * unit.Pa )
     << " " << Qv / ( pico * unit.sec * pico * unit.sec * unit.Pa / pow( Angstro * unit.m, 3 ) )
     << " " << H0 / ( kilo * unit.J / mol ) << endl;
}



void NosePoincareAndersen::WriteNPAZ( const Unit& unit, ostream& to )
{
  to << "@NPAZ\n"
     << kT / ( unit.K() )
     << " " << Pext / ( unit.atm() )
     << " " << s
     << " " << ps / ( pico * unit.sec * kilo * unit.J / mol )
     << " " << Qs / ( pico * unit.sec * pico * unit.sec * kilo * unit.J / mol )
     << " " << pv / ( pico * unit.sec * unit.Pa )
     << " " << Qv / ( pico * unit.sec * pico * unit.sec * unit.Pa / pow( Angstro * unit.m, 3 ) )
     << " " << H0 / ( kilo * unit.J / mol ) << endl;
}



double
NosePoincareAndersen::RecalcH0( double ek, double ep, int dof, double volume )
{
  return ek + ep + Eks() + Eps( dof ) + Ekv() + Epv( volume );
}



//Sturgeon's NPH
void NosePoincareAndersen::ReadSNPH( const Unit& unit, FILE* file )
{
  char buf[1000];
  fgets( buf, sizeof ( buf ), file );
  sscanf( buf, "%lf %lf %lf %lf\n", &Pext, &pv, &Qv, &H0);
  H0 *= ( kilo * unit.J / mol );
  Pext *= ( unit.atm() );
  pv *= ( pico * unit.sec ) * ( unit.Pa );
  Qv *= ( pico * unit.sec ) * ( pico * unit.sec ) * ( unit.Pa ) / pow( Angstro * unit.m, 3 );
  s = 1.0;
  ps = 0.0;
  Qs = 1e99;
  //Replica交換の際に使うので、破綻しないように0でない値を入れておく。
  kT = 1.0;
  disableTempControl = true;
}



void NosePoincareAndersen::WriteSNPH( const Unit& unit, ostream& to )
{
  to << "@SNPH\n"
     << " " << Pext / ( unit.atm() )
     << " " << pv / ( pico * unit.sec * unit.Pa )
     << " " << Qv / ( pico * unit.sec * pico * unit.sec * unit.Pa / pow( Angstro * unit.m, 3 ) )
     << " " << H0 / ( kilo * unit.J / mol ) << "\n";
}



void NosePoincareAndersen::Write( const Unit& unit, ostream& to )
{
  if ( disableTempControl ){
    WriteSNPH( unit, to );
  }
  else{
    WriteNOPA( unit, to );
  }
}
