#include <iostream>
#include <cmath>
#include <cstdio>
#include "Interaction/Interaction.hpp"
#include "debug.hpp"
#include "Vector3.hpp"


Interaction::~Interaction()
{
  mesg("~Interaction\n");
}



LJ::LJ()
{
  set( 1.0L, 1.0L );
}



LJ::LJ( double e, double s )
{
  set( e, s );
}



void LJ::set( double e, double s )
{
  eps = e;
  sig = s;
}




void
LJ::Force( const Vector3& d, double& fo, double& potential ) const 
{
  double rri = 1.0 / (d.x*d.x + d.y*d.y + d.z*d.z);
  double drs = sig * sig * rri;
  double dr6 = drs * drs * drs;
  double dr12 = dr6 * dr6;
  fo = -4.0 * eps * (-12.0 * dr12 + 6.0 * dr6) * rri;
  potential = 4.0 * eps * (dr12-dr6);                
  //cerr << "LJSIG: " << sig << " < " << sqrt(1.0/rri) << " > " << potential << endl;
}



Coulomb::Coulomb()
{
  set2( 0 );
}



Coulomb::Coulomb( double c )
{
  set2( c );
}



//速度重視のため二乗値だけを設定し、chargeは設定しない。使用注意
void Coulomb::set2( double c )
{
  //charge = c;
  sqrch  = c;
}




void Coulomb::Force( const Vector3& d, double& fo, double& potential ) const
{
  double rri = 1.0 / (d.x*d.x + d.y*d.y + d.z*d.z);
  double radiusi = sqrt(rri);
  fo = sqrch*radiusi*rri;
  potential = sqrch*radiusi;                
  //cerr << "CSIG: " << rri << " = " << potential << endl;
}






LJC::LJC( const IntrParams& ljc1, const IntrParams& ljc2 )
{
  //Lorentz-Berthelot Rule
  set2( ljc1.epssqrt * ljc2.epssqrt,
        ljc1.sighalf + ljc2.sighalf,
        ljc1.charge * ljc2.charge );
}




LJC::LJC( double e, double s, double c )
{
  set2( e, s, c );
}



void LJC::set2( double e, double s, double c )
{
  lj.set( e, s );
  coulomb.set2( c );
}




void
LJC::Force( const Vector3& d, double& fo, double& potential ) const 
{
  double rri = 1.0 / (d.x*d.x + d.y*d.y + d.z*d.z);
  double radiusi = sqrt(rri);
  double drs = lj.sig * lj.sig * rri;
  double dr6 = drs * drs * drs;
  double dr12 = dr6 * dr6;
  fo = -4.0 * lj.eps * (-12.0 * dr12 + 6.0 * dr6) * rri + coulomb.sqrch*radiusi*rri;
  potential = 4.0 * lj.eps * (dr12-dr6) + coulomb.sqrch*radiusi;                
  //cerr << "LJCSIG: " << lj.sig << " < " << drs << " > " << potential << endl;
}



IntrParams::IntrParams( double e, double s, double c )
{
  set( e, s, c );
}

#define BUFSIZE 1024
void
IntrParams::read( FILE* const file, const Unit& u )
{
  char buf[BUFSIZE];
  double e, s, c;
  fgets(buf,BUFSIZE,file);
  sscanf(buf,"%lf %lf %lf\n",&e, &s, &c);
  //cerr << "E"<< e <<":S" << s <<":C" << c << endl;
  e *= SI::kilo * u.J / mol;
  s *= Angstro * u.m;
  c *= u.coulomb;
  set( e, s, c );
}



void IntrParams::set( double e, double s, double c )
{
  eps = e;
  epssqrt = sqrt( eps );
  sig = s;
  sighalf = sig / 2;
  charge = c;
}




IntrParams::IntrParams( FILE* const file, const Unit& unit )
{
  read( file, unit );
}



void
IntrParams::Write( const Unit& u, ostream& to ) const
{
  double e,s,c;
  e = eps / ( SI::kilo * u.J / mol );
  s = sig / ( Angstro * u.m );
  c = charge / ( u.coulomb );
  to << e
     << " " << s
     << " " << c << '\n';
}

