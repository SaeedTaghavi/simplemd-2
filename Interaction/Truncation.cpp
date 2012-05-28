#include <iostream>
#include "Interaction/Truncation.hpp"
#include "debug.hpp"
#include <cstdio>

using namespace SI;

Truncation::Truncation()
{
  //mesg("Truncation(): must be initialized explicitly.\n");
  //exit(1);
  inner = outer = rq = 0;
}



Truncation::Truncation( const Unit& u, FILE* file )
{
  char buf[1024];
  double in, out;
  fgets(buf,sizeof(buf),file);
  sscanf(buf,"%lf\n",&in);
  fgets(buf,sizeof(buf),file);
  sscanf(buf,"%lf\n",&out);
  Set( in * ( Angstro * u.m ), 
       out * ( Angstro * u.m ) );
}


void 
Truncation::Set( double in, double out )
{
  double dr;
  inner = in;
  outer = out;
  dr = inner - outer;
  rq = 1.0/(dr*dr*dr*dr*dr);
}



void
Truncation::ReadRCOA( const Unit& u, FILE* file )
{
  char buf[1024];
  double width, out;
  fgets(buf,sizeof(buf),file);
  sscanf(buf,"%lf %lf\n",&out, &width);
  out       *= ( Angstro * u.m );
  width     *= ( Angstro * u.m );
  double in = out - width;
  Set( in, out );
}



Truncation::~Truncation()
  {
    mesg("~Truncation()\n");
  }



void Truncation::WriteRCOA( const Unit& u, ostream& to )
{
  double in = inner / ( Angstro * u.m );
  double out = outer / ( Angstro * u.m );
  to << "@RCOA\n"
     << out
     << " " <<  out-in << '\n';
}


void Truncation::Smooth( double r, double& sfe, double& sff ) const
{
  double rrl,rrs,rrs2,rrl2;
  rrl = r - outer;
  rrs = r - inner;
  rrl2 = rrl*rrl;
  rrs2 = rrs*rrs;
  sfe = rrl2*rrl*rq*(10.0*rrs2-5.0*rrl*rrs+rrl2);
  sff = 30.0*rq*rrl2*rrs2/r;
}
  
  
  
