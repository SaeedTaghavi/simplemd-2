#include <cstdio>
#include "SingleMol/MonatomicMol.hpp"
#include "Vector3.hpp"
#include "debug.hpp"
#include "Plugin/CollectorPlugin.hpp"
//Copy constructor
MonatomicMol::MonatomicMol( const MonatomicMol& m )
{
  *this = m;
}



//Normal Constructor
MonatomicMol::MonatomicMol( double m ) : mass(m), massi(1.0/m)
{
  velocity.Set(0,0,0);
  center.coord.Set(0,0,0);
  center.force.Set(0,0,0);
  order = -1;
}



//Destructor
MonatomicMol::~MonatomicMol()
{
  mesg("~MonatomicMol\n");
}



const Vector3&
MonatomicMol::Position() const
{
  return center.coord;
}



void MonatomicMol::Translate( const Vector3& offset )
{
  center.coord.x += offset.x;
  center.coord.y += offset.y;
  center.coord.z += offset.z;
}



void MonatomicMol::ScaleVelocity( const Vector3& r )
{
  velocity.x *= r.x;
  velocity.y *= r.y;
  velocity.z *= r.z;
}



void MonatomicMol::ScaleVelocity( double r )
{
  velocity.x *= r;
  velocity.y *= r;
  velocity.z *= r;
}



void MonatomicMol::AddVelocity( const Vector3& v )
{
  velocity.x += v.x;
  velocity.y += v.y;
  velocity.z += v.z;
}



void MonatomicMol::ScalePosition( const Vector3& r )
{
  center.coord.x *= r.x;
  center.coord.y *= r.y;
  center.coord.z *= r.z;
}



void MonatomicMol::Preforce()
{
  center.force = Vector3(0,0,0);
}





void
MonatomicMol::Postforce( //double vesseld1, double volume, 
                          double mass )
{
}



const Vector3&
MonatomicMol::GetVelocity() const
{
  return velocity;
}



double
MonatomicMol::GetEk() const
{
  return 0.5*mass*velocity.square();
}



void
MonatomicMol::BoxCoordinate( const Box& box )
{
  box.BoxCoordinate( center.coord );
}



void
MonatomicMol::ProgressPosition( double dt )
{
  center.coord.x += velocity.x * dt;
  center.coord.y += velocity.y * dt;
  center.coord.z += velocity.z * dt;
}



void
MonatomicMol::ProgressMomentum( double dt )
{
  velocity.x += center.force.x * dt * massi;
  velocity.y += center.force.y * dt * massi;
  velocity.z += center.force.z * dt * massi;
}



void
MonatomicMol::WriteATG5( const Unit& u, ostream& to )
{
  to << center.coord.x / ( Angstro * u.m )
     << " " << center.coord.y / ( Angstro * u.m )
     << " " << center.coord.z / ( Angstro * u.m )
    //convert from internal unit to A/ps
     << " " << velocity.x / ( Angstro * u.m / (SI::pico * u.sec) )
     << " " << velocity.y / ( Angstro * u.m / (SI::pico * u.sec) )
     << " " << velocity.z / ( Angstro * u.m / (SI::pico * u.sec) )
    //  unit is kJ/mol/A
     << " " << center.force.x / ( SI::kilo * u.J / (mol * Angstro * u.m ) )
     << " " << center.force.y / ( SI::kilo * u.J / (mol * Angstro * u.m ) )
     << " " << center.force.z / ( SI::kilo * u.J / (mol * Angstro * u.m ) )
     << " " << order << endl;
}



void
MonatomicMol::WriteAR3A( const Unit& u, ostream& to )
{
  to << center.coord.x / ( Angstro * u.m )
     << " " << center.coord.y / ( Angstro * u.m )
     << " " << center.coord.z / ( Angstro * u.m ) << endl;
}



SingleMolEntity::~SingleMolEntity()
{
  mesg("~SingleMolEntity\n");
}



void
MonatomicMol::ReadATG5( int o, const Unit& u, FILE* file )
{
  char buf[1000];
  order = o;
  
  fgets(buf,sizeof(buf),file);
  //split
  center.coord.x = atof(strtok(buf," \t")) * Angstro * u.m;
  center.coord.y = atof(strtok(NULL," \t")) * Angstro * u.m;
  center.coord.z = atof(strtok(NULL," \t")) * Angstro * u.m;
  velocity.x = atof(strtok(NULL," \t")) * ( Angstro * u.m / (SI::pico * u.sec) );
  velocity.y = atof(strtok(NULL," \t")) * ( Angstro * u.m / (SI::pico * u.sec) );
  velocity.z = atof(strtok(NULL," \t")) * ( Angstro * u.m / (SI::pico * u.sec) );
  center.force.x    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / (mol * Angstro * u.m ) );
  center.force.y    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / (mol * Angstro * u.m ) );
  center.force.z    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / (mol * Angstro * u.m ) );
}



void
MonatomicMol::ReadAR3A( int o, const Unit& u, FILE* file )
{
  char buf[1000];
  char *b;
  order = o;
  double xx[5], yy[5], zz[5];
  
  fgets(buf,sizeof(buf),file);
  
  b = buf;
  for( int k=0; k<6; k++ ){
    char *tok;
    double x,y,z;
    tok = strtok(b," \t");
    if(tok==NULL)
      return;
    b = NULL;
    x = atof(tok);
    y = atof(strtok(NULL," \t"));
    z = atof(strtok(NULL," \t"));
    x *= Angstro * u.m;
    y *= Angstro * u.m;
    z *= Angstro * u.m;
    if ( k == 0 ){
      center.coord.Set(x,y,z);
    }
    else{
      xx[k-1] = x;
      yy[k-1] = y;
      zz[k-1] = z;
    }    
  }
  velocity.x = 0;
  velocity.y = 0;
  velocity.z = 0;
}



void MonatomicMol::AddForce( const Vector3& f )
{
  center.force.x += f.x;
  center.force.y += f.y;
  center.force.z += f.z;
}



double*
MonatomicMol::unserialize( double* const p )
{
  const Vector3 v( p[0], p[1], p[2] );
  SetPosition( v );
  return p+3;
}



double*
MonatomicMol::serialize( double* p ) const
{
  const Vector3& center = Position();
  p[0] = center.x;
  p[1] = center.y;
  p[2] = center.z;
  return p+3;
}



double*
MonatomicMol::serializeforce( double* xi ) const
{
  const Vector3& force = GetForce();
  xi[0] = force.x * massi;
  xi[1] = force.y * massi;
  xi[2] = force.z * massi;
  return xi+3;
}
