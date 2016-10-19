#include <cassert>
#include <iostream>
#include <cstdio>
#include "MolProperty.hpp"
#include "debug.hpp"

Flexible::Flexible( unsigned int nsite, const IntrParamsArray& i, int DoF ) : MolProperty( string("") )
{
  intr = i;
  mass = vector<double>(nsite, 1.0);
  massi = vector<double>(nsite, 1.0);
  //totalNum = totalnum;
  numSite  = nsite;
  dof = DoF;
}

  

Flexible::Flexible( unsigned int nsite, const IntrParamsArray& i, vector<double>& mass_, vector<string>& name, string id, int DoF ) : 
  numSite( nsite ), intr( i ), mass( mass_ ), MolProperty( id ), atomName(name), dof(DoF)
{
  assert( nsite == mass.size() );
  massi.resize(nsite);
  for( unsigned int j=0; j< nsite; j++ ){
    massi[j] = 1 / mass[j];
  }
}

  

Flexible::Flexible() : MolProperty( string("") )
{
  intr = IntrParamsArray(0);
  mass  = vector<double>(0);
  massi = vector<double>(0);
  //totalNum = 0;
  numSite  = 0;
  dof = 0;
}



Flexible::~Flexible()
{
  mesg("~Flexible\n");
}




#define BUFSIZE 1024
//Read from file, put in the cells

//Compatibility for sample input data
void Flexible::ReadMonatomic( FILE* file )
{
  char buf[BUFSIZE];
  double m;
  fgets(buf,BUFSIZE,file);
  sscanf(buf,"%lf\n",&m);
  mass.resize(1);
  massi.resize(1);
  mass[0] = m;
  massi[0] = 1 / m;
  dof = 3;
  //set numbe of sites in a molecule
  numSite = 1;
  // and resize intr
  intr = IntrParamsArray( 1 );
  Unit unit;
  assert( 0 );
  //read interaction parameters from file and put into intr
  fgets( buf, sizeof(buf), file );
  //sscanf( buf, "%d", &totalNum );
  assert( 0 );
}




const IntrParamsArray&
Flexible::GetIntrParams() const
{
  return intr;
}



const vector<double>&
Flexible::GetMassArray() const
{
  return mass;
}



double
Flexible::GetMass() const
{
  double m=0;
  for(int i=0;i<numSite;i++)
    m += mass[i];
  return m;
}



const vector<double>&
Flexible::GetMassI() const
{
  return massi;
}



void
Flexible::Write( const Unit& unit, ostream& to ) const
{
  WriteDEFP( unit, to );
}






void
Flexible::WriteDEFP( const Unit& u, ostream& to ) const
{
  to << "@DEFP\n"
     << id08 << "\n";
  int nsite = GetNumSite();
  to << nsite << '\n';
  for( int i=0; i<nsite; i++ ){
    to << mass[i] / ( u.g / mol )
       << " " << atomName[i] << '\n';
  }
  for( int i=0; i<nsite; i++ ){
    intr[i].Write( u, to );
  }
}



Rigid::Rigid( int nsite, const IntrParamsArray& i, vector<double>& m, vector<Vector3>& s, vector<string>& name, string id, int DoF )
{
  intr     = i;
  mass     = m;
  massi.resize(nsite);
  //totalNum = totalnum;
  site     = s;
  id08     = id;
  numSite  = nsite;
  atomName = name;
  dof      = DoF;
  total_mass = 0;
  for( int i=0; i<nsite; i++ ){
      total_mass += mass[i];
      massi[i] = 1 / mass[i];
  }
  //Inertia
  inertia.Set(0,0,0);
  for( int i=0; i<nsite; i++ ){
    Vector3& s = site[i];
    inertia.x += mass[i] * (s.y*s.y + s.z*s.z);
    inertia.y += mass[i] * (s.z*s.z + s.x*s.x);
    inertia.z += mass[i] * (s.x*s.x + s.y*s.y);
  }
  //perfect
}

  

Rigid::Rigid()
{
  intr = IntrParamsArray(0);
  mass = vector<double>(0);
  massi = vector<double>(0);
  //totalNum = 0;
  numSite  = 0;
  site = vector<Vector3>( numSite );
}



Rigid::~Rigid()
{
  mesg("~Rigid\n");
}




#define BUFSIZE 1024
//Read from file, put in the cells

//Compatibility for sample input data
void Rigid::ReadMonatomic( FILE* file )
{
  char buf[BUFSIZE];
  double m;
  fgets(buf,BUFSIZE,file);
  sscanf(buf,"%lf\n",&m);
  total_mass = m;
  //set numbe of sites in a molecule
  numSite = 1;
  // and resize intr
  Vector3 com = Vector3(0,0,0);
  site = vector<Vector3>( numSite, com );
  intr = IntrParamsArray( 1 );
  Unit unit;
  assert( 0 );
  //read interaction parameters from file and put into intr
  //!!!intr[0] = IntrParamsHandle( new LJ( file, unit ) );
  fgets( buf, sizeof(buf), file );
  //sscanf( buf, "%d", &totalNum );
  assert(0);
}




const IntrParamsArray&
Rigid::GetIntrParams() const
{
  return intr;
}



double
Rigid::GetMass() const
{
  return total_mass;
}



const Vector3&
Rigid::GetInertia() const
{
  return inertia;
}



void
Rigid::SetInertia( const Vector3& v )
{
  inertia = v;
}




void
Rigid::Write( const Unit& unit, ostream& to ) const
{
  WriteDEFR( unit, to );
}



void
Rigid::WriteDEFR( const Unit& u, ostream& to ) const
{
  to << "@DEFR\n"
     << id08 << '\n';
  int nsite = GetNumSite();
  to << nsite << '\n';
  for( int i=0; i<nsite; i++ ){
    to << site[i].x / ( Angstro * u.m )
       << " " << site[i].y / ( Angstro * u.m )
       << " " << site[i].z / ( Angstro * u.m )
       << " " << mass[i] / ( u.g / mol )
       << " " << atomName[i] << '\n';
  }
  for( int i=0; i<nsite; i++ ){
    intr[i].Write( u, to );
  }
}
