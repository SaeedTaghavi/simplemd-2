#include <cstring>
#include "FileIO.hpp"

/*
 *File I/O is not an essential feature of the object.
 *So it should not be included in the object.
 */



FlexibleHandle
ReadDEFR( string id08, const Unit& u, FILE* file )
{
  char buf[1024];
  
  //define new rigid molecule
  //number of sites
  fgets( buf, sizeof( buf ), file );
  int nsite = atoi( buf );
  //mass
  vector<double> mass;
  vector<Vector3> site;
  vector<string> atomName;
  for( int i=0; i<nsite; i++ ){
    fgets( buf, sizeof( buf ), file );
    double x,y,z,m;
    char b2[1024];
    sscanf( buf, "%lf %lf %lf %lf %s\n", &x, &y, &z, &m, b2 );
    x *= Angstro * u.m;
    y *= Angstro * u.m;
    z *= Angstro * u.m;
    m *= u.g / mol;
    site.push_back( Vector3( x,y,z ) );
    mass.push_back( m );
    atomName.push_back( string( b2 ) );
  }
  //interaction parameters
  IntrParamsArray ia;
  for( int i=0; i<nsite; i++ ){
    IntrParams ip( file, u );
    ia.push_back( ip );
  }
  int dof  = 6;
  return FlexibleHandle( new Rigid( nsite, ia, mass, site, atomName, id08, dof ) );
}



FlexibleHandle
ReadDEFP( string id08, const Unit& u, FILE* file )
{
  char buf[1024];
  
  //define new molecule
  //number of sites
  fgets( buf, sizeof( buf ), file );
  int nsite = atoi( buf );
  //mass
  vector<double> mass;
  vector<Vector3> site;
  vector<string> atomName;
  for( int i=0; i<nsite; i++ ){
    fgets( buf, sizeof( buf ), file );
    double m;
    char b2[1024];
    sscanf( buf, "%lf %s\n", &m, b2 );
    m *= u.g / mol;
    mass.push_back( m );
    atomName.push_back( string( b2 ) );
  }
  //interaction parameters
  IntrParamsArray ia;
  for( int i=0; i<nsite; i++ ){
    IntrParams ip( file, u );
    ia.push_back( ip );
  }
  int dof  = 3 * nsite;
  return FlexibleHandle( new Flexible( nsite, ia, mass, atomName, id08, dof ) );
}




void
WriteMDVW( const RigidBody2& r, const Unit& unit, ostream& to )
{
  Rigid* p = dynamic_cast<Rigid*> ( r.prop.get() );
  assert( p != 0 );
  int nsite = r.atom.size();
  for( int s=0; s<nsite; s++ ){
    to << p->GetAtomName(s)
       << " " << r.atom[s].center.coord.x + r.com.center.coord.x
       << " " << r.atom[s].center.coord.y + r.com.center.coord.y
       << " " << r.atom[s].center.coord.z + r.com.center.coord.z << endl;
  }
}


void
WriteMDVWHeader( const MolPropertyHandle& mph, const RigidBodies2& r, const BoxHandle& box, const Unit& unit, ostream& to )
{
  int nmol = r.Size();
  to << "-center 0 0 0" << endl;
  to << "-fold" << endl;
  to << "-length '("
     << box->x << ", " << box->y << ", " << box->z << ")'" << endl;
  int totalsites = 0;
  for(int i=0; i<nmol; i++){
    totalsites += r.mols[i].atom.size();
  }
  to << totalsites << endl;
}



void
WriteMDVW( RigidBodies2& r, const Unit& unit, ostream& to )
{
  int nmol = r.Size();
  sort( r.mols.begin(), r.mols.end() );
  for(int i=0; i<nmol; i++){
    r.mols[i].prepareintra();
    WriteMDVW( r.mols[i], unit, to );
  }
}



void
Read( PolyatomicMol& m, int o, FILE* file, const Unit& u )
{
  char buf[1000];
  char *b;
  double xx[6],yy[6],zz[6];
  m.order = o;
  
  for( int i=0; i<m.numAtom; i++ ){
    for(int j=0; j<6; j++ ){
      xx[j] = 0;
      yy[j] = 0;
      zz[j] = 0;
    }
    
    //read a line
    fgets(buf,sizeof(buf),file);
    
    //split
    b = buf;
    for( int k=0; k<6; k++ ){
      char* tok = strtok(b," \t");
      if(tok==NULL)
        break;
      b = NULL;
      xx[k] = atof(tok) * Angstro * u.m;
      yy[k] = atof(strtok(NULL," \t")) * Angstro * u.m;
      zz[k] = atof(strtok(NULL," \t")) * Angstro * u.m;
    }
    //put the values
    Vector3 position(xx[0],yy[0],zz[0]);
    m.atom[i].SetPosition( position );
    Vector3 velocity(0,0,0);
    m.atom[i].SetVelocity(velocity);
  }
  //report();
}



void WriteDTPS( double dt, const Unit& unit, ostream& to )
{
  to << "@DTPS" << endl
     << dt / ( pico * unit.sec ) << endl;
};



void
Write( map<string, FlexibleHandle>& moldict, const Unit& unit, ostream& to )
{
  typedef map<string, FlexibleHandle>::const_iterator CI;
  for( CI p = moldict.begin(); p!=moldict.end(); ++p ){
    Flexible& fh = *(p->second);
    fh.Write( unit, to );
  }
}



