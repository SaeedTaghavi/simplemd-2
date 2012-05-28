#include <cstdio>
#include "Mols/RigidBodies2.hpp"

int main( int argc, char* argv[] )
{
  Unit unit;
  MolPropertyHandle prop( MolPropertyHandle( new Rigid ) );
  char buf[1024];
  int nmol;

  while( NULL != fgets( buf, sizeof(buf), stdin ) ){
    string tag = string( buf );
    tag.resize(5);
    if ( tag == "@WTG6" ){
      RigidBodies2 mol( prop );
      fgets( buf, sizeof( buf ), stdin );
      nmol = atoi( buf );
      mol.ReadWTG6( nmol, unit, stdin );
      mol.WriteNX4A( unit, cout );
    }
    else{
      cout << string( buf );
    }
  }
}

    
