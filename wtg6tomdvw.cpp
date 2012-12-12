#include <cstdio>
#include <map>
#include "Mols/RigidBodies2.hpp"
#include "FileIO.hpp"




int main( int argc, char* argv[] )
{
  Unit unit;
  //MolPropertyHandle prop( MolPropertyHandle( new Rigid ) );
  char buf[1024];
  map<string, FlexibleHandle> moldict;
  string current_id08;
  BoxHandle box;
  const int isfixed = 0;

  while( NULL != fgets( buf, sizeof(buf), stdin ) ){
    string tag = string( buf );
    tag.resize(5);
    if ( tag == "@WTG6" || tag == "@NX4A" ){
      MolPropertyHandle ph = moldict[current_id08];
      RigidBodies2 mol( ph, isfixed );
      fgets( buf, sizeof( buf ), stdin );
      int nmol = atoi( buf );
      if ( tag == "@WTG6" ) {
	mol.ReadWTG6( nmol, unit, stdin );
      }
      else if ( tag == "@NX4A" ){
	mol.ReadNX4A( nmol, *box, unit, stdin );
      }	
      WriteMDVWHeader( ph, mol, box, unit, cout );
      WriteMDVW( mol, unit, cout );
    }
    else if ( tag == "@BOX3" ) {
      const int isP[3] = {1,1,1};
      SemiperiodicBox* sbox = new SemiperiodicBox( isP );
      sbox->ReadBOX3( unit, stdin );
      box = BoxHandle( sbox );
    }
    else if ( tag == "@DEFR" || tag == "@DEFP" ) {
      // define a molecule, put in the moldict
      //new ID08
      fgets( buf, sizeof( buf ), stdin );
      current_id08 = string( buf );
      current_id08.resize(8, ' ');
      //read molecular structure
      FlexibleHandle mp;
      if ( tag == "@DEFR" ){
        mp = ReadDEFR( current_id08, unit, stdin );
      }
      else {
        mp = ReadDEFP( current_id08, unit, stdin );
      }
      moldict.insert(
        map<string, FlexibleHandle>::value_type(
          current_id08,
          mp ));
    }
    else if ( tag == "@ID08" ) {
      fgets( buf, sizeof(buf), stdin );
      current_id08 = string( buf );
      current_id08.resize( 8, ' ' );
    }
    /*
    else{
      cout << string( buf );
    }
    */
  }
}
