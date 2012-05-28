#include <fstream>
#include <string>
#include "FileIO.hpp"
#include "Cell/Cell.hpp"
#include "Cell/GridCell.hpp"
#include "SingleMol/SingleMol.hpp"
#include "SingleMol/MonatomicMol.hpp"
#include "SingleMol/PolyatomicMol.hpp"
#include "SingleMol/RigidBody.hpp"
#include "SingleMol/RigidBody2.hpp"
#include "Mols/MonatomicMols.hpp"
#include "Mols/PolyatomicMols.hpp"
#include "Mols/RigidBodies.hpp"
#include "debug.hpp"
#include "System/System.hpp"
#include "MolCollection.hpp"
#include "PairTool.hpp"

//Compatibility w/ Genesis6
MolCollection*
PairTool::Read( FILE* input )
{
  Cell* cell = 0;
  char buf[1024];
  string current_id08;
  int ncompo = 1;
  int coordcount = 0;
  int ndivx=0, ndivy=0, ndivz=0;
  MolCollection* molecules = 0;

  //cerr << "Simple::Read()" << endl;
  while( NULL != fgets( buf, sizeof(buf), input ) ){

    string tag = string( buf );
    tag.resize( 5 );
    //cerr << tag << endl;
    if ( tag == "@BOX3" ) {
      const int isP[3] = {1,1,1};
      box = new SemiperiodicBox( isP );
      box->ReadBOX3( unit, input );
    }
    else if ( tag == "@VOXN" ) {
      // Cell division (by fixed number)
      fgets( buf, sizeof( buf ), input );
      sscanf( buf, "%d %d %d", &ndivx, &ndivy, &ndivz );
    }
    else if ( tag == "@RCOA" ) {
      // Truncation
      rc.ReadRCOA( unit, input );
    }
    else if ( tag == "@NCMP" )  {
      // Number of components
      fgets( buf, sizeof( buf ), input );
      ncompo = atoi( buf );
    }
    else if ( tag == "@DEFR" || tag == "@DEFP" ) {
      // define a molecule, put in the moldict
      //new ID08
      fgets( buf, sizeof( buf ), input );
      current_id08 = string( buf );
      current_id08.resize(8, ' ');
      //read molecular structure and interaction
      FlexibleHandle mp;
      if ( tag == "@DEFR" ){
	mp = ReadDEFR( current_id08, unit, input );
      }
      else {
	mp = ReadDEFP( current_id08, unit, input );
      }
      if ( mp )
	moldict.insert(
		       map<string, FlexibleHandle>::value_type(
							       current_id08,
							       mp ));
    }
    else if ( tag == "@ID08" ) {
      fgets( buf, sizeof(buf), input );
      current_id08 = string( buf );
      current_id08.resize( 8, ' ' );
    }
    else if ( tag == "@AR3A" 
	      || tag == "@NX4A" 
	      //|| tag == "@NX4B" 
	      || tag == "@ATG5" 
	      //|| tag == "@WTG5" 
	      || tag == "@WTG6" ) {
      coordcount ++;
      //Get the number of molecules from the input
      fgets( buf, sizeof(buf), input );
      int nmol = atoi( buf );
      //Make a copy from the moldict
      FlexibleHandle ph = moldict[current_id08];
      //ph->SetSize( nmol );
      MolsHandle mh;
      if ( tag == "@AR3A" ){
	MonatomicMols* mm = new MonatomicMols( ph );
	mh = MolsHandle( mm );
	mm->ReadAR3A( nmol, *box, unit, input );
      }
      else if ( tag == "@NX4B" ){
	RigidBodies* mm = new RigidBodies( ph );
	mh = MolsHandle( mm );
	mm->ReadNX4A( nmol, *box, unit, input );
      }
      else if ( tag == "@NX4A" ){
	RigidBodies* mm = new RigidBodies( ph );
	mh = MolsHandle( mm );
	mm->ReadNX4A( nmol, *box, unit, input );
      }
      else if ( tag == "@WTG5" ){
	RigidBodies* mm = new RigidBodies( ph );
	mh = MolsHandle( mm );
	mm->ReadWTG5( nmol, unit, input );
      }
      else if ( tag == "@ATG5" ){
	mh = MolsHandle( new MonatomicMols( ph, nmol, unit, input ) );
      }
      if ( cell == 0 ){
	if ( box != NULL ){
	  if ( ndivx && ndivy && ndivz ){
	    cell = new GridCell( *box, ncompo, ndivx, ndivy, ndivz );
	  }
	  else{
	    //Bookkeeping method is not so fast.
	    // and is temporarily unavailable for monatomic mols.
	    //cell = new BookPeriodicCell( box, ncompo );
	    //cell = new PeriodicCell( box, ncompo );
	    Vector3 offset(0,0,0);
	    cell = new RelocatableCell( *box, ncompo, offset );
	  }
	}
	else{
	  cell = new SimpleCell( ncompo );
	  assert(0);
	}
      }
      if ( molecules == 0 ){
	molecules = new MolCollection2( cell );
      }
      molecules->push_back( ph, mh );
    }
    else{
      if ( buf[0] == '@' ){
	cerr << buf;
      }
    }
    if ( coordcount == ncompo ){
      break;
    }
  }
  return molecules;
}
