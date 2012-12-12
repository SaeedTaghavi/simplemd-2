/*
 * main.cppをベースに、Quenchを作成する。ほとんど変更しなくていいはずなので、最終的にはmain.cppに融合してもいいと思う。
 */


#include <fstream>
#include <string>
#include "FileIO.hpp"
#include "System/System.hpp"
#include "System/QuenchSystem.hpp"
#include "Cell/Cell.hpp"
#include "Cell/PeriodicCell.hpp"
#include "Cell/GridCell.hpp"
#include "SingleMol/SingleMol.hpp"
#include "SingleMol/MonatomicMol.hpp"
#include "SingleMol/PolyatomicMol.hpp"
#include "SingleMol/RigidBody.hpp"
#include "SingleMol/RigidBody2.hpp"
#include "Mols/MonatomicMols.hpp"
#include "Mols/PolyatomicMols.hpp"
#include "Mols/RigidBodies.hpp"
#include "Integrator/NosePoincare.hpp"
#include "Integrator/NosePoincareAndersen.hpp"
#include "debug.hpp"

using namespace std;
using namespace boost;
using namespace SI;

#define BUFSIZE 1024



#ifdef DEBUG
int count;
#endif



class Sample : public Main
{
  /*
  double vf;
  PotVir pv;
  //single (Andersen's) box
  // Molecule container ( separated from the vessel )
  // Common cut-off function
  int dof;
  */

private:

  // time interval
  double dt;

   // is 1 if the cell stretches only in z-axis.
  int    zonly;

  // total steps for MD ( unused for Replica Exchange )
  int    nstep;

  // absolute time from initial run.
  double absoluteTime;

  // cell and molecules
  MolCollection2* molecules;

  // dictionary for molecular name and properties
  map<string, FlexibleHandle> moldict;

  // simulation cell size
  Box* box;

  // Standard unit <=> internal unit
  Unit* unit;

  // Interaction truncation
  Truncation rc;

  // The molecular dynamics algorithm
  SimpleMDHandle simplemd;

  // filename ( for replicas )
  string infilename;

  // output streams. Assigned to files in case of RXMD
  ostream* fout;
  ostream* ferr;

  // coeff for quenching
  double coeff;

  //
  void GetOptions( int argc, char* argv[] )
  {
    char c;
    int errflg=0;
    extern int optopt,optind;

    //nstep = 4;
    //return;
    while ((c = getopt(argc, argv, ":c:")) != -1)
      switch (c) {
        case 'c':
          coeff = atof( optarg );
          // TReset Plugin to reset the temperature every step
          break;
        case ':':        /* -t without arguments */
          fprintf( stderr, "Option -%c requires an argument\n",optopt);
          errflg++;
          break;
        case '?':
          fprintf( stderr, "Unrecognized option: - %c\n",optopt);
          errflg++;
      }
        if(optind+1!=argc)
          errflg++;
    if (errflg) {
      fprintf(stderr, "usage: %s [-t temp] nstep\n",argv[0]);
      exit (2);
    }
    nstep   = atoi(argv[argc-1]);
  }



  
public:
  Sample( int argc, char* argv[] ) : fout( &cout ), ferr( &cerr )
  {
    //FILE* file;
    unit    = new Unit();
    coeff = 1;

    GetOptions( argc, argv );
    molecules = 0;
    dt = 0;
    absoluteTime = 0;
    zonly = 0;
    box = 0;

    DefaultDict();
    int forceloaded = Read( stdin );

    //Plugins are passed to simplemd.
    //They are invoked in initialize(), running(), terminate() of simplemd.
    ProcessPluginArray pl;

    simplemd = SimpleMDHandle( new SimpleQuench2( 
                                 molecules,
                                 unit,
                                 rc,
                                 moldict,
                                 pl,
                                 coeff,
                                 nstep,
                                 fout,
                                 ferr
                                 ) );
    /*/
      simplemd = SimpleMDHandle( new SimpleQuench( 
                                 molecules,
                                 unit,
                                 rc,
                                 moldict,
                                 pl,
                                 coeff,
                                 fout,
                                 ferr
                                 ) );
    */
  }

  virtual ~Sample()
  {
    delete unit;
    delete box;

    //These are passed to simplemd and deleted there.
    //if ( nosepoincare ) delete nosepoincare;
    //if ( nosepoincareandersen ) delete nosepoincareandersen;
    //delete molecules;

    //delete ferr;
    //delete fout;
    cerr << "~Sample" << endl;
  }
  

  
  void
    DefaultDict()
  {
    //Sample code to set a new molecule in the dictionary.
    //Unit Atomic LJ molecule
    const Unit& u = *unit;
    int nsite = 1;
    //int nmol  = 0;
    double eps = 1.0 * kilo * u.J / mol;
    double sig = 1.0 * Angstro * u.m;
    int dof = 3;
    double m = 1.0 * u.g / mol;
    vector<string> atomName( 1, "ReducedAr" );
    vector<double> mass = vector<double>( nsite, m );
    IntrParams ip( eps, sig, 0 );
    IntrParamsArray  ia = IntrParamsArray( nsite, ip );
    string id08("REDUCEAR");
    moldict.insert(
      map<string, FlexibleHandle>::value_type(
        string( id08 ),
        FlexibleHandle( new Flexible( nsite, ia, mass, atomName, id08, dof ))));
  }
  
  
  
  //Compatibility w/ Genesis6
  int Read( FILE* input )
  {
    Cell* cell = 0;
    char buf[1024];
    string current_id08;
    int ncompo = 1;
    int coordcount = -1;
    int ndivx=0, ndivy=0, ndivz=0;
    int forceloaded = 1;
    int innerloop = 0;
    const int isfixed = 0;
    //cerr << "Simple::Read()" << endl;
    while( NULL != fgets( buf, sizeof(buf), input ) ){
      string tag = string( buf );
      tag.resize( 5 );
      //cerr << tag << endl;
      if ( tag == "@DTPS" ){
        fgets(buf,BUFSIZE, input);
        dt = atof( buf );
        // convert to internal unit
        dt *= pico * unit->sec;
      }
      else if ( tag == "@BOX3" ) {
        const int isP[3] = {1,1,1};
        box = new SemiperiodicBox( isP );
        box->ReadBOX3( *unit, input );
      }
      else if ( tag == "@VOXN" ) {
        // Cell division (by fixed number)
        fgets( buf, sizeof( buf ), input );
        sscanf( buf, "%d %d %d", &ndivx, &ndivy, &ndivz );
      }
      else if ( tag == "@RCOA" ) {
        // Truncation
        rc.ReadRCOA( *unit, input );
      }
      else if ( tag == "@NCMP" )  {
        // Number of components
        fgets( buf, sizeof( buf ), input );
        ncompo = atoi( buf );
      }
      else if ( tag == "@TABS" )  {
        // Absolute time past initial
        fgets( buf, sizeof( buf ), input );
        absoluteTime = atof( buf ) * ( pico * unit->sec );
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
	  mp = ReadDEFR( current_id08, *unit, input );
	}
	else {
	  mp = ReadDEFP( current_id08, *unit, input );
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
      else if ( tag == "@INNR" ) {
	fgets( buf, sizeof(buf), input );
	innerloop = atoi( buf );
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
          MonatomicMols* mm = new MonatomicMols( ph, isfixed );
          mh = MolsHandle( mm );
          mm->ReadAR3A( nmol, *box, *unit, input );
          forceloaded = 0;
        }
        else if ( tag == "@NX4B" ){
          RigidBodies* mm = new RigidBodies( ph, isfixed );
          mh = MolsHandle( mm );
          mm->ReadNX4A( nmol, *box, *unit, input );
          forceloaded = 0;
        }
        else if ( tag == "@NX4A" ){
          RigidBodies* mm = new RigidBodies( ph, isfixed );
          mh = MolsHandle( mm );
          mm->ReadNX4A( nmol, *box, *unit, input );
          forceloaded = 0;
        }
        else if ( tag == "@WTG5" ){
          RigidBodies* mm = new RigidBodies( ph, isfixed );
          mh = MolsHandle( mm );
          mm->ReadWTG5( nmol, *unit, input );
        }
        else if ( tag == "@ATG5" ){
          mh = MolsHandle( new MonatomicMols( ph, nmol, *unit, input, isfixed ) );
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
    }

    return forceloaded;
  }

  
  
  void initialize()
    {
      simplemd->initialize();
    }
  
  
  bool running()
  {
    return simplemd->running();
  }
  
  
  
  void terminate()
  {
    simplemd->terminate();
  }
  
  
  
};



int main(int argc,char *argv[])
{
  Sample *system = new Sample( argc, argv );
  system->initialize();
  while( system->running() );
  system->terminate();
  delete system;
}



