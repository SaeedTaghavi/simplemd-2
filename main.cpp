/*
 * nve + variable volume = nph
 * 完全にep+ek+vrが安定しないのは、近似のせいか。
 * ともかく、とりあえず運動エネルギーと体積変化PVの相殺は確認した。

changes from nph2m.cpp:
  System.cpp is separated.
changes from nph2k.cpp:
  Compatibility with Genesis6 for REDUCEAR
changes from nph2i.cpp:
  Interaction is generalized.
changes to nph2i.cpp
  vector and shared_ptr is utilized.
  splitted according to classes.
changes from nph2e.cc:
  use boost and simplify
changes from nph2d.cc:
  introduce STL
changes from nph2.cc:
  molecular set in introduced in place of MonatomicMolPosition

 */


extern "C" {
#include "FakeMPI.h"
}
#include <fstream>
#include <string>
#include <cstring>
#include "Plugin/Plugin.hpp"
#include "Plugin/SnapShot.hpp"
#include "Plugin/SnapShotOneComp.hpp"
#include "Plugin/Log.hpp"
#include "Plugin/Interval.hpp"
#include "Plugin/Terminate.hpp"
#include "Plugin/TReset.hpp"
#include "Plugin/StopDrift.hpp"
#include "Plugin/Replica.hpp"
#include "Plugin/CPReplica.hpp"
#include "Plugin/LoadBalancer.hpp"
#include "Plugin/LoadBalancer2.hpp"
#include "Plugin/CollectorPlugin.hpp"
#include "Plugin/Zumbrella.hpp"
#include "Plugin/FixPlugin.hpp"
#include "Plugin/FixComponentPlugin.hpp"
#include "FileIO.hpp"
#include "System/System.hpp"
#include "System/SimpleMD.hpp"
#include "System/NosePoincareMD.hpp"
#include "System/NosePoincareAndersenMD.hpp"
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
#include "Mols/RigidBodies2.hpp"
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
  MolCollection* molecules;

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

  // Container of the variables for Nose-Poincare Integrator
  NosePoincare* nosepoincare;

  // Container of the variables for Nose-Poincare-Andersen Integrator
  NosePoincareAndersen* nosepoincareandersen;

  //for replica exchange
  LoadBalancer2* balancer;

  //Plugin list for simplemd
  ProcessPluginArray plugins;

  // filename ( for replicas )
  string infilename;

  // output streams. Assigned to files in case of RXMD
  // this MUST be a pointer instead of reference
  ostream* fout;
  ostream* ferr;
  bool redirected;

  //
  sFakeMPI* fmpi;

  //コマンドラインオプションは、エラー出力にもそのまま表示するべき。
  void GetOptions( int argc, char* argv[] )
  {
    double bathTemp;

    char c;
    int errflg=0;
    extern char *optarg;
    extern int optopt,optind;

    //nstep = 4;
    //return;
    while ((c = getopt(argc, argv, ":t:")) != -1)
      switch (c) {
        case 't':
          bathTemp = atof( optarg ) * unit->K();
          // TReset Plugin to reset the temperature every step
          plugins.push_back( ProcessPluginHandle( new TReset( bathTemp ) ) );
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
    fout->precision(17);
    ferr->precision(17);
    unit    = new Unit();
    box     = 0;
    GetOptions( argc, argv );
    molecules = 0;
    dt = 0;
    nosepoincare = 0;
    nosepoincareandersen = 0;//new NosePoincareAndersen( 1, 0, 1e2 );
    absoluteTime = 0;
    zonly = 0;
    balancer = 0; // 並列処理の負荷均衡。場合によってはReplica交換法
    fmpi = 0;
    redirected = 0;

    DefaultDict();
    int forceloaded = Read( stdin );
    if ( balancer ){
      //
      //複数のノード用の個別設定ファイルを読みこむ。
      //ここにくる前にすでにプロセスはforkしている。
      //
      FILE* input = fopen( infilename.c_str(), "r" );
      assert( input != NULL );
      bool result = Read( input );
      forceloaded = forceloaded || result;
      fclose( input );
    }
    else {
      //Simple plugin to terminate infinite loop of simplemd.
      //LoadBalancer automatically terminates.
      plugins.push_back( ProcessPluginHandle( new Terminate( nstep ) ) );
    }
    //Plugins are passed to simplemd.
    //They are invoked in initialize(), running(), terminate() of simplemd.
    if ( nosepoincare ){
      cerr << "NosePoincare\n";
      simplemd = SimpleMDHandle( new NosePoincareMD( 
                                                     absoluteTime,
                                                     dt,
                                                     molecules,
                                                     unit,
                                                     rc,
                                                     moldict,
                                                     plugins,
                                                     nosepoincare,
                                                     forceloaded,
						     fout,
						     ferr
						     ) );
    }
    else if ( nosepoincareandersen ){
      if ( zonly ){
        simplemd = SimpleMDHandle( new NosePoincareAndersenZMD(
                                                                absoluteTime,
                                                                dt,
                                                                molecules,
                                                                unit,
                                                                rc,
                                                                moldict,
								plugins,
                                                                nosepoincareandersen,
                                                                forceloaded,
								fout, ferr
								) );
      }
      else{
        simplemd = SimpleMDHandle( new NosePoincareAndersenMD(
                                                               absoluteTime,
                                                               dt,
                                                               molecules,
                                                               unit,
                                                               rc,
                                                               moldict,
							       plugins,
                                                               nosepoincareandersen,
                                                               forceloaded,
							       fout,
							       ferr
							       ) );
        
      }
    }
    else{
      simplemd = SimpleMDHandle( new SimpleMD(
                                               absoluteTime,
                                               dt,
                                               molecules,
                                               unit,
                                               rc,
                                               moldict,
                                               plugins,
					       fout,
					       ferr
					       ) );
                                               
    }
  }



  virtual ~Sample()
  {
    delete unit;
    delete box;

    if ( redirected ){
      delete ferr;
      delete fout;
    }
    if ( fmpi ){
      fakempi_finalize( fmpi );
    }
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
  //継続情報に含めない設定は、コマンドラインで指定すべき。
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
      else if ( tag == "@TSTO" ){
        fgets( buf, sizeof(buf), input );
        int interval = atoi( buf );
        // Interval plugin containing StopDrift plugin.
        plugins.push_back( ProcessPluginHandle( new Interval( interval, new StopDrift() ) ) );
      }
      else if ( tag == "@NOPO" ){
	assert( nosepoincare == 0 );
        nosepoincare = new NosePoincare( 0, 0, 0 );
        // read variables for nose-poincare integrator
        nosepoincare->ReadNOPO( *unit, input );
      }
      else if ( tag == "@NOPA" ){
	assert( nosepoincareandersen == 0 );
	//cerr << "@NOPA" << endl;
        nosepoincareandersen = new NosePoincareAndersen( 0, 0, 0 );
        // read variables for nose-poincare-andersen integrator
        nosepoincareandersen->ReadNOPA( *unit, input );
      }
      else if ( tag == "@NPAZ" ){
	assert( nosepoincareandersen == 0 );
        nosepoincareandersen = new NosePoincareAndersen( 0, 0, 0 );
        // read variables for nose-poincare-andersen-z integrator
        nosepoincareandersen->ReadNOPA( *unit, input );
        zonly = 1;
      }
      else if ( tag == "@SNPH" ){
	assert( nosepoincareandersen == 0 );
	//Sturgeon's NPH
        nosepoincareandersen = new NosePoincareAndersen( 0, 0, 0 );
        nosepoincareandersen->ReadSNPH( *unit, input );
      }
      else if ( tag == "@BOX3" ) {
        const int isP[3] = {1,1,1};
        box = new SemiperiodicBox( isP );
        box->ReadBOX3( *unit, input );
      }
      else if ( tag == "@BXLA" ) {
        const int isP[3] = {1,1,1};
        box = new SemiperiodicBox( isP );
        box->ReadBXLA( *unit, input );
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
        moldict.insert(
                       map<string, FlexibleHandle>::value_type(
                                                                  current_id08,
                                                                  mp ));
      }
      else if ( tag == "@SNPO" ){
	fgets( buf, sizeof(buf), input );
	int interval = atoi( buf );
        // interval plugin containing snapshot plugin
        plugins.push_back( ProcessPluginHandle( new Interval( interval, new SnapShot() ) ) );
      }
      else if ( tag == "@SNP1" ){
	fgets( buf, sizeof(buf), input );
	int interval = atoi( buf );
        fgets( buf, sizeof(buf), input );
        string id08 = string( buf );
        id08.resize( 8, ' ' );
        // interval plugin containing snapshot plugin
        plugins.push_back( ProcessPluginHandle( new Interval( interval, new SnapShotOneComp( id08 ) ) ) );
      }
      else if ( tag == "@NLOG" )  {
	fgets( buf, sizeof(buf), input );
	int interval = atoi( buf );
        // interval plugin containing logging plugin
        plugins.push_back( ProcessPluginHandle( new Interval( interval, new Log() ) ) );
      }
      else if ( tag == "@ID08" ) {
        fgets( buf, sizeof(buf), input );
        current_id08 = string( buf );
        current_id08.resize( 8, ' ' );
      }
      else if ( tag == "@ZUMB" )  {
        // Z-axis umbrella        
	fgets( buf, sizeof( buf ), input );
	double distance, fconst;
        sscanf( buf, "%lf %lf", &distance, &fconst );
	distance *= Angstro * unit->m;
	fconst   *= kilo * unit->J / mol / (Angstro * unit->m) / (Angstro * unit->m);
	//cerr << distance << "#" << fconst << endl;
        plugins.push_back( ProcessPluginHandle( new Zumbrella( distance, fconst, current_id08 ) ) );
      }
/*
      else if ( tag == "@FIXB" )  {
				//@FIXB is insufficiently implemented @FIXC (not taking constraints into account for temperature calculations.)
	//FixPluginを使うと、水2成分のうち片方を固定することができない。成分番号で指定するか、
				fgets( buf, sizeof( buf ), input );
				int isfix = atoi(buf);
				if ( isfix ){
	  			plugins.push_back( ProcessPluginHandle( new FixPlugin( current_id08 ) ) );
				}
      }
*/
      else if ( tag == "@FIXc" )  {
				//@FIXc is insufficiently implemented @FIXC (not taking constraints into account for temperature calculations.)
				//現在の成分のみを固定。
				fgets( buf, sizeof( buf ), input );
				int isfix = atoi(buf);
				if ( isfix ){
					//next component will be fixed.
	  			plugins.push_back( ProcessPluginHandle( new FixComponentPlugin( coordcount+1 ) ) );
				}
      }
      else if ( tag == "@INNR" ) {
	fgets( buf, sizeof(buf), input );
	innerloop = atoi( buf );
      }
      /*
      else if ( tag == "@NRXG" ) {
        assert( balancer != 0 );

	fgets( buf, sizeof(buf), input );
	balancer->ntrial = atoi( buf );
        }*/
      else if ( tag == "@MPIN" ) {
	fgets( buf, sizeof(buf), input );
	int nprocs = atoi( buf );
	fmpi = ::fakempi_new( nprocs );
	for( int i=0; i<nprocs; i++ ){
	  fgets( buf, sizeof(buf), input );
	  if ( i == fmpi->myrank ){
            //remove cr
            buf[strlen(buf)-1]='\0';
	    infilename = string( buf );
            string outfilename = infilename + ".out";
            string errfilename = infilename + ".err";
            fout = new ofstream( outfilename.c_str(), ios::out );
            ferr = new ofstream( errfilename.c_str(), ios::out );
	    fout->precision(17);
	    ferr->precision(17);
	    redirected = 1;
	  }
	}
	assert( innerloop != 0 );
        int count     = nstep / innerloop;
	//balancer = new LoadBalancer( fmpi, innerloop, count, replica );
	// child plugin will be inserted later. 2006-6-21
	//balancer = new LoadBalancer( fmpi, innerloop, count );
	balancer = new LoadBalancer2( fmpi, innerloop, count );
        // loadbalancer plugin containing replica exchange plugin
	plugins.push_back( ProcessPluginHandle( balancer ) );
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
          mm->ReadAR3A( nmol, *box, *unit, input );
          forceloaded = 0;
        }
        else if ( tag == "@NX4B" ){
          RigidBodies* mm = new RigidBodies( ph );
          mh = MolsHandle( mm );
          mm->ReadNX4A( nmol, *box, *unit, input );
          forceloaded = 0;
        }
        else if ( tag == "@NX4A" ){
          RigidBodies2* mm = new RigidBodies2( ph );
          mh = MolsHandle( mm );
          mm->ReadNX4A( nmol, *box, *unit, input );
          forceloaded = 0;
        }
        else if ( tag == "@WTG5" ){
          RigidBodies* mm = new RigidBodies( ph );
          mh = MolsHandle( mm );
          mm->ReadWTG5( nmol, *unit, input );
        }
        else if ( tag == "@WTG6" ){
          RigidBodies2* mm = new RigidBodies2( ph );
          mh = MolsHandle( mm );
          mm->ReadWTG6( nmol, *unit, input );
        }
        else if ( tag == "@ATG5" ){
          mh = MolsHandle( new MonatomicMols( ph, nmol, *unit, input ) );
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
          molecules = new MolCollection( cell );
        }
        molecules->push_back( ph, mh );
      }
      else{
        if ( buf[0] == '@' ){
          cerr << buf;
        }
      }
    }


    //interconnection between plugins
    if ( balancer ){
      if ( nosepoincare ){
	int ntrial = 1;
	ProcessPlugin* replica;
	replica = new Replica( fmpi, 5678, ntrial, unit, nosepoincare );
	balancer->set( replica );
	//cerr << "NOPO" << endl;
      }
      else if ( nosepoincareandersen ){
	int ntrial = 1;
	ProcessPlugin* replica;
	replica = new CPReplica( fmpi, 5678, ntrial, unit, nosepoincareandersen );
	balancer->set( replica );
	//cerr << "NOPA" << endl;
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



