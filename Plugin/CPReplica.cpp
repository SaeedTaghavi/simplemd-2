#include <fstream>
#include "Plugin/CPReplica.hpp"
#include "System/SimpleMD.hpp"
#include "Plugin/CollectorPlugin.hpp"

CPReplica::CPReplica( sFakeMPI* f, int seed, int nt, Unit* const u,
		      NosePoincareAndersen* npa ) : fmpi( f ), ntrial( nt ), nosepoincareandersen( npa ), fout( 0 ), ferr( 0 ), alter( false ), unit( u )
{
  if ( fmpi->myrank == 0 ){
    random = RandomD01Handle( new RandomD01( seed ) );
    node_in_charge = vector<int>( fmpi->nprocs, 0 );
  }
  replica = vector<ExchangePT>( fmpi->nprocs, ExchangePT( random ));
}



void
CPReplica::initialize( SimpleMD* sys )
{
  system = sys;
  if ( fout == 0 ) fout = system->GetFout();
  if ( ferr == 0 ) ferr = system->GetFerr();

  if ( fmpi->myrank == 0 ){
    for( int i=0; i< fmpi->nprocs; i++ ){
      node_in_charge[i] = i;
    }
  }
}



bool
CPReplica::running()
{
  //replica exchange here
  //collect energy and beta
  const PotVir& pv = system->GetPv();
  MolCollection* const mol = system->GetMolCollection();
  double volume = mol->Volume();
  double Energy = pv.GetEp();
  double kT     = nosepoincareandersen->GetkT();
  double lastkT = kT;
  double pext   = nosepoincareandersen->GetPext();
  double lastpext = pext;
  double kTs[fmpi->nprocs];
  double Pexts[fmpi->nprocs];
  double Energies[fmpi->nprocs];
  double Volumes[fmpi->nprocs];
  fakempi_gather( fmpi, 8, &kT, (char*)kTs );
  fakempi_gather( fmpi, 8, &Energy,  (char*)Energies ); 
  fakempi_gather( fmpi, 8, &volume,  (char*)Volumes ); 
  fakempi_gather( fmpi, 8, &pext, (char*)Pexts );
  //exchange 
  int nprocs = fmpi->nprocs;
  if ( fmpi->myrank == 0 ){
    for( int i=0; i<nprocs; i++ ){
      vector<double> vars;
      vector<double> consts;
      vars.push_back( kTs[i] );
      vars.push_back( Pexts[i] );
      consts.push_back( Energies[i] );
      consts.push_back( Volumes[i] );
      replica[i].Set( consts, vars );
      //cerr << kTs[i] << " kt " << i << endl;
    }
    //今のところ、隣としか交換しないのだが、これでは何かと問題多い。
    //やはり、ランダムな交換をゆるし、交換頻度を上げた方がいい。
    /*
    int i = 0;
    if ( alter ){
      i = 1;
    }
    alter = ! alter;
    */
    for( int rep=0; rep < 2*nprocs; rep++ ){
      int i = random->get() * nprocs;
      int j = random->get() * nprocs;
      if ( i != j ){
        int p1 = node_in_charge[i];
        int p2 = node_in_charge[j];
        cerr << p1 << "x" << p2 << " trial.\n";
        bool accept = replica[p1].Trial( replica[p2] );
        if ( accept ){
          replica[p1].exchange( replica[p2] );
          node_in_charge[i] = p2;
          node_in_charge[j] = p1;
          *ferr << p1 << "x" << p2 << "exchanged.\n";
        }
        //多分これはループ外でいいはず。
        /*
          for( int i=0; i<nprocs; i++ ){
          kTs[i] = replica[i].GetkT();
          Enthalpys[i]  = replica[i].GetEnergy();
          }
        */
      }
    }
    for( int i=0; i<nprocs; i++ ){
      const vector<double>& vars = replica[i].GetVars();
      kTs[i] = vars[0];
      Pexts[i] = vars[1];
      cerr << vars[0] << " " << vars[1] << endl;
    }
  }
  fakempi_scatter( fmpi, 8, &pext, (char*)Pexts );
  fakempi_scatter( fmpi, 8, &kT, (char*)kTs );
  if ( lastpext != pext ){
    nosepoincareandersen->SetPext( pext );
    nosepoincareandersen->SetPv( 0 );
    *ferr << lastpext << "->" << pext << " p changed.\n";
  }
  if ( lastkT != kT ){
    *ferr << lastkT << "->" << kT << " kT changed.\n";
    *fout << "@MDRT" << endl;
    *fout << kT / unit->K() << endl;
    nosepoincareandersen->SetkT( kT );

    //double r = sqrt( kT / lastkT );
    //Vector3 ratio( r, r, r );
    //system->GetMolCollection()->ScaleVelocity( ratio );
    
    //Reset S and its derivetives
    double r = sqrt( kT / lastkT );
    double s = nosepoincareandersen->S();
    MolCollection* mol = system->GetMolCollection();
    ScaleVelocity2Plugin p( r/s );
    p.Execute0( mol );
    //mol->ScaleVelocity( ratio );
    nosepoincareandersen->SetS( 1 );
    //nosepoincareandersen->SetH0( 0 );
    nosepoincareandersen->SetPs( 0 );
  }
  if ( lastpext != pext || kT != lastkT ){
    const PotVir& pv = system->GetPv();
    nosepoincareandersen->SetH0( nosepoincareandersen->RecalcH0( system->GetEk(), pv.GetEp(), mol->GetDoF(), mol->Volume() ) );
    cerr << "Reset H0 to " << nosepoincareandersen->GetH0() << endl;
  }
  return true;
}


