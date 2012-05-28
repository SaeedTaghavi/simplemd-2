#include <fstream>
#include "Plugin/Replica.hpp"
#include "System/SimpleMD.hpp"
#include "Plugin/CollectorPlugin.hpp"

Replica::Replica( sFakeMPI* f, int seed, int nt, Unit* const u, NosePoincare* np ) : fmpi( f ), ntrial( nt ), nosepoincare( np ), fout( 0 ), ferr( 0 ), alter( false ), unit( u )
{
  if ( fmpi->myrank == 0 ){
    random = RandomD01Handle( new RandomD01( seed ) );
    node_in_charge = vector<int>( fmpi->nprocs, 0 );
  }
  replica = vector<Exchange>( fmpi->nprocs, Exchange( random ));
}



void
Replica::initialize( SimpleMD* sys )
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
Replica::running()
{
  //replica exchange here
  //collect energy and beta
  double kT = nosepoincare->GetkT();
  const PotVir& pv = system->GetPv();
  double U  = pv.GetEp();
  double lastkT = kT;
  double kTs[fmpi->nprocs];
  double Us[fmpi->nprocs];
  fakempi_gather( fmpi, 8, &kT, (char*)kTs );
  fakempi_gather( fmpi, 8, &U,  (char*)Us );
  //exchange 
  int nprocs = fmpi->nprocs;
  if ( fmpi->myrank == 0 ){
    for( int i=0; i<nprocs; i++ ){
      vector<double> vars;
      vector<double> consts;
      vars.push_back( kTs[i] );
      consts.push_back( Us[i] );
      replica[i].Set( consts, vars );
    }
    //今のところ、隣としか交換しないのだが、これでは何かと問題多い。
    //やはり、ランダムな交換をゆるし、交換頻度を上げた方がいい。
    /*
    int i = 0;
    if ( alter ){
      i = 1;
    }
    alter = ! alter;
    for( ; i<nprocs-1; i+=2 ){
    */
    for( int rep=0; rep < 2*nprocs; rep++ ){
      int i = random->get() * nprocs;
      int j = random->get() * nprocs;
      if ( i != j ){
        int p1 = node_in_charge[i];
        int p2 = node_in_charge[j];
        //int p2 = node_in_charge[i+1];
        cerr << p1 << "x" << p2 << " trial.\n";
        bool accept = replica[p1].Trial( replica[p2] );
        if ( accept ){
          replica[p1].exchange( replica[p2] );
          node_in_charge[i] = p2;
          node_in_charge[j] = p1;
          //node_in_charge[i+1] = p1;
          *ferr << p1 << "x" << p2 << "exchanged.\n";
        }
        //多分これはループ外でいいはず。
        /*
          for( int i=0; i<nprocs; i++ ){
          kTs[i] = replica[i].GetkT();
          Us[i]  = replica[i].GetEnergy();
          }
        */
      }
    }
    for( int i=0; i<nprocs; i++ ){
      const vector<double>& vars = replica[i].GetVars();
      kTs[i] = vars[0];
    }
  }
  fakempi_scatter( fmpi, 8, &kT, (char*)kTs );
  if ( lastkT != kT ){
    *ferr << lastkT << "->" << kT << " kT changed.\n";
    *fout << "@MDRT" << endl;
    *fout << kT / unit->K() << endl;
    nosepoincare->SetkT( kT );

    //double r = sqrt( kT / lastkT );
    //Vector3 ratio( r, r, r );
    //system->GetMolCollection()->ScaleVelocity( ratio );

    //Reset S and its derivetives
    double r = sqrt( kT / lastkT );
    double s = nosepoincare->S();
    Vector3 ratio( r/s, r/s, r/s );
    MolCollection* mol = system->GetMolCollection();
    ScaleVelocityPlugin p( ratio );
    p.Execute0( mol );
    //mol->ScaleVelocity( ratio );
    nosepoincare->SetS( 1 );
    //nosepoincare->SetH0( 0 );
    nosepoincare->SetPs( 0 );
    nosepoincare->ResetH0( system->GetEk(), pv.GetEp(), mol->GetDoF() );
  }
  return true;
}



