#include <cstdio>
#include "MolCollection.hpp"
#include "debug.hpp"
#include "Cell/Cell.hpp"
#include "Plugin/CollectorPlugin.hpp"
#include "Plugin/PairProcessPlugin.hpp"




MolCollection::MolCollection( Cell* c )
{
    cell = c;
    //dof    = 0;
    //prop   = new MolPropertyArray( ncompo );
    combi  = Combination();
    mesg("MolCollection(Cell*)\n");
    //cerr << "MolCollection" << endl;
}




void
MolCollection::push_back( const FlexibleHandle& ph, const MolsHandle& mh )
{
  int compo = combi.push_back( ph );
  cell->Set( compo-1, mh );
}



MolCollection::~MolCollection()
{
  mesg("~MolCollection\n");
  delete cell;
}



int MolCollection::GetDoF()
{
  int dof = 0;
  const FlexibleArray& prop = combi.GetProperty();
  int ncompo = prop.size();
  for(int compo=0; compo<ncompo; compo++ ){
    int size = cell->GetSize( compo );
    dof += prop[compo]->GetDoF() * size;
    printf("%d dof %d nmol\n", dof, size);
  }
  return dof;
}



int
MolCollection::Force( const Truncation& rc, PotVir &pv )
{
  //mesg("force2\n");
  //cell->Force( rc, pv );!!!
  return cell->Force( combi, rc, pv );
}



int
MolCollection::PairProcess( PairProcessPlugin& p, const Truncation& rc ) const
{
  p.Monitor1( *cell, combi, rc );
  return 1;
}



void 
MolCollection::Preforce()
{
  cell->Preforce();
}



void 
MolCollection::Postforce()
{
  cell->Postforce();
}




/*
void 
MolCollection::ScaleVelocity( const Vector3& r )
{
  ScaleVelocityPlugin p( r );
  p.Hook1( cell );
  //cell->ScaleVelocity( r );
}
*/



void 
MolCollection::Expand( const Vector3& r )
{
  Vector3 ri( 1/r.x, 1/r.y, 1/r.z );
  ScaleVelocityPlugin p( ri );
  //cell->ScaleVelocity( ri );
  p.HookL2( cell );
  //cell->ScalePosition( r );

  //座標の拡張
  cell->Rescale( r );
}



double 
MolCollection::GetEk()
{
  return cell->GetEk();
}



void
MolCollection::GetEkt( Vector3& ekt ) const
{
  cell->GetEkt( ekt );
}



void 
MolCollection::Write( const Unit& unit, double dt, ostream& to )
{
  const FlexibleArray& prop = combi.GetProperty();
  //Write number of components
  int ncompo = cell->GetNumCompo();
  to << "@NCMP\n"
     << ncompo << '\n';
  //Write grid (if any)
  cell->Write( unit, to );
  //Write id and coordinates
  //実際には、cellの中に分子がはいっているので、Cellに出力させる。
  //単一のMolsArrayに、Gridの全分子を集約し、ソートしてから出力する。
  for(int compo=0; compo<ncompo; compo++){
    const FlexibleHandle& fh = prop[compo];
    to << "@ID08\n"
       << fh->id08 << '\n';
    MolsHandle mols = cell->PullAll( compo );
    mols->Write( unit, to );
  }
}



void 
MolCollection::SnapShot( const Unit& unit, double dt, ostream& to )
{
  const FlexibleArray& prop = combi.GetProperty();
  //Write number of components
  int ncompo = cell->GetNumCompo();
  to << "@NCMP\n"
     << ncompo << '\n';
  cell->Write( unit, to );
  for(int compo=0; compo<ncompo; compo++){
    const FlexibleHandle& fh = prop[compo];
    to << "@ID08\n"
       << fh->id08 << '\n';
    MolsHandle mols = cell->CopyAll( compo );
    mols->SnapShot( unit, to );
  }
}



void MolCollection::ProgressPosition( double dt )
{
  cell->ProgressPosition( dt );
}



void MolCollection::ProgressMomentum( double dt )
{
  cell->ProgressMomentum( dt );
}



void MolCollection::StopDrift()
{
  Vector3 momentum;
  cell->TotalMomentum( momentum );
  double mass  = cell->GetMass();
  double coeff = -1 / mass;
  momentum.x  *= coeff;
  momentum.y  *= coeff;
  momentum.z  *= coeff;
  cout << momentum.x
       << " " << momentum.y
       << " " << momentum.z << " drift\n";

  //Add vecoloty to all molecules via plugin.
  AddVelocityPlugin p( momentum );
  p.HookL2( cell );
}



void
MolCollection::PluginHookL2( CollectorPlugin& plugin )
{
  plugin.HookL2( cell );
}




void MolCollection2::unserialize( double* const p )
{
  cell->unserialize( p );
}



void MolCollection2::serialize( double* p ) const
{
  cell->serialize( p );
}



void MolCollection2::serializeforce( double* xi ) const
{
  cell->serializeforce( xi );
}



int MolCollection2::qdof() const
{
  return cell->qdof();
}



