#include "Plugin/Zumbrella.hpp"
#include "System/SimpleMD.hpp"
#include "Cell/Cell.hpp"
#include "Mols/Mols.hpp"
#include "Mols/MonatomicMols.hpp"
#include "Mols/RigidBodies.hpp"
#include "MolProperty.hpp"
#include "MolCollection.hpp"
#include "Vector3.hpp"

/*
 *This plugin is called from Zumbrella.
 *うんどうにこうそくをあたえるので、おんどのけいさんにえいきょうするはず。いまはおんどのほせいはおこなわない。
 */

void
ZumbrellaHelper::Execute2( Mols* mols )
{
  MolPropertyHandle property = mols->GetProperty();
  if ( property->id08 == id08 ){
    MonatomicMols* mmh = dynamic_cast<MonatomicMols*> ( mols );
    assert( mmh );
    int nmol = mols->Size();
    assert( nmol == 1 );
    //単原子分子で、しかもid08が一致する場合に限り、下位を呼びだす。
    mols->Execute( *this );
  }
  //それ以外の場合は下位に降りない。
}



void
ZumbrellaHelper::Execute3( SingleMolEntity* mol )
{
  //Z方向には外力を加え、XY方向は力を消す。
  MonatomicMol* mmh = dynamic_cast<MonatomicMol*> ( mol );
  //This test is unnecessary.
  if ( mmh ){
    const Vector3& force=mmh->GetForce();
    //Add external force here ????
    const Vector3& coord=mmh->Position();
    //cerr << force.z << ":" << fconst << ":" << coord.z << ":" << distance << endl;
    mmh->center.force.Set(0,0,force.z - fconst * ( coord.z - distance ) );
  }
}
Zumbrella::Zumbrella( double distance_, double fconst_, string id08_ ) : distance( distance_ ), fconst( fconst_ ), id08( id08_ )
{
}


void
Zumbrella::forcepatch()
{
  //ここで、CollectorPluginを初期化し、使いすてる。
  ZumbrellaHelper zumb( distance, fconst, id08 );
  zumb.Execute0( system->GetMolCollection() );
}


void
Zumbrella::initialize( SimpleMD* sys )
{
  system = sys;
}

