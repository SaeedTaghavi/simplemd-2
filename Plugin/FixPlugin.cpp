#include "Plugin/CollectorPlugin.hpp"
#include "System/SimpleMD.hpp"
#include "Plugin/FixPlugin.hpp"
#include "Cell/Cell.hpp"
#include "Mols/Mols.hpp"
#include "Mols/MonatomicMols.hpp"
#include "Mols/RigidBodies.hpp"
#include "MolProperty.hpp"
#include "MolCollection.hpp"
#include "Vector3.hpp"

/*
 *ぶんしにくわわるちからを0にすることで、ぶんしのへいしんをとめるぷらぐいん。とうめんたんげんしぶんしのみたいおう。
 *固定した分子の間の、力の計算は行われるので、計算が速くはならない。
 */

void
FixHelper::HookL1( Mols* mols )
{
  MolPropertyHandle property = mols->GetProperty();
  if ( property->id08 == id08 ){
    MonatomicMols* mmh = dynamic_cast<MonatomicMols*> ( mols );
    assert( mmh );
    //id08が一致する場合に限り、下位を呼びだす。
    mols->PluginHookL0( *this );
  }
  //それ以外の場合は下位に降りない。
}



void
FixHelper::HookL0( SingleMolEntity* mol )
{
  //Z方向には外力を加え、XY方向は力を消す。
  MonatomicMol* mmh = dynamic_cast<MonatomicMol*> ( mol );
  assert(mmh);
  mmh->center.force.Set(0,0,0);
}



void
FixPlugin::forcepatch()
{
  //ここで、CollectorPluginを初期化し、使いすてる。
  FixHelper fixhelper( id08 );
  fixhelper.HookL3( system->GetMolCollection() );
}


void
FixPlugin::initialize( SimpleMD* sys )
{
  system = sys;
}

