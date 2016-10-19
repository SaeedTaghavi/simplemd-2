#include <cassert>
#include "Plugin/CollectorPlugin.hpp"
#include "System/SimpleMD.hpp"
#include "Plugin/FixPlugin.hpp"
#include "Cell/Cell.hpp"
#include "Mols/Mols.hpp"
#include "Mols/MonatomicMols.hpp"
#include "Mols/RigidBodies.hpp"
#include "Mols/RigidBodies2.hpp"
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
    mols->PluginHookL0( *this );
  }
  //cout << "FixHelper::HookL1" << endl;
}



void
FixHelper::HookL0( SingleMolEntity* mol )
{
  MonatomicMol* mmh = dynamic_cast<MonatomicMol*> ( mol );
  if ( mmh ) {
    mmh->center.force.Set(0,0,0);
  }
  else{
    RigidBody2* mmh = dynamic_cast<RigidBody2*> ( mol );
    assert(mmh);
    //Hook position is moved from AFTER PostForce to BEFORE PostForce.
    //That is, not the collective force and torque but the force on the atom must be changed.
    //mmh->com.center.force.Set(0,0,0);
    //mmh->torque.Set(0,0,0);
    int nsite = mmh->atom.size();
    for(int site=0;site<nsite;site++ ){
      mmh->atom[site].center.force.Set(0,0,0);
    }
  }
  //  cout << "FixHelper::HookL0" << endl;
}



void
FixPlugin::ForceHook()
{
  //ここで、CollectorPluginを初期化し、使いすてる。
  FixHelper fixhelper( id08 );
  //cout << ":" << id08 << ":" << endl;
  fixhelper.HookL3( system->GetMolCollection() );
}


void
FixPlugin::initialize( SimpleMD* sys )
{
  system = sys;
}

