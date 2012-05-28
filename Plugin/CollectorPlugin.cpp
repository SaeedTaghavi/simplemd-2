#include "Plugin/CollectorPlugin.hpp"
#include "Cell/Cell.hpp"
#include "Mols/Mols.hpp"
#include "Mols/MonatomicMols.hpp"
#include "Mols/RigidBodies.hpp"
#include "MolCollection.hpp"
#include "Vector3.hpp"



void CollectorPlugin::HookL3( MolCollection* coll )
{
  coll->PluginHookL2( *this );
}



void CollectorPlugin::HookL2( Cell* cell )
{
  cell->PluginHookL1( *this );
}



void CollectorPlugin::HookL1( Mols* mols )
{
  mols->PluginHookL0( *this );
}


void Test::HookL0( SingleMolEntity* mol )
{
  cerr << mol->GetOrder() << " TEST:ORDER" << endl;
}

  
void
QDoFPlugin::HookL1( Mols* mols )
{
  MonatomicMols* mmh = dynamic_cast<MonatomicMols*> ( mols );
  if ( mmh ){
    int    nmol = mols->Size();
    dof += nmol * 3;
  }
  else{
    RigidBodies* mmh = dynamic_cast<RigidBodies*> ( mols );
    assert( mmh );
    int    nmol = mols->Size();
    dof += nmol * 7;
  }
}



void
AddVelocityPlugin::HookL0( SingleMolEntity* mol )
{
  mol->AddVelocity( velocity );
}



void
ScaleVelocityPlugin::HookL0( SingleMolEntity* mol )
{
  mol->ScaleVelocity( ratio );
}



void
ScaleVelocity2Plugin::HookL0( SingleMolEntity* mol )
{
  mol->ScaleVelocity( ratio );
}




void
ScalePositionPlugin::HookL0( SingleMolEntity* mol )
{
  mol->ScalePosition( ratio );
}



/*
void
SerializePlugin::Hook( SingleMolEntity* mol )
{
  MonatomicMols* mmh = dynamic_cast<MonatomicMol*> ( mol );
  if ( mmh ){
    int    nmol = mols->Size();
    dof += nmol * 3;
  }
  else{
    RigidBodies* mmh = dynamic_cast<RigidBody*> ( mol );
    assert( mmh );
    int    nmol = mols->Size();
    dof += nmol * 7;
  }
}
*/
