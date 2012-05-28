#include "Plugin/CollectorPlugin.hpp"
#include "Cell/Cell.hpp"
#include "Mols/Mols.hpp"
#include "Mols/MonatomicMols.hpp"
#include "Mols/RigidBodies.hpp"
#include "MolCollection.hpp"
#include "Vector3.hpp"



void CollectorPlugin::Execute0( MolCollection* coll )
{
  coll->Execute( *this );
}



void CollectorPlugin::Execute1( Cell* cell )
{
  cell->Execute( *this );
}



void CollectorPlugin::Execute2( Mols* mols )
{
  mols->Execute( *this );
}


void Test::Execute3( SingleMolEntity* mol )
{
  cerr << mol->GetOrder() << " TEST:ORDER" << endl;
}

  
void
QDoFPlugin::Execute2( Mols* mols )
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
AddVelocityPlugin::Execute3( SingleMolEntity* mol )
{
  mol->AddVelocity( velocity );
}



void
ScaleVelocityPlugin::Execute3( SingleMolEntity* mol )
{
  mol->ScaleVelocity( ratio );
}



void
ScaleVelocity2Plugin::Execute3( SingleMolEntity* mol )
{
  mol->ScaleVelocity( ratio );
}




void
ScalePositionPlugin::Execute3( SingleMolEntity* mol )
{
  mol->ScalePosition( ratio );
}



/*
void
SerializePlugin::Execute( SingleMolEntity* mol )
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
