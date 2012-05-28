#include "Plugin/TReset.hpp"
#include "System/SimpleMD.hpp"
#include "Plugin/CollectorPlugin.hpp"

TReset::TReset( double temp ) : bathTemp( temp )
{
}



void
TReset::initialize( SimpleMD* sys )
{
  system = sys;
}



bool
TReset::running()
{
  MolCollection* mol = system->GetMolCollection();
  double ek = mol->GetEk();
  double r  = sqrt( bathTemp / system->Temperature( ek ) );
  Vector3 ratio( r, r, r );
  ScaleVelocityPlugin p( ratio );
  p.Execute0( mol );
  //mol->ScaleVelocity( ratio );
  return true;
}



void
TReset::terminate()
{
  //do nothing.
}
