#include "Plugin/Terminate.hpp"
#include "System/SimpleMD.hpp"

Terminate::Terminate( int s ) : steps( s )
{
}



void
Terminate::initialize( SimpleMD* sys )
{
  system = sys;
}



bool
Terminate::running()
{
  return ( system->GetStep() +1 < steps );
}
