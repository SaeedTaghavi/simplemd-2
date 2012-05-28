#include "Plugin/Interval.hpp"
#include "System/SimpleMD.hpp"

Interval::Interval( int intv, ProcessPlugin* p ) : interval( intv ), plugin( p )
{
}



void
Interval::initialize( SimpleMD* sys )
{
  system = sys;
  plugin->initialize(sys);
}



bool
Interval::running()
{
  if ( system->GetStep() % interval == 0 ){
    return plugin->running();
  }
  return true;
}



void
Interval::terminate()
{
  plugin->terminate();
}
