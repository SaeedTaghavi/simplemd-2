#include "Plugin/StopDrift.hpp"
#include "System/SimpleMD.hpp"

StopDrift::StopDrift()
{
}



void
StopDrift::initialize( SimpleMD* sys )
{
  system = sys;
}



bool
StopDrift::running()
{
  system->GetMolCollection()->StopDrift();
  return true;
}
