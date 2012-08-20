#include "Plugin/SnapShot.hpp"
#include "System/SimpleMD.hpp"







bool
SnapShot::running()
{
  system->SnapShot(os, mode);
  return true;
}



