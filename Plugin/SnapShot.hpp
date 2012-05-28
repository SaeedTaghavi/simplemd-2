#ifndef SNAPSHOT_HPP
#define SNAPSHOT_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;

class SnapShot : public ProcessPlugin
{
  SimpleMD* system;

  ostream& os;
public:
  explicit SnapShot( ostream& out=cerr ) : os(out) {};
  void initialize( SimpleMD* sys ){system = sys;};
  bool running();
  ~SnapShot()
  {
    cerr << "~Plugin::SnapShot" << endl;
  }
};


#endif
