#ifndef SNAPSHOT_HPP
#define SNAPSHOT_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;

class SnapShot : public ProcessPlugin
{
  SimpleMD* system;

  ostream& os;
  int mode;
public:
  explicit SnapShot( ostream& out=cerr, int mode_=0 ) : os(out), mode(mode_) {};
  void initialize( SimpleMD* sys ){system = sys;};
  bool running();
  ~SnapShot()
  {
    cerr << "~Plugin::SnapShot" << endl;
  }
};


#endif
