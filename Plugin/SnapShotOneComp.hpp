#ifndef SNAPSHOTONECOMP_HPP
#define SNAPSHOTONECOMP_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;

class SnapShotOneComp : public ProcessPlugin
{
  SimpleMD* system;
  string id08;
  ostream* fout;
  ostream* ferr;
  int format;
public:
  explicit SnapShotOneComp( ostream* out, ostream* err, const string& id08, int fmt );
  explicit SnapShotOneComp( const string& id08, int fmt );
  void initialize( SimpleMD* sys );
  bool running();
  ~SnapShotOneComp()
  {
    cerr << "~Plugin::SnapShotOneComp" << endl;
  }
};


#endif
