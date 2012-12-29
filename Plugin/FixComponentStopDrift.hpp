#ifndef FIXCOMPONENTSTOPDRIFT_HPP
#define FIXCOMPONENTSTOPDRIFT_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;

class MolCollection;
class Cell;
class Mols;
//class MonatomicMols;
class SingleMolEntity;

class FixComponentStopDriftHelper : public CollectorPlugin{
	int compo;
public:
  FixComponentStopDriftHelper( int c ) : compo(c) {};
  void HookL2( Cell* cell );
  void HookL0( SingleMolEntity* mol );
};


class FixComponentStopDriftPlugin : public ProcessPlugin
{
  SimpleMD* system;
  int compo;
public:
  explicit FixComponentStopDriftPlugin( int c ): compo( c ) {};
  void initialize( SimpleMD* sys );
  void ForceHook();
  ~FixComponentStopDriftPlugin()
  {
    cerr << "~Plugin::FixComponentStopDrift" << endl;
  }
};


#endif
