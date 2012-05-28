#ifndef FIXCOMPONENTPLUGIN_HPP
#define FIXCOMPONENTPLUGIN_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;

class MolCollection;
class Cell;
class Mols;
//class MonatomicMols;
class SingleMolEntity;

class FixComponentHelper : public CollectorPlugin{
	int compo;
public:
  FixComponentHelper( int c ) : compo(c) {};
  void HookL2( Cell* cell );
  void HookL0( SingleMolEntity* mol );
};


class FixComponentPlugin : public ProcessPlugin
{
  SimpleMD* system;
  int compo;
public:
  explicit FixComponentPlugin( int c ): compo( c ) {};
  void initialize( SimpleMD* sys );
  void ForceHook();
  ~FixComponentPlugin()
  {
    cerr << "~Plugin::FixComponentPlugin" << endl;
  }
};


#endif
