#ifndef FIXPLUGIN_HPP
#define FIXPLUGIN_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;

class MolCollection;
class Cell;
class Mols;
//class MonatomicMols;
class SingleMolEntity;

class FixHelper : public CollectorPlugin{
  string id08;
public:
  FixHelper( string id08_ ) : id08( id08_ ) {};
  void Execute2( Mols* mols );
  void Execute3( SingleMolEntity* mol );
};


class FixPlugin : public ProcessPlugin
{
  SimpleMD* system;
  string    id08;
public:
  explicit FixPlugin( string id08_ ): id08( id08_ ) {};
  void initialize( SimpleMD* sys );
  void forcepatch();
  ~FixPlugin()
  {
    cerr << "~Plugin::FixPlugin" << endl;
  }
};


#endif
