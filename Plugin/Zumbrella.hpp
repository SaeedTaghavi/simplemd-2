#ifndef ZUMBRELLA_HPP
#define ZUMBRELLA_HPP


#include "Vector3.hpp"
#include "Plugin/CollectorPlugin.hpp"

class MolCollection;
class Cell;
class Mols;
//class MonatomicMols;
class SingleMolEntity;

class ZumbrellaHelper : public CollectorPlugin{
  double distance;
  double fconst;
  string id08;
public:
  ZumbrellaHelper( double distance_, double fconst_, string id08_ ) :
    distance( distance_ ), fconst( fconst_ ), id08( id08_ ) {};
  void Execute2( Mols* mols );
  void Execute3( SingleMolEntity* mol );
};


#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;

class Zumbrella : public ProcessPlugin
{
  SimpleMD* system;
  double    distance;
  double    fconst;
  string    id08;
public:
  explicit Zumbrella( double distance_, double fconst_, string id08_ );
  void initialize( SimpleMD* sys );
  void forcepatch();
  ~Zumbrella()
  {
    cerr << "~Plugin::Zumbrella" << endl;
  }
};


#endif
