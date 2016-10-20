#ifndef MAGNETICFIELDPLUGIN_HPP
#define MAGNETICFIELDPLUGIN_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;

//class MolCollection;
//class Cell;
//class Mols;
//class MonatomicMols;
class SingleMolEntity;

class MagneticFieldHelper : public CollectorPlugin{
  const Vector3& B;
  MolPropertyHandle property;
public:
  MagneticFieldHelper( const Vector3& B_ ) : B( B_ ) {};
  void HookL1( Mols* mols );
  void HookL0( SingleMolEntity* mol );
};


class MagneticFieldPlugin : public ProcessPlugin
{
  SimpleMD* system;
  Vector3 B; //in Tesla
public:
  explicit MagneticFieldPlugin( Vector3 B_ ): B( B_ ) {};
  void initialize( SimpleMD* sys ){ system = sys; };
  void ForceHook();
  ~MagneticFieldPlugin()
  {
    cerr << "~Plugin::MagneticFieldPlugin" << endl;
  }
};


#endif
