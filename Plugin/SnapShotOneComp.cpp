#include "Plugin/CollectorPlugin.hpp"
#include "Plugin/SnapShotOneComp.hpp"
#include "System/SimpleMD.hpp"
#include "Const.hpp"

class SnapShotOneCompHelper : public CollectorPlugin{
  string id08;
  SimpleMD& system;
  int format;
public:
  SnapShotOneCompHelper( SimpleMD& system_, string id08_, int fmt ) : id08( id08_ ), system(system_), format(fmt) {};
  void HookL1( Mols* mols );
  void HookL0( SingleMolEntity* mol );
};



void
SnapShotOneCompHelper::HookL1( Mols* mols )
{
  MolPropertyHandle property = mols->GetProperty();
  if ( property->id08 == id08 ){
    ostream* ferr = system.GetFerr();
    //*ferr << "COMPO:" << id08 << endl;
    if ( format == TRAJ_COORDONLY ){
      mols->SnapShot( system.GetUnit(), *ferr );
    }
    else if ( format == TRAJ_FULL ){
      mols->Write( system.GetUnit(), *ferr );
    }
  }
  //それ以外の場合は下位に降りない。
}


void
SnapShotOneCompHelper::HookL0( SingleMolEntity* mol )
{
}



SnapShotOneComp::SnapShotOneComp( ostream* out, ostream* err, const string& id08_, int fmt)  : fout( out ), ferr( err ), id08(id08_), system(0), format(fmt)
{
}


SnapShotOneComp::SnapShotOneComp( const string& id08_, int fmt ) : fout( 0 ), ferr( 0 ), id08(id08_), format(fmt)
{
}





void
SnapShotOneComp::initialize( SimpleMD* sys )
{
  system = sys;
  if ( fout == 0 ) fout = system->GetFout();
  if ( ferr == 0 ) ferr = system->GetFerr();
}



bool
SnapShotOneComp::running()
{
  //ここで、CollectorPluginを初期化し、使いすてる。
  SnapShotOneCompHelper helper( *system, id08, format );
  helper.HookL3( system->GetMolCollection() );
  return true;
}
