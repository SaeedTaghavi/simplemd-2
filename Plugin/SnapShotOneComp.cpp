#include "Plugin/CollectorPlugin.hpp"
#include "Plugin/SnapShotOneComp.hpp"
#include "System/SimpleMD.hpp"

class SnapShotOneCompHelper : public CollectorPlugin{
  string id08;
  SimpleMD& system;
public:
  SnapShotOneCompHelper( SimpleMD& system_, string id08_ ) : id08( id08_ ), system(system_) {};
  void Execute2( Mols* mols );
  void Execute3( SingleMolEntity* mol );
};



void
SnapShotOneCompHelper::Execute2( Mols* mols )
{
  MolPropertyHandle property = mols->GetProperty();
  if ( property->id08 == id08 ){
    ostream* ferr = system.GetFerr();
    *ferr << "COMPO:" << id08 << endl;
    mols->SnapShot( system.GetUnit(), *ferr );
  }
  //それ以外の場合は下位に降りない。
}


void
SnapShotOneCompHelper::Execute3( SingleMolEntity* mol )
{
}



SnapShotOneComp::SnapShotOneComp( ostream* out, ostream* err, const string& id08_)  : fout( out ), ferr( err ), id08(id08_), system(0)
{
}


SnapShotOneComp::SnapShotOneComp( const string& id08_ ) : fout( 0 ), ferr( 0 ), id08(id08_)
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
  SnapShotOneCompHelper helper( *system, id08 );
  helper.Execute0( system->GetMolCollection() );
  return true;
}
