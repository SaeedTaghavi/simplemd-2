#include "Plugin/Log.hpp"
#include "System/SimpleMD.hpp"

Log::Log( ostream* out, ostream* err ) : fout( out ), ferr( err ), tag( "@LOGE" )
{
}



Log::Log() : fout( 0 ), ferr( 0 ), tag( "@LOGE" )
{
}



void
Log::initialize( SimpleMD* sys )
{
  system = sys;
  if ( fout == 0 ) fout = system->GetFout();
  if ( ferr == 0 ) ferr = system->GetFerr();
}



bool
Log::running()
{
  *ferr << tag << endl;
  *ferr << system->log() << endl;
  return true;
}



void
Log::terminate()
{
  //*fout << "@NLOG"  << endl
  //	<< interval << endl;
}
