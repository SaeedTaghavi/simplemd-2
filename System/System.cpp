#include "System/System.hpp"
#include "System/SimpleMD.hpp"
#include "debug.hpp"
#include "FileIO.hpp"
#include "Const.hpp"
#include <sstream>
#include "Plugin/PairProcessPlugin.hpp"

//Write final data for continuation
void 
SimpleMD::write(){
  WriteDTPS( dt, *unit, *fout );
  *fout << "@TABS" << endl 
     << absoluteTime / ( pico * unit->sec ) << endl;
  rc.WriteRCOA( *unit, *fout );
  Write( moldict, *unit, *fout );
  molecules->Write( *unit, dt, *fout );
}



//Write snapshot (under construction)
void 
SimpleMD::SnapShot( ostream& out, int fmt ){
  out << "@TABS" << endl
         << absoluteTime / ( pico * unit->sec ) << endl;
  if ( fmt == TRAJ_COORDONLY ){
    molecules->SnapShot( *unit, dt, out );
  }
  else if ( fmt == TRAJ_FULL ){
    molecules->Write( *unit, dt, out );
  }
}



//Calculate pressure ( in internal pressure unit )
double
SimpleMD::Pressure( double idealpv )
{
  /* See Ueda P.162 (9.19) */
  double vrsum = pv.Tr();
  return ( idealpv + vrsum/3.0 ) / molecules->Volume();
}



//Calculate temperature ( in internal energy unit )
double 
SimpleMD::Temperature( double ek )
{
  // kT = 2 Ek / DoF
  // ek excludes fixed molecules
  return ek*2.0 / dof;
}



// Calculate force
void
SimpleMD::Force()
{
  molecules->Preforce();
  molecules->Force( rc, pv );
  molecules->Postforce();
  //2010-2-9 added for external field plugins
  int nplugin = plugins.size();
  for( int i=0; i<nplugin; i++ ){
    plugins[i]->ForceHook();
  }
}



//One-line log for monitor.
string 
SimpleMD::log()
{
  ostringstream logging;
  const Unit& u = *unit;
  double ek = GetEk(); //excludes fixed molecules
  Vector3 ekt;
  GetEkt( ekt );       //excludes fixed molecules
  double ektrans = ekt.x + ekt.y + ekt.z;
  double idealpv = 2 * ektrans / 3;
  logging.precision(17);
  logging << step
     << " " << absoluteTime / ( pico * u.sec )
     << " " << Temperature( ek ) / u.K()
     << " " << Pressure( idealpv ) / u.atm()
     << " " << molecules->Volume() / pow( Angstro * u.m, 3 )
     << " " << pv.GetEp() / ( kilo * u.J / mol )
     << " " << ek / ( kilo * u.J / mol )
     << " " << ( pv.GetEp()+ek ) / ( kilo * u.J / mol );
  return logging.str();
}    


  
//Simplest integrator
double 
SimpleMD::Symplectic2nd()
{
  molecules->ProgressMomentum( dt / 2 );
  molecules->ProgressPosition( dt );
  Force();
  molecules->ProgressMomentum( dt / 2 );
    
  return dt;
}



double 
SimpleMD::GetEk()
{
  return molecules->GetEk();
}



void 
SimpleMD::GetEkt( Vector3& ekt )
{
  molecules->GetEkt( ekt );
}




SimpleMD::SimpleMD(
                    double absoluteTime0,
                    double interval,
                    MolCollection* mol,
                    Unit* const unit0,
                    Truncation rc0,
                    map<string, FlexibleHandle>& dict,
                    const ProcessPluginArray& pl,
		    ostream* fout0,
		    ostream* ferr0
  ) : fout( fout0 ), ferr( ferr0 ), plugins( pl ), unit( unit0 ), dt( interval ), step( 0 ), absoluteTime( absoluteTime0 ),
      molecules(mol),rc(rc0),moldict(dict),input(stdin),dof(molecules->GetDoF())
{
  fout->precision(17);
  ferr->precision(17);
  //*fout << absoluteTime << "END\n";
}



SimpleMD::~SimpleMD()
{
  mesg("~SimpleMD\n");
  delete molecules;
}
   
  

void
SimpleMD::initialize()
{
  int nplugin = plugins.size();
  //*fout << nplugin << " nplugin\n";
  for( int i=0; i<nplugin; i++ ){
    plugins[i]->initialize( this );
  }
}



double
SimpleMD::OneStep()
{
  return Symplectic2nd();
}



bool
SimpleMD::running()
{
  double dt = OneStep();
  absoluteTime += dt;
  //cerr << dt << "SimpleMD::running()" << endl;

  bool cont = true;
  int nplugin = plugins.size();
  for( int i=0; i<nplugin; i++ ){
    if ( ! plugins[i]->running() ){
      cont = false;
    }
  }
  step ++;
  return cont;
}
  
  
  
void 
SimpleMD::terminate()
{
  write();
  int nplugin = plugins.size();
  for( int i=0; i<nplugin; i++ ){
    plugins[i]->terminate();
  }
}



double 
Simple4th::Symplectic4th()
{
  const double cr2 = cbrt(2);
  const double v0 = 1.0 / ( 2 - cr2 );
  const double v1 = - cr2 / ( 2 - cr2 );
  const double v2 = v0 / 2;
  const double v3 = ( 1 - cr2 ) * v2;
  //
  //New Ueda 
  //
  molecules->ProgressMomentum( dt * v2 );
  molecules->ProgressPosition( dt * v0 );
  Force();
  molecules->ProgressMomentum( dt * v3 );
  molecules->ProgressPosition( dt * v1 );
  Force();
  molecules->ProgressMomentum( dt * v3 );
  molecules->ProgressPosition( dt * v0 );
  Force();
  molecules->ProgressMomentum( dt * v2 );
    
  //ek = molecules->GetEk();
  return dt;
}



Simple4th::Simple4th(
                      double absoluteTime0,
                      double interval,
                      MolCollection* mol,
                      Unit* const unit0,
                      Truncation rc0,
                      map<string, FlexibleHandle>& dict,
                      const ProcessPluginArray& pl,
			   ostream* fout0,
			   ostream* ferr0

  ) : SimpleMD(
                absoluteTime0,
                interval,
                mol,
                unit0,
                rc0,
                dict,
                pl,
		fout0,
		ferr0
		)
{
}



double 
Simple4th::OneStep()
{
  return Symplectic4th();
}






