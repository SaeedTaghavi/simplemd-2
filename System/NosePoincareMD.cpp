#include <sstream>
#include "System/NosePoincareMD.hpp"
#include "debug.hpp"
#include "Plugin/CollectorPlugin.hpp"
#include "Integrator/NosePoincare.hpp"


NosePoincareMD::NosePoincareMD(
                                double absoluteTime0,
                                double interval,
                                MolCollection* mol,
                                Unit* const unit0,
                                Truncation rc0,
                                map<string, FlexibleHandle>& dict,
                                const ProcessPluginArray& pl,
                                NosePoincare* nosepoincare0,
                                bool forceloaded,
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
  nosepoincare = nosepoincare0;
  assert( nosepoincare != 0 );

  //運動量をスケールしておく。
  double s = nosepoincare->S();
  ScaleVelocity2Plugin p( s );
  p.HookL3( molecules );
  //molecules->ScaleVelocity( ratio );
  
  //力、ポテンシャル、ヴィリアルを計算する。
  Force();
  //ただし、一部の成分が、座標だけの情報からスタートする場合は、
  //力は0にしておく。(あるいは力を半分にしてソフトスタートする?)
  if ( ! forceloaded ){
    //力は最初は0にしておく。
    molecules->Preforce();
  }
  double H0 = nosepoincare->Hzero();
  if ( H0 == 0 ){
    //double ekss = GetEk();
    //H0 =  ekss + pv.GetEp() + nosepoincare->Ek() + nosepoincare->Ep( dof );
    //nosepoincare->SetH0( H0 );
    nosepoincare->ResetH0( GetEk(), pv.GetEp(), dof );
  }
}




NosePoincareMD::~NosePoincareMD()
{
  mesg("~NosePoincareMD\n");
  delete nosepoincare;
}


  
void 
NosePoincareMD::write()
{
  SimpleMD::write();
  nosepoincare->WriteNOPO( *unit, *fout );
}



double 
NosePoincareMD::GetEk()
{
  double ek = molecules->GetEk();
  double s = nosepoincare->S();
  ek /= (s*s);
  return ek;
}



void 
NosePoincareMD::GetEkt( Vector3& ekt )
{
  molecules->GetEkt( ekt );
  double s = nosepoincare->S();
  double ssi = 1/(s*s);
  ekt.x *= ssi;
  ekt.y *= ssi;
  ekt.z *= ssi;
}



string
NosePoincareMD::log()
{
  ostringstream logging;
  const Unit& u = *unit;
  double ek = GetEk();
  double s = nosepoincare->S();
  logging.precision(17);
  logging << SimpleMD::log()
          << " " << s 
          << " " << ( nosepoincare->Ep( dof ) ) / ( kilo * u.J / mol )
          << " " << ( nosepoincare->Ek() ) / ( kilo * u.J / mol )
          << " " << s * ( pv.GetEp() + ek + nosepoincare->Ep( dof ) + nosepoincare->Ek() - nosepoincare->Hzero() ) / ( kilo * u.J / mol )
	  << " " << nosepoincare->GetkT() / u.K();
  return logging.str();
}    



double
NosePoincareMD::Symplectic2nd()
{
  //Nose-Poincare Ueda P.214
  nosepoincare->ProgressL3( dt / 2 );
  
  molecules->ProgressMomentum( nosepoincare->S() * dt / 2 );
  nosepoincare->ProgressL2( dt / 2, pv.GetEp() );
  
  double ek = molecules->GetEk();
  molecules->ProgressPosition( dt / nosepoincare->S() );
  //dof includes the molecular degrees of freedom only.
  nosepoincare->ProgressL1( dt, ek, dof );
  
  Force();
  molecules->ProgressMomentum( nosepoincare->S() * dt / 2 );
  nosepoincare->ProgressL2( dt / 2, pv.GetEp() );
  
  nosepoincare->ProgressL3( dt / 2 );
  
  return dt;
}



void 
NosePoincareMD::terminate()
{
  //運動量をスケールしておく。
  double s = nosepoincare->S();
  ScaleVelocity2Plugin p( 1.0/s );
  p.HookL3( molecules );
  //molecules->ScaleVelocity( ratio );
  SimpleMD::terminate();
}
