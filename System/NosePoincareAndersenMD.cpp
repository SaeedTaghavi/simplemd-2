#include <sstream>
#include "debug.hpp"
#include "System/NosePoincareAndersenMD.hpp"
#include "Plugin/CollectorPlugin.hpp"
#include "Integrator/NosePoincareAndersen.hpp"

NosePoincareAndersenMD::NosePoincareAndersenMD(
                                                double absoluteTime0,
                                                double interval,
                                                MolCollection* mol,
                                                Unit* const unit0,
                                                Truncation rc0,
                                                map<string, FlexibleHandle>& dict,
                                                const ProcessPluginArray& pl,
                                                NosePoincareAndersen* nosepoincareandersen0,
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
  nosepoincareandersen = nosepoincareandersen0;
  assert( nosepoincareandersen != 0 );

  NosePoincareAndersen& npa = *nosepoincareandersen;
  //運動量をスケールしておく。
  double ek = GetEk();
  double s = npa.S();
  ScaleVelocity2Plugin p( s );
  p.HookL3( molecules );
  
  //力、ポテンシャル、ヴィリアルを計算する。
  //ポテンシャルはsに影響されていないことを確認。2006-5-22
  Force();
  //ただし、一部の成分が、座標だけの情報からスタートする場合は、
  //力は0にしておく。(あるいは力を半分にしてソフトスタートする?)
  if ( ! forceloaded ){
    //力は最初は0にしておく。
    molecules->Preforce();
  }
  double H0 = npa.Hzero();//ファイルから読みこんだ値。
  ek = GetEk();
  double instantH0 = npa.RecalcH0( ek, pv.GetEp(), dof, molecules->Volume() );
  if ( H0 == 0 ){
    npa.SetH0( instantH0 );
  }
  else{
    //指定されたH0が、初期配置から計算されるH0に一致しない場合は、分子の速度をスケール
    //する
    double deltaH0 = H0 - instantH0;
    double ss = sqrt( ( deltaH0 + ek ) / ek );
    ScaleVelocity2Plugin p( ss );
    p.HookL3( molecules );
    npa.SetH0( H0 );
  }    
}



 
NosePoincareAndersenMD::~NosePoincareAndersenMD()
{
  mesg("~NosePoincareAndersenMD\n");
  delete nosepoincareandersen;
}


  
void 
NosePoincareAndersenMD::write()
{
  SimpleMD::write();
  nosepoincareandersen->Write( *unit, *fout );
}



double 
NosePoincareAndersenMD::GetEk()
{
  double ek = molecules->GetEk();
  double s = nosepoincareandersen->S();
  ek /= (s*s);
  return ek;
}



void 
NosePoincareAndersenMD::GetEkt( Vector3& ekt )
{
  molecules->GetEkt( ekt );
  double s = nosepoincareandersen->S();
  ekt.x /= (s*s);
  ekt.y /= (s*s);
  ekt.z /= (s*s);
}



string
NosePoincareAndersenMD::log()
{
  ostringstream logging;
  const Unit& u = *unit;
  NosePoincareAndersen& npa = *nosepoincareandersen;
  double s = nosepoincareandersen->S();
  double ek = GetEk();
  logging.precision(17);
  logging << SimpleMD::log()
          << " " << s 
          << " " << ( npa.Eps( dof ) ) / ( kilo * u.J / mol )
          << " " << ( npa.Eks() ) / ( kilo * u.J / mol )
          << " " << ( npa.Epv( molecules->Volume() ) ) / ( kilo * u.J / mol )
          << " " << ( npa.Ekv() ) / ( kilo * u.J / mol )
          << " " << s * ( pv.GetEp() + ek + npa.Eps( dof ) + npa.Eks() + npa.Epv( molecules->Volume() ) + npa.Ekv() - npa.Hzero() ) / ( kilo * u.J / mol )
	  << " " << npa.GetkT()  / u.K()
	  << " " << npa.GetPext() / u.atm();
  return logging.str();
}    



double
NosePoincareAndersenMD::Symplectic2nd()
{
  NosePoincareAndersen& npa = *nosepoincareandersen;
  //cerr << pv.ep << " first ep\n";
  //Nose-Poincare-Andersen by Sturgeon
  double s= npa.S();
  molecules->ProgressMomentum( s * dt / 2 );
  
  double ekss = GetEk();
  Vector3 ekt;
  GetEkt( ekt );
  double ektrans = ekt.x + ekt.y + ekt.z;
  double idealpv = 2 * ektrans / 3;
  double pressure = Pressure( idealpv );
  npa.ProgressPv( dt / 2, pressure );
  npa.ProgressPs( dt / 2, ekss, pv.GetEp(), dof, molecules->Volume() );
  //9f-1
  molecules->ProgressPosition( dt / ( 2*s ) );
  
  double oldv = molecules->Volume();
  double newv = oldv + npa.ProgressSV( dt / 2 );
  //box->Set( newv );Expandが体積も変更する。
  //cout << box->x << " new box\n";
  double r = pow( newv / oldv, 1.0/3.0 );
  Vector3 ratio( r,r,r );
  molecules->Expand( ratio );
  
  //9f-2
  s = npa.S();
  molecules->ProgressPosition( dt / ( 2*s ) );
  
  Force();
  ekss = GetEk();
  npa.ProgressPs( dt / 2, ekss, pv.GetEp(), dof, molecules->Volume() );
  
  GetEkt( ekt );
  ektrans = ekt.x + ekt.y + ekt.z;
  idealpv = 2 * ektrans / 3;
  pressure = Pressure( idealpv );
  npa.ProgressPv( dt / 2, pressure );
  
  //v13 = pow( box->Volume(), 1.0/3.0 );
  molecules->ProgressMomentum( s * dt / 2 );
  //cerr << pv.GetEp() << " last ep\n";
  return dt;
}



NosePoincareAndersenZMD::NosePoincareAndersenZMD(
                                                double absoluteTime0,
                                                double interval,
                                                MolCollection* mol,
                                                Unit* const unit0,
                                                Truncation rc0,
                                                map<string, FlexibleHandle>& dict,
                                                  const ProcessPluginArray& pl,
                                                NosePoincareAndersen* nosepoincareandersen0,
						  bool forceloaded,
			   ostream* fout0,
			   ostream* ferr0

  ) : NosePoincareAndersenMD(
                              absoluteTime0,
                              interval,
                              mol,
                              unit0,
                              rc0,
                              dict,
                              pl,
                              nosepoincareandersen0,
                              forceloaded,
			      fout0,
			      ferr0
 )
{
}



double
NosePoincareAndersenZMD::Symplectic2nd()
{
  NosePoincareAndersen& npa = *nosepoincareandersen;
  //cerr << pv.GetEp() << " first ep\n";
  //Nose-Poincare-Andersen by Sturgeon
  double s= npa.S();
  molecules->ProgressMomentum( s * dt / 2 );
  
  double ekss = GetEk();
  Vector3 ekt;
  GetEkt( ekt );
  //double ektrans = ekt.x + ekt.y + ekt.z;
  //double idealpv = 2 * ektrans / 3;
  double ektrans = ekt.z;
  double idealpv = 2 * ektrans;
  double pressure = Pressure( idealpv );
  npa.ProgressPv( dt / 2, pressure );
  npa.ProgressPs( dt / 2, ekss, pv.GetEp(), dof, molecules->Volume() );
  //9f-1
  molecules->ProgressPosition( dt / ( 2*s ) );
  
  double oldv = molecules->Volume();
  double newv = oldv + npa.ProgressSV( dt / 2 );
  //cout << box->x << " new box\n";
  Vector3 ratio( 1.0, 1.0, newv / oldv );
  //box->Scale( ratio );boxの大きさはExpandで変更されるはず。
  molecules->Expand( ratio );
  
  //9f-2
  s = npa.S();
  molecules->ProgressPosition( dt / ( 2*s ) );
  
  Force();
  ekss = GetEk();
  npa.ProgressPs( dt / 2, ekss, pv.GetEp(), dof, molecules->Volume() );
  
  GetEkt( ekt );
  //ektrans = ekt.x + ekt.y + ekt.z;
  //idealpv = 2 * ektrans / 3;
  ektrans = ekt.z;
  idealpv = 2 * ektrans;
  pressure = Pressure( idealpv );
  npa.ProgressPv( dt / 2, pressure );
  
  //v13 = pow( box->Volume(), 1.0/3.0 );
  molecules->ProgressMomentum( s * dt / 2 );
  //cerr << pv.GetEp() << " last ep\n";
  return dt;
}
 


void 
NosePoincareAndersenZMD::write()
{
  SimpleMD::write();
  nosepoincareandersen->WriteNPAZ( *unit, *fout );
}



void 
NosePoincareAndersenMD::terminate()
{
  //運動量をスケールしておく。
  //この方法だと、ランの途中で速度をwriteできない。writeそのものに、sでスケールする処理を組込まなければいけない。さらに、現在のコードでは集約とソートが必要になるので非常に中途出力処理が重い。order付きの座標出力を作るべき。
  double s = nosepoincareandersen->S();
  ScaleVelocity2Plugin p( 1.0/s );
  p.HookL3( molecules );
  //molecules->ScaleVelocity( ratio );
  SimpleMD::terminate();
}
