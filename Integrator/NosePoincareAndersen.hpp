#ifndef NOSEPOINCAREANDERSEN_HPP
#define NOSEPOINCAREANDERSEN_HPP
#include <iostream>
#include "Unit.hpp"

using namespace std;

/*
 *Symplectic Integratorの補助関数。
 *
 *System.hppに入れてしまってもよさそうだが、分離しておいたほうが取り扱いに便利。
 *
 *
 */


//とりあえず継承しない。壊してしまいそうなので。

class NosePoincareAndersen
{
private:
  double s;
  double ps;
  double pv;
  double Qs;
  double Qv;
  double kT;
  double H0;
  double Pext;
  bool   disableTempControl;
public:
  NosePoincareAndersen( double s_ini, double ps_ini, double Qs_ini ) : s(s_ini), ps(ps_ini), Qs(Qs_ini), disableTempControl( false )
  {
  }
  double GetkT() const
  {
    return kT;
  }
  double GetPext() const
  {
    return Pext;
  }
  double GetH0() const
  {
    return H0;
  }
  void SetkT( double kT_ini )
  {
    kT = kT_ini;
  }
  void SetH0( double H0_ini )
  {
    H0 = H0_ini;
  }
  void SetPext( double Pext_ini )
  {
    Pext = Pext_ini;
  }
  void SetPs( double x )
  {
    ps = x;
  }
  void SetS( double x )
  {
    s = x;
  }
  void SetPv( double x )
  {
    pv = x;
  }
  double S() const {
    return s;
  }
  void ProgressPv( double tau, double pressure );
  void ProgressPs( double tau, double ekss, double U, int dof, double volume );
  double ProgressSV( double tau );
  double Eps( int dof );
  double Eks();
  double Epv( double volume );
  double Ekv();
  double Hzero(){ return H0; }
  void ReadNOPA( const Unit& unit, FILE* file );
  void WriteNOPA( const Unit& unit, ostream& to );
  void ReadSNPH( const Unit& unit, FILE* file );
  void WriteSNPH( const Unit& unit, ostream& to );
  void WriteNPAZ( const Unit& unit, ostream& to );
  void Write( const Unit& unit, ostream& to );
  double RecalcH0( double ek, double ep, int dof, double volume );
};


#endif
