#ifndef NOSEPOINCARE_HPP
#define NOSEPOINCARE_HPP
#include <iostream>
#include "Unit.hpp"

using namespace std;

/*
 *Symplectic Integratorの補助関数。
 *
 *System.hppに入れてしまってもよさそうだが、分離しておいたほうが取り扱いに便利。
 */



class NosePoincare
{
private:
  double s;
  double ps;
  double Q;
  double kT;
  double H0;
public:
  NosePoincare( double s_ini, double ps_ini, double Q_ini ) : s(s_ini), ps(ps_ini), Q(Q_ini)
  {
    H0 = 0;
  }
  void SetH0( double H0_ini )
  {
    H0 = H0_ini;
  }
  void SetS( double x ){ s = x; }
  void SetPs( double x ){ ps = x; }
  void SetkT( double x ){ kT = x; }
  double S() const { return s; }
  void ProgressL3( double tau ); 
  void ProgressL2( double tau, double U ); 
  void ProgressL1( double tau, double Ek, int dof ); 
  double Ep( int dof );
  double Ek();
  double Hzero(){ return H0; }
  void ReadNOPO( const Unit& unit, FILE* file );
  void WriteNOPO( const Unit& unit, ostream& to );
  double GetkT() const
  {
    return kT;
  }
  void ResetH0( double ek, double ep, int dof );
};



#endif
