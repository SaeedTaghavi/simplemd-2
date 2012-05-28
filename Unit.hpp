#ifndef UNIT_HPP
#define UNIT_HPP
#include <cmath>
#include <iostream>

namespace SI {
  //aux unit
  static const double kilo  = 1e+3;
  static const double milli = 1e-3;
  static const double micro = 1e-6;
  static const double nano  = 1e-9;
  static const double pico  = 1e-12;
  static const double femto = 1e-15;
  static const double mol   = 6.022045e+23;
  static const double Angstro = 1e-10; //as an aux unit
  static const double kB    = 1.380662e-23;
};


//Conversion between conventional and internal units 
static const double NA = SI::mol;
static const double PA=101326;
static const double PIright =  3.14159265358979324;
static const double PItanaka = 3.1415926535897932;
static const double PI = PItanaka;
static const double CAright = 4.186;
//American Calorie
static const double CAtanaka = 4.184;
static const double CA = CAtanaka;
static const double Kelvin = 1.0 / ( NA * SI::kB );
static const double EE = 1.6021892e-19;
static const double EPSright = 0.885418782e-11;
static const double EPStanaka = (EE*EE*NA*1e7/(4*PI*CA *332.17752));
static const double EPS = EPStanaka;
static const double ANG2M = 1e-8;
static const double COEFFright = (EE*EE*NA/(4*PI*EPS*ANG2M * 1000));
static const double COEFFtanaka = (CAtanaka*332.17752);
static const double COEFF = COEFFtanaka;

//For compatibility with nph2b.cc
/*
class Unit
{
public:
  virtual ~Unit(){};
  virtual double From_kJ_mol( double energy ) const;
  virtual double To_kJ_mol( double energy ) const;
  virtual double To_K( double energy ) const;
  virtual double From_atm( double pressure ) const;
  virtual double To_atm( double pressure ) const;
  virtual double From_ps( double time ) const;
  virtual double To_ps( double time ) const;
  virtual double From_Angstrom( double length ) const;
  virtual double To_Angstrom( double length ) const;
  virtual double From_g_mol( double mass ) const;
  virtual double To_g_mol( double mass ) const;
  virtual double From_Coulomb( double charge ) const;
  virtual double To_Coulomb( double charge ) const;
};
*/

using namespace std;
using namespace SI;

//For compatibility with Genesis6
class Unit
{
  /*
   *内部単位系は今のところ以下の通り
   *長さ: Angstrom
   *エネルギー：dJ/mol
   *質量：g/mol
   *時間：ps
   *電荷：単位電荷e

   内部数値を、出力する場合にはこんな風にする。
   printf("%f\n", energy / (kilo*J/mol) );

   読みこむ場合は逆にかけ算になる。
   */
public:
  double sec, m, J, g, radian, Pa, coulomb;
  Unit()
    {
      //conversion coeff
      sec   = 1e12;
      m     = 1e10;
      J     = 0.1 * NA;
      g     = NA;
      radian= 1;
      Pa    = J / (m*m*m);
      coulomb = sqrt( COEFF * (kilo*J/mol) );
    }
  //custom units
  double K() const
    {
      return J / mol / Kelvin;
    }
  double atm() const
    {
      return Pa * PA;
    }

  virtual ~Unit()
  {
    cerr << "~Unit" << endl;
  }
};

#endif
