#include <cmath>
#include "Unit.hpp"



//For compatibility with nph2b.cc
double Unit::From_kJ_mol( double energy ) const
{
  return energy;
}



double Unit::To_kJ_mol( double energy ) const
{
  return energy;
}



double Unit::To_K( double energy ) const
{
  return energy;
}



//internal pressure unit is Pa
double Unit::From_atm( double pressure ) const
{
  return pressure;
}



double Unit::To_atm( double pressure ) const
{
  return pressure;
}



//internal time unit is pico seconds
double Unit::From_ps( double time ) const
{
  return time;
}



double Unit::To_ps( double time ) const
{
  return time;
}



//internal length unit is Angstrom
double Unit::From_Angstrom( double length ) const
{
  return length;
}



double Unit::To_Angstrom( double length ) const
{
  return length;
}



double Unit::From_g_mol( double mass ) const
{
  return mass;
}



double Unit::To_g_mol( double mass ) const
{
  return mass;
}



double Unit::From_Coulomb( double charge ) const
{
  return charge;
}



double Unit::To_Coulomb( double charge ) const
{
  return charge;
}




//internal energy unit is dJ/mol
double StandardUnit::From_kJ_mol( double energy ) const
{
  return energy * 100;
}



double StandardUnit::To_kJ_mol( double energy ) const
{
  return energy / 100.0;
}



double StandardUnit::To_K( double energy ) const
{
  return To_kJ_mol( energy ) * 1000 * Kelvin;
}



//internal pressure unit is dJ/mol / A^3
double StandardUnit::From_atm( double pressure ) const
{
  return pressure*PA * NB / (10 * 1e+30);
}



double StandardUnit::To_atm( double pressure ) const
{
  return pressure/(PA*NB) * ( 10 * 1e+30);
}



//internal time unit is pico seconds
double StandardUnit::From_ps( double time ) const
{
  return time;
}



double StandardUnit::To_ps( double time ) const
{
  return time;
}



//internal length unit is Angstrom
double StandardUnit::From_Angstrom( double length ) const
{
  return length;
}



double StandardUnit::To_Angstrom( double length ) const
{
  return length;
}



double StandardUnit::From_g_mol( double mass ) const
{
  return mass;
}



double StandardUnit::To_g_mol( double mass ) const
{
  return mass;
}



/*
 * Energy = Charge^2 / length
 */
double StandardUnit::From_Coulomb( double charge ) const
{
  return charge * sqrt( From_kJ_mol( COEFF ) );
}



double StandardUnit::To_Coulomb( double charge ) const
{
  return charge / sqrt( From_kJ_mol( COEFF ) );
}



