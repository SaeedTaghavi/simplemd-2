#ifndef POTVIR_HPP
#define POTVIR_HPP

//Washed.

#include "Vector3.hpp"
/*class for collecting potential energy and virial*/
class PotVir
{
  double ep;
  Matrix33 vr;
public:
  PotVir();
  PotVir( double ep0, const Matrix33& vr0 );
  void Clear();
  void Add( const PotVir& pv );
  void Set( double ep0, const Matrix33& vr0 );
  double Tr() const;
  double GetEp() const { return ep; }
};
#endif

