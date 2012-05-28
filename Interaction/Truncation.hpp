#ifndef TRUNCATION__HPP
#define TRUNCATION__HPP
#include "Unit.hpp"
#include <iostream>

//WASHED.


using namespace std;

/*
 *いわゆるカットオフを規定するクラス。現在(2006-6-3)のコードでは、カッ
 *トオフ関数は1つしか設定できないが、本来は成分の組合せ
 *(Combination.hppで規定)ごとに、異なるCutoff長さが適用されてよいはず。

 *なお、SimpleMDでは、カットオフ関数は原子対の間に作用するのではなく、
 *分子対(通常は重心間)に作用する。
 */ 

/* Class of cutoff function */
class Truncation
{
  double outer, inner;
  double rq;
public:
  Truncation();
  Truncation( const Unit& unit, FILE* file );
  ~Truncation();
  void Smooth( double r, double& sfe, double& sff ) const;
  void WriteRCOA( const Unit& unit, ostream& to );
  void ReadRCOA( const Unit& unit, FILE* file );
  void Set( double in, double out );
  double GetOuter() const { return outer; }
  double GetInner() const { return inner; }
  double GetRQ() const    { return rq; }
};
#endif
