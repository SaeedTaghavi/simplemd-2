#ifndef COMBINATION_HPP
#define COMBINATION_HPP
#include <vector>
#include <boost/shared_ptr.hpp>
#include "Unit.hpp"
#include "Vector3.hpp"
#include "MolProperty.hpp"
#include "Interaction/Interaction.hpp"


/*
 * Intersiteは、2つの多原子分子の間の原子間相互作用の組み合わ
 * せをまとめたもの。
 * 
 * N原子分子とM原子分子の間は最大NM種類の相互作用がありうるが、実際に
 * は、相互作用しない対もある。そこで、実際に相互作用しうる対のみを
 * InteractionArrayに格納しておく。
 */

class Intersite
{
public:
  InteractionArray intr;
  int nsite1;
  int nsite2;
  vector<int> site1;
  vector<int> site2;

  void Set( const FlexibleHandle& p1,
	    const FlexibleHandle& p2 );
  void Set( int s1, const IntrParams& i1,
	    int s2, const IntrParams& i2 );
};


typedef shared_ptr<Intersite> IntersiteHandle;


/*
 *溶液における、成分ごとの相互作用を規定するクラス。
 */

class Combination
{
  FlexibleArray prop;
public:
  vector<IntersiteHandle> ima;
  vector<int> compo1;
  vector<int> compo2;

  //"Get" functions
  const FlexibleArray& GetProperty() const { return prop; }

  //"Set" functions
  int push_back( const FlexibleHandle& ph );
};

#endif
