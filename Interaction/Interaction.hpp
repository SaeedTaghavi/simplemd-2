#ifndef INTERACTION_HPP
#define INTERACTION_HPP
#include <vector>
//#include <boost/shared_ptr.hpp>
#include "Unit.hpp"
#include "Vector3.hpp"

using namespace std;
//using namespace boost;

/*
 *原子間2体相互作用を規定するクラス。まだもう少し種類を増やす必要があるだろう。
 */



//Virtual base class must have at least one virtual function!
//Otherwise dynamic_cast fails. ( on g++ )
class Interaction
{
public:
  virtual ~Interaction();
  virtual void Force( const Vector3& d, double& fo, double& potential ) const = 0;  
};



class LJ : public Interaction
{
public:
  double eps, sig;
  LJ();
  LJ( double e, double s );
  void read( FILE* const file, const Unit& unit );
  void set( double e, double s );
  void Force( const Vector3& d, double& fo, double& potential ) const ;  
};



class Coulomb : public Interaction
{
public:
  double sqrch;
  Coulomb();
  Coulomb( double qq );
  void read( FILE* const file, const Unit& unit );
  void set2( double sqrcharge );
  void Force( const Vector3& d, double& fo, double& potential ) const ;  
};



class IntrParams
{
public:
  double eps,sig,charge;
  double epssqrt, sighalf;
  IntrParams()
  {
    set(0,0,0);
  }
  IntrParams( FILE* const file, const Unit& unit );
  IntrParams( double e, double s, double c );
  void read( FILE* const file, const Unit& unit );
  void set( double e, double s, double charge );
  void Write( const Unit& unit, ostream& to ) const;
};



class LJC : public Interaction
{
public:
  LJ lj;
  Coulomb coulomb;
  LJC();
  LJC( const IntrParams&, const IntrParams& );
  LJC( FILE* const file, const Unit& unit );
  LJC( double e, double s, double c );
  void set2( double e, double s, double sqrcharge );
  void Force( const Vector3& d, double& fo, double& potential ) const;  
};



typedef std::shared_ptr<Interaction> InteractionHandle;
typedef vector<InteractionHandle> InteractionArray;
typedef vector<IntrParams> IntrParamsArray;





typedef std::shared_ptr<LJ> LJHandle;
typedef std::shared_ptr<LJC> LJCHandle;
typedef std::shared_ptr<Coulomb> CoulombHandle;

#endif
