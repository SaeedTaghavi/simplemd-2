#ifndef MONATOMICMOL_HPP
#define MONATOMICMOL_HPP
#include "SingleMol/SingleMol.hpp"
//#include "MolProperty.hpp"
//#include "Interaction/Truncation.hpp"
//#include "Interaction/PotVir.hpp"
//#include "Interaction/Combination.hpp"

class CollectorPlugin;
struct SiteOfAction;
class Vector3;
class PotVir;
class Unit;

/*
 *一分子の特性を規定するクラス。単原子分子、多原子柔軟分子、多原子剛体分子を準備してある。
 *
 *直線型剛体分子も分離したほうがいいかも。
 *柔軟分子は未完成。まだ使えない。
 */



/*Position and force of a single molecule*/
class MonatomicMol : public SingleMolEntity
{
  //for symplectic integrator only
  Vector3 velocity;
  
  int       order;
protected:
    //mass, copied from prop in the constructor
  double mass, massi;

public:
  //Temporarily published
  SiteOfAction center;
  //MonatomicMol(){};
  MonatomicMol( const MonatomicMol & );
  MonatomicMol( double m );
  ~MonatomicMol();

  //"Get" functions
  const Vector3& Position() const;
  const Vector3& GetVelocity() const;
  int GetOrder() const {return order;}
  const Vector3& GetForce() const { return center.force; }
  double GetEk() const;

  //"Set" functions
  void SetPosition( const Vector3& v ) { center.coord = v; }
  void AddForce( const Vector3& f );
  void SetVelocity( const Vector3& v ) { velocity = v; }
  void SetMass( double m ) { mass = m; massi = 1.0 / m; }

  void Preforce();
  void Postforce( double mass );
  void BoxCoordinate( const Box& box );
  void ProgressPosition( double dt );
  void ProgressMomentum( double dt );

  //I/O (should be elsewhere)
  void WriteAR3A( const Unit& unit, ostream& to );
  void WriteATG5( const Unit& unit, ostream& to );
  void ReadATG5( int o, const Unit& u, FILE* file );
  void ReadAR3A( int o, const Unit& u, FILE* file );

  //for sorting
  int operator<( const MonatomicMol& that ) const { return (order < that.order);
  }

  //Serialize/Unserialize configuration.
  double* unserialize( double* const p );
  double* serialize( double* const p ) const;
  double* serializeforce( double* const xi ) const;

  //features for plugin
  void AddVelocity( const Vector3& velocity );
  void Translate( const Vector3& offset );
  void ScaleVelocity( const Vector3& r );
  void ScaleVelocity( double r );
  void ScalePosition( const Vector3& r );

private:
  //crude copy is allowed.
  //MonatomicMol( const MonatomicMol & );
  //MonatomicMol& operator=( const MonatomicMol & );
};


typedef std::shared_ptr<MonatomicMol>   MonatomicMolHandle;


/*
template<class Mol, class Agent, class Arg> class execute
{
public:
  execute( Mol m, Arg r )
  {
    Agent x( m, r );  //execute constructor
  }
};
*/

class SingleMolAgent;

#endif
