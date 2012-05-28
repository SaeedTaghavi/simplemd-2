#ifndef RIGIDBODY_HPP
#define RIGIDBODY_HPP
#include "SingleMol/PolyatomicMol.hpp"
#include "MolProperty.hpp"

class Rigid;
class Truncation;
class Intersite;
class ListItem;


/*
 *剛体分子をQuaternionで表現するクラス。こちらはすでに役目を終え、使わない予定だが、RigidBody2が継承しているので消せない。
 */


/*Position and force of a polyatomic molecule*/
class RigidBody : public PolyatomicMol
{
protected:
    //inertia, copied from prop in the constructor
    Vector3 in;
public:
  //Center-of-Mass
  MonatomicMol com;
  //position (in addition to Polyatomic sites)
  //SiteOfAction center;
  Quaternion  q;
  Quaternion dq;
  Vector3    w;
  Matrix33  rot;
  //vector<Vector3> intra; //intramolecular coordinates relative to com.
  //derivatives
  Vector3   torque;
  MolPropertyHandle prop;
  //end
  RigidBody( int nsite, const MolPropertyHandle& p );
  virtual ~RigidBody();

  //"Get" function
  double GetEk() const;

  const Vector3& Position() const;
  void Preforce();
  void Postforce( double mass );

  void CollectForce();
  void BoxCoordinate( const Box& box );
  void ProgressPosition( double dt );

  //I/O
  void WriteWTG5( const Unit& unit, ostream& to );
  void WriteNX4A( const Unit& u, ostream& to );
  void ReadWTG5( int o, const Unit& u, FILE* file );
  void ReadNX4A( int o, const Unit& u, FILE* file );

  //for sorting
  int operator<( const RigidBody& that ) const { return (order < that.order);}

  //local extension
  double GetSquareVelocity( double mass ) const;
  const Vector3& GetMomentOfInertia() const { return in; }
  void ProgressMomentum( double dt );
  void prepareintra();

  //features for plugin
  void ScaleVelocity( const Vector3& r );
  void ScaleVelocity( double r );
  void ScalePosition( const Vector3& r );
  void Translate( const Vector3& );
  void AddVelocity( const Vector3& velocity );

  //Serialize/Unserialize configuration.
  double* unserialize( double* const p );
  double* serialize( double* const p ) const;
  double* serializeforce( double* const xi ) const;

protected:
  //local function
  void resetsiteforce();
  void preparerotationmatrix();

  //RigidBody( const RigidBody & );
  //RigidBody& operator=( const RigidBody & );
};



int
force_offset(
        RigidBody& r1, Rigid& p1,
        RigidBody& r2, Rigid& p2,
        const Vector3 &offset, 
        const Truncation& rc,
        PotVir &pv
        );
int
force_pbc(
        RigidBody& r1, Rigid& p1,
        RigidBody& r2, Rigid& p2,
        const Box& box, 
        const Truncation& rc,
        PotVir &pv
        );
int
force_simple(
        RigidBody& r1, Rigid& p1,
        RigidBody& r2, Rigid& p2,
        const Truncation& rc,
        PotVir &pv
        );


int
force_offset2(
             RigidBody& r1,
             RigidBody& r2,
             const Intersite& im,
              const Vector3 &offset, 
             const Truncation& rc,
             PotVir &pv
             );
int
force_pbc2(
          RigidBody& r1,
          RigidBody& r2,
          const Intersite& im,
           const Box& box, 
          const Truncation& rc,
          PotVir &pv
          );
int
force_simple2(
             RigidBody& r1,
             RigidBody& r2,
             const Intersite& im,
              const Truncation& rc,
             PotVir &pv
             );


/*
bool
force_core_core(
		RigidBody& r1,
		RigidBody& r2,
		const Intersite& im,
		const Vector3& d0,
		PotVir &pv,
		double sfe,
		double sff
		);
*/


double
Potential( const SingleMolEntity& r1, const SingleMolEntity& r2, const Intersite& im, const ListItem& li );


#endif
