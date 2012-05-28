#ifndef RIGIDBODIES_HPP
#define RIGIDBODIES_HPP

#include "MolProperty.hpp"
#include "Vector3.hpp"
#include "Interaction/Truncation.hpp"
#include "Interaction/PotVir.hpp"
#include "SingleMol/RigidBody.hpp"
#include "Mols/Mols.hpp"

using namespace std;
using namespace boost;

/*
 *RigidBodies is a set of polyatomic molecules.
 *NOTE: the 0th site of a molecule must be the center of mass position.
 *In some molecule type, you may have to add extra site for com 
 * even if com does not interact.
 */
class RigidBodies : public Mols
{
protected:
private:
  MolPropertyHandle prop;
public:
  vector<RigidBody> mols;
  //additional properties
  //end
  RigidBodies( const MolPropertyHandle& );
  ~RigidBodies();
  void BoxCoordinate( const Box& box );
  MolsHandle Emmigrate( const Box& box );
  //force between the same molecular set in the single periodic simulation cell.
  int Force_PBC( const Intersite& im, 
                 const Box& box,
                 const Truncation& rc,
                 PotVir &pv );
  //force between the same molecular set in the single non-periodic simulation cell.
  int Force( const Intersite& im,
	      const Truncation& rc,
	      PotVir &pv );
  //force between the same molecular set in a pair of simulation cell.
  int Force_Simple( const MolsHandle& m,
		     const Intersite& im,
		     const Truncation& rc,
		     PotVir &pv
		     );
  int Force_Offset( const MolsHandle& m,
		     const Intersite& im,
		     const Vector3 &offset,
		     const Truncation& rc,
		     PotVir &pv
		     );
  int Force_OffsetPBC( const MolsHandle& m,
		       const Intersite& im,
		       const Vector3 &offset,
		       const Box& box,
		       const Truncation& rc,
		       PotVir &pv );
  int Force_PBC( const MolsHandle& m,
		  const Intersite& im,
		  const Box& box,
		  const Truncation& rc,
		  PotVir &pv
		  );
  int Size() const;
  MolPropertyHandle GetProperty() const;
  SingleMolHandle Peek( int i ) const;
  const SingleMolEntity& Molecule( int i ) const;
  int  Push( SingleMolHandle mol );
  void Postforce();
  void SetProperty( const MolPropertyHandle& p );
  void Preforce();
  double GetEk() const;
  void GetEkt( Vector3& ekt ) const;
  void ProgressPosition( double dt );
  void ProgressMomentum( double dt );
  void Translate( const Vector3& offset );
  MolsHandle EmptyClone() const;
  //void Report() const;
  //void ScaleVelocity( const Vector3& r );
  //void ScalePosition( const Vector3& r );
  const Vector3& Position( int i ) const;
  int Force_General( const MolsHandle& m, const Intersite& im, const TruncPair& truncpair, PotVir &pv );
  void TotalMomentum( Vector3& momentum ) const;
  //void AddVelocity( const Vector3& velocity );
  //void Concat( const MolsHandle& src );
  void PluginHookL0( CollectorPlugin& plugin );

  //Serialize/Unserialize configurations of all components.
  virtual double* unserialize( double* const p );
  virtual double* serialize( double* const p ) const;
  virtual double* serializeforce( double* const xi ) const;

  //I/O
  void Write( const Unit& unit, ostream& to );
  void SnapShot( const Unit& unit, ostream& to );
  void WriteNX4A( const Unit& unit, ostream& to );
  void ReadWTG5( int n, const Unit& unit, FILE* file );
  void ReadNX4A( int n, const Box& box, const Unit& unit, FILE* file );

private:
  void WriteWTG5( const Unit& unit, ostream& to );
  RigidBody* peek1( int i ) const;
  //molecule-type-specific force subroutines
  int force( 
	       RigidBodies& m2,
	       const Vector3 &offset, 
	       const Truncation& rc,
	       PotVir &pv
	       );
  RigidBody*   pull( int i );
  int  push1( const RigidBody& mol );
  void size( int n );
  //void force2accel( double dt, const Vessel& vessel );
  void collectforce();
};


#include "Mols/MonatomicMols.hpp"
#include "Interaction/ListVector.hpp"


int
force_simple(
             RigidBodies& m1,
             RigidBodies& m2,
	     const Intersite& im,
             const Truncation& rc,
             PotVir &pv
             );



int
force_offset(
             RigidBodies& m1,
             RigidBodies& m2,
	     const Intersite& im,
             const Vector3 &offset, 
             const Truncation& rc,
             PotVir &pv
             );





int
force_pbc(
	  RigidBodies& m1,
	  RigidBodies& m2,
	  const Intersite& im,
	  const Box& box, 
	  const Truncation& rc,
	  PotVir &pv );





/*heteromolecular interactions*/

int
force_common_2nd(
  RigidBodies& m1,
  MonatomicMols& m2,
  const Intersite& im,
  const TruncPair& lv,
  PotVir &pv
  );




int
force_common_3rd(
  RigidBodies& m1,
  MonatomicMols& m2,
  const Intersite& im,
  const TruncPair& lv,
  PotVir &pv
  );





int
force_offset(
             RigidBodies& m1,
             MonatomicMols& m2,
	     const Intersite& im,
             const Vector3 &offset, 
             const Truncation& rc,
             PotVir &pv
             );


int
force_offsetpbc(
		RigidBodies& m1,
		MonatomicMols& m2,
		const Intersite& im,
		const Vector3 &offset, 
		const Box& box,
		const Truncation& rc,
		PotVir &pv
		);


int
force_simple(
             RigidBodies& m1,
             MonatomicMols& m2,
	     const Intersite& im,
             const Truncation& rc,
             PotVir &pv
             );



int
force_pbc(
             RigidBodies& m1,
             MonatomicMols& m2,
	     const Intersite& im,
             const Box& box,
             const Truncation& rc,
             PotVir &pv
             );



#endif
