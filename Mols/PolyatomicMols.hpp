#ifndef POLYATOMICMOLS_H
#define POLYATOMICMOLS_H

#include <iostream>
#include "SingleMol/PolyatomicMol.hpp"
#include "Mols/Mols.hpp"

using namespace std;
//using namespace boost;

/*
 *PolyatomicMols is a set of polyatomic molecules.
 *It is not necessary rigid, i.e. one can make spring-bead molecule 
 * with this class when spring interaction is given.
 *It is a natural extension of MonatomicMols3.
 */
class PolyatomicMols : public Mols
{
protected:
private:
  MolPropertyHandle prop;
  int isFixed;
public:
  vector<PolyatomicMol> mols;

  PolyatomicMols( MolPropertyHandle, int isfixed );
  ~PolyatomicMols();
  void BoxCoordinate( const Box& box );
  MolsHandle Emmigrate( const Box &box );
  //force between the same molecular set in the single periodic simulation cell.
  int Force_PBC( const Intersite& im,
                 const Box& box,
                 const Truncation& rc,
                 PotVir &pv );
  //force between the same molecular set in the single non-periodic simulation cell.
  int Force( const Intersite& im, const Truncation& rc, PotVir &pv );
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
		    PotVir &pv
		    );
  int Force_PBC( const MolsHandle& m,
                 const Intersite& im,
                 const Box& box,
                 const Truncation& rc,
                 PotVir &pv
    );
  int Force_General( const MolsHandle& m, const Intersite& im, const TruncPair& truncpair, PotVir &pv ){assert(0);return 0;}
  int Size() const;
  MolPropertyHandle GetProperty() const;
  SingleMolHandle Peek( int i ) const;
  const SingleMolEntity& Molecule( int i ) const;
  int  Push( SingleMolHandle mol );
  void Postforce( //const Vessel&
                  );
  void SetProperty( const MolPropertyHandle& p );
  void Preforce();
  //double GetEk() const;
  void GetEkt( Vector3& ekt ) const;
  void Translate( const Vector3& offset );
  MolsHandle EmptyClone() const;
  //void Report() const;
  void Init( int n, int nsite );
  void ProgressMomentum( double dt );
  void ProgressPosition( double dt );
  void Write( const Unit& unit, ostream& to ){}
  const Vector3& Position( int i ) const;
  void TotalMomentum( Vector3& momentum ) const;
  //void AddVelocity( const Vector3& );
  void PluginHookL0( CollectorPlugin& plugin );
  //void Concat( const MolsHandle& src );
private:
  PolyatomicMol* peek1( int i ) const;
  //molecule-type-specific force subroutines
  int force( 
	       PolyatomicMols& m2,
	       const Vector3 &offset, 
	       const Truncation& rc,
	       PotVir &pv
	       );
  PolyatomicMol*   pull( int i );
  int  push1( const PolyatomicMol& mol );
  //void ScaleVelocity( const Vector3& r );
  //void ScalePosition( const Vector3& r );
  void size( int n );
};


int
force_offset(
	PolyatomicMols& m1,
	PolyatomicMols& m2,
	const Vector3 &offset, 
	const Truncation& rc,
	PotVir &pv
	);




int
force_pbc(
             PolyatomicMols& m1,
             PolyatomicMols& m2,
             const Box& box, 
             const Truncation& rc,
             PotVir &pv
             );




#endif
