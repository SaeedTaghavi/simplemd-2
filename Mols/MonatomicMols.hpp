#ifndef MONATOMICMOLS2_H
#define MONATOMICMOLS2_H

#include <iostream>
#include "MolProperty.hpp"
#include "Vector3.hpp"
#include "SingleMol/MonatomicMol.hpp"
#include "Mols/Mols.hpp"
#include "Interaction/Combination.hpp"
#include "Interaction/ListVector.hpp"
#include "Interaction/Truncation.hpp"
#include "Interaction/PotVir.hpp"


class CollectorPlugin;


using namespace std;
//using namespace boost;

class MonatomicMols : public Mols
{
private:
  MolPropertyHandle prop;
  int isFixed;
public:
  vector<MonatomicMol> mols;
  //vector<Vector3> f;
  MonatomicMols( int isfixed );
  MonatomicMols( MolPropertyHandle, int isfixed );
  MonatomicMols( MolPropertyHandle, int nmol, const Unit& unit, FILE* input, int isfixed );
  virtual ~MonatomicMols();
  void BoxCoordinate( const Box& box );
  MolsHandle Emmigrate( const Box& box );
  int Force_PBC( const Intersite& im,
	      const Box& box, const Truncation& rc, PotVir &pv );
  int Force( const Intersite& im, const Truncation& rc, PotVir &pv );
  int Force_Simple( const MolsHandle& m,
		     const Intersite& im,
		     const Truncation& rc,
		     PotVir &pv );
  int Force_Offset( const MolsHandle& m,
		     const Intersite& im,
		     const Vector3 &offset,
		     const Truncation& rc,
		     PotVir &pv );
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
		 PotVir &pv );
  int Force_General( const MolsHandle& m, const Intersite& im, const TruncPair& truncpair, PotVir &pv );
  virtual double GetEk() const;
  virtual int IsFixed() const;
  void GetEkt( Vector3& ekt ) const;
  void Preforce();
  SingleMolHandle Peek( int i ) const;
  const SingleMolEntity& Molecule( int i ) const;
  int  Push( SingleMolHandle mol );
  //void ScaleVelocity( const Vector3& r );
  //void ScalePosition( const Vector3& r );
  void Translate( const Vector3& offset );
  MolsHandle EmptyClone() const;
  void SetProperty( const MolPropertyHandle& p );
  MolPropertyHandle GetProperty() const;
  int Size() const;
  void Postforce();
  //void Report() const {}
  void ProgressPosition( double dt );
  void ProgressMomentum( double dt );
  void Write( const Unit& unit, ostream& to );
  void SnapShot( const Unit& unit, ostream& to );
  const Vector3& Position( int i ) const;
  void ReadATG5( int nmol, const Unit& unit, FILE* input );
  void ReadAR3A( int n, const Box& box, const Unit& unit, FILE* file );
  void TotalMomentum( Vector3& momentum ) const;
  //void AddVelocity( const Vector3& );
  //void Concat( const MolsHandle& src );
  void PluginHookL0( CollectorPlugin& plugin );
    
  //たぶんMols.cppに一括してしまえると思う。
  //Serialize/Unserialize configurations of all components.
  virtual double* unserialize( double* const p );
  virtual double* serialize( double* const p ) const;
  virtual double* serializeforce( double* const xi ) const;

private:
  void WriteATG5( const Unit& unit, ostream& to );
  void WriteAR3A( const Unit& unit, ostream& to );
  int force( 
	       MonatomicMols& m2,
	       const Vector3 &offset, 
	       const Truncation& rc,
	       PotVir &pv
	       );
  void init( int isfixed );
  MonatomicMol* pull( int i );
  void size( int n );
  int  push1( const MonatomicMol& mol );
  MonatomicMol* peek1( int i ) const;
};


int
force_offset(
      MonatomicMols& m1,
      MonatomicMols& m2,
      const Vector3 &offset, 
      const Truncation& rc,
      PotVir &pv
      );



int
force_pbc(
      MonatomicMols& m1,
      MonatomicMols& m2,
      const Box& box, 
      const Truncation& rc,
      PotVir &pv
      );



/*
template<class Mols, class Agent, class Arg> class mols_execute
{
public:
  mols_execute( const MolsHandle& m, Arg r )
  {
    int nmol = ((Mols*)m.get())->mols.size();
    for(int i=0; i<nmol; i++){
      Agent x( ((Mols*)m.get())->mols[i], r );
    }
  }
};
*/


#endif
