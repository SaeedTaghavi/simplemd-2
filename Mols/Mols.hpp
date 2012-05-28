#ifndef MOLS_HPP
#define MOLS_HPP

//#include <iostream>
#include "MolProperty.hpp"
#include "Vector3.hpp"
#include "Interaction/Truncation.hpp"
#include "Interaction/PotVir.hpp"
#include "SingleMol/SingleMol.hpp"
#include "Interaction/Combination.hpp"
#include "Interaction/ListVector.hpp"

using namespace std;
using namespace boost;

/*
 *同種分子の集合体を規定するクラス。Layer 1

 *力計算の計算速度が問題になる場合は、より下位のクラスであるSingleMolに落とさず、このクラスで直接計算したほうがいい場合もある。
 */


class Mols;
class TruncPair;
class CollectorPlugin;
class PairProcessPlugin;

typedef shared_ptr<Mols> MolsHandle;
typedef vector<MolsHandle> MolsArray;

class Mols
{
public:
  virtual ~Mols(){}

  //"Get" functions
  virtual double GetEk() const = 0;

  virtual void BoxCoordinate( const Box& box ) = 0;
  virtual MolsHandle Emmigrate( const Box& box ) = 0;
  //force between the same molecular set in the single periodic simulation cell.
  virtual int Force_PBC( const Intersite& im,
                         const Box& box,
                         const Truncation& rc,
                         PotVir &pv ) = 0;
  //force between the same molecular set in the single non-periodic simulation cell.
  virtual int Force( const Intersite& im,
		      const Truncation& rc,
		      PotVir &pv ) = 0;
  //force between the same molecular set in a pair of simulation cell.
  virtual int Force_Simple( const MolsHandle& m,
			     const Intersite& im,
			     const Truncation& rc,
			     PotVir &pv
			     ) = 0;
  virtual int Force_Offset( const MolsHandle& m,
			     const Intersite& im,
			     const Vector3 &offset,
			     const Truncation& rc,
			     PotVir &pv
			     ) = 0;
  virtual int Force_OffsetPBC( const MolsHandle& m,
			       const Intersite& im,
			       const Vector3 &offset,
			       const Box& box,
			       const Truncation& rc,
			       PotVir &pv
			       ) = 0;
  virtual int Force_PBC( const MolsHandle& m,
			  const Intersite& im,
			  const Box& box,
			  const Truncation& rc,
			  PotVir &pv
			  ) = 0;
  virtual int Force_General( const MolsHandle& m, const Intersite& im, const TruncPair& truncpair, PotVir &pv ) = 0;
  virtual void Preforce() = 0;
  virtual int Size() const= 0;
  //will be obsolete
  virtual SingleMolHandle Peek( int i ) const = 0;
  virtual const SingleMolEntity& Molecule( int i ) const = 0;
  virtual MolPropertyHandle GetProperty() const = 0;
  virtual int  Push( SingleMolHandle mol ) = 0;
  virtual MolsHandle EmptyClone() const = 0;
  virtual void Postforce() = 0;
  virtual void SetProperty( const MolPropertyHandle& p ) = 0;
  //virtual void Report() const = 0;
  virtual void ProgressMomentum( double dt ) = 0;
  virtual const Vector3& Position( int i ) const = 0;
  virtual void Translate( const Vector3& offset ) = 0;
  virtual void PluginHookL0( CollectorPlugin& plugin ) = 0;

  //could be plugin
  virtual void TotalMomentum( Vector3& momentum ) const = 0;
  virtual void GetEkt( Vector3& ekt ) const = 0;

  //I/O
  virtual void Write( const Unit& unit, ostream& to ) = 0;
  virtual void SnapShot( const Unit& unit, ostream& to ) = 0;

  //features replaced by plugin
  //virtual void AddVelocity( const Vector3& ) = 0;
  //virtual void ScaleVelocity( const Vector3& r ) = 0;
  //virtual void ScalePosition( const Vector3& r ) = 0;

  //functions for quenching
  virtual double* unserialize( double* const p ) = 0;
  virtual double* serialize( double* const p ) const = 0;
  virtual double* serializeforce( double* const p ) const = 0;

  //implemented common functions
  virtual void Concat( const MolsHandle& src ) ;
};


/*
template<class Mols, class Agent, class Arg> class mols_execute{
public:
  mols_execute( const Mols& m, const Arg& arg )
  {
    int nmol = m.Size();
    for( int i=0; i < nmol; i++ ){
      Agent a( m.mols[i], arg );
    }
  }
};
*/
  

//public
int PairProcess( PairProcessPlugin& p, const Mols& m1, const Mols& m2, const Intersite& im, const TruncPair& truncpair );


int pairprocess( PairProcessPlugin& p, const Mols& m1, const Intersite& im, const Truncation& rc );

int pairprocess_simple( PairProcessPlugin& p, const Mols& m1, const Mols& m2, const Intersite& im, const Truncation& rc );

int pairprocess_offsetpbc( PairProcessPlugin& p, const Mols& m1, const Mols& m2, const Intersite& im, const Vector3& offset, const Box& box, const Truncation& rc );

int pairprocess_pbc( PairProcessPlugin& p, const Mols& m1, const Mols& m2, const Intersite& im, const Box& box, const Truncation& rc );

int pairprocess_pbc( PairProcessPlugin& p, const Mols& m1, const Intersite& im, const Box& box, const Truncation& rc );

int pairprocess_offset( PairProcessPlugin& p, const Mols& m1, const Mols& m2, const Intersite& im, const Vector3& offset, const Truncation& rc );


#endif
