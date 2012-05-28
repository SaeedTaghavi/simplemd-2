#ifndef CELL_HPP
#define CELL_HPP
#include "Vector3.hpp"
#include "Mols/Mols.hpp"
#include "Interaction/Combination.hpp"
#include "Interaction/PotVir.hpp"
#include "Interaction/Truncation.hpp"

class CollectorPlugin;
class PairProcessPlugin;
/*
 *Cellは、多成分の分子群を入れる容器。それは具体的な形を持つ場合もある
 *し、周期境界な容器の場合もあるし、オープンである場合もある。容器を分
 *割する場合はこのオブジェクトを拡張する。(GridCell.hpp参照)
 */

//cell is a box containing molecules
//Virtual base class
class Cell
{
public:
  virtual ~Cell();

  //"Get" functions
  virtual double GetEk() const = 0;
  virtual void GetEkt( Vector3& ekt ) const = 0;
  virtual int GetSize( const int compo ) const = 0;
  virtual double GetMass() const = 0;
  virtual int GetNumCompo() const = 0;
  virtual void TotalMomentum( Vector3& momentum ) const = 0;
  virtual double Volume() = 0;

  //"Set" functions
  virtual void Set( const int compo, const MolsHandle &h ) = 0;

  //virtual void ScaleVelocity( const Vector3& r ) = 0;
  //pre-process of force calculation (initialization)
  virtual void Preforce() = 0;
  //force calculation
  virtual int  Force( const Combination& combi, const Truncation& rc, PotVir &pv ) = 0;
  //post-process of force calculation (collect statistics)
  virtual void Postforce() = 0;
  virtual MolsHandle PullAll( const int compo ) = 0;
  virtual MolsHandle CopyAll( const int compo ) = 0;
  virtual void ProgressPosition( double ) = 0;
  virtual void ProgressMomentum( double ) = 0;
  virtual void Write( const Unit& unit, ostream& to ) = 0;
  //virtual void AddVelocity( const Vector3& velocity ) = 0;
  //functions for quenching
  virtual double* unserialize( double* const p ) = 0;
  virtual double* serialize( double* p ) const = 0;
  virtual double* serializeforce( double* xi ) const = 0;
  virtual int  qdof() const = 0;
  virtual void PluginHookL1( CollectorPlugin& plugin ) = 0;
  virtual int  PairProcess( PairProcessPlugin& p, const Combination& combi, const Truncation& rc ) const = 0;
  virtual void Rescale( const Vector3& r ) = 0;
};



//simple cell is an open box containing molecules
class SimpleCell : public Cell
{
private:
protected:
public:
  //groups of the molecules
  MolsArray mols;
  explicit SimpleCell();
  explicit SimpleCell( int ncompo );
  virtual ~SimpleCell();
  void SetNumCompo( int compo );
  int GetNumCompo() const;
  void Preforce();
  void Set( const int compo, const MolsHandle &h );
  void Push( int compo, const SingleMolHandle &m );
  void Postforce();
  int  Force( const Combination& combi, const Truncation& rc, PotVir &pv );
  //void ScaleVelocity( const Vector3& r );
  double GetEk() const;
  void GetEkt( Vector3& ekt ) const;
  int GetSize( const int compo ) const;
  double GetMass() const;
  void ProgressPosition( double );
  void ProgressMomentum( double );
  virtual void Write( const Unit& unit, ostream& to ){}
  virtual void SnapShot( const Unit& unit, const int compo, ostream& to );
  MolsHandle PullAll( const int compo );
  MolsHandle CopyAll( const int compo );
  void TotalMomentum( Vector3& momentum ) const;
  //void AddVelocity( const Vector3& velocity );
  //functions for quenching
  //Serialize/Unserialize configurations of all components.
  double* unserialize( double* const p );
  double* serialize( double* const p ) const;
  double* serializeforce( double* const xi ) const;
  //Serialize/Unserialize configurations of single component.
  double* unserialize( int compo, double* const p );
  double* serialize( int compo, double* const p ) const;
  double* serializeforce( int compo, double* const xi ) const;
  int  qdof() const;

  void PluginHookL1( CollectorPlugin& plugin );
  int  PairProcess( PairProcessPlugin& p, const Combination& combi, const Truncation& rc ) const;
  virtual double Volume(){ return 1; }
  void Rescale( const Vector3& r );
private:
  //void reportsize();
};






  


// relocatablecell is a closed box with an origin containing molecules
class RelocatableCell : public SimpleCell
{
protected:
  Vector3 origin;

  //Handleで渡すということは、知らないうちにどこかで書換えられる可能性がある。バグの温床になる。
  //BoxHandle boundary;
  //Boxは多形を使いたい。しかし、その操作は内部だけで行いたい。(外部でいつのまにか書換えられることは避けたい。)
  //とりあえず、ポインタとして実装し、その確保と開放はすべてクラス内で行うことにする。
  Box* boundary;
  //total number of molecules in the cell
  //updated at migration.
  int     residents;
public:
  RelocatableCell();
  RelocatableCell( const Box& b, int ncompo, const Vector3 &o );
  ~RelocatableCell();
  int IsEmpty() const {return residents == 0; }
  //void SetBox( const Box& b, const Vector3 &o );
  //extract the molecules outside the boundary.
  //emmigrated molecules are removed from the list.
  const Vector3& GetOrigin() const;
  MolsArray*  Emmigrate();
  void PushAbs( int compo, const SingleMolHandle &m );
  int Force_Intercell( const Combination& combi, const Vector3 &o, const RelocatableCell& c, const Truncation& rc, PotVir &pv );
  int  Force( const Combination& combi, const Truncation& rc, PotVir &pv );
  MolsHandle PullAll( const int compo );
  void Preforce();
  void Rescale( const Vector3& r );
  virtual double Volume(){ return boundary->Volume(); }
  void Write( const Unit& unit, ostream& to ){ boundary->WriteBOX3( unit, to ); }
  int PairProcess( PairProcessPlugin& p, const Combination& combi, const Truncation& rc ) const;
  int PairProcess_Intercell( PairProcessPlugin& p, const Combination& combi, const Vector3 &o, const RelocatableCell& c, const Truncation& rc ) const;

};


typedef shared_ptr<RelocatableCell> RelocatableCellHandle;

/*
template<class Cell, class Agent, class Arg> class cell_execute
{
 public:
  cell_execute( Cell cell, Arg arg )
    {
      int ncompo = cell.mols.size();
      for( int i=0;i<ncompo;i++)
	Agent x( cell.mols[i], arg );
    }
};
*/

#endif
