#ifndef GRIDCELL_HPP
#define GRIDCELL_HPP
#include "Cell/Cell.hpp"

/*
 *GridCellは、等間隔に分割されたCell。座標から直ちに所属する小セルを特
 *定できる。どの小セルが隣接しているか、そこにどの分子が属しているかが
 *すぐわかるので、分子間相互作用をN^2通り調べなくてもよくなる。大きな
 *システムでは速度が著しく向上する。
 */

/*GridCell consists of multiple cells*/
class GridCell : public Cell
{
private:
  Box* boundary;
  int     nx, ny, nz, nc;
  vector<RelocatableCellHandle> cells;
  //relative cell offsets for homomolecules
  int    (*delta0)[13];
  int ndelta0;
  Vector3 offset0[13];
  //relative cell offsets for heteromolecules
  //it is not currently used
  //int    (*delta1)[26];
  //Vector3 offset1[26];
  //boxes are controlled in another class. They are smart links.
  //BoxHandle boundary;
  //BoxHandle smallbox;
public:
  explicit GridCell();
  explicit GridCell( const Box& box, int ncompo, int _nx, int _ny, int _nz );
  virtual ~GridCell();
  int GetNumCompo() const;
  void Set( const int compo, const MolsHandle& mols );
  //void ScaleVelocity( const Vector3& r );
  //void ScalePosition( const Vector3& r );
  void Preforce();
  void Postforce();
  int Force( const Combination& combi, const Truncation& rc, PotVir &pv );
  double GetEk() const;
  void GetEkt( Vector3& ekt ) const;
  int GetSize( const int compo ) const;
  double GetMass() const;
  void ProgressPosition( double );
  void ProgressMomentum( double );
  void Write( const Unit& unit, ostream& to );
  MolsHandle PullAll( const int compo );
  MolsHandle CopyAll( const int compo );
  void TotalMomentum( Vector3& momentum ) const;
  //void AddVelocity( const Vector3& velocity );
  double* unserialize( double* const p ){ assert(0); return p; }
  double* serialize( double* p ) const;
  double* serializeforce( double* xi ) const { assert(0); return xi; }
  int  qdof() const { assert(0); return 0; }
  void Execute( CollectorPlugin& plugin );
  int  PairProcess( PairProcessPlugin& p, const Combination& combi, const Truncation& rc ) const;
  void Rescale( const Vector3& r );
private:
  //Conversion from global coordinate to cell coordinate
  int whichcell( const Vector3& coord );
  void relocate();
  void SetBox( const Box& box );
  virtual double Volume(){ return boundary->Volume(); }
  void ScaleRelativeCellVectors( const Vector3& r );
};



#endif
