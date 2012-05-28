#ifndef MOLCOLLECTION_HPP
#define MOLCOLLECTION_HPP
#include "Vector3.hpp"
#include "Cell/Cell.hpp"
#include "Interaction/Interaction.hpp"
#include "Interaction/Combination.hpp"

class PairProcessPlugin;

/*
 *MolCollectionは、セル(分子群と、分子を入れる容器をあわせたもの)と、
 *分子間相互作用を束ねたもの。一般的な相互作用パラメータを与えられると、
 *それをもとに相互作用行列(どの分子とどの分子がどんな相互作用をするか)
 *を自動的に構築する。
 *
 *このレイヤは冗長かもしれない。
 */
class MolCollection
{
protected:
  Cell*   cell;
private:
  //MolPropertyArray* prop;
   /*
   *Interaction Combination.
   */
  Combination combi;
public:
  explicit MolCollection( Cell* c );
  ~MolCollection();

  //"Get" functions
  double GetEk();
  void GetEkt( Vector3& ekt ) const;
  int GetDoF();
  double Volume(){ return cell->Volume(); }
  const Combination& GetCombination() const { return combi; }
  const Cell& GetCell() const { return *cell; }

  //"Set" functions
  void push_back( const FlexibleHandle& ph, const MolsHandle& mh );
  void Set( Cell* cell );

  void Preforce();
  int Force( const Truncation& rc, PotVir &pv );
  void Postforce();
  //void ScaleVelocity( const Vector3& r );
  void Expand( const Vector3& r );
  void ProgressPosition( double dt );
  void ProgressMomentum( double dt );
  void Write( const Unit& unit, double dt, ostream& to );
  void SnapShot( const Unit& unit, double dt, ostream& to );
  void StopDrift();
  void Execute( CollectorPlugin& plugin );
  int  PairProcess( PairProcessPlugin& p, const Truncation& rc ) const;
  //for test
  //void AddVelocity( const Vector3& v ){ cell->AddVelocity( v ); }

private:
};



class MolCollection2 : public MolCollection
{
public:
  explicit MolCollection2( Cell* c ) : MolCollection( c ){}
  void unserialize( double* const p );
  void serialize( double* p ) const;
  void serializeforce( double* xi ) const;
  int  qdof() const;
};



#endif
