#ifndef PAIRPROCESSPLUGIN_GUARD
#define PAIRPROCESSPLUGIN_GUARD


#include <cassert>
#include "MolProperty.hpp"

/*
 * 分子対ごとの処理を行う汎用プラグイン。
 */

struct ListItem;
class SingleMolEntity;
class Mols;
class Intersite;
class TruncPair;
class MolCollection;
class Cell;
class Combination;
class Truncation;

class PairProcessPlugin{
  int combi;
public:
  //"Get" function
  virtual int GetCombi(){ return combi; } // returns current combi

  //"Set" function
  virtual void SetCombi( int c ){ combi = c; }

  virtual void Monitor0( const MolCollection& coll, const Truncation& rc );
  virtual void Monitor1( const Cell& cell, const Combination& combi, const Truncation& rc );
  virtual void Monitor2( const Mols& m1, const Mols& m2, const Intersite& im, const TruncPair& truncpair );
  virtual void Monitor( const SingleMolEntity& mol1, const SingleMolEntity& mol2, const Intersite& im, const ListItem& listitem ) = 0;
};



class PairProcessTest: public PairProcessPlugin{
public:
  void Monitor( const SingleMolEntity& mol1, const SingleMolEntity& mol2, const Intersite& im, const ListItem& listitem );
};



class WaterHB: public PairProcessPlugin{
public:
  //1成分決め打ち。第0成分の分子数をノード数として出力する。
  //内部で、セル間相互作用をどんな順序で計算するかは既定されていないので、その場で出力するのはまずい。出力データは、一旦全部vectorに蓄積すべき。!!!
  void Monitor1( const Cell& cell, const Combination& combi, const Truncation& rc );
  void Monitor2( const Mols& m1, const Mols& m2, const Intersite& im, const TruncPair& truncpair );
  void Monitor( const SingleMolEntity& mol1, const SingleMolEntity& mol2, const Intersite& im, const ListItem& listitem ){ assert(0); }
};




#endif
