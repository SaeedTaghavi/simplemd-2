#ifndef PERIODICCELL_HPP
#define PERIODICCELL_HPP
#include "Cell/Cell.hpp"

/*
 *単純な周期境界セル。RelocatableCellに吸収されたので今はもう使わない。
 */


class PeriodicCell : public SimpleCell
{
protected:
  //the boundary box is shared, i.e. it expands externally.
  BoxHandle boundary;
public:
  PeriodicCell( const BoxHandle& box, int ncompo );
  ~PeriodicCell();
  void Preforce();
  int  Force( const Combination& combi, const Truncation& rc, PotVir &pv );
private:
  void relocate();
};



class BookPeriodicCell : public PeriodicCell
{
private:
  vector<SimplePair> simplepair;
public:
  BookPeriodicCell( const BoxHandle& box, int ncompo );
  int  Force( const Combination& combi, const Truncation& rc, PotVir &pv );
};

#endif
