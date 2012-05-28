#ifndef LISTVECTOR_H
#define LISTVECTOR_H
#include "Vector3.hpp"

//WASHED.

using namespace std;

class Truncation;
class Mols;
/*
 *対リストを規定するクラス。対リスト方式は、Vector機でなければそれほど
 *性能的なメリットがないし、並列化の妨げになる可能性もある。他の方法も
 *使えるようにコーディングするのが望ましい。
 */


struct SimpleListItem {
  int pair1;
  int pair2;
};


struct ListItem {
  int pair1;
  int pair2;
  Vector3 d;
  double sfe, sff;
};




/*
 * もしBookkeepingをするのなら、TruncPairのもう少し簡易なバージョンを作る。
 * あるリストから、距離閾値を与えて別のリストを生成するようなメソッドを準備する。
 */


class TruncPair;
class Mols;

/*
 *Just record pairs (not record relative vectors)
 * for Bookkeeping method.
 */
class SimplePair{
  vector<SimpleListItem> listitems;
public:
  SimplePair(){}
  int size() const { return listitems.size(); }
  int list_offset( const Mols& m1, const Mols& m2, const Vector3& offset, double rc );
  int list_offsetpbc( const Mols& m1,
		      const Mols& m2,
		      const Vector3& offset,
		      const Box& box,
		      double rc );
  int list_pbc( const Mols& m1, const Mols& m2, const Box& box, double rc );
  int list_simple( const Mols& m1, const Mols& m2, double rc );
  int list_pbc( const Mols& m, const Box& box, double rc );
  int list_simple( const Mols& m, double rc );
  TruncPair* Compress( const Mols& m1, const Mols& m2, const Box& box, const Truncation& rc );
  const vector<SimpleListItem>& GetListItems() const { return listitems; }
};


class TruncPair{
  vector<ListItem> listitems;
public:
  TruncPair(){}
  int size() const { return listitems.size(); }
  int list_offset( const Mols& m1, const Mols& m2, const Vector3& offset, const Truncation& rc );
  int list_offsetpbc( const Mols& m1,
		      const Mols& m2,
		      const Vector3& offset,
		      const Box& box,
		      const Truncation& rc );
  int list_pbc( const Mols& m1, const Mols& m2, const Box& box, const Truncation& rc );
  int list_simple( const Mols& m1, const Mols& m2, const Truncation& rc );
  int list_pbc( const Mols& m, const Box& box, const Truncation& rc );
  int list_simple( const Mols& m, const Truncation& rc );
  void push_back( ListItem li ){ listitems.push_back( li ); }
  const vector<ListItem>& GetListItems() const { return listitems; }
};

#endif
