#ifndef LISTVECTOR_H
#define LISTVECTOR_H
#include "Vector3.hpp"

//WASHED.

using namespace std;

class Truncation;
class Mols;
/*
 *�Хꥹ�Ȥ��ꤹ�륯�饹���Хꥹ�������ϡ�Vector���Ǥʤ���Ф���ۤ�
 *��ǽŪ�ʥ��åȤ��ʤ��������󲽤�˸���ˤʤ��ǽ���⤢�롣¾����ˡ��
 *�Ȥ���褦�˥����ǥ��󥰤���Τ�˾�ޤ�����
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
 * �⤷Bookkeeping�򤹤�Τʤ顢TruncPair�Τ⤦�����ʰפʥС��������롣
 * ����ꥹ�Ȥ��顢��Υ���ͤ�Ϳ�����̤Υꥹ�Ȥ���������褦�ʥ᥽�åɤ�������롣
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
