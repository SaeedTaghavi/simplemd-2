#include "Mols/Mols.hpp"
#include "Plugin/PairProcessPlugin.hpp"



void
Mols::Concat( const MolsHandle& src )
{
  int nmol = src->Size();
  for( int i=0; i < nmol; i++ ){
    Push( src->Peek( i ) );
  }
}



int
pairprocess0(
  PairProcessPlugin& p,
  const Mols& m1,
  const Mols& m2,
  const Intersite& im,
  const TruncPair& truncpair
  )
{
  p.Monitor2( m1, m2, im, truncpair );
  return 1;
}



int
PairProcess(
  PairProcessPlugin& p,
  const Mols& m1,
  const Mols& m2,
  const Intersite& im,
  const TruncPair& truncpair
  )
{
  int truncpairsize = truncpair.size();
  for(int kk=0; kk<truncpairsize; kk++){
    const ListItem& li = truncpair.GetListItems()[kk];
    int i = li.pair1;
    int j = li.pair2;
    const SingleMolEntity& mol1 = m1.Molecule(i);
    const SingleMolEntity& mol2 = m2.Molecule(j);
    p.Monitor( mol1, mol2, im, li );
  }
  return 1;
}




int pairprocess( PairProcessPlugin& p, const Mols& m1, const Intersite& im, const Truncation& rc )
{
  TruncPair truncpair;
  if ( truncpair.list_simple( m1, rc ) ){
    return pairprocess0( p, m1, m1, im, truncpair );
  }
  return 0;
}



int pairprocess_simple( PairProcessPlugin& p, const Mols& m1, const Mols& m2, const Intersite& im, const Truncation& rc )
{
  TruncPair truncpair;
  if ( truncpair.list_simple( m1, m2, rc ) ){
    return pairprocess0( p, m1, m2, im, truncpair );
  }
  return 0;
}



int pairprocess_offsetpbc( PairProcessPlugin& p, const Mols& m1, const Mols& m2, const Intersite& im, const Vector3& offset, const Box& box, const Truncation& rc )
{
  TruncPair truncpair;
  if ( truncpair.list_offsetpbc( m1, m2, offset, box, rc ) ){
    return pairprocess0( p, m1, m2, im, truncpair );
  }
  return 0;
}



int pairprocess_pbc( PairProcessPlugin& p, const Mols& m1, const Mols& m2, const Intersite& im, const Box& box, const Truncation& rc )
{
  TruncPair truncpair;
  if ( truncpair.list_pbc( m1, m2, box, rc ) ){
    return pairprocess0( p, m1, m2, im, truncpair );
  }
  return 0;
}



int pairprocess_pbc( PairProcessPlugin& p, const Mols& m1, const Intersite& im, const Box& box, const Truncation& rc )
{
  TruncPair truncpair;
  if ( truncpair.list_pbc( m1, box, rc ) ){
    return pairprocess0( p, m1, m1, im, truncpair );
  }
  return 0;
}



int pairprocess_offset( PairProcessPlugin& p, const Mols& m1, const Mols& m2, const Intersite& im, const Vector3& offset, const Truncation& rc )
{
  TruncPair truncpair;
  if ( truncpair.list_offset( m1, m2, offset, rc ) ){
    return pairprocess0( p, m1, m2, im, truncpair );
  }
  return 0;
}



