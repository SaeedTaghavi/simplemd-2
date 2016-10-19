#include <cassert>
#include "Plugin/PairProcessPlugin.hpp"
#include "Interaction/ListVector.hpp"
#include "SingleMol/SingleMol.hpp"
#include "Mols/Mols.hpp"
#include "Mols/RigidBodies.hpp"
#include "MolCollection.hpp"
#include "Cell/Cell.hpp"
#include "SingleMol/RigidBody.hpp"



void PairProcessPlugin::Monitor0( const MolCollection& coll, const Truncation& rc )
{
  coll.PairProcess( *this, rc );
}



void PairProcessPlugin::Monitor1( const Cell& cell, const Combination& combi, const Truncation& rc )
{
  cell.PairProcess( *this, combi, rc );
}



void PairProcessPlugin::Monitor2( const Mols& m1, const Mols& m2, const Intersite& im, const TruncPair& truncpair )
{
  ::PairProcess( *this, m1, m2, im, truncpair );
}



void PairProcessTest::Monitor( const SingleMolEntity& mol1, const SingleMolEntity& mol2, const Intersite& im, const ListItem& listitem )
{
  int o1 = mol1.GetOrder();
  int o2 = mol2.GetOrder();
  //const Vector3& delta = listitem.d;
  cout << o1 << " " << o2 << " pair\n";
}



void WaterHB::Monitor2( const Mols& m1, const Mols& m2, const Intersite& im, const TruncPair& truncpair )
{
  //Get a set of molecules
  const RigidBodies& r1 = dynamic_cast<const RigidBodies&> ( m1 );
  const RigidBodies& r2 = dynamic_cast<const RigidBodies&> ( m2 );
  //end if they are not rigid molecules
  if ( &r1 == 0 || &r2 == 0 )
    return;
  //Get the property( atom names ) of the molecule set
  MolProperty* const ph1 = r1.GetProperty().get();
  MolProperty* const ph2 = r2.GetProperty().get();
  Rigid* const  p1 = dynamic_cast<Rigid*>( ph1 );
  Rigid* const  p2 = dynamic_cast<Rigid*>( ph2 );
  //Must be same type of water models.
  if ( ! ( p1->id08 == "TIP4P   " || 
           p1->id08 == "NVDE____" || 
           p1->id08 == "SPC_E___" ||
           p1->id08 == "TIP4P2K5" ) ) return;
  if ( p1 != p2 ) return;
  //Get the interaction list vector ( of a subcell )
  int truncpairsize = truncpair.size();
  for(int kk=0; kk<truncpairsize; kk++){
    const ListItem& li = truncpair.GetListItems()[kk];
    int i = li.pair1;
    int j = li.pair2;
    //Get the single molecule
    const RigidBody& mol1 = dynamic_cast<const RigidBody&>( m1.Molecule(i) );
    const RigidBody& mol2 = dynamic_cast<const RigidBody&>( m2.Molecule(j) );
    int order1 = mol1.GetOrder();
    int order2 = mol2.GetOrder();
    //結合するのは、相互作用する対とは限らない。
    //int npair = im.intr.size();
    
    //Find the nearest O-H pair
    double dmin = 2.5*2.5;
    int hmol = -1;
    int omol = -1;
    int nsite1 = p1->GetNumSite();
    int nsite2 = p2->GetNumSite();
    for( int s1=0; s1<nsite1; s1++ ){
      for( int s2=0; s2<nsite2; s2++ ){
        string atom1 = p1->GetAtomName( s1 );
        string atom2 = p2->GetAtomName( s2 );
        if ( ( atom1 == "H" && atom2 == "O" ) ||
             ( atom1 == "O" && atom2 == "H" ) ){
          // Get the position of the atom in the internal coordinate
          const SiteOfAction& soa1 = mol1.atom[s1].center;
          const SiteOfAction& soa2 = mol2.atom[s2].center;
          Vector3 delta;
          // intermolecular vector is prepared in li.
          delta.x  = li.d.x + soa1.coord.x - soa2.coord.x;
          delta.y  = li.d.y + soa1.coord.y - soa2.coord.y;
          delta.z  = li.d.z + soa1.coord.z - soa2.coord.z;
          double sq = delta.square();
          if ( sq < dmin ){
            dmin = sq;
            if ( atom1 == "H" ){
              hmol = order1;
              omol = order2;
            }
            else{
              hmol = order2;
              omol = order1;
            }
          }
        }
      }
    }
    if ( 0 <= hmol ){
      cout << hmol << " " << omol << endl;
    }
  }
}




void WaterHB::Monitor1( const Cell& cell, const Combination& combi, const Truncation& rc )
{
  cout << "@NGPH" << endl;
  //第0成分の分子数
  cout << cell.GetSize( 0 ) << endl;
  cell.PairProcess( *this, combi, rc );
  cout << "-1 -1" << endl;
}



