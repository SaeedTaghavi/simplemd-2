#include <fstream>
#include <string>
#include "FileIO.hpp"
#include "Cell/Cell.hpp"
#include "Cell/GridCell.hpp"
#include "SingleMol/SingleMol.hpp"
#include "SingleMol/MonatomicMol.hpp"
#include "SingleMol/PolyatomicMol.hpp"
#include "SingleMol/RigidBody.hpp"
#include "SingleMol/RigidBody2.hpp"
#include "Mols/MonatomicMols.hpp"
#include "Mols/PolyatomicMols.hpp"
#include "Mols/RigidBodies.hpp"
#include "debug.hpp"
#include "System/System.hpp"
#include "MolCollection.hpp"
#include "Plugin/PairProcessPlugin.hpp"
#include "PairTool.hpp"

using namespace std;
using namespace boost;
using namespace SI;


//private implimentation
class PairPotential: public PairProcessPlugin{
  int com;
  int nmol;
  const Unit& unit;
  //targetted combi
public:
  PairPotential( int c, int n, const Unit& u ) : com(c), nmol(n), unit(u) {}
  void Monitor0( const MolCollection& coll, const Truncation& rc );
  void Monitor2( const Mols& m1, const Mols& m2, const Intersite& im, const TruncPair& truncpair );
  void Monitor( const SingleMolEntity& mol1, const SingleMolEntity& mol2, const Intersite& im, const ListItem& listitem );
};



void PairPotential::Monitor2( const Mols& m1, const Mols& m2, const Intersite& im, const TruncPair& truncpair )
{
  //barrier
  if ( com == GetCombi() ){
    ::PairProcess( *this, m1, m2, im, truncpair );
  }
}



void PairPotential::Monitor( const SingleMolEntity& mol1, const SingleMolEntity& mol2, const Intersite& im, const ListItem& listitem )
{
  using namespace SI;
  //::potential resides in RigidBody.cpp
  double pot = ::Potential( mol1, mol2, im, listitem );
  int o1 = mol1.GetOrder();
  int o2 = mol2.GetOrder();
  pot /= ( unit.J / mol );
  double doo = sqrt( listitem.d.square() );
  cout << o1 << " " << o2 << " " << pot << " " << doo << endl;
}



void PairPotential::Monitor0( const MolCollection& coll, const Truncation& rc )
{
  cout << "@PAIR" << endl << nmol << endl;
  coll.PairProcess( *this, rc );
  cout << "-1 -1 0 0" << endl;
}












class PairPotentialTool : public PairTool
{

public:
  PairPotentialTool() : PairTool(){}
  virtual ~PairPotentialTool()
  {
    cerr << "~PairPotentialTool" << endl;
  }
  virtual void initialize(){};
  virtual void terminate(){};
  bool running()
  {
    MolCollection* molecules = Read( stdin );
    if ( 0 == molecules ){
      return false;
    }
    molecules->Preforce();
    PairPotential pairp( 0, molecules->GetCell().GetSize(0), unit );
    pairp.Monitor0( *molecules, rc );
    delete molecules;
    return true;
  }
};




int main(int argc,char *argv[])
{
  PairPotentialTool ppt;
  ppt.initialize();
  while( ppt.running() );
  ppt.terminate();
}



