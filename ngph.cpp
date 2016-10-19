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
//using namespace boost;
//using namespace SI;

class GraphTool : public PairTool
{

public:
  GraphTool() : PairTool(){}
  virtual ~GraphTool()
  {
    cerr << "~GraphTool" << endl;
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
    WaterHB waterhb;
    waterhb.Monitor0( *molecules, rc );
    delete molecules;
    return true;
  }
};



int main(int argc,char *argv[])
{
  GraphTool gt;
  gt.initialize();
  while( gt.running() );
  gt.terminate();
}
