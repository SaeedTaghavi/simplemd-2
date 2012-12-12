#include "Plugin/CollectorPlugin.hpp"
#include "System/SimpleMD.hpp"
#include "Plugin/FixComponentPlugin.hpp"
#include "Cell/Cell.hpp"
#include "Cell/GridCell.hpp"
#include "Mols/Mols.hpp"
#include "Mols/MonatomicMols.hpp"
#include "Mols/RigidBodies.hpp"
#include "Mols/RigidBodies2.hpp"
#include "MolProperty.hpp"
#include "MolCollection.hpp"
#include "Vector3.hpp"

/*
 *ぶんしにくわわるちからを0にすることで、ぶんしのへいしんをとめるぷらぐいん。とうめんたんげんしぶんしのみたいおう。
 *固定した分子の間の、力の計算は行われるので、計算が速くはならない。
 */

void
FixComponentHelper::HookL2( Cell* cell )
{
	//bypass the dig-in mechanism
	SimpleCell* sch = dynamic_cast<SimpleCell*> ( cell );
	if ( sch ){
	  this->HookL1( sch->mols[compo].get() );
	}
	else{
		GridCell* gch = dynamic_cast<GridCell*> ( cell );
		assert (gch );
		int nc = gch->GetNumCells();
	  for( int c=0; c<nc; c++ ){
			const RelocatableCellHandle rc = gch->GetCell(c);
		//			/should be called recursively, but failed.
		  this->HookL1( rc->mols[compo].get() );
		}
  }
	//	cout << "FixComponentHelper::HookL2 - " << compo << endl;
}


void
FixComponentHelper::HookL0( SingleMolEntity* mol )
{
  MonatomicMol* mmh = dynamic_cast<MonatomicMol*> ( mol );
  if ( mmh ) {
    mmh->center.force.Set(0,0,0);
	}
	else{
	  RigidBody2* mmh = dynamic_cast<RigidBody2*> ( mol );
	  assert(mmh);
		mmh->com.center.force.Set(0,0,0);
		mmh->torque.Set(0,0,0);
	}
}



void
FixComponentPlugin::ForceHook()
{
  //ここで、CollectorPluginを初期化し、使いすてる。
  FixComponentHelper fixcomponenthelper( compo );
  //cout << ":" << id08 << ":" << endl;
  fixcomponenthelper.HookL3( system->GetMolCollection() );
}


void
FixComponentPlugin::initialize( SimpleMD* sys )
{
  system = sys;
}

