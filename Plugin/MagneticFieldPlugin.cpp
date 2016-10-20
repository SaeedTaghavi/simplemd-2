#include "Plugin/CollectorPlugin.hpp"
#include "System/SimpleMD.hpp"
#include "Plugin/MagneticfieldPlugin.hpp"
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
 *Lorenzian force by the magnetic field. Requested by Prof. Nomura.
 */

void
MagneticFieldHelper::HookL1( Mols* mols )
{
  //the molecule does not know his own physical properties.
  //They are recorded in different data structure IntrParams.
  //So they must be made accessible here.

  //get the property and keep it in the private variable.
  property = mols->GetProperty();
  //go deeper
  mols->PluginHookL0( *this );
}


Vector3
cross_product( const Vector3& a, const Vector3& b )
{
  Vector3 prod;
  prod.x = a.y*b.z - a.z*b.y;
  prod.y = a.z*b.x - a.x*b.z;
  prod.z = a.x*b.y - a.y*b.x;
  return prod; //it will be copied.
}



void
MagneticFieldHelper::HookL0( SingleMolEntity* mol )
{
  //here calculate the Lorentz force
  //change the behaviour based on the molecular type
  MonatomicMol* mmh = dynamic_cast<MonatomicMol*> ( mol );
  if ( mmh ) {
    assert(0);
  }
  else{
    RigidBody2* mmh = dynamic_cast<RigidBody2*> ( mol );
    assert(mmh);
    Rigid* rigid = dynamic_cast<Rigid*> ( property.get() );
    assert(rigid);
    cerr << "Hey we are here!" << endl;
    cerr << "B " << B.print() << endl;
    //calculate the velocities of the atoms.
    //1. group velocity  A / ps
    //Q1 does it really contain the velocity at the moment?
    const Vector3& vmol = mmh->com.GetVelocity();
    cerr << "velocity " << vmol.print() << endl;
    //2. angular velocity  rad / ps
    Vector3& av = mmh->w;
    int nsite = mmh->atom.size();
    for( int site=0; site<nsite; site++ ){
      //3. intramolecular positions  A
      //We assume it is calculated in force() in advance properly.
      Vector3 intra( mmh->atom[site].center.coord.x,
		     mmh->atom[site].center.coord.y,
		     mmh->atom[site].center.coord.z );
      //4. atomic velocity A / ps
      Vector3 vatom = cross_product( av, intra );
      vatom.x += vmol.x;
      vatom.y += vmol.y;
      vatom.z += vmol.z;
      //5. Lorentz force to atoms
      //Unit of B is unknown.
      const IntrParamsArray& intr = rigid->GetIntrParams();
      double q = intr[site].charge;
      cerr << "vatom" << site << " " << vatom.print() << endl;
      Vector3 F = cross_product(vatom, B);
      cerr << "F" << site << " " << F.print() << endl;
      cerr << "q " << q << endl;
      mmh->atom[site].center.force.x += q * F.x;
      mmh->atom[site].center.force.y += q * F.y;
      mmh->atom[site].center.force.z += q * F.z;
    }
  }
}



void
MagneticFieldPlugin::ForceHook()
{
  //ここで、CollectorPluginを初期化し、使いすてる。
  MagneticFieldHelper magneticfieldhelper( B );
  magneticfieldhelper.HookL3( system->GetMolCollection() );
}



