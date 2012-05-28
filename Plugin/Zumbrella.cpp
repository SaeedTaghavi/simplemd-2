#include "Plugin/Zumbrella.hpp"
#include "System/SimpleMD.hpp"
#include "Cell/Cell.hpp"
#include "Mols/Mols.hpp"
#include "Mols/MonatomicMols.hpp"
#include "Mols/RigidBodies.hpp"
#include "MolProperty.hpp"
#include "MolCollection.hpp"
#include "Vector3.hpp"

/*
 *This plugin is called from Zumbrella.
 *����ɤ��ˤ��������򤢤�����Τǡ�����ɤΤ�������ˤ������礦����Ϥ������ޤϤ���ɤΤۤ����Ϥ����ʤ�ʤ���
 */

void
ZumbrellaHelper::Execute2( Mols* mols )
{
  MolPropertyHandle property = mols->GetProperty();
  if ( property->id08 == id08 ){
    MonatomicMols* mmh = dynamic_cast<MonatomicMols*> ( mols );
    assert( mmh );
    int nmol = mols->Size();
    assert( nmol == 1 );
    //ñ����ʬ�Ҥǡ�������id08�����פ�����˸¤ꡢ���̤�ƤӤ�����
    mols->Execute( *this );
  }
  //����ʳ��ξ��ϲ��̤˹ߤ�ʤ���
}



void
ZumbrellaHelper::Execute3( SingleMolEntity* mol )
{
  //Z�����ˤϳ��Ϥ�ä���XY�������Ϥ�ä���
  MonatomicMol* mmh = dynamic_cast<MonatomicMol*> ( mol );
  //This test is unnecessary.
  if ( mmh ){
    const Vector3& force=mmh->GetForce();
    //Add external force here ????
    const Vector3& coord=mmh->Position();
    //cerr << force.z << ":" << fconst << ":" << coord.z << ":" << distance << endl;
    mmh->center.force.Set(0,0,force.z - fconst * ( coord.z - distance ) );
  }
}
Zumbrella::Zumbrella( double distance_, double fconst_, string id08_ ) : distance( distance_ ), fconst( fconst_ ), id08( id08_ )
{
}


void
Zumbrella::forcepatch()
{
  //�����ǡ�CollectorPlugin�����������Ȥ����Ƥ롣
  ZumbrellaHelper zumb( distance, fconst, id08 );
  zumb.Execute0( system->GetMolCollection() );
}


void
Zumbrella::initialize( SimpleMD* sys )
{
  system = sys;
}

