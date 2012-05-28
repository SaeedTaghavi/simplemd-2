#include "Cell.hpp"
#include "debug.hpp"
#include "Mols/MonatomicMols.hpp"
#include "Mols/RigidBodies.hpp"
#include "Mols/RigidBodies2.hpp"
#include "Plugin/CollectorPlugin.hpp"
#include "Plugin/PairProcessPlugin.hpp"

Cell::~Cell(){
  mesg("~Cell\n");
}




void SimpleCell::PluginHookL1( CollectorPlugin& plugin )
{
  int ncompo = mols.size();
  for( int i=0;i<ncompo;i++)
    plugin.HookL1( mols[i].get() );
}



SimpleCell::SimpleCell()
{
  SetNumCompo( 0 );
}



SimpleCell::SimpleCell( int ncompo )
{
  SetNumCompo( ncompo );
}



SimpleCell::~SimpleCell()
{
  mesg("~SimpleCell\n");
}



void SimpleCell::SetNumCompo( int compo )
{
  mols = MolsArray( compo );
}



int SimpleCell::GetNumCompo() const
{
  return mols.size();
}



/*
void SimpleCell::ScaleVelocity( const Vector3& r )
{
  int ncompo = mols.size();
  for( int i=0;i<ncompo;i++)
    mols[i]->ScaleVelocity( r );
}
*/








void SimpleCell::Set( const int compo, const MolsHandle &h )
{
  mols[compo] = h;
}



void SimpleCell::Push( int compo, const SingleMolHandle &m )
{
  //mesg("push6\n");
  mols[compo]->Push( m );
}






void SimpleCell::Preforce()
{
  int ncompo = mols.size();
  for( int i=0;i<ncompo;i++){
    mols[i]->Preforce();
  }
}



/*
void SimpleCell::reportsize()
{
  int ncompo = mols.size();
  for( int i=0;i<ncompo;i++){
    mesg2("size: %d\n", mols[i]->size() );
    mols[i]->Report();
  }
}
*/



void SimpleCell::Postforce()
{
  int ncompo = mols.size();
  for( int i=0;i<ncompo;i++)
    mols[i]->Postforce( //vessel
                        );
}



double SimpleCell::GetEk() const
{
  double ek = 0;
  int ncompo = mols.size();
  for( int i=0; i<ncompo; i++ ){
    ek += mols[i]->GetEk();
  }
  return ek;
}



void SimpleCell::GetEkt( Vector3& ekt ) const
{
  ekt.Set(0,0,0);
  int ncompo = mols.size();
  for( int i=0; i<ncompo; i++ ){
    Vector3 ek;
    mols[i]->GetEkt( ek );
    ekt.x += ek.x;
    ekt.y += ek.y;
    ekt.z += ek.z;
  }
}



int SimpleCell::Force( const Combination& combi, const Truncation& rc, PotVir &pv )
{
  //mesg2("force5 %d\n", ncompo);
  pv.Clear();
  int executed = 0;
  int ncombi = combi.ima.size();
  for( int co=0; co<ncombi; co++ ){
    int compo1 = combi.compo1[co];
    int compo2 = combi.compo2[co];
    PotVir p;
    Intersite* im = combi.ima[co].get();
    int result;
    if ( compo1 == compo2 ){
      //force between homomolecules
      result = mols[compo1]->Force( *im, rc, p );
    }
    else{
      //force between heteromolecules
      result = mols[compo1]->Force_Simple( mols[compo2], *im, rc, p );
    }
    if ( result ){
      executed++;
      pv.Add( p );
    }
  }
  return executed;
}



int SimpleCell::PairProcess( PairProcessPlugin& p, const Combination& combi, const Truncation& rc ) const
{
  int executed = 0;
  int ncombi = combi.ima.size();
  for( int co=0; co<ncombi; co++ ){
    int compo1 = combi.compo1[co];
    int compo2 = combi.compo2[co];
    Intersite* im = combi.ima[co].get();
    p.SetCombi( co );
    int result;
    if ( compo1 == compo2 ){
      //force between homomolecules
      result = pairprocess( p, *(mols[compo1]), *im, rc );
    }
    else{
      //force between heteromolecules
      result = pairprocess_simple( p, *(mols[compo1]), *(mols[compo2]), *im, rc );
    }
    if ( result ){
      executed++;
    }
  }
  return executed;
}



int 
SimpleCell::GetSize( const int compo ) const
{
  return mols[compo]->Size();
}
















/*
double*
unserialize0( MolsHandle& mh, double* const p )
{
  MonatomicMols* mmh = dynamic_cast<MonatomicMols*> ( mh.get() );
  if ( mmh ){
    int    nmol = mh->Size();
    for( int i=0; i < nmol; i++ ){
      MonatomicMol& m = mmh->mols[i];
      p = ::unserialize( m, p );
    }
    return p;
  }
  else{
    RigidBodies* mmh = dynamic_cast<RigidBodies*> ( mh.get() );
    assert( mmh );
    int    nmol = mh->Size();
    for( int i=0; i < nmol; i++ ){
      RigidBody& m = mmh->mols[i];
      p = ::unserialize( m, p );
    }
    return p;
  }
}



double*
serialize0( const MolsHandle& mh, double* p )
{
  MolPropertyHandle h = mh->GetProperty();
  MonatomicMols* mmh = dynamic_cast<MonatomicMols*> ( mh.get() );
  if ( mmh ){
    Flexible* e = dynamic_cast<Flexible*> ( h.get() );
    assert( e );
    int    nmol = mh->Size();
    for( int i=0; i < nmol; i++ ){
      MonatomicMol& m = mmh->mols[i];
      p = ::serialize( m, p );
    }
    return p;
  }
  else{
    //本当はRigidBodies2にしたい。
    RigidBodies* mmh = dynamic_cast<RigidBodies*> ( mh.get() );
    assert( mmh );
    Rigid* e = dynamic_cast<Rigid*> ( h.get() );
    assert( e );
    int    nmol = mh->Size();
    for( int i=0; i < nmol; i++ ){
      RigidBody& m = mmh->mols[i];
      p = ::serialize( m, p );
    }
    return p;
  }
}



double*
serializeforce0( const MolsHandle& mh, double* p )
{
  MolPropertyHandle h = mh->GetProperty();
  MonatomicMols* mmh = dynamic_cast<MonatomicMols*> ( mh.get() );
  if ( mmh ){
    Flexible* e = dynamic_cast<Flexible*> ( h.get() );
    assert ( e );
    vector<double> massi = e->GetMassI();

    int    nmol = mh->Size();
    for( int i=0; i < nmol; i++ ){
      MonatomicMol& m = mmh->mols[i];
      p = ::serializeforce( m, massi[0], p );
    }
    return p;
  }
  else{
    RigidBodies* mmh = dynamic_cast<RigidBodies*> ( mh.get() );
    assert( mmh );
    Rigid* e = dynamic_cast<Rigid*> ( h.get() );
    assert( e );
    double massi = 1.0 / e->GetMass();
    
    int    nmol = mh->Size();
    for( int i=0; i < nmol; i++ ){
      RigidBody& m = mmh->mols[i];
      p = ::serializeforce( m, massi, p );
    }
    return p;
  }
}
*/


int
qdof( const MolsHandle& mh )
{
  MonatomicMols* mmh = dynamic_cast<MonatomicMols*> ( mh.get() );
  if ( mmh ){
    int    nmol = mh->Size();
    return nmol * 3;
  }
  else{
    RigidBodies* mmh = dynamic_cast<RigidBodies*> ( mh.get() );
    assert( mmh );
    int    nmol = mh->Size();
    return nmol * 7;
  }
}






//double* unserialize0( MolsHandle& mh, double* p );
//double* serializeforce0( const MolsHandle& mh, double* p );



//set coordinates
int
SimpleCell::qdof() const
{
  int dof = 0;
  int ncompo = mols.size();
  for( int i=0; i<ncompo; i++ ){
    dof += ::qdof( mols[i] );
  }
  return dof;
}



//set coordinates
double*
SimpleCell::unserialize( double* const p )
{
  double* ptr = p;
  int ncompo = mols.size();
  for( int i=0; i<ncompo; i++ ){
    ptr = unserialize( i, ptr );
  }
  return ptr;
}



//get coordinates
double*
SimpleCell::serialize( double* xi ) const
{
  double* ptr = xi;
  int ncompo = mols.size();
  for( int i=0; i<ncompo; i++ ){
    ptr = serialize( i, ptr );
  }
  return ptr;
}



//get force vector
double*
SimpleCell::serializeforce( double* xi ) const
{
  double* ptr = xi;
  int ncompo = mols.size();
  for( int i=0; i<ncompo; i++ ){
    ptr = serializeforce( i, ptr );
  }
  return ptr;
}



//set coordinates
double*
SimpleCell::unserialize( int compo, double* const p )
{
  double* ptr = p;
  ptr = mols[compo]->unserialize( p );
  return ptr;
}



//get coordinates
double*
SimpleCell::serialize( int compo, double* xi ) const
{
  double* ptr = xi;
  ptr = mols[compo]->serialize( xi );
  return ptr;
}



//get force vector
double*
SimpleCell::serializeforce( int compo, double* xi ) const
{
  double* ptr = xi;
  ptr = mols[compo]->serializeforce( xi );
  return ptr;
}



void SimpleCell::ProgressPosition( double dt )
{
  int ncompo = mols.size();
  for( int i=0;i<ncompo;i++){
    MonatomicMols* m = dynamic_cast<MonatomicMols*>( mols[i].get() );
    if ( m ){
      m->ProgressPosition( dt );
    }
    else{
      RigidBodies* m = dynamic_cast<RigidBodies*>( mols[i].get() );
      if ( m ){
        m->ProgressPosition( dt );
      }
      else{
        RigidBodies2* m = dynamic_cast<RigidBodies2*>( mols[i].get() );
        assert( m );
        m->ProgressPosition( dt );
      }
    }
  }
}



void SimpleCell::ProgressMomentum( double dt )
{
  int ncompo = mols.size();
  for( int i=0;i<ncompo;i++){
    mols[i]->ProgressMomentum( dt );
  }
}




MolsHandle
SimpleCell::PullAll( const int compo )
{
  Box b(0,0,0);
  return mols[compo]->Emmigrate( b );
}



MolsHandle
SimpleCell::CopyAll( const int compo )
{
  MolsHandle mh = mols[compo]->EmptyClone();
  //ありあわせの関数だけで作ってみる。
  mh->Concat( mols[compo] );
  return mh;
}



void SimpleCell::TotalMomentum( Vector3& momentum ) const
{
  int ncompo = mols.size();
  momentum.Set(0,0,0);
  for( int i=0;i<ncompo;i++){
    Vector3 m;
    mols[i]->TotalMomentum( m );
    momentum.x += m.x;
    momentum.y += m.y;
    momentum.z += m.z;
  }
}



/*
void SimpleCell::AddVelocity( const Vector3& velocity )
{
  int ncompo = mols.size();
  for( int i=0;i<ncompo;i++){
    mols[i]->AddVelocity( velocity );
  }
}
*/



double SimpleCell::GetMass() const
{
  int ncompo = mols.size();
  double mass = 0;
  for( int i=0;i<ncompo;i++){
    mass += mols[i]->Size() * mols[i]->GetProperty()->GetMass();
  }
  return mass;
}

  

void
SimpleCell::SnapShot( const Unit& unit, const int compo, ostream& to )
{
  ///file入出力は枝葉末節の部分なので、極力objectに組みこまない。
  //to << mols[compo]->Container( Textify( mols[compo], unit ) );
}



void
SimpleCell::Rescale( const Vector3& r )
{
  ScalePositionPlugin p( r );
  p.HookL2( this );
}



RelocatableCell::RelocatableCell() : SimpleCell(), boundary( new Box(0,0,0) ), residents(0)
{
}



RelocatableCell::RelocatableCell( const Box& b, int ncompo, const Vector3 &o ) : SimpleCell( ncompo ), origin( o ), boundary( b.clone() ), residents(0)
{
}



RelocatableCell::~RelocatableCell()
{
  delete boundary;
  mesg("~RelocatableCell\n");
}



/*
void RelocatableCell::SetBox( const Box& b, const Vector3 &o )
{
  boundary = b;
  origin = o;
}
*/



//extract the molecules outside the boundary.
//emmigrated molecules are removed from the list.
MolsArray* RelocatableCell::Emmigrate(){
  int ncompo = mols.size();
  MolsArray* emmigrants = new MolsArray (ncompo);
  //reportsize();
  mesg("Emmigrate\n");
  //recalculate residents
  residents = 0;
  for( int i=0;i<ncompo;i++){
    MolsHandle& m = emmigrants->at(i);
    m = mols[i]->Emmigrate( *boundary );
    //printf("Top order+: %d\n", m->order[0]);
    /*{
      int nmol = m->size();
      for( int j=0; j<nmol; j++ ){
      printf("%d ", m->order[j]);
      }
      printf("\n");
      }*/
    /*絶対座標で返す。ただし、外boxの内にある保証はない。*/
    m->Translate( origin );
    residents += mols[i]->Size();
  }
  mesg("Emmigrated.\n");
  return emmigrants;
}




int RelocatableCell::Force( const Combination& combi, const Truncation& rc, PotVir &pv )
{
  SemiperiodicBox* sbox = dynamic_cast<SemiperiodicBox*>( boundary );
  if ( 0 == sbox ){
    return SimpleCell::Force( combi, rc, pv );
  }
  pv.Clear();
  int executed = 0;
  int ncombi = combi.ima.size();
  for( int co=0; co<ncombi; co++ ){
    int compo1 = combi.compo1[co];
    int compo2 = combi.compo2[co];
    PotVir p;
    Intersite* im = combi.ima[co].get();
    int result;
    if ( compo1 == compo2 ){
      result = mols[compo1]->Force_PBC( *im, *boundary, rc, p );
    }
    else{
      result = mols[compo1]->Force_PBC( mols[compo2], *im, *boundary, rc, p );
    }
    if ( result ){
      executed ++;
      //printf("%d %d %f\n", compo1, compo2, p.ep );
      pv.Add(p);
    }
  }
  return executed;
}



int RelocatableCell::Force_Intercell( const Combination& combi, const Vector3 &o, const RelocatableCell& c, const Truncation& rc, PotVir &pv )
{
  if ( IsEmpty() )return 0;
  if ( c.IsEmpty() )return 0;
  pv.Clear();
  int executed = 0;
  int ncombi = combi.ima.size();
  SemiperiodicBox* sbox = dynamic_cast<SemiperiodicBox*>( boundary );
  if ( sbox ){
    for( int co=0; co<ncombi; co++ ){
      PotVir p;
      int compo1 = combi.compo1[co];
      int compo2 = combi.compo2[co];
      Intersite* im = combi.ima[co].get();
      if ( mols[compo1]->Force_OffsetPBC( c.mols[compo2], *im, o, *boundary, rc, p ) ){
        executed ++;
	pv.Add(p);
      }
      if ( compo1 != compo2 ){
        /*
         *When two components are different, the following two differs.
         *interaction between the compo 1 in cell 1 and compo2 in cell 2
         *interaction between the compo 1 in cell 2 and compo2 in cell 1
         */
        if ( mols[compo2]->Force_OffsetPBC( c.mols[compo1], *im, o, *boundary, rc, p ) ){
          executed ++;
	  pv.Add(p);
        }
      }
    }
  }
  else{ // !sbox
    for( int co=0; co<ncombi; co++ ){
      PotVir p;
      int compo1 = combi.compo1[co];
      int compo2 = combi.compo2[co];
      Intersite* im = combi.ima[co].get();
      if ( mols[compo1]->Force_Offset( c.mols[compo2], *im, o, rc, p ) ){
        executed ++;
	pv.Add(p);
      }
      if ( compo1 != compo2 ){
        /*
         *When two components are different, the following two differs.
         *interaction between the compo 1 in cell 1 and compo2 in cell 2
         *interaction between the compo 1 in cell 2 and compo2 in cell 1
         */
        if ( mols[compo2]->Force_Offset( c.mols[compo1], *im, o, rc, p ) ){
          executed ++;
	  pv.Add(p);
        }
      }
    }
  }
  return executed;
}



int RelocatableCell::PairProcess( PairProcessPlugin& p, const Combination& combi, const Truncation& rc ) const
{
  SemiperiodicBox* sbox = dynamic_cast<SemiperiodicBox*>( boundary );
  if ( 0 == sbox ){
    return SimpleCell::PairProcess( p, combi, rc );
  }
  int ncombi = combi.ima.size();
  for( int co=0; co<ncombi; co++ ){
    int compo1 = combi.compo1[co];
    int compo2 = combi.compo2[co];
    Intersite* im = combi.ima[co].get();
    p.SetCombi( co );
    if ( compo1 == compo2 ){
      pairprocess_pbc( p, *(mols[compo1]), *im, *boundary, rc );
    }
    else{
      pairprocess_pbc( p, *(mols[compo1]), *(mols[compo2]), *im, *boundary, rc );
    }
  }
  return 1;
}



int RelocatableCell::PairProcess_Intercell( PairProcessPlugin& p, const Combination& combi, const Vector3 &o, const RelocatableCell& c, const Truncation& rc ) const
{
  if ( IsEmpty() )return 0;
  if ( c.IsEmpty() )return 0;
  int executed = 0;
  int ncombi = combi.ima.size();
  SemiperiodicBox* sbox = dynamic_cast<SemiperiodicBox*>( boundary );
  if ( sbox ){
    for( int co=0; co<ncombi; co++ ){
      int compo1 = combi.compo1[co];
      int compo2 = combi.compo2[co];
      Intersite* im = combi.ima[co].get();
      if ( pairprocess_offsetpbc( p, *(mols[compo1]), *(c.mols[compo2]), *im, o, *boundary, rc ) ){
        executed ++;
      }
      if ( compo1 != compo2 ){
        /*
         *When two components are different, the following two differs.
         *interaction between the compo 1 in cell 1 and compo2 in cell 2
         *interaction between the compo 1 in cell 2 and compo2 in cell 1
         */
        if ( pairprocess_offsetpbc( p, *(mols[compo2]), *(c.mols[compo1]), *im, o, *boundary, rc ) ){
          executed ++;
        }
      }
    }
  }
  else{
    for( int co=0; co<ncombi; co++ ){
      int compo1 = combi.compo1[co];
      int compo2 = combi.compo2[co];
      Intersite* im = combi.ima[co].get();
      if ( pairprocess_offset( p, *(mols[compo1]), *(c.mols[compo2]), *im, o, rc ) ){
        executed ++;
      }
      if ( compo1 != compo2 ){
        /*
         *When two components are different, the following two differs.
         *interaction between the compo 1 in cell 1 and compo2 in cell 2
         *interaction between the compo 1 in cell 2 and compo2 in cell 1
         */
        if ( pairprocess_offset( p, *(mols[compo2]), *(c.mols[compo1]), *im, o, rc ) ){
          executed ++;
        }
      }
    }
  }
  return executed;
}



const Vector3&
RelocatableCell::GetOrigin() const
{
  return origin;
}


/*
 *Push a molecule, which is in absolute coorditane system.
 */
void
RelocatableCell::PushAbs( int compo, const SingleMolHandle &m )
{
  //mesg("push6\n");
  Vector3 orig,rel;
  orig = m->Position();
  Vector3 offset;
  offset.x = -origin.x;
  offset.y = -origin.y;
  offset.z = -origin.z;
  m->Translate( offset );
  rel = m->Position();
  /*
  fprintf(stderr, "%d %f,%f,%f = %f,%f,%f + %f,%f,%f\n",
          m->GetOrder(), orig.x, orig.y, orig.z,
          origin.x, origin.y, origin.z,
          rel.x, rel.y, rel.z );
  */
  //cerr << m->GetOrder() << endl;
  mols[compo]->Push( m );
  residents ++;
}




MolsHandle
RelocatableCell::PullAll( const int compo )
{
  MolsHandle m = SimpleCell::PullAll( compo );
  m->Translate( origin );
  return m;
}



void RelocatableCell::Preforce()
{
  int ncompo = mols.size();
  for( int i=0;i<ncompo;i++ ){
    mols[i]->BoxCoordinate( *boundary );
  }
  SimpleCell::Preforce();
}



void
RelocatableCell::Rescale( const Vector3& r )
{
  boundary->Scale( r );
  origin.Scale( r );

  ScalePositionPlugin p( r );
  p.HookL2( this );
}
