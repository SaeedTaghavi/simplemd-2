#include <cassert>
#include <valarray>
#include "Mols/RigidBodies.hpp"
#include "debug.hpp"
#include "Plugin/CollectorPlugin.hpp"

void
RigidBodies::SetProperty( const MolPropertyHandle& p )
{
  prop = p;
  Rigid* r = dynamic_cast<Rigid*> ( p.get() );
  assert( r != 0 );
}



MolPropertyHandle
RigidBodies::GetProperty() const
{
  return prop;
}



int
RigidBodies::Size() const
{
  return mols.size();
}



RigidBodies::RigidBodies( const MolPropertyHandle& p )
{
    SetProperty( p );
}



MolsHandle
RigidBodies::EmptyClone() const
{
    MolsHandle h = MolsHandle( new RigidBodies( prop ) );
  return h;
}



RigidBodies::~RigidBodies()
{
  mesg("~RigidBodies\n");
  //Arrays are freed in ~MonatomicMols3()

}



int RigidBodies::push1( const RigidBody &mol )
{
  mols.push_back( mol );
  return Size();
}



int RigidBodies::Push( SingleMolHandle mol )
{
  RigidBody* m = dynamic_cast<RigidBody*> ( mol.get() );
  assert( m != 0 );
  return push1( *m );
}



/*
 *Peek the specified element of RigidBodies.
 *Returns a ptr to RigidBody, which must be deleted later.
 */
RigidBody*
RigidBodies::peek1( int which ) const
{
  RigidBody* mol = new RigidBody( mols[which] );
  return mol;
}



//simple reference to a molecule
const SingleMolEntity&
RigidBodies::Molecule( int i ) const
{
  return mols[i];
}



SingleMolHandle
RigidBodies::Peek( int which ) const
{
  return SingleMolHandle(peek1(which));
}



/*
 *Pull out the specified element of RigidBodies.
 *The element is removed.
 *Returns RigidBody, which must be deleted later.
 */
RigidBody*
RigidBodies::pull( int which )
{
  RigidBody* mol = peek1( which );
  if ( mol == 0 ) exit(1);
  mols.erase( mols.begin() + which );
  return mol;
}



MolsHandle
RigidBodies::Emmigrate( const Box &box )
{
  RigidBodies* emmigrants = new RigidBodies( prop );
  mesg("RigidBodies3::emmigrate\n");
  int i=0;
  while( i < Size() ){
    Vector3 com = mols[i].Position();
    if ( box.Contains( com ) ){
      i++;
    }
    else{
      // copy the molecule temporarily.
      RigidBody* mol = pull( i );
      emmigrants->push1( *mol );
      delete mol;
    }
  }
  return MolsHandle( emmigrants );
}



double RigidBodies::GetEk() const
{
  int nmol  = Size();
  double ek = 0;
  for( int i=0; i<nmol; i++ ){
    ek += mols[i].GetEk();
  }
  return ek;
}



void RigidBodies::GetEkt( Vector3& ekt ) const
{
  MolPropertyHandle h = GetProperty();
  Rigid*       e = dynamic_cast<Rigid*> ( h.get() );
  assert( e != 0 );
  double mass = e->GetMass();
  int nmol  = Size();
  ekt.Set(0,0,0);
  for( int i=0; i<nmol; i++ ){
    Vector3 velocity = mols[i].com.GetVelocity();
    ekt.x += velocity.x * velocity.x;
    ekt.y += velocity.y * velocity.y;
    ekt.z += velocity.z * velocity.z;
    //ek += velocity.square();
  }
  ekt.x *= mass * 0.5;
  ekt.y *= mass * 0.5;
  ekt.z *= mass * 0.5;
}



void
RigidBodies::TotalMomentum( Vector3& momentum ) const
{
  MolPropertyHandle h = GetProperty();
  Rigid*       e = dynamic_cast<Rigid*> ( h.get() );
  assert( e != 0 );
  double mass = e->GetMass();
  int nmol  = Size();

  momentum.Set(0,0,0);
  for( int i=0; i<nmol; i++ ){
    Vector3 velocity = mols[i].com.GetVelocity();
    momentum.x += velocity.x;
    momentum.y += velocity.y;
    momentum.z += velocity.z;
  }
  momentum.x *= mass;
  momentum.y *= mass;
  momentum.z *= mass;
}




/*
void RigidBodies::ScaleVelocity( const Vector3& r )
{
  int nmol = Size();
  for( int i=0; i < nmol; i++ ){
    mols[i].ScaleVelocity( r );
  }
}
*/



/*
void RigidBodies::ScalePosition( const Vector3& r )
{
  int nmol = Size();
  for( int i=0; i < nmol; i++ ){
    mols[i].ScalePosition( r );
  }
}
*/





void RigidBodies::Preforce()
{
  int nmol = Size();
  for( int i=0; i<nmol; i++ ) {
    mols[i].Preforce();
  }
}




void RigidBodies::BoxCoordinate( const Box& box )
{
  int nmol = Size();
  for( int i=0; i<nmol; i++ ){
    mols[i].BoxCoordinate( box );
  }
}





void RigidBodies::Postforce()
{
    MolPropertyHandle h = GetProperty();
    Rigid* e = dynamic_cast<Rigid*> ( h.get() );
    assert( e != 0 );
    double mass   = e->GetMass();

    int nmol = Size();
    for( int i=0; i<nmol; i++ ){
        mols[i].Postforce( mass );
    }
}








void 
RigidBodies::Translate( const Vector3& offset )
{
    int nmol = Size();
    for( int j=0; j<nmol; j++ ){
        mols[j].Translate( offset );
    }
}


/*
void RigidBodies::Report() const
{
}
*/



void 
RigidBodies::ProgressMomentum( double dt )
{
  int nmol = Size();
  for( int j=0; j<nmol; j++ ){
    mols[j].ProgressMomentum( dt );
  }
}



void 
RigidBodies::ProgressPosition( double dt )
{
  int nmol = Size();
  for( int j=0; j<nmol; j++ ){
    mols[j].ProgressPosition( dt );
  }
}










/*
 *試作1：リストベクトルの事前形成。速度にはほとんど貢献しない。
 */











/*
int
force_common_1st(
  RigidBodies& m1,
  RigidBodies& m2,
  const Intersite& im,
  const ListVector& truncpair,
  PotVir &pv
  )
{
  epsum = 0;
  vrsum.Clear();
  int executed = 0;
  int npair = truncpair.size();
  for(int k=0; k<npair; k++){
    const ListItem& li = truncpair.listitems[k];
    int i = li.pair1;
    int j = li.pair2;
    RigidBody& r1 = m1.mols[i];
    RigidBody& r2 = m2.mols[j];
    PotVir pv0;
    if ( force_core_core( r1, r2, im, li.d, pv0, li.sfe, li.sff ) ){
      executed ++;
      epsum += pv0.ep;
      for(int i=0;i<3;i++){
	for(int j=0;j<3;j++){
	  vrsum.v[i][j] += pv0.vr.v[i][j];
	}
      }
    }
  }
  return executed;
}
*/


/*
 *試作2：サイト対のループを外に出す。
 */
int
force_common_2nd(
  RigidBodies& m1,
  RigidBodies& m2,
  const Intersite& im,
  const TruncPair& truncpair,
  PotVir &pv
  )
{
  //exit(0);
  double epsum = 0;
  Matrix33 vrsum;
  //pv.Clear();
  int truncpairsize = truncpair.size();
  int npair = im.intr.size();
  //valarray<double> ep(0.0, truncpairsize);
  vector<double> ep(truncpairsize, 0.0);
  vector<Vector3> ff(truncpairsize, Vector3(0,0,0));
  for( int k=0; k<npair; k++ ){
    int s1 = im.site1[k];
    int s2 = im.site2[k];
    for(int kk=0; kk<truncpairsize; kk++){
      //const ListItem& li = truncpair.listitems[kk];
      const ListItem& li = truncpair.GetListItems()[kk];
      int i = li.pair1;
      int j = li.pair2;
      RigidBody& r1 = m1.mols[i];
      RigidBody& r2 = m2.mols[j];
      SiteOfAction& soa1 = r1.atom[s1].center;
      SiteOfAction& soa2 = r2.atom[s2].center;
      Vector3 delta;
      double fo;
      double potential;
      delta.x  = li.d.x + soa1.coord.x - soa2.coord.x;
      delta.y  = li.d.y + soa1.coord.y - soa2.coord.y;
      delta.z  = li.d.z + soa1.coord.z - soa2.coord.z;
      im.intr[k]->Force( delta, fo, potential );
      ep[kk] += potential;
      //printf("%f %f %f %f %f %f\n", delta.x, delta.y, delta.z, fo, potential,ep[kk] );
      fo *= li.sfe;
      double ffx = delta.x*fo;
      double ffy = delta.y*fo;
      double ffz = delta.z*fo;
      soa1.force.x += ffx;
      soa1.force.y += ffy;
      soa1.force.z += ffz;
      soa2.force.x -= ffx;
      soa2.force.y -= ffy;
      soa2.force.z -= ffz;
      ff[kk].x += ffx;
      ff[kk].y += ffy;
      ff[kk].z += ffz;
    }
  }
  for(int kk=0; kk<truncpairsize; kk++){
    const ListItem& li = truncpair.GetListItems()[kk];
    int i = li.pair1;
    int j = li.pair2;
    RigidBody& r1 = m1.mols[i];
    RigidBody& r2 = m2.mols[j];
    double comfx, comfy, comfz;
    //force induced by interaction truncation
    comfx = -li.sff * ep[kk] * li.d.x;
    comfy = -li.sff * ep[kk] * li.d.y;
    comfz = -li.sff * ep[kk] * li.d.z;
    r1.com.center.force.x += comfx;
    r1.com.center.force.y += comfy;
    r1.com.center.force.z += comfz;
    r2.com.center.force.x -= comfx;
    r2.com.center.force.y -= comfy;
    r2.com.center.force.z -= comfz;
    epsum         += li.sfe * ep[kk];
    comfx += ff[kk].x;
    comfy += ff[kk].y;
    comfz += ff[kk].z;
    vrsum.v[0][0] += li.d.x * comfx;
    vrsum.v[0][1] += li.d.x * comfy;
    vrsum.v[0][2] += li.d.x * comfz;
    vrsum.v[1][0] += li.d.y * comfx;
    vrsum.v[1][1] += li.d.y * comfy;
    vrsum.v[1][2] += li.d.y * comfz;
    vrsum.v[2][0] += li.d.z * comfx;
    vrsum.v[2][1] += li.d.z * comfy;
    vrsum.v[2][2] += li.d.z * comfz;
  }
  //cout << epsum << "  epsum\n";
  pv.Set( epsum, vrsum );
  return 1;
}





/*
 *試作3：さらに、相互作用ごとにループを分ける。
 *polymorphismを使わない分速くなると思ったのだが、遅くなってしまった。
 */
int
force_common_3rd(
  RigidBodies& m1,
  RigidBodies& m2,
  const Intersite& im,
  const TruncPair& truncpair,
  PotVir &pv
  )
{
  double epsum = 0;
  Matrix33 vrsum;
  //pv.Clear();
  int truncpairsize = truncpair.size();
  int npair = im.intr.size();
  valarray<double> ep(0.0, truncpairsize);
  vector<Vector3> ff(truncpairsize, Vector3(0,0,0));
  for( int k=0; k<npair; k++ ){
    int s1 = im.site1[k];
    int s2 = im.site2[k];
    {
      Coulomb* intr = dynamic_cast<Coulomb*>( im.intr[k].get() );
      if ( intr ){
        double sqrch = intr->sqrch;
        for(int kk=0; kk<truncpairsize; kk++){
          //const ListItem& li = truncpair.listitems[kk];
	  const ListItem& li = truncpair.GetListItems()[kk];
          int i = li.pair1;
          int j = li.pair2;
          RigidBody& r1 = m1.mols[i];
          RigidBody& r2 = m2.mols[j];
          SiteOfAction& soa1 = r1.atom[s1].center;
          SiteOfAction& soa2 = r2.atom[s2].center;
          double dx,dy,dz;
          double fo;
          double potential;
          dx  = li.d.x + soa1.coord.x - soa2.coord.x;
          dy  = li.d.y + soa1.coord.y - soa2.coord.y;
          dz  = li.d.z + soa1.coord.z - soa2.coord.z;
          //intr->Force( d, fo, potential );
          double rri = 1.0 / (dx*dx + dy*dy + dz*dz);
          double radiusi = sqrt(rri);
          potential = sqrch *radiusi;                
          fo = potential * rri;
          ep[kk] += potential;
          double ffx = dx*fo*li.sfe;
          double ffy = dy*fo*li.sfe;
          double ffz = dz*fo*li.sfe;
          soa1.force.x += ffx;
          soa1.force.y += ffy;
          soa1.force.z += ffz;
          soa2.force.x -= ffx;
          soa2.force.y -= ffy;
          soa2.force.z -= ffz;
          ff[kk].x += ffx;
          ff[kk].y += ffy;
          ff[kk].z += ffz;
        }
      }
    }
//  }
//  for( int k=0; k<npair; k++ ){
//    int s1 = im.site1[k];
//    int s2 = im.site2[k];
    {
      LJ* intr = dynamic_cast<LJ*>( im.intr[k].get() );
      if ( intr ){
        double sig2 = intr->sig * intr->sig;
        double eps4 = 4.0 * intr->eps;
        for(int kk=0; kk<truncpairsize; kk++){
	  const ListItem& li = truncpair.GetListItems()[kk];
          //const ListItem& li = truncpair.listitems[kk];
          int i = li.pair1;
          int j = li.pair2;
          RigidBody& r1 = m1.mols[i];
          RigidBody& r2 = m2.mols[j];
          SiteOfAction& soa1 = r1.atom[s1].center;
          SiteOfAction& soa2 = r2.atom[s2].center;
          //Vector3 d;
          double dx,dy,dz;
          double fo;
          double potential;
          dx  = li.d.x + soa1.coord.x - soa2.coord.x;
          dy  = li.d.y + soa1.coord.y - soa2.coord.y;
          dz  = li.d.z + soa1.coord.z - soa2.coord.z;
          //intr->Force( d, fo, potential );
          double rri = 1.0 / (dx*dx + dy*dy + dz*dz);
          double drs = sig2 * rri;
          double dr6 = drs * drs * drs;
          double dr12 = dr6 * dr6;
          fo = -eps4 * (-12.0 * dr12 + 6.0 * dr6) * rri;
          potential = eps4 * (dr12-dr6);                
          ep[kk] += potential;
          double ffx = dx*fo*li.sfe;
          double ffy = dy*fo*li.sfe;
          double ffz = dz*fo*li.sfe;
          soa1.force.x += ffx;
          soa1.force.y += ffy;
          soa1.force.z += ffz;
          soa2.force.x -= ffx;
          soa2.force.y -= ffy;
          soa2.force.z -= ffz;
          ff[kk].x += ffx;
          ff[kk].y += ffy;
          ff[kk].z += ffz;
        }
      }
    }
//  }
//  for( int k=0; k<npair; k++ ){
//    int s1 = im.site1[k];
//    int s2 = im.site2[k];
    {
      LJC* intr = dynamic_cast<LJC*>( im.intr[k].get() );
      if ( intr ){
        double sig2 = intr->lj.sig * intr->lj.sig;
        double eps4 = 4.0 * intr->lj.eps;
        for(int kk=0; kk<truncpairsize; kk++){
	  const ListItem& li = truncpair.GetListItems()[kk];
          //const ListItem& li = truncpair.listitems[kk];
          int i = li.pair1;
          int j = li.pair2;
          RigidBody& r1 = m1.mols[i];
          RigidBody& r2 = m2.mols[j];
          SiteOfAction& soa1 = r1.atom[s1].center;
          SiteOfAction& soa2 = r2.atom[s2].center;
          //Vector3 d;
          double dx,dy,dz;
          double fo;
          double potential;
          dx  = li.d.x + soa1.coord.x - soa2.coord.x;
          dy  = li.d.y + soa1.coord.y - soa2.coord.y;
          dz  = li.d.z + soa1.coord.z - soa2.coord.z;
          //intr->Force( d, fo, potential );
          double rri = 1.0 / (dx*dx + dy*dy + dz*dz);
          double radiusi = sqrt(rri);
          double drs = sig2 * rri;
          double dr6 = drs * drs * drs;
          double dr12 = dr6 * dr6;
          double csr = intr->coulomb.sqrch*radiusi;
          fo =  rri*( -eps4 * (-12.0 * dr12 + 6.0 * dr6) + csr);
          potential = eps4 * (dr12-dr6) + csr;
          ep[kk] += potential;
          double ffx = dx*fo*li.sfe;
          double ffy = dy*fo*li.sfe;
          double ffz = dz*fo*li.sfe;
          soa1.force.x += ffx;
          soa1.force.y += ffy;
          soa1.force.z += ffz;
          soa2.force.x -= ffx;
          soa2.force.y -= ffy;
          soa2.force.z -= ffz;
          ff[kk].x += ffx;
          ff[kk].y += ffy;
          ff[kk].z += ffz;
        }
      }
    }
  }
  for(int kk=0; kk<truncpairsize; kk++){
    //const ListItem& li = truncpair.listitems[kk];
    const ListItem& li = truncpair.GetListItems()[kk];
    int i = li.pair1;
    int j = li.pair2;
    RigidBody& r1 = m1.mols[i];
    RigidBody& r2 = m2.mols[j];
    double comfx, comfy, comfz;
    //force induced by interaction truncation
    comfx = -li.sff * ep[kk] * li.d.x;
    comfy = -li.sff * ep[kk] * li.d.y;
    comfz = -li.sff * ep[kk] * li.d.z;
    r1.com.center.force.x += comfx;
    r1.com.center.force.y += comfy;
    r1.com.center.force.z += comfz;
    r2.com.center.force.x -= comfx;
    r2.com.center.force.y -= comfy;
    r2.com.center.force.z -= comfz;
    epsum         += li.sfe * ep[kk];
    comfx += ff[kk].x;
    comfy += ff[kk].y;
    comfz += ff[kk].z;
    vrsum.v[0][0] += li.d.x * comfx;
    vrsum.v[0][1] += li.d.x * comfy;
    vrsum.v[0][2] += li.d.x * comfz;
    vrsum.v[1][0] += li.d.y * comfx;
    vrsum.v[1][1] += li.d.y * comfy;
    vrsum.v[1][2] += li.d.y * comfz;
    vrsum.v[2][0] += li.d.z * comfx;
    vrsum.v[2][1] += li.d.z * comfy;
    vrsum.v[2][2] += li.d.z * comfz;
  }
  pv.Set( epsum, vrsum );
  return 1;
}





int
force_offset(
             RigidBodies& m1,
             RigidBodies& m2,
	     const Intersite& im,
             const Vector3 &offset, 
             const Truncation& rc,
             PotVir &pv
             )
{
  TruncPair truncpair;
  if ( truncpair.list_offset( m1, m2, offset, rc ) ){
    //force_common_1st( m1, m2, im, truncpair, pv );
    return force_common_2nd( m1, m2, im, truncpair, pv );
    //force_common_3rd( m1, m2, im, truncpair, pv );
  }
  return 0;
}



int
force_offsetpbc(
             RigidBodies& m1,
             RigidBodies& m2,
	     const Intersite& im,
             const Vector3 &offset, 
	     const Box& box,
             const Truncation& rc,
             PotVir &pv
             )
{
  TruncPair truncpair;
  if ( truncpair.list_offsetpbc( m1, m2, offset, box, rc ) ){
    //force_common_1st( m1, m2, im, truncpair, pv );
    return force_common_2nd( m1, m2, im, truncpair, pv );
    //force_common_3rd( m1, m2, im, truncpair, pv );
  }
  return 0;
}



int
force_simple(
             RigidBodies& m1,
             RigidBodies& m2,
	     const Intersite& im,
             const Truncation& rc,
             PotVir &pv
             )
{
  TruncPair truncpair;
  if ( truncpair.list_simple( m1, m2, rc ) ){
    //force_common_1st( m1, m2, im, truncpair, pv );
    return force_common_2nd( m1, m2, im, truncpair, pv );
    //force_common_3rd( m1, m2, im, truncpair, pv );
  }
  return 0;
}



int
force_pbc(
             RigidBodies& m1,
             RigidBodies& m2,
	     const Intersite& im,
             const Box& box,
             const Truncation& rc,
             PotVir &pv
             )
{
  TruncPair truncpair;
  if ( truncpair.list_pbc( m1, m2, box, rc ) ){
    //force_common_1st( m1, m2, im, truncpair, pv );
    return force_common_2nd( m1, m2, im, truncpair, pv );
    //force_common_3rd( m1, m2, im, truncpair, pv );
  }
  return 0;
}





//force w/o PBC
int RigidBodies::Force( const Intersite& im, const Truncation& rc, PotVir &pv )
{
  TruncPair truncpair;
  if ( truncpair.list_simple( *this, rc ) ){
    //force_common_1st( *this, *this, im, truncpair, pv );
    return force_common_2nd( *this, *this, im, truncpair, pv );
    //force_common_3rd( *this, *this, im, truncpair, pv );
  }
  return 0;
}



//force w/ PBC
int RigidBodies::Force_PBC( const Intersite& im, const Box& box, const Truncation& rc, PotVir &pv )
{
  TruncPair truncpair;
  if ( truncpair.list_pbc( *this, box, rc ) ){
    //force_common_1st( *this, *this, im, truncpair, pv );
    return force_common_2nd( *this, *this, im, truncpair, pv );
    //force_common_3rd( *this, *this, im, truncpair, pv );
  }
  return 0;
}



void
RigidBodies::Write( const Unit& unit, ostream& to )
{
  WriteWTG5( unit, to );
}



void
RigidBodies::WriteWTG5( const Unit& unit, ostream& to )
{
  int nmol = Size();
  to << "@WTG5\n"
     << nmol << '\n';
  sort( mols.begin(), mols.end() );
  for(int i=0; i<nmol; i++){
    mols[i].WriteWTG5( unit, to );
  }
}



void
RigidBodies::SnapShot( const Unit& unit, ostream& to )
{
  WriteNX4A( unit, to );
}



void
RigidBodies::WriteNX4A( const Unit& unit, ostream& to )
{
  int nmol = Size();
  to << "@NX4A\n"
     << nmol << '\n';
  sort( mols.begin(), mols.end() );
  for(int i=0; i<nmol; i++){
    mols[i].WriteNX4A( unit, to );
  }
}





/* Heteromolecular interactions */


/* list vectorの処理は、Position()関数を使えば共通化できる(が遅くなるかも)*/







/*
 *試作2：サイト対のループを外に出す。
 */
int
force_common_2nd(
  RigidBodies& m1,
  MonatomicMols& m2,
  const Intersite& im,
  const TruncPair& truncpair,
  PotVir &pv
  )
{
  double epsum = 0;
  Matrix33 vrsum;
  //pv.Clear();
  int truncpairsize = truncpair.size();
  int npair = im.intr.size();
  valarray<double> ep(0.0, truncpairsize);
  vector<Vector3> ff(truncpairsize, Vector3(0,0,0));
  for( int k=0; k<npair; k++ ){
    int s1 = im.site1[k];
    //int s2 = im.site2[k];
    for(int kk=0; kk<truncpairsize; kk++){
      const ListItem& li = truncpair.GetListItems()[kk];
      //const ListItem& li = truncpair.listitems[kk];
      int i = li.pair1;
      int j = li.pair2;
      RigidBody& r1 = m1.mols[i];
      MonatomicMol& r2 = m2.mols[j];
      SiteOfAction& soa1 = r1.atom[s1].center;
      SiteOfAction& soa2 = r2.center;
      Vector3 delta;
      double fo;
      double potential;
      delta.x  = li.d.x + soa1.coord.x;
      delta.y  = li.d.y + soa1.coord.y;
      delta.z  = li.d.z + soa1.coord.z;
      im.intr[k]->Force( delta, fo, potential );
      //if ( potential > 10000 )
	//printf("%d %d %f %f %f %f %f\n", r1.order,r2.order, delta.x, delta.y, delta.z, fo, potential );
      ep[kk] += potential;
      fo *= li.sfe;
      double ffx = delta.x*fo;
      double ffy = delta.y*fo;
      double ffz = delta.z*fo;
      soa1.force.x += ffx;
      soa1.force.y += ffy;
      soa1.force.z += ffz;
      soa2.force.x -= ffx;
      soa2.force.y -= ffy;
      soa2.force.z -= ffz;
      ff[kk].x += ffx;
      ff[kk].y += ffy;
      ff[kk].z += ffz;
    }
  }
  for(int kk=0; kk<truncpairsize; kk++){
    const ListItem& li = truncpair.GetListItems()[kk];
    //const ListItem& li = truncpair.listitems[kk];
    int i = li.pair1;
    int j = li.pair2;
    RigidBody& r1 = m1.mols[i];
    MonatomicMol& r2 = m2.mols[j];
    double comfx, comfy, comfz;
    //force induced by interaction truncation
    comfx = -li.sff * ep[kk] * li.d.x;
    comfy = -li.sff * ep[kk] * li.d.y;
    comfz = -li.sff * ep[kk] * li.d.z;
    r1.com.center.force.x += comfx;
    r1.com.center.force.y += comfy;
    r1.com.center.force.z += comfz;
    r2.center.force.x -= comfx;
    r2.center.force.y -= comfy;
    r2.center.force.z -= comfz;
    epsum         += li.sfe * ep[kk];
    comfx += ff[kk].x;
    comfy += ff[kk].y;
    comfz += ff[kk].z;
    vrsum.v[0][0] += li.d.x * comfx;
    vrsum.v[0][1] += li.d.x * comfy;
    vrsum.v[0][2] += li.d.x * comfz;
    vrsum.v[1][0] += li.d.y * comfx;
    vrsum.v[1][1] += li.d.y * comfy;
    vrsum.v[1][2] += li.d.y * comfz;
    vrsum.v[2][0] += li.d.z * comfx;
    vrsum.v[2][1] += li.d.z * comfy;
    vrsum.v[2][2] += li.d.z * comfz;
  }
  pv.Set( epsum, vrsum );
  return 1;
}





/*
 *試作3：さらに、相互作用ごとにループを分ける。
 *polymorphismを使わない分速くなると思ったのだが、遅くなってしまった。
 */
int
force_common_3rd(
  RigidBodies& m1,
  MonatomicMols& m2,
  const Intersite& im,
  const TruncPair& truncpair,
  PotVir &pv
  )
{
  double epsum = 0;
  Matrix33 vrsum;
  //pv.Clear();
  int truncpairsize = truncpair.size();
  int npair = im.intr.size();
  valarray<double> ep(0.0, truncpairsize);
  vector<Vector3> ff(truncpairsize, Vector3(0,0,0));
  for( int k=0; k<npair; k++ ){
    int s1 = im.site1[k];
    //int s2 = im.site2[k];
    {
      Coulomb* intr = dynamic_cast<Coulomb*>( im.intr[k].get() );
      if ( intr ){
        double sqrch = intr->sqrch;
        for(int kk=0; kk<truncpairsize; kk++){
	  const ListItem& li = truncpair.GetListItems()[kk];
          //const ListItem& li = truncpair.listitems[kk];
          int i = li.pair1;
          int j = li.pair2;
          RigidBody& r1 = m1.mols[i];
          MonatomicMol& r2 = m2.mols[j];
          SiteOfAction& soa1 = r1.atom[s1].center;
          SiteOfAction& soa2 = r2.center;
          double dx,dy,dz;
          double fo;
          double potential;
          dx  = li.d.x + soa1.coord.x - soa2.coord.x;
          dy  = li.d.y + soa1.coord.y - soa2.coord.y;
          dz  = li.d.z + soa1.coord.z - soa2.coord.z;
          //intr->Force( d, fo, potential );
          double rri = 1.0 / (dx*dx + dy*dy + dz*dz);
          double radiusi = sqrt(rri);
          potential = sqrch *radiusi;                
          fo = potential * rri;
          ep[kk] += potential;
          double ffx = dx*fo*li.sfe;
          double ffy = dy*fo*li.sfe;
          double ffz = dz*fo*li.sfe;
          soa1.force.x += ffx;
          soa1.force.y += ffy;
          soa1.force.z += ffz;
          soa2.force.x -= ffx;
          soa2.force.y -= ffy;
          soa2.force.z -= ffz;
          ff[kk].x += ffx;
          ff[kk].y += ffy;
          ff[kk].z += ffz;
        }
      }
    }
//  }
//  for( int k=0; k<npair; k++ ){
//    int s1 = im.site1[k];
//    int s2 = im.site2[k];
    {
      LJ* intr = dynamic_cast<LJ*>( im.intr[k].get() );
      if ( intr ){
        double sig2 = intr->sig * intr->sig;
        double eps4 = 4.0 * intr->eps;
        for(int kk=0; kk<truncpairsize; kk++){
	  const ListItem& li = truncpair.GetListItems()[kk];
          //const ListItem& li = truncpair.listitems[kk];
          int i = li.pair1;
          int j = li.pair2;
          RigidBody& r1 = m1.mols[i];
          MonatomicMol& r2 = m2.mols[j];
          SiteOfAction& soa1 = r1.atom[s1].center;
          SiteOfAction& soa2 = r2.center;
          //Vector3 d;
          double dx,dy,dz;
          double fo;
          double potential;
          dx  = li.d.x + soa1.coord.x - soa2.coord.x;
          dy  = li.d.y + soa1.coord.y - soa2.coord.y;
          dz  = li.d.z + soa1.coord.z - soa2.coord.z;
          //intr->Force( d, fo, potential );
          double rri = 1.0 / (dx*dx + dy*dy + dz*dz);
          double drs = sig2 * rri;
          double dr6 = drs * drs * drs;
          double dr12 = dr6 * dr6;
          fo = -eps4 * (-12.0 * dr12 + 6.0 * dr6) * rri;
          potential = eps4 * (dr12-dr6);                
          ep[kk] += potential;
          double ffx = dx*fo*li.sfe;
          double ffy = dy*fo*li.sfe;
          double ffz = dz*fo*li.sfe;
          soa1.force.x += ffx;
          soa1.force.y += ffy;
          soa1.force.z += ffz;
          soa2.force.x -= ffx;
          soa2.force.y -= ffy;
          soa2.force.z -= ffz;
          ff[kk].x += ffx;
          ff[kk].y += ffy;
          ff[kk].z += ffz;
        }
      }
    }
//  }
//  for( int k=0; k<npair; k++ ){
//    int s1 = im.site1[k];
//    int s2 = im.site2[k];
    {
      LJC* intr = dynamic_cast<LJC*>( im.intr[k].get() );
      if ( intr ){
        double sig2 = intr->lj.sig * intr->lj.sig;
        double eps4 = 4.0 * intr->lj.eps;
        for(int kk=0; kk<truncpairsize; kk++){
	  const ListItem& li = truncpair.GetListItems()[kk];
          //const ListItem& li = truncpair.listitems[kk];
          int i = li.pair1;
          int j = li.pair2;
          RigidBody& r1 = m1.mols[i];
          MonatomicMol& r2 = m2.mols[j];
          SiteOfAction& soa1 = r1.atom[s1].center;
          SiteOfAction& soa2 = r2.center;
          //Vector3 d;
          double dx,dy,dz;
          double fo;
          double potential;
          dx  = li.d.x + soa1.coord.x - soa2.coord.x;
          dy  = li.d.y + soa1.coord.y - soa2.coord.y;
          dz  = li.d.z + soa1.coord.z - soa2.coord.z;
          //intr->Force( d, fo, potential );
          double rri = 1.0 / (dx*dx + dy*dy + dz*dz);
          double radiusi = sqrt(rri);
          double drs = sig2 * rri;
          double dr6 = drs * drs * drs;
          double dr12 = dr6 * dr6;
          double csr = intr->coulomb.sqrch*radiusi;
          fo =  rri*( -eps4 * (-12.0 * dr12 + 6.0 * dr6) + csr);
          potential = eps4 * (dr12-dr6) + csr;
          ep[kk] += potential;
          double ffx = dx*fo*li.sfe;
          double ffy = dy*fo*li.sfe;
          double ffz = dz*fo*li.sfe;
          soa1.force.x += ffx;
          soa1.force.y += ffy;
          soa1.force.z += ffz;
          soa2.force.x -= ffx;
          soa2.force.y -= ffy;
          soa2.force.z -= ffz;
          ff[kk].x += ffx;
          ff[kk].y += ffy;
          ff[kk].z += ffz;
        }
      }
    }
  }
  for(int kk=0; kk<truncpairsize; kk++){
    const ListItem& li = truncpair.GetListItems()[kk];
    //const ListItem& li = truncpair.listitems[kk];
    int i = li.pair1;
    int j = li.pair2;
    RigidBody& r1 = m1.mols[i];
    MonatomicMol& r2 = m2.mols[j];
    double comfx, comfy, comfz;
    //force induced by interaction truncation
    comfx = -li.sff * ep[kk] * li.d.x;
    comfy = -li.sff * ep[kk] * li.d.y;
    comfz = -li.sff * ep[kk] * li.d.z;
    r1.com.center.force.x += comfx;
    r1.com.center.force.y += comfy;
    r1.com.center.force.z += comfz;
    r2.center.force.x -= comfx;
    r2.center.force.y -= comfy;
    r2.center.force.z -= comfz;
    epsum         += li.sfe * ep[kk];
    comfx += ff[kk].x;
    comfy += ff[kk].y;
    comfz += ff[kk].z;
    vrsum.v[0][0] += li.d.x * comfx;
    vrsum.v[0][1] += li.d.x * comfy;
    vrsum.v[0][2] += li.d.x * comfz;
    vrsum.v[1][0] += li.d.y * comfx;
    vrsum.v[1][1] += li.d.y * comfy;
    vrsum.v[1][2] += li.d.y * comfz;
    vrsum.v[2][0] += li.d.z * comfx;
    vrsum.v[2][1] += li.d.z * comfy;
    vrsum.v[2][2] += li.d.z * comfz;
  }
  pv.Set( epsum, vrsum );
  return 1;
}





int
force_offset(
             RigidBodies& m1,
             MonatomicMols& m2,
	     const Intersite& im,
             const Vector3 &offset, 
             const Truncation& rc,
             PotVir &pv
             )
{
  TruncPair truncpair;
  if ( truncpair.list_offset( m1, m2, offset, rc ) ){
    //force_common_1st( m1, m2, im, truncpair, pv );
    return force_common_2nd( m1, m2, im, truncpair, pv );
    //force_common_3rd( m1, m2, im, truncpair, pv );
  }
  return 0;
}



int
force_offsetpbc(
             RigidBodies& m1,
             MonatomicMols& m2,
	     const Intersite& im,
             const Vector3 &offset, 
	     const Box& box,
             const Truncation& rc,
             PotVir &pv
             )
{
  TruncPair truncpair;
  if ( truncpair.list_offsetpbc( m1, m2, offset, box, rc ) ){
    //force_common_1st( m1, m2, im, truncpair, pv );
    return force_common_2nd( m1, m2, im, truncpair, pv );
    //force_common_3rd( m1, m2, im, truncpair, pv );
  }
  return 0;
}



int
force_simple(
             RigidBodies& m1,
             MonatomicMols& m2,
	     const Intersite& im,
             const Truncation& rc,
             PotVir &pv
             )
{
  TruncPair truncpair;
  if ( truncpair.list_simple( m1, m2, rc ) ){
    //force_common_1st( m1, m2, im, truncpair, pv );
    return force_common_2nd( m1, m2, im, truncpair, pv );
    //force_common_3rd( m1, m2, im, truncpair, pv );
  }
  return 0;
}



int
force_pbc(
             RigidBodies& m1,
             MonatomicMols& m2,
	     const Intersite& im,
             const Box& box,
             const Truncation& rc,
             PotVir &pv
             )
{
  TruncPair truncpair;
  if ( truncpair.list_pbc( m1, m2, box, rc ) ){
    //force_common_1st( m1, m2, im, truncpair, pv );
    return force_common_2nd( m1, m2, im, truncpair, pv );
    //force_common_3rd( m1, m2, im, truncpair, pv );
  }
  return 0;
}





//select force functions dynamically
int RigidBodies::Force_Simple(
			       const MolsHandle& m,
			       const Intersite& im,
			       const Truncation& rc,
			       PotVir &pv
			       )
{
  RigidBodies* mm = dynamic_cast<RigidBodies*> ( m.get() );
  if ( mm != 0 ){
    return ::force_simple( *this, *mm, im, rc, pv );
  }
  else{
    MonatomicMols* mm2 = dynamic_cast<MonatomicMols*> ( m.get() );
    assert ( mm2 != 0 );
    return ::force_simple( *this, *mm2, im, rc, pv );
  }
  return 0;
}




//select force functions dynamically
int RigidBodies::Force_Offset(
			       const MolsHandle& m,
			       const Intersite& im,
			       const Vector3 &offset,
			       const Truncation& rc,
			       PotVir &pv
			       )
{
  RigidBodies* mm = dynamic_cast<RigidBodies*> ( m.get() );
  if ( mm != 0 ){
    return ::force_offset( *this, *mm, im, offset, rc, pv );
  }
  else{
    MonatomicMols* mm2 = dynamic_cast<MonatomicMols*> ( m.get() );
    assert ( mm2 != 0 );
    return ::force_offset( *this, *mm2, im, offset, rc, pv );
  }
  return 0;
}



//select force functions dynamically
int RigidBodies::Force_OffsetPBC(
			       const MolsHandle& m,
			       const Intersite& im,
			       const Vector3 &offset,
			       const Box& box,
			       const Truncation& rc,
			       PotVir &pv
			       )
{
  RigidBodies* mm = dynamic_cast<RigidBodies*> ( m.get() );
  if ( mm != 0 ){
    return ::force_offsetpbc( *this, *mm, im, offset, box, rc, pv );
  }
  else{
    MonatomicMols* mm2 = dynamic_cast<MonatomicMols*> ( m.get() );
    assert ( mm2 != 0 );
    return ::force_offsetpbc( *this, *mm2, im, offset, box, rc, pv );
  }
  return 0;
}




//select force functions dynamically
int RigidBodies::Force_PBC( 
			    const MolsHandle& m,
			    const Intersite& im,
			    const Box& box,
			    const Truncation& rc,
			    PotVir &pv
			    )
{
  RigidBodies* mm = dynamic_cast<RigidBodies*> ( m.get() );
  if ( mm != 0 ){
    return ::force_pbc( *this, *mm, im, box, rc, pv );
  }
  else{
    MonatomicMols* mm2 = dynamic_cast<MonatomicMols*> ( m.get() );
    assert ( mm2 != 0 );
    return ::force_pbc( *this, *mm2, im, box, rc, pv );
  }
  return 0;
}




//select force functions dynamically
int RigidBodies::Force_General( 
			    const MolsHandle& m,
			    const Intersite& im,
			    const TruncPair& truncpair,
			    PotVir &pv
			    )
{
  RigidBodies* mm = dynamic_cast<RigidBodies*> ( m.get() );
  if ( mm != 0 ){
    return ::force_common_2nd( *this, *mm, im, truncpair, pv );
  }
  else{
    MonatomicMols* mm2 = dynamic_cast<MonatomicMols*> ( m.get() );
    assert ( mm2 != 0 );
    return ::force_common_2nd( *this, *mm2, im, truncpair, pv );
  }
  return 0;
}



const Vector3&
RigidBodies::Position( int i ) const
{
  return mols[i].Position();
}



void
RigidBodies::ReadWTG5( int n, const Unit& unit, FILE* file )
{
  MolPropertyHandle pp = GetProperty();
  Rigid* p = dynamic_cast<Rigid*> ( pp.get() );
  assert( p != 0 );
  for( int i=0; i<n; i++ ){
    RigidBody* m = new RigidBody( p->GetNumSite(), pp );
    m->ReadWTG5( i, unit, file );
    Push( SingleMolHandle( m ) );
  }
}



void
RigidBodies::ReadNX4A( int n, const Box& box, const Unit& unit, FILE* file )
{
  MolPropertyHandle pp = GetProperty();
  Rigid* p = dynamic_cast<Rigid*> ( pp.get() );
  assert( p != 0 );
  for( int i=0; i<n; i++ ){
    RigidBody* m = new RigidBody( p->GetNumSite(), pp );
    m->ReadNX4A( i, unit, file );
    if ( &box != 0 )
      box.BoxCoordinate( m->com.center.coord );
    Push( SingleMolHandle( m ) );
  }
}



/*
void
RigidBodies::AddVelocity( const Vector3& velocity )
{
  int nmol  = Size();
  for( int i=0; i<nmol; i++ ){
    mols[i].AddVelocity( velocity );
  }
}
*/






void RigidBodies::PluginHookL0( CollectorPlugin& plugin )
{
  int nmol = mols.size();
  for(int i=0; i<nmol; i++)
    plugin.HookL0( &mols[i] );
}



//functions for quenching
double*
RigidBodies::unserialize( double* const p )
{
  double* ptr = p;
  int    nmol = Size();
  for( int i=0; i < nmol; i++ ){
    ptr = mols[i].unserialize( ptr );
  }
  return ptr;
}



double*
RigidBodies::serialize( double* const p ) const
{
  double* ptr = p;
  int    nmol = Size();
  for( int i=0; i < nmol; i++ ){
    ptr = mols[i].serialize( ptr );
  }
  return ptr;
}



double*
RigidBodies::serializeforce( double* const p ) const
{
  double* ptr = p;
  int    nmol = Size();
  for( int i=0; i < nmol; i++ ){
    ptr = mols[i].serializeforce( ptr );
  }
  return ptr;
}


