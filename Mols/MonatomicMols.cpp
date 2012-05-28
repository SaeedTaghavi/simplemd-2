#include "Mols/MonatomicMols.hpp"
#include "debug.hpp"
#include "Mols/RigidBodies.hpp"
#include "Mols/RigidBodies2.hpp"
#include "Plugin/CollectorPlugin.hpp"

/*Rewritten from MonatomicMols.cpp by regarding it as a set of MonatomicMol.
 *It is, however, slower than the original one.
 *because of bad ordering of molecular coordinate on the memory, maybe.
 *Original one is as fast as nph2b.cc, the reference c++ code.
 *Dispose it, and polyatomic molecules should be written based on MonatomicMols.cpp,
 * not from this code.
 * NOTE: When using Fujitsu's compiler for hp2500, templates process is too slow.
 */

void
MonatomicMols::SetProperty( const MolPropertyHandle& p )
{
  prop = p;
}



MolPropertyHandle
MonatomicMols::GetProperty() const
{
  return prop;
}



int
MonatomicMols::Size() const
{
  return mols.size();
}




MonatomicMols::MonatomicMols()
{
  init();
}



MonatomicMols::MonatomicMols( MolPropertyHandle p )
{
  init();
  SetProperty( p );
}



MonatomicMols::MonatomicMols( MolPropertyHandle p, int nmol, const Unit& unit, FILE* input )
{
  init();
  SetProperty( p );
  ReadATG5( nmol, unit, input );
}
  




MolsHandle
MonatomicMols::EmptyClone() const
{
  MolsHandle h = MolsHandle( new MonatomicMols( prop ) );
  return h;
}



void
MonatomicMols::init()
{
}



MonatomicMols::~MonatomicMols()
{
  mesg("~MonatomicMols\n");
}



int MonatomicMols::push1( const MonatomicMol &mol )
{
  mols.push_back( mol );
  return mols.size();
}



int MonatomicMols::Push( SingleMolHandle mol )
{
  MonatomicMol* m = dynamic_cast<MonatomicMol*> ( mol.get() );
  assert( m != 0 );
  return push1( *m );
}



/*
 *Get the reference to the specified element of MonatomicMols.
 */
MonatomicMol*
MonatomicMols::peek1( int i ) const
{
    return new MonatomicMol( mols[i] );
}



//simple reference to a molecule
const SingleMolEntity&
MonatomicMols::Molecule( int i ) const
{
  return mols[i];
}



SingleMolHandle
MonatomicMols::Peek( int i ) const
{
    return SingleMolHandle( peek1( i ) );
}




/*
 *Pull out the specified element of MonatomicMols.
 *The element is removed.
 *Returns MonatomicMol, which must be deleted later.
 *NOTICE: force is not preserved.
 */
MonatomicMol*
MonatomicMols::pull( int i )
{
    MonatomicMol* mol = peek1(i);
    int nmol = mols.size();
    nmol --;
  //move the last one
  mols[i] = mols[nmol];
  //truncate the last
  mols.pop_back();
  return mol;
}



MolsHandle
MonatomicMols::Emmigrate( const Box& box )
{
  MonatomicMols* emmigrants = new MonatomicMols;
  unsigned int i=0;
  while( i<mols.size() ){
    MonatomicMol& m = mols[i];
    if ( box.Contains( m.Position() ) ){
      i++;
    }
    else{
      // copy the molecule temporarily.
      MonatomicMol* mol = pull( i );
      emmigrants->push1( *mol );
      delete mol;
    }
  }
  return MolsHandle( emmigrants );
}




double MonatomicMols::GetEk() const
{
  int nmol = Size();
  double ek = 0;
  for( int i=0; i<nmol; i++ ){
    ek += mols[i].GetEk();
  }
  return ek;
}



void MonatomicMols::GetEkt( Vector3& ekt ) const
{
  MolPropertyHandle h = GetProperty();
  Flexible*       e = dynamic_cast<Flexible*> ( h.get() );
  assert( e != 0 );
  assert( e->GetNumSite() == 1 );
  const vector<double>& massarray = e->GetMassArray();
  double mass = massarray[0];
  int    nmol = Size();
  ekt.Set(0,0,0);
  for( int i=0; i<nmol; i++ ){
    Vector3 velocity = mols[i].GetVelocity();
    ekt.x += velocity.x * velocity.x;
    ekt.y += velocity.y * velocity.y;
    ekt.z += velocity.z * velocity.z;
  }
  ekt.x *= mass * 0.5;
  ekt.y *= mass * 0.5;
  ekt.z *= mass * 0.5;
}



void
MonatomicMols::TotalMomentum( Vector3& momentum ) const
{
  MolPropertyHandle h = GetProperty();
  Flexible*       e = dynamic_cast<Flexible*> ( h.get() );
  assert( e != 0 );
  assert( e->GetNumSite() == 1 );
  const vector<double>& massarray = e->GetMassArray();
  double mass = massarray[0];
  int    nmol = Size();

  momentum.Set(0,0,0);
  for( int i=0; i<nmol; i++ ){
    Vector3 velocity = mols[i].GetVelocity();
    momentum.x += velocity.x;
    momentum.y += velocity.y;
    momentum.z += velocity.z;
  }
  momentum.x *= mass;
  momentum.y *= mass;
  momentum.z *= mass;
}



/*
void
MonatomicMols::AddVelocity( const Vector3& velocity )
{
  int    nmol = Size();
  for( int i=0; i<nmol; i++ ){
    mols[i].AddVelocity( velocity );
  }
}
*/



/*
void MonatomicMols::ScaleVelocity( const Vector3& r )
{
  int    nmol = mols.size();
  for( int i=0; i < nmol; i++ ){
    mols[i].ScaleVelocity( r );
  }
}
*/


/*
void MonatomicMols::ScalePosition( const Vector3& r )
{
  int    nmol = mols.size();
  for( int i=0; i < nmol; i++ ){
    mols[i].ScalePosition( r );
  }
}
*/


void MonatomicMols::Preforce()
{
  int    nmol = mols.size();
  for( int i=0; i< nmol; i++ ){
    mols[i].Preforce();
  }
}





void MonatomicMols::BoxCoordinate( const Box& box )
{
  int    nmol = mols.size();
  for( int i=0; i<nmol; i++ ){
    MonatomicMol& m = mols[i];
    //Vector3 position = m.Position();
    //fprintf(stderr, "%f %f %f\n", position.x, position.y, position.z );
    Vector3 pos( m.Position() );
    box.BoxCoordinate( pos );
    m.SetPosition( pos );
  }
}


int MonatomicMols::Force_PBC( const Intersite& im, const Box &box, const Truncation& rc, PotVir &pv )
{
  double ep = 0;
  Matrix33 vr;
  MolPropertyHandle h  = GetProperty();
  Flexible*       e = dynamic_cast<Flexible*> ( h.get() );
  if ( e ){
    assert( e->GetNumSite() == 1 );
    const IntrParams& ip = e->GetIntrParams()[0];
    const double sig = ip.sig;
    const double eps = ip.eps;
    const int    nmol = mols.size();
    const double inner = rc.GetInner();
    const double outer = rc.GetOuter();
    const double rq    = rc.GetRQ();
    for( int i=0; i<nmol-1;i++ ){
      MonatomicMol& mi = mols[i];
      for( int j=i+1; j<nmol; j++ ){
	MonatomicMol& mj = mols[j];
        Vector3 d;
        const Vector3& vi = mi.Position();
        const Vector3& vj = mj.Position();
	d.x = vi.x - vj.x;
	d.y = vi.y - vj.y;
	d.z = vi.z - vj.z;
        box.RelativePosition( d );
        const double& dx = d.x;
        const double& dy = d.y;
        const double& dz = d.z;
	
	double drr = (dx*dx+dy*dy+dz*dz);
	if( drr < outer*outer ){
	  double sfe,sff,r;
	  if( drr < inner*inner){
	    sfe = 1.0;
	    sff = 0.0;
	  }
	  else {
	    double rrl,rrs,rrs2,rrl2;
	    r = sqrt(drr);
	    rrl = r - outer;
	    rrs = r - inner;
	    rrl2 = rrl*rrl;
	    rrs2 = rrs*rrs;
	    sfe = rrl2*rrl*rq*(10.0*rrs2-5.0*rrl*rrs+rrl2);
	    sff = 30.0*rq*rrl2*rrs2/r;
	  }
	  double drs = sig * sig / drr;
	  double dr6 = drs * drs * drs;
	  double dr12 = dr6 * dr6;
	  double fo = -4.0 * eps * (-12.0 * dr12 + 6.0 * dr6)/drr;
	  double fullpot = 4.0 * eps * (dr12-dr6);
	  ep += sfe*fullpot;
	  double ffx = dx*(fo*sfe - sff*fullpot);
	  double ffy = dy*(fo*sfe - sff*fullpot);
	  double ffz = dz*(fo*sfe - sff*fullpot);
	  mi.center.force.x += ffx;
	  mi.center.force.y += ffy;
	  mi.center.force.z += ffz;
	  mj.center.force.x -= ffx;
	  mj.center.force.y -= ffy;
	  mj.center.force.z -= ffz;
	  vr.v[0][0] += dx * ffx;
	  vr.v[0][1] += dx * ffy;
	  vr.v[0][2] += dx * ffz;
	  vr.v[1][0] += dy * ffx;
	  vr.v[1][1] += dy * ffy;
	  vr.v[1][2] += dy * ffz;
	  vr.v[2][0] += dz * ffx;
	  vr.v[2][1] += dz * ffy;
	  vr.v[2][2] += dz * ffz;
#ifdef DEBUG
	  //printf("# %d %d %f\n", order[i], order[j], r);
	  ::count++;
#endif
	}
      }
    }
    pv.Set(ep,vr);
#ifdef DEBUG
    fprintf(stderr, "%d\n", ::count);
#endif
    return 1;
  }
  return 0;
}

/*
 *OBSOLETEの方が3倍ぐらい速い
 *結局、listvectorの生成コストがかなり大きいし、単原子分子の場合、相対位置ベクトルだけですべて計算できてしまうので、それを局所変数として確保できると、速度にかなり効く。その意味では、Rigidでも速くなるかどうかはわからない。(Rigidの場合、ListVectorの重要度は相対的に低い。)
 */

//force w/o PBC
int MonatomicMols::Force( const Intersite& im, const Truncation& rc, PotVir &pv )
{
  double ep = 0;
  Matrix33 vr;
  
  MolPropertyHandle h  = GetProperty();
  Flexible*       e  = dynamic_cast<Flexible*> ( h.get() );
  if ( e ){
    if ( e->GetNumSite() != 1 )
      exit(1);
    const IntrParams& ih = e->GetIntrParams()[0];
    double sig = ih.sig;
    double eps = ih.eps;
    int    nmol = mols.size();
    const double inner = rc.GetInner();
    const double outer = rc.GetOuter();
    //const double rq    = rc.GetRQ();
    for( int i=0; i<nmol-1;i++ ){
      MonatomicMol& mi = mols[i];
      for( int j=i+1; j<nmol; j++ ){
	MonatomicMol& mj = mols[j];
        const Vector3& vi = mi.Position();
        const Vector3& vj = mj.Position();
	double dx = vi.x - vj.x;
	double dy = vi.y - vj.y;
	double dz = vi.z - vj.z;
	
	double drr = (dx*dx+dy*dy+dz*dz);
	if( drr < outer*outer ){
	  double sfe,sff;
	  if( drr < inner*inner){
	    sfe = 1.0;
	    sff = 0.0;
	  }
	  else {
	    rc.Smooth( sqrt( drr ), sfe, sff );
	    /*ここをSmoothに置きかえると、1秒/5.3秒も速くなった。
	    double rrl,rrs,rrs2,rrl2,r;
	    r = sqrt(drr);
	    rrl = r - outer;
	    rrs = r - inner;
	    rrl2 = rrl*rrl;
	    rrs2 = rrs*rrs;
	    sfe = rrl2*rrl*rq*(10.0*rrs2-5.0*rrl*rrs+rrl2);
	    sff = 30.0*rq*rrl2*rrs2/r;
	    */
	  }
	  double drs = sig * sig / drr;
	  double dr6 = drs * drs * drs;
	  double dr12 = dr6 * dr6;
	  double fo = -4.0 * eps * (-12.0 * dr12 + 6.0 * dr6)/drr;
	  double fullpot = 4.0 * eps * (dr12-dr6);
	  ep += sfe*fullpot;
	  double ffx = dx*(fo*sfe - sff*fullpot);
	  double ffy = dy*(fo*sfe - sff*fullpot);
	  double ffz = dz*(fo*sfe - sff*fullpot);
	  mi.center.force.x += ffx;
	  mi.center.force.y += ffy;
	  mi.center.force.z += ffz;
	  mj.center.force.x -= ffx;
	  mj.center.force.y -= ffy;
	  mj.center.force.z -= ffz;
	  vr.v[0][0] += dx * ffx;
	  vr.v[0][1] += dx * ffy;
	  vr.v[0][2] += dx * ffz;
	  vr.v[1][0] += dy * ffx;
	  vr.v[1][1] += dy * ffy;
	  vr.v[1][2] += dy * ffz;
	  vr.v[2][0] += dz * ffx;
	  vr.v[2][1] += dz * ffy;
	  vr.v[2][2] += dz * ffz;
#ifdef DEBUG
	  //printf("# %d %d %f\n", order[i], order[j], r);
	  ::count++;
#endif
	}
      }
    }
    pv.Set(ep,vr);
#ifdef DEBUG
    fprintf(stderr, "%d\n", ::count);
#endif
    return 1;
  }
  return 0;
}



void MonatomicMols::Postforce( //const Vessel& vessel
                                )
{
  MolPropertyHandle h = GetProperty();
  Flexible*       e = dynamic_cast<Flexible*> ( h.get() );
  assert( e != 0 );
  assert( e->GetNumSite() == 1 );
  //MolProperty* e = h.get();
  const vector<double>& massarray = e->GetMassArray();
  double mass   = massarray[0];
  //double volume = vessel.volume;
  //double d1     = vessel.pc.d1();

  int    nmol = mols.size();
  for( int i=0; i< nmol; i++ ){
      mols[i].Postforce( //d1, volume, 
                         mass );
  }
}



//select force functions dynamically
int MonatomicMols::Force_Simple( 
				  const MolsHandle& m,
				  const Intersite& im,
				  const Truncation& rc,
				  PotVir &pv
				  )
{
  Vector3 offset(0,0,0);
  return Force_Offset( m,im,offset,rc,pv);
}



int
force_offset_simpler(
	MonatomicMols& m1,
	MonatomicMols& m2,
	const Intersite& im,
	const Vector3 &offset, 
	const Truncation& rc,
	PotVir &pv
	);

//select force functions dynamically
int MonatomicMols::Force_Offset( 
			   const MolsHandle& m,
			   const Intersite& im,
			   const Vector3 &offset,
			   const Truncation& rc,
			   PotVir &pv
			   )
{
  MonatomicMols* mm = dynamic_cast<MonatomicMols*> ( m.get() );
  if ( mm != 0 ){
    return ::force_offset( *this, *mm, offset, rc, pv );
    //return ::force_offset_simpler( *this, *mm, im, offset, rc, pv );
  }
  else{
    RigidBodies* mm2 = dynamic_cast<RigidBodies*> ( m.get() );
    assert ( mm2 != 0 );
    Vector3 revert;
    revert.x = -offset.x;
    revert.y = -offset.y;
    revert.z = -offset.z;
    return ::force_offset( *mm2, *this, im, revert, rc, pv );
 }

  return 0;
}



int
force_offsetpbc(
		MonatomicMols& m1,
		MonatomicMols& m2,
		const Vector3 &offset, 
		const Box& box,
		const Truncation& rc,
		PotVir &pv
		);


//select force functions dynamically
int MonatomicMols::Force_OffsetPBC( 
			   const MolsHandle& m,
			   const Intersite& im,
			   const Vector3 &offset,
			   const Box& box,
			   const Truncation& rc,
			   PotVir &pv
			   )
{
  MonatomicMols* mm = dynamic_cast<MonatomicMols*> ( m.get() );
  if ( mm != 0 ){
    return ::force_offsetpbc( *this, *mm, offset, box, rc, pv );
    //return ::force_offset_simpler( *this, *mm, im, offset, rc, pv );
  }
  else{
    RigidBodies2* mm2 = dynamic_cast<RigidBodies2*> ( m.get() );
    assert ( mm2 != 0 );
    Vector3 revert;
    revert.x = -offset.x;
    revert.y = -offset.y;
    revert.z = -offset.z;
    return ::force_offsetpbc( *mm2, *this, im, revert, box, rc, pv );
 }
  return 0;
}



//select force functions dynamically
int MonatomicMols::Force_PBC( 
                            const MolsHandle& m,
			    const Intersite& im, 
                            const Box& box,
                            const Truncation& rc,
                            PotVir &pv
                            )
{
  MonatomicMols* mm = dynamic_cast<MonatomicMols*> ( m.get() );
  if ( mm != 0 ){
    return ::force_pbc( *this, *mm, box, rc, pv );
  }
  else{
    RigidBodies* mm2 = dynamic_cast<RigidBodies*> ( m.get() );
    assert ( mm2 != 0 );
    return ::force_pbc( *mm2, *this, im, box, rc, pv );
  }
  return 0;
}



int MonatomicMols::Force_General( 
                            const MolsHandle& m,
			    const Intersite& im, 
                            const TruncPair& truncpair,
                            PotVir &pv
                            )
{
  assert( 0 );
  MonatomicMols* mm = dynamic_cast<MonatomicMols*> ( m.get() );
  if ( mm != 0 ){
    //return ::force_common_2nd( *this, *mm, truncpair, pv );
    //MonatomicMol同士の相互作用はlistvectorを使っていないので、そもそもbookkeepingにするのが難しい。保留
  }
  else{
    RigidBodies* mm2 = dynamic_cast<RigidBodies*> ( m.get() );
    assert ( mm2 != 0 );
    //MonatomicMolとRigidBodyの場合、主客反転しなければいけないが、その際にtruncpairの相対vectorも反転する必要がある。これも保留。
    //return ::force_common_2nd( *mm2, *this, im, truncpair, pv );
  }
  return 0;
}




void 
MonatomicMols::Translate( const Vector3& offset )
{
  int    nmol = mols.size();
  for( int j=0; j<nmol; j++ ){
    mols[j].Translate( offset );
  }
}




void
MonatomicMols::ProgressPosition( double dt )
{
  int    nmol = mols.size();
  for( int j=0; j<nmol; j++ ){
    mols[j].ProgressPosition( dt );
  }
}



void
MonatomicMols::ProgressMomentum( double dt )
{
  int    nmol = mols.size();
  for( int j=0; j<nmol; j++ ){
    mols[j].ProgressMomentum( dt );
  }
}



int
force_offset(
	MonatomicMols& m1,
	MonatomicMols& m2,
	const Vector3 &offset, 
	const Truncation& rc,
	PotVir &pv
	)
{
  double ep = 0;
  Matrix33 vr;
  
  int nmol  = m1.Size();
  int nmol2 = m2.Size();
  MolPropertyHandle h1  = m1.GetProperty();
  Flexible*       e1  = dynamic_cast<Flexible*> ( h1.get() );
  MolPropertyHandle h2  = m2.GetProperty();
  Flexible*       e2  = dynamic_cast<Flexible*> ( h2.get() );
  if ( e1 && e2 ){
    assert( e1->GetNumSite() == 1 );
    assert( e2->GetNumSite() == 1 );
    const IntrParams& i1 = e1->GetIntrParams()[0];
    const IntrParams& i2 = e2->GetIntrParams()[0];
    double sig = i1.sighalf + i2.sighalf;
    double eps = i1.epssqrt * i2.epssqrt;
    const double inner = rc.GetInner();
    const double outer = rc.GetOuter();
    //const double rq    = rc.GetRQ();
    for( int i=0; i<nmol; i++ ){
      MonatomicMol& mi = m1.mols[i];
      const Vector3& coordi = mi.Position();
      //Vector3& coordi = mi.center.coord;
      double dxi = coordi.x;
      double dyi = coordi.y;
      double dzi = coordi.z;
      Vector3& mif = m1.mols[i].center.force;
      for( int j=0; j<nmol2; j++ ){
        MonatomicMol& mj = m2.mols[j];
        //Vector3& coordj = mj.center.coord;
        const Vector3& coordj = mj.Position();
        Vector3& mjf = m2.mols[j].center.force;
        double dx = dxi - coordj.x + offset.x;
        double dy = dyi - coordj.y + offset.y;
        double dz = dzi - coordj.z + offset.z;
        
        double drr = dx*dx + dy*dy + dz*dz;
        if( drr < outer*outer ){
          double sfe,sff;
          if( drr < inner*inner){
            sfe = 1.0;
            sff = 0.0;
          }
          else {
	    rc.Smooth( sqrt( drr ), sfe, sff );
          }
          double drs = sig * sig / drr;
          double dr6 = drs * drs * drs;
          double dr12 = dr6 * dr6;
          double f = -4.0 * eps * (-12.0 * dr12 + 6.0 * dr6)/drr;
          double fullpot = 4.0 * eps * (dr12-dr6);
          ep += sfe*fullpot;
          double ffx = dx*(f*sfe - sff*fullpot);
          double ffy = dy*(f*sfe - sff*fullpot);
          double ffz = dz*(f*sfe - sff*fullpot);
          mif.x += ffx;
          mif.y += ffy;
          mif.z += ffz;
          mjf.x -= ffx;
          mjf.y -= ffy;
          mjf.z -= ffz;
          vr.v[0][0] += dx * ffx;
          vr.v[0][1] += dx * ffy;
          vr.v[0][2] += dx * ffz;
          vr.v[1][0] += dy * ffx;
          vr.v[1][1] += dy * ffy;
          vr.v[1][2] += dy * ffz;
          vr.v[2][0] += dz * ffx;
          vr.v[2][1] += dz * ffy;
          vr.v[2][2] += dz * ffz;
#ifdef DEBUG
          //printf("# %d %d %f %f %f %f\n", m1.order[i], m2.order[j], drr, dx, dy, dz );
          ::count++;
#endif
        }
      }
    }
    pv.Set( ep, vr );
#ifdef DEBUG
    fprintf(stderr, "%d\n", ::count);
#endif
    return 1;
  }
  return 0;
}



int
force_offsetpbc(
	MonatomicMols& m1,
	MonatomicMols& m2,
	const Vector3 &offset, 
	const Box& box,
	const Truncation& rc,
	PotVir &pv
	)
{
  double ep = 0;
  Matrix33 vr;
  
  int nmol  = m1.Size();
  int nmol2 = m2.Size();
  MolPropertyHandle h1  = m1.GetProperty();
  Flexible*       e1  = dynamic_cast<Flexible*> ( h1.get() );
  MolPropertyHandle h2  = m2.GetProperty();
  Flexible*       e2  = dynamic_cast<Flexible*> ( h2.get() );
  if ( e1 && e2 ){
    assert( e1->GetNumSite() == 1 );
    assert( e2->GetNumSite() == 1 );
    const IntrParams& i1 = e1->GetIntrParams()[0];
    const IntrParams& i2 = e2->GetIntrParams()[0];
    double sig = i1.sighalf + i2.sighalf;
    double eps = i1.epssqrt * i2.epssqrt;
    const double inner = rc.GetInner();
    const double outer = rc.GetOuter();
    //const double rq    = rc.GetRQ();
    for( int i=0; i<nmol; i++ ){
      MonatomicMol& mi = m1.mols[i];
      const Vector3& coordi = mi.Position();
      double dxi = coordi.x;
      double dyi = coordi.y;
      double dzi = coordi.z;
      Vector3& mif = m1.mols[i].center.force;
      for( int j=0; j<nmol2; j++ ){
        MonatomicMol& mj = m2.mols[j];
        const Vector3& coordj = mj.Position();
        //Vector3& coordj = mj.center.coord;
        Vector3& mjf = m2.mols[j].center.force;
        Vector3 d;
        d.x = dxi - coordj.x + offset.x;
        d.y = dyi - coordj.y + offset.y;
        d.z = dzi - coordj.z + offset.z;
        box.RelativePosition( d );
        const double& dx = d.x;
        const double& dy = d.y;
        const double& dz = d.z;

        double drr = dx*dx + dy*dy + dz*dz;
        if( drr < outer*outer ){
          double sfe,sff;
          if( drr < inner*inner){
            sfe = 1.0;
            sff = 0.0;
          }
          else {
	    rc.Smooth( sqrt( drr ), sfe, sff );
          }
          double drs = sig * sig / drr;
          double dr6 = drs * drs * drs;
          double dr12 = dr6 * dr6;
          double f = -4.0 * eps * (-12.0 * dr12 + 6.0 * dr6)/drr;
          double fullpot = 4.0 * eps * (dr12-dr6);
          ep += sfe*fullpot;
          double ffx = dx*(f*sfe - sff*fullpot);
          double ffy = dy*(f*sfe - sff*fullpot);
          double ffz = dz*(f*sfe - sff*fullpot);
          mif.x += ffx;
          mif.y += ffy;
          mif.z += ffz;
          mjf.x -= ffx;
          mjf.y -= ffy;
          mjf.z -= ffz;
          vr.v[0][0] += dx * ffx;
          vr.v[0][1] += dx * ffy;
          vr.v[0][2] += dx * ffz;
          vr.v[1][0] += dy * ffx;
          vr.v[1][1] += dy * ffy;
          vr.v[1][2] += dy * ffz;
          vr.v[2][0] += dz * ffx;
          vr.v[2][1] += dz * ffy;
          vr.v[2][2] += dz * ffz;
#ifdef DEBUG
          //printf("# %d %d %f %f %f %f\n", m1.order[i], m2.order[j], drr, dx, dy, dz );
          ::count++;
#endif
        }
      }
    }
    pv.Set(ep,vr);
#ifdef DEBUG
    fprintf(stderr, "%d\n", ::count);
#endif
    return 1;
  }
  return 0;
}



int
force_offset_simpler(
	MonatomicMols& m1,
	MonatomicMols& m2,
	const Intersite& im,
	const Vector3 &offset, 
	const Truncation& rc,
	PotVir &pv
	)
{
  double ep = 0;
  Matrix33 vr;
  
  int nmol  = m1.Size();
  int nmol2 = m2.Size();
  MolPropertyHandle h1  = m1.GetProperty();
  Flexible*       e1  = dynamic_cast<Flexible*> ( h1.get() );
  MolPropertyHandle h2  = m2.GetProperty();
  Flexible*       e2  = dynamic_cast<Flexible*> ( h2.get() );
  if ( e1 && e2 ){
    assert( e1->GetNumSite() == 1 );
    assert( e2->GetNumSite() == 1 );
    //const IntrParams& i1 = e1->GetIntrParams()[0];
    //const IntrParams& i2 = e2->GetIntrParams()[0];
    //double sig = i1.sighalf + i2.sighalf;
    //double eps = i1.epssqrt * i2.epssqrt;
    const double inner = rc.GetInner();
    const double outer = rc.GetOuter();
    //const double rq    = rc.GetRQ();
    for( int i=0; i<nmol; i++ ){
      MonatomicMol& mi = m1.mols[i];
      const Vector3& coordi = mi.Position();
      double dxi = coordi.x;
      double dyi = coordi.y;
      double dzi = coordi.z;
      Vector3& mif = m1.mols[i].center.force;
      for( int j=0; j<nmol2; j++ ){
        MonatomicMol& mj = m2.mols[j];
        const Vector3& coordj = mj.Position();
        Vector3& mjf = m2.mols[j].center.force;
	Vector3 delta(dxi - coordj.x + offset.x,
		      dyi - coordj.y + offset.y,
		      dzi - coordj.z + offset.z );
        
        double drr = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
        if( drr < outer*outer ){
          double sfe,sff;
          if( drr < inner*inner){
            sfe = 1.0;
            sff = 0.0;
          }
          else {
	    rc.Smooth( sqrt( drr ), sfe, sff );
	    /*ところがここをSmoothにおきかえると0.5sec/5.3sec遅くなる。(gcc -O2)何故?*/
	    //iccだと変化しない。
            /*double rrl,rrs,rrs2,rrl2,r;
            r = sqrt(drr);
            rrl = r - outer;
            rrs = r - inner;
            rrl2 = rrl*rrl;
            rrs2 = rrs*rrs;
            sfe = rrl2*rrl*rq*(10.0*rrs2-5.0*rrl*rrs+rrl2);
            sff = 30.0*rc.rq*rrl2*rrs2/r;*/
          }
	  double fo, potential;
	  //Forceを使うと、+2sec/3.3sec遅くなる。
	  //deltaをオブジェクトにしたせいでレジスター化できないせいか。
	  //でも、Rigid-Monatomic相互作用ではこの書きかたをせざるをえない。
	  im.intr[0]->Force( delta, fo, potential );
          ep += sfe*potential;
          double ffx = delta.x*(fo*sfe - sff*potential);
          double ffy = delta.y*(fo*sfe - sff*potential);
          double ffz = delta.z*(fo*sfe - sff*potential);
          mif.x += ffx;
          mif.y += ffy;
          mif.z += ffz;
          mjf.x -= ffx;
          mjf.y -= ffy;
          mjf.z -= ffz;
          vr.v[0][0] += delta.x * ffx;
          vr.v[0][1] += delta.x * ffy;
          vr.v[0][2] += delta.x * ffz;
          vr.v[1][0] += delta.y * ffx;
          vr.v[1][1] += delta.y * ffy;
          vr.v[1][2] += delta.y * ffz;
          vr.v[2][0] += delta.z * ffx;
          vr.v[2][1] += delta.z * ffy;
          vr.v[2][2] += delta.z * ffz;
#ifdef DEBUG
          //printf("# %d %d %f %f %f %f\n", m1.order[i], m2.order[j], drr, dx, dy, dz );
          ::count++;
#endif
        }
      }
    }
    pv.Set(ep,vr);
#ifdef DEBUG
    fprintf(stderr, "%d\n", ::count);
#endif
    return 1;
  }
  return 0;
}



int
force_pbc(
      MonatomicMols& m1,
      MonatomicMols& m2,
      const Box& box, 
      const Truncation& rc,
      PotVir &pv
      )
{
  double ep = 0;
  Matrix33 vr;
  
  int nmol  = m1.Size();
  int nmol2 = m2.Size();
  MolPropertyHandle h1  = m1.GetProperty();
  Flexible*       e1  = dynamic_cast<Flexible*> ( h1.get() );
  MolPropertyHandle h2  = m2.GetProperty();
  Flexible*       e2  = dynamic_cast<Flexible*> ( h2.get() );
  if ( e1 && e2 ){
    assert( e1->GetNumSite() == 1 );
    assert( e2->GetNumSite() == 1 );
    const IntrParams& i1 = e1->GetIntrParams()[0];
    const IntrParams& i2 = e2->GetIntrParams()[0];
    double sig = i1.sighalf + i2.sighalf;
    double eps = i1.epssqrt * i2.epssqrt;
    const double inner = rc.GetInner();
    const double outer = rc.GetOuter();
    const double rq    = rc.GetRQ();
    for( int i=0; i<nmol; i++ ){
      MonatomicMol& mi = m1.mols[i];
      const Vector3& coordi = mi.Position();
      double dxi = coordi.x;
      double dyi = coordi.y;
      double dzi = coordi.z;
      Vector3& mif = m1.mols[i].center.force;
      for( int j=0; j<nmol2; j++ ){
        MonatomicMol& mj = m2.mols[j];
        const Vector3& coordj = mj.Position();
        Vector3& mjf = m2.mols[j].center.force;
        Vector3 d;
        d.x = dxi - coordj.x;
        d.y = dyi - coordj.y;
        d.z = dzi - coordj.z;
        box.RelativePosition( d );
        const double& dx = d.x;
        const double& dy = d.y;
        const double& dz = d.z;
        
        double drr = dx*dx + dy*dy + dz*dz;
        if( drr < outer*outer ){
          double sfe,sff,r;
          if( drr < inner*inner){
            sfe = 1.0;
            sff = 0.0;
          }
          else {
            double rrl,rrs,rrs2,rrl2;
            r = sqrt(drr);
            rrl = r - outer;
            rrs = r - inner;
            rrl2 = rrl*rrl;
            rrs2 = rrs*rrs;
            sfe = rrl2*rrl*rq*(10.0*rrs2-5.0*rrl*rrs+rrl2);
            sff = 30.0*rq*rrl2*rrs2/r;
          }
          double drs = sig * sig / drr;
          double dr6 = drs * drs * drs;
          double dr12 = dr6 * dr6;
          double f = -4.0 * eps * (-12.0 * dr12 + 6.0 * dr6)/drr;
          double fullpot = 4.0 * eps * (dr12-dr6);
          ep += sfe*fullpot;
          double ffx = dx*(f*sfe - sff*fullpot);
          double ffy = dy*(f*sfe - sff*fullpot);
          double ffz = dz*(f*sfe - sff*fullpot);
          mif.x += ffx;
          mif.y += ffy;
          mif.z += ffz;
          mjf.x -= ffx;
          mjf.y -= ffy;
          mjf.z -= ffz;
          vr.v[0][0] += dx * ffx;
          vr.v[0][1] += dx * ffy;
          vr.v[0][2] += dx * ffz;
          vr.v[1][0] += dy * ffx;
          vr.v[1][1] += dy * ffy;
          vr.v[1][2] += dy * ffz;
          vr.v[2][0] += dz * ffx;
          vr.v[2][1] += dz * ffy;
          vr.v[2][2] += dz * ffz;
#ifdef DEBUG
          //printf("# %d %d %f %f %f %f\n", m1.order[i], m2.order[j], drr, dx, dy, dz );
          ::count++;
#endif
        }
      }
    }
    pv.Set(ep,vr);
#ifdef DEBUG
    fprintf(stderr, "%d\n", ::count);
#endif
    return 1;
  }
return 0;
}




void
MonatomicMols::Write( const Unit& unit, ostream& to )
{
  WriteATG5( unit, to );
}



void
MonatomicMols::WriteATG5( const Unit& unit, ostream& to )
{
  int nmol = mols.size();
  to << "@ATG5\n"
     << nmol << '\n';
  sort( mols.begin(), mols.end() );
  for(int i=0; i<nmol; i++){
    mols[i].WriteATG5( unit, to );
  }
}



void
MonatomicMols::SnapShot( const Unit& unit, ostream& to )
{
  WriteAR3A( unit, to );
}



void
MonatomicMols::WriteAR3A( const Unit& unit, ostream& to )
{
  int nmol = mols.size();
  to << "@AR3A\n"
     << nmol << '\n';
  sort( mols.begin(), mols.end() );
  for(int i=0; i<nmol; i++){
    mols[i].WriteAR3A( unit, to );
  }
}






void
MonatomicMols::ReadATG5( int n, const Unit& unit, FILE* file )
{
  for( int i=0; i<n; i++ ){
    MonatomicMol* m = new MonatomicMol( prop->GetMass() );
    m->ReadATG5( i, unit, file );
    Push( SingleMolHandle( m ) );
  }
}



void
MonatomicMols::ReadAR3A( int n, const Box& box, const Unit& unit, FILE* file )
{
  for( int i=0; i<n; i++ ){
    MonatomicMol* m = new MonatomicMol( prop->GetMass() );
    m->ReadAR3A( i, unit, file );
    if ( &box != 0 )
      BoxCoordinate( box );
    Push( SingleMolHandle( m ) );
  }
}



const Vector3&
MonatomicMols::Position( int i ) const
{
  return mols[i].Position();
}


void MonatomicMols::PluginHookL0( CollectorPlugin& plugin )
{
  int nmol = mols.size();
  for(int i=0; i<nmol; i++)
    plugin.HookL0( &mols[i] );
}



//functions for quenching
double*
MonatomicMols::unserialize( double* const p )
{
  double* ptr = p;
  int    nmol = Size();
  for( int i=0; i < nmol; i++ ){
    ptr = mols[i].unserialize( ptr );
  }
  return ptr;
}



double*
MonatomicMols::serialize( double* const p ) const
{
  double* ptr = p;
  int    nmol = Size();
  for( int i=0; i < nmol; i++ ){
    ptr = mols[i].serialize( ptr );
  }
  return ptr;
}



double*
MonatomicMols::serializeforce( double* const p ) const
{
  double* ptr = p;
  int    nmol = Size();
  for( int i=0; i < nmol; i++ ){
    ptr = mols[i].serializeforce( ptr );
  }
  return ptr;
}


