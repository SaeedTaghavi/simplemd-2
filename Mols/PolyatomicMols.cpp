#include <cassert>
#include "Mols/PolyatomicMols.hpp"
#include "debug.hpp"


/*
PolyatomicMols::PolyatomicMols()
{
  dx = 0;
}
*/


PolyatomicMols::PolyatomicMols( MolPropertyHandle p )
{
  SetProperty( p );
}







void
PolyatomicMols::SetProperty( const MolPropertyHandle& p )
{
  prop = p;
}



MolPropertyHandle
PolyatomicMols::GetProperty() const
{
  return prop;
}



int
PolyatomicMols::Size() const
{
  return mols.size();
}




MolsHandle
PolyatomicMols::EmptyClone() const
{
  MolsHandle h = MolsHandle( new PolyatomicMols( prop ) );
  return h;
}



PolyatomicMols::~PolyatomicMols()
{
  mesg("~PolyatomicMols\n");
  //Arrays are freed in ~MonatomicMols3()

}



int PolyatomicMols::push1( const PolyatomicMol &mol )
{
  mols.push_back( mol );
  return mols.size();
}



int PolyatomicMols::Push( SingleMolHandle mol )
{
  PolyatomicMol* m = dynamic_cast<PolyatomicMol*> ( mol.get() );
  assert( m != 0 );
  return push1( *m );
}



/*
 *Peek the specified element of PolyatomicMols.
 *Returns a ptr to PolyatomicMol, which must be deleted later.
 */
PolyatomicMol*
PolyatomicMols::peek1( int which ) const
{
  PolyatomicMol* mol = new PolyatomicMol( mols[which] );
  return mol;
}



//simple reference to a molecule
const SingleMolEntity&
PolyatomicMols::Molecule( int i ) const
{
  return mols[i];
}



SingleMolHandle
PolyatomicMols::Peek( int which ) const
{
  return SingleMolHandle(peek1(which));
}



/*
 *Pull out the specified element of PolyatomicMols.
 *The element is removed.
 *Returns PolyatomicMol, which must be deleted later.
 */
PolyatomicMol*
PolyatomicMols::pull( int which )
{
  PolyatomicMol* mol = peek1( which );
  if ( mol == 0 ) exit(1);
  mols.erase( mols.begin() + which );
  return mol;
}



MolsHandle
PolyatomicMols::Emmigrate( const Box &box )
{
  PolyatomicMols* emmigrants = new PolyatomicMols( prop );
  mesg("PolyatomicMols3::emmigrate\n");
  unsigned int i=0;
  while( i<mols.size() ){
    Vector3 c = mols[i].Position();
    if ( box.Contains( c ) ){
      i++;
    }
    else{
      // copy the molecule temporarily.
      PolyatomicMol* mol = pull( i );
      emmigrants->push1( *mol );
      delete mol;
    }
  }
  return MolsHandle( emmigrants );
}



/*
double PolyatomicMols::GetEk() const
{
  double ek = 0;
  int nmol  = mols.size();
  for( int i=0; i<nmol; i++ ){
    ek += mols[i].GetEk();
  }
  return ek;
}
*/



void
PolyatomicMols::TotalMomentum( Vector3& momentum ) const
{
  assert(0);
}



/*
void
PolyatomicMols::AddVelocity( const Vector3& velocity )
{
  int nmol = Size();
  for( int i=0; i<nmol; i++ ){
    mols[i].AddVelocity( velocity );
  }
}
*/



/*
void PolyatomicMols::ScaleVelocity( const Vector3& r )
{
  int nmol = Size();
  for( int i=0; i < nmol; i++ ){
    mols[i].ScaleVelocity( r );
  }
}
*/



/*
void PolyatomicMols::ScalePosition( const Vector3& r )
{
  int nmol = mols.size();
  for( int i=0; i < nmol; i++ ){
    mols[i].ScalePosition( r );
  }
}
*/


void PolyatomicMols::Preforce()
{
  int nmol = mols.size();
  for( int i=0; i<nmol; i++ ){
    mols[i].Preforce();
  }
}





void PolyatomicMols::BoxCoordinate( const Box & box )
{
  int nmol = mols.size();
  for( int i=0; i<nmol; i++ ){
    mols[i].BoxCoordinate( box );
  }
}



int PolyatomicMols::Force_PBC( const Intersite& im, const Box &box, const Truncation& rc, PotVir &pv )
{
  double ep = 0;
  Matrix33 vr;
  double bx = box.x / 2;
  double by = box.y / 2;
  double bz = box.z / 2;
  MolPropertyHandle h  = GetProperty();
  Flexible*       p = dynamic_cast<Flexible*> ( h.get() );
  if ( p ){
    int nsite = p->GetNumSite();
    int nmol  = mols.size();

    const IntrParamsArray& ia = p->GetIntrParams();
    //printf("NMOL: %d\n", nmol);
    for( int si=0; si<nsite; si++ ){
      for( int sj=0; sj<nsite; sj++ ){
        const IntrParams& lji = ia[si];
        const IntrParams& ljj = ia[sj];        
	{
	  double sig = lji.sighalf + ljj.sighalf;
	  double eps = lji.epssqrt * ljj.epssqrt;
	  for( int i=0; i<nmol-1;i++ ){
	    PolyatomicMol& mi = mols[i];
	    for( int j=i+1; j<nmol; j++ ){
	      PolyatomicMol& mj = mols[j];
	      double dx = mi.atom[si].center.coord.x - mj.atom[sj].center.coord.x;
	      double dy = mi.atom[si].center.coord.y - mj.atom[sj].center.coord.y;
	      double dz = mi.atom[si].center.coord.z - mj.atom[sj].center.coord.z;
	      if ( dx < -bx ) 
		dx += box.x;
	      else if ( bx < dx )
		dx -= box.x;
	      if ( dy < -by ) 
		dy += box.y;
	      else if ( by < dy )
		dy -= box.y;
	      if ( dz < -bz ) 
		dz += box.z;
	      else if ( bz < dz )
		dz -= box.z;
	      
	      double drr = (dx*dx+dy*dy+dz*dz);
	      if( drr < rc.outer*rc.outer ){
		double sfe,sff,r;
		if( drr < rc.inner*rc.inner){
		  sfe = 1.0;
		  sff = 0.0;
		}
		else {
		  double rrl,rrs,rrs2,rrl2;
		  r = sqrt(drr);
		  rrl = r - rc.outer;
		  rrs = r - rc.inner;
		  rrl2 = rrl*rrl;
		  rrs2 = rrs*rrs;
		  sfe = rrl2*rrl*rc.rq*(10.0*rrs2-5.0*rrl*rrs+rrl2);
		  sff = 30.0*rc.rq*rrl2*rrs2/r;
		}
		double drs = sig * sig / drr;
		double dr6 = drs * drs * drs;
		double dr12 = dr6 * dr6;
		double fo = 4.0 * eps * (-12.0 * dr12 + 6.0 * dr6)/drr;
		double fullpot = 4.0 * eps * (dr12-dr6);
		ep += sfe*fullpot;
		double ffx = dx*(fo*sfe + sff*fullpot);
		double ffy = dy*(fo*sfe + sff*fullpot);
		double ffz = dz*(fo*sfe + sff*fullpot);
		mi.atom[si].center.force.x -= ffx;
		mi.atom[si].center.force.y -= ffy;
		mi.atom[si].center.force.z -= ffz;
		mj.atom[sj].center.force.x += ffx;
		mj.atom[sj].center.force.y += ffy;
		mj.atom[sj].center.force.z += ffz;
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
	}
      }
    }
    pv.ep = ep;
    pv.vr = vr;
#ifdef DEBUG
    fprintf(stderr, "%d\n", ::count);
#endif
    return 1;
  }
  return 0;
}



//force w/o PBC
int PolyatomicMols::Force( const Intersite& im, const Truncation& rc, PotVir &pv )
{
  double ep = 0;
  Matrix33 vr;
  MolPropertyHandle h  = GetProperty();
  Flexible*       p = dynamic_cast<Flexible*> ( h.get() );
  if ( p ){
    int nsite = p->GetNumSite();
    const IntrParamsArray& ia = p->GetIntrParams();
    //printf("NMOL: %d\n", nmol);
    for( int si=0; si<nsite; si++ ){
      for( int sj=0; sj<nsite; sj++ ){
        const IntrParams& lji = ia[si];
        const IntrParams& ljj = ia[sj];        
        {
          double sig = lji.sighalf + ljj.sighalf;
          double eps = lji.epssqrt * ljj.epssqrt;
          int nmol = mols.size();
	  for( int i=0; i<nmol-1;i++ ){
	    PolyatomicMol& mi = mols[i];
	    for( int j=i+1; j<nmol; j++ ){
	      PolyatomicMol& mj = mols[j];
	      double dx = mi.atom[si].center.coord.x - mj.atom[sj].center.coord.x;
	      double dy = mi.atom[si].center.coord.y - mj.atom[sj].center.coord.y;
	      double dz = mi.atom[si].center.coord.z - mj.atom[sj].center.coord.z;
	      double drr = (dx*dx+dy*dy+dz*dz);
	      if( drr < rc.outer*rc.outer ){
		double sfe,sff,r;
		if( drr < rc.inner*rc.inner){
		  sfe = 1.0;
		  sff = 0.0;
		}
		else {
		  double rrl,rrs,rrs2,rrl2;
		  r = sqrt(drr);
		  rrl = r - rc.outer;
		  rrs = r - rc.inner;
		  rrl2 = rrl*rrl;
		  rrs2 = rrs*rrs;
		  sfe = rrl2*rrl*rc.rq*(10.0*rrs2-5.0*rrl*rrs+rrl2);
		  sff = 30.0*rc.rq*rrl2*rrs2/r;
		}
		double drs = sig * sig / drr;
		double dr6 = drs * drs * drs;
		double dr12 = dr6 * dr6;
		double fo = 4.0 * eps * (-12.0 * dr12 + 6.0 * dr6)/drr;
		double fullpot = 4.0 * eps * (dr12-dr6);
		ep += sfe*fullpot;
		double ffx = dx*(fo*sfe + sff*fullpot);
		double ffy = dy*(fo*sfe + sff*fullpot);
		double ffz = dz*(fo*sfe + sff*fullpot);
		mi.atom[si].center.force.x -= ffx;
		mi.atom[si].center.force.y -= ffy;
		mi.atom[si].center.force.z -= ffz;
		mj.atom[sj].center.force.x += ffx;
		mj.atom[sj].center.force.y += ffy;
		mj.atom[sj].center.force.z += ffz;
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
	}
      }
    }
    pv.ep = ep;
    pv.vr = vr;
#ifdef DEBUG
    fprintf(stderr, "%d\n", ::count);
#endif
    return 1;
  }
  return 0;
}



void PolyatomicMols::Postforce( //const Vessel& vessel 
                                )
{
  MolPropertyHandle h = GetProperty();
  Flexible*       e = dynamic_cast<Flexible*> ( h.get() );
  if ( e == 0 ) exit(1);
  const vector<double>& massarray = e->GetMassArray();
  //double volume = vessel.volume;
  //double d1     = vessel.pc.d1();

  //change unit, subtract force from volume changes, etc..
  int nmol = mols.size();
  for( int i=0; i<nmol; i++ ) {
      mols[i].Postforce( //d1, volume, 
                         massarray );
  }
}





//select force functions dynamically
int PolyatomicMols::Force_Simple( 
			   const MolsHandle& m,
			   const Intersite& im,
			   const Truncation& rc,
			   PotVir &pv
			   )
{
  Vector3 offset(0,0,0);
  return Force_Offset( m,im,offset,rc,pv);
}



//select force functions dynamically
int PolyatomicMols::Force_Offset( 
			   const MolsHandle& m,
			   const Intersite& im,
			   const Vector3 &offset,
			   const Truncation& rc,
			   PotVir &pv
			   )
{
  PolyatomicMols* mm = dynamic_cast<PolyatomicMols*> ( m.get() );
  if ( mm != 0 ){
    return ::force_offset( *this, *mm, offset, rc, pv );
  }
  return 0;
}



int
force_offsetpbc(
		PolyatomicMols& m1,
		PolyatomicMols& m2,
		const Vector3 &offset, 
		const Box& box,
		int* isPeriodic,
		const Truncation& rc,
		PotVir &pv
		);




//select force functions dynamically
int PolyatomicMols::Force_OffsetPBC( 
			   const MolsHandle& m,
			   const Intersite& im,
			   const Vector3 &offset,
			   const Box& box,
			   int* isPeriodic,
			   const Truncation& rc,
			   PotVir &pv
			   )
{
  PolyatomicMols* mm = dynamic_cast<PolyatomicMols*> ( m.get() );
  if ( mm != 0 ){
    return ::force_offsetpbc( *this, *mm, offset, box, isPeriodic, rc, pv );
  }
  return 0;
}



//select force functions dynamically
int PolyatomicMols::Force_PBC( 
                            const MolsHandle& m,
			    const Intersite& im,
                            const Box& box,
                            const Truncation& rc,
                            PotVir &pv
                            )
{
  PolyatomicMols* mm = dynamic_cast<PolyatomicMols*> ( m.get() );
  if ( mm != 0 ){
    return ::force_pbc( *this, *mm, box, rc, pv );
  }
  return 0;
}



void 
PolyatomicMols::Translate( const Vector3& offset )
{
  int nmol = mols.size();
  for( int i=0; i<nmol; i++ ){
    mols[i].Translate( offset );
  }  
}




void PolyatomicMols::ProgressMomentum( double dt )
{
  MolPropertyHandle h = GetProperty();
  Flexible* e = dynamic_cast<Flexible*> ( h.get() );
  assert( e != 0 );
  vector<double> massi = e->GetMassI();
  
  int    nmol = mols.size();
  for( int j=0; j<nmol; j++ ){
    mols[j].ProgressMomentum2( dt, massi );
  }
}



int
force_offset(
	PolyatomicMols& m1,
	PolyatomicMols& m2,
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
  Flexible*       p1  = dynamic_cast<Flexible*> ( h1.get() );
  MolPropertyHandle h2  = m2.GetProperty();
  Flexible*       p2  = dynamic_cast<Flexible*> ( h2.get() );
  if ( p1 && p2 ){
    int nsite1 = p1->GetNumSite();
    int nsite2 = p2->GetNumSite();
    const IntrParamsArray& ia1 = p1->GetIntrParams();
    const IntrParamsArray& ia2 = p2->GetIntrParams();
    for( int s1=0; s1<nsite1; s1++ ){
      for( int s2=0; s2<nsite2; s2++ ){
        const IntrParams& lji = ia1[s1];
        const IntrParams& ljj = ia2[s2];        
        {
          double sig = lji.sighalf + ljj.sighalf;
          double eps = lji.epssqrt * ljj.epssqrt;
          for( int i=0; i<nmol; i++ ){
	    PolyatomicMol& mi = m1.mols[i];
	    for( int j=0; j<nmol2; j++ ){
	      PolyatomicMol& mj = m2.mols[j];
	      double dx = mi.atom[s1].center.coord.x - mj.atom[s2].center.coord.x + offset.x;
	      double dy = mi.atom[s1].center.coord.y - mj.atom[s2].center.coord.y + offset.y;
	      double dz = mi.atom[s1].center.coord.z - mj.atom[s2].center.coord.z + offset.z;

	      double drr = dx*dx + dy*dy + dz*dz;
	      if( drr < rc.outer*rc.outer ){
		double sfe,sff,r;
		if( drr < rc.inner*rc.inner){
		  sfe = 1.0;
		  sff = 0.0;
		}
		else {
		  double rrl,rrs,rrs2,rrl2;
		  r = sqrt(drr);
		  rrl = r - rc.outer;
		  rrs = r - rc.inner;
		  rrl2 = rrl*rrl;
		  rrs2 = rrs*rrs;
		  sfe = rrl2*rrl*rc.rq*(10.0*rrs2-5.0*rrl*rrs+rrl2);
		  sff = 30.0*rc.rq*rrl2*rrs2/r;
		}
		double drs = sig * sig / drr;
		double dr6 = drs * drs * drs;
		double dr12 = dr6 * dr6;
		double fo = 4.0 * eps * (-12.0 * dr12 + 6.0 * dr6)/drr;
		double fullpot = 4.0 * eps * (dr12-dr6);
		ep += sfe*fullpot;
		double ffx = dx*(fo*sfe + sff*fullpot);
		double ffy = dy*(fo*sfe + sff*fullpot);
		double ffz = dz*(fo*sfe + sff*fullpot);
		mi.atom[s1].center.force.x -= ffx;
		mi.atom[s1].center.force.y -= ffy;
		mi.atom[s1].center.force.z -= ffz;
		mj.atom[s2].center.force.x += ffx;
		mj.atom[s2].center.force.y += ffy;
		mj.atom[s2].center.force.z += ffz;
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
	}
      }
    }
    pv.ep = ep;
    pv.vr = vr;
#ifdef DEBUG
    fprintf(stderr, "%d\n", ::count);
#endif
    return 1;
  }
  return 0;
}



int
force_offsetpbc(
	PolyatomicMols& m1,
	PolyatomicMols& m2,
	const Vector3 &offset, 
	const Box& box,
	int* isPeriodic,
	const Truncation& rc,
	PotVir &pv
	)
{
  assert(0);//Cutoff is wrong!!
  double ep = 0;
  Matrix33 vr;
  double bx = box.x / 2;
  double by = box.y / 2;
  double bz = box.z / 2;
  
  int nmol  = m1.Size();
  int nmol2 = m2.Size();
  MolPropertyHandle h1  = m1.GetProperty();
  Flexible*       p1  = dynamic_cast<Flexible*> ( h1.get() );
  MolPropertyHandle h2  = m2.GetProperty();
  Flexible*       p2  = dynamic_cast<Flexible*> ( h2.get() );
  if ( p1 && p2 ){
    int nsite1 = p1->GetNumSite();
    int nsite2 = p2->GetNumSite();
    const IntrParamsArray& ia1 = p1->GetIntrParams();
    const IntrParamsArray& ia2 = p2->GetIntrParams();
    for( int s1=0; s1<nsite1; s1++ ){
      for( int s2=0; s2<nsite2; s2++ ){
        const IntrParams& lji = ia1[s1];
        const IntrParams& ljj = ia2[s2];        
        {
          double sig = lji.sighalf + ljj.sighalf;
          double eps = lji.epssqrt * ljj.epssqrt;
          for( int i=0; i<nmol; i++ ){
	    PolyatomicMol& mi = m1.mols[i];
	    for( int j=0; j<nmol2; j++ ){
	      PolyatomicMol& mj = m2.mols[j];
	      double dx = mi.atom[s1].center.coord.x - mj.atom[s2].center.coord.x + offset.x;
	      double dy = mi.atom[s1].center.coord.y - mj.atom[s2].center.coord.y + offset.y;
	      double dz = mi.atom[s1].center.coord.z - mj.atom[s2].center.coord.z + offset.z;
	      if ( isPeriodic[0] ){
		if ( dx < -bx ) 
		  dx += box.x;
		else if ( bx < dx )
		  dx -= box.x;
	      }
	      if ( isPeriodic[1] ){
		if ( dy < -by ) 
		  dy += box.y;
		else if ( by < dy )
		  dy -= box.y;
	      }
	      if ( isPeriodic[2] ){
		if ( dz < -bz ) 
		  dz += box.z;
		else if ( bz < dz )
		  dz -= box.z;
	      }

	      double drr = dx*dx + dy*dy + dz*dz;
	      if( drr < rc.outer*rc.outer ){
		double sfe,sff,r;
		if( drr < rc.inner*rc.inner){
		  sfe = 1.0;
		  sff = 0.0;
		}
		else {
		  double rrl,rrs,rrs2,rrl2;
		  r = sqrt(drr);
		  rrl = r - rc.outer;
		  rrs = r - rc.inner;
		  rrl2 = rrl*rrl;
		  rrs2 = rrs*rrs;
		  sfe = rrl2*rrl*rc.rq*(10.0*rrs2-5.0*rrl*rrs+rrl2);
		  sff = 30.0*rc.rq*rrl2*rrs2/r;
		}
		double drs = sig * sig / drr;
		double dr6 = drs * drs * drs;
		double dr12 = dr6 * dr6;
		double fo = 4.0 * eps * (-12.0 * dr12 + 6.0 * dr6)/drr;
		double fullpot = 4.0 * eps * (dr12-dr6);
		ep += sfe*fullpot;
		double ffx = dx*(fo*sfe + sff*fullpot);
		double ffy = dy*(fo*sfe + sff*fullpot);
		double ffz = dz*(fo*sfe + sff*fullpot);
		mi.atom[s1].center.force.x -= ffx;
		mi.atom[s1].center.force.y -= ffy;
		mi.atom[s1].center.force.z -= ffz;
		mj.atom[s2].center.force.x += ffx;
		mj.atom[s2].center.force.y += ffy;
		mj.atom[s2].center.force.z += ffz;
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
	}
      }
    }
    pv.ep = ep;
    pv.vr = vr;
#ifdef DEBUG
    fprintf(stderr, "%d\n", ::count);
#endif
    return 1;
  }
  return 0;
}




int
force_pbc(
      PolyatomicMols& m1,
      PolyatomicMols& m2,
      const Box& box, 
      const Truncation& rc,
      PotVir &pv
      )
{
  double ep = 0;
  Matrix33 vr;
  double bx = box.x / 2;
  double by = box.y / 2;
  double bz = box.z / 2;
  
  int nmol  = m1.Size();
  int nmol2 = m2.Size();
  MolPropertyHandle h1  = m1.GetProperty();
  Flexible*       p1  = dynamic_cast<Flexible*> ( h1.get() );
  MolPropertyHandle h2  = m2.GetProperty();
  Flexible*       p2  = dynamic_cast<Flexible*> ( h2.get() );
  if ( p1 && p2 ){
    int nsite1 = p1->GetNumSite();
    int nsite2 = p2->GetNumSite();
    const IntrParamsArray& ia1 = p1->GetIntrParams();
    const IntrParamsArray& ia2 = p2->GetIntrParams();
    for( int s1=0; s1<nsite1; s1++ ){
      for( int s2=0; s2<nsite2; s2++ ){
        const IntrParams& lji = ia1[s1];
        const IntrParams& ljj = ia2[s2];        
        {
          double sig = lji.sighalf + ljj.sighalf;
          double eps = lji.epssqrt * ljj.epssqrt;
          for( int i=0; i<nmol; i++ ){
            PolyatomicMol& mi = m1.mols[i];
            for( int j=0; j<nmol2; j++ ){
              PolyatomicMol& mj = m2.mols[j];
              double dx = mi.atom[s1].center.coord.x - mj.atom[s2].center.coord.x;
              double dy = mi.atom[s1].center.coord.y - mj.atom[s2].center.coord.y;
              double dz = mi.atom[s1].center.coord.z - mj.atom[s2].center.coord.z;
              if ( dx < -bx ) 
                dx += box.x;
              else if ( bx < dx )
                dx -= box.x;
              if ( dy < -by ) 
                dy += box.y;
              else if ( by < dy )
                dy -= box.y;
              if ( dz < -bz ) 
                dz += box.z;
              else if ( bz < dz )
                dz -= box.z;
              
              double drr = dx*dx + dy*dy + dz*dz;
              if( drr < rc.outer*rc.outer ){
                double sfe,sff,r;
                if( drr < rc.inner*rc.inner){
                  sfe = 1.0;
                  sff = 0.0;
                }
                else {
                  double rrl,rrs,rrs2,rrl2;
                  r = sqrt(drr);
                  rrl = r - rc.outer;
                  rrs = r - rc.inner;
                  rrl2 = rrl*rrl;
                  rrs2 = rrs*rrs;
                  sfe = rrl2*rrl*rc.rq*(10.0*rrs2-5.0*rrl*rrs+rrl2);
                  sff = 30.0*rc.rq*rrl2*rrs2/r;
                }
                double drs = sig * sig / drr;
                double dr6 = drs * drs * drs;
                double dr12 = dr6 * dr6;
                double fo = 4.0 * eps * (-12.0 * dr12 + 6.0 * dr6)/drr;
                double fullpot = 4.0 * eps * (dr12-dr6);
                ep += sfe*fullpot;
                double ffx = dx*(fo*sfe + sff*fullpot);
                double ffy = dy*(fo*sfe + sff*fullpot);
                double ffz = dz*(fo*sfe + sff*fullpot);
                mi.atom[s1].center.force.x -= ffx;
                mi.atom[s1].center.force.y -= ffy;
                mi.atom[s1].center.force.z -= ffz;
                mj.atom[s2].center.force.x += ffx;
                mj.atom[s2].center.force.y += ffy;
                mj.atom[s2].center.force.z += ffz;
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
        }
      }
    }
    pv.ep = ep;
    pv.vr = vr;
#ifdef DEBUG
    fprintf(stderr, "%d\n", ::count);
#endif
    return 1;
  }
  return 0;
}




void PolyatomicMols::Report() const
{
}



const Vector3&
PolyatomicMols::Position( int i ) const
{
  return mols[i].Position();
}



void PolyatomicMols::PluginHookL0( CollectorPlugin& plugin )
{
  int nmol = mols.size();
  for(int i=0; i<nmol; i++)
    plugin.HookL0( mols[i] );
}
