#include "Interaction/ListVector.hpp"
#include "Interaction/Truncation.hpp"
#include "Mols/Mols.hpp"

int SimplePair::list_offset(
  const Mols& m1,
  const Mols& m2,
  const Vector3&     offset,
  double rc
  )
{
  const int nmol1 = m1.Size();
  const int nmol2 = m2.Size();
  for( int i=0; i<nmol1;i++ ){
    const Vector3& v1 = m1.Position(i);
    for( int j=0; j<nmol2; j++ ){
      const Vector3& v2 = m2.Position(j);
      Vector3 d0;
      d0.x = v1.x - v2.x + offset.x;
      d0.y = v1.y - v2.y + offset.y;
      d0.z = v1.z - v2.z + offset.z;
      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < rc*rc ){
	SimpleListItem li;
	li.pair1 = i;
	li.pair2 = j;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}



int SimplePair::list_offsetpbc(
  const Mols& m1,
  const Mols& m2,
  const Vector3&     offset,
  const Box& box,
  double rc
  )
{
  const int nmol1 = m1.Size();
  const int nmol2 = m2.Size();
  for( int i=0; i<nmol1;i++ ){
    const Vector3& v1 = m1.Position(i);
    for( int j=0; j<nmol2; j++ ){
      const Vector3& v2 = m2.Position(j);
      Vector3 d0;
      d0.x = v1.x - v2.x + offset.x;
      d0.y = v1.y - v2.y + offset.y;
      d0.z = v1.z - v2.z + offset.z;
      box.RelativePosition( d0 );

      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < rc*rc ){
	SimpleListItem li;
	li.pair1 = i;
	li.pair2 = j;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}



int SimplePair::list_pbc(
  const Mols& m1,
  const Mols& m2,
  const Box&     box,
  double rc
  )
{
  const int nmol1 = m1.Size();
  const int nmol2 = m2.Size();
  for( int i=0; i<nmol1;i++ ){
    const Vector3& v1 = m1.Position(i);
    for( int j=0; j<nmol2; j++ ){
      const Vector3& v2 = m2.Position(j);
      Vector3 d0;
      d0.x = v1.x - v2.x;
      d0.y = v1.y - v2.y;
      d0.z = v1.z - v2.z;
      box.RelativePosition( d0 );

      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < rc*rc ){
	SimpleListItem li;
	li.pair1 = i;
	li.pair2 = j;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}



int SimplePair::list_pbc(
  const Mols& m,
  const Box&     box,
  double rc
  )
{
  const int nmol = m.Size();
  for( int i=0; i<nmol;i++ ){
    const Vector3& v1 = m.Position(i);
    for( int j=i+1; j<nmol; j++ ){
      const Vector3& v2 = m.Position(j);
      Vector3 d0;
      d0.x = v1.x - v2.x;
      d0.y = v1.y - v2.y;
      d0.z = v1.z - v2.z;              
      box.RelativePosition( d0 );

      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < rc*rc ){
	SimpleListItem li;
	li.pair1 = i;
	li.pair2 = j;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}



int SimplePair::list_simple(
  const Mols& m1,
  const Mols& m2,
  double rc
  )
{
  const int nmol1 = m1.Size();
  const int nmol2 = m2.Size();
  for( int i=0; i<nmol1;i++ ){
    const Vector3& v1 = m1.Position( i );
    for( int j=0; j<nmol2; j++ ){
      const Vector3& v2 = m2.Position( j );
      Vector3 d0;
      //d0.x = r1.com.center.coord.x - r2.com.center.coord.x;
      //d0.y = r1.com.center.coord.y - r2.com.center.coord.y;
      //d0.z = r1.com.center.coord.z - r2.com.center.coord.z;
      d0.x = v1.x - v2.x;
      d0.y = v1.y - v2.y;
      d0.z = v1.z - v2.z;
      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < rc*rc ){
	SimpleListItem li;
	li.pair1 = i;
	li.pair2 = j;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}



int SimplePair::list_simple(
  const Mols& m,
  double rc
  )
{
  const int nmol = m.Size();
  for( int i=0; i<nmol;i++ ){
    const Vector3& v1 = m.Position( i );
    for( int j=i+1; j<nmol; j++ ){
      const Vector3& v2 = m.Position( j );
      Vector3 d0;
      d0.x = v1.x - v2.x;
      d0.y = v1.y - v2.y;
      d0.z = v1.z - v2.z;
      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < rc*rc ){
	SimpleListItem li;
	li.pair1 = i;
	li.pair2 = j;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}




int TruncPair::list_offset(
  const Mols& m1,
  const Mols& m2,
  const Vector3&     offset,
  const Truncation& rc
  )
{
  const int nmol1 = m1.Size();
  const int nmol2 = m2.Size();
  const double inner = rc.GetInner();
  const double outer = rc.GetOuter();
  for( int i=0; i<nmol1;i++ ){
    const Vector3& v1 = m1.Position(i);
    for( int j=0; j<nmol2; j++ ){
      const Vector3& v2 = m2.Position(j);
      Vector3 d0;
      d0.x = v1.x - v2.x + offset.x;
      d0.y = v1.y - v2.y + offset.y;
      d0.z = v1.z - v2.z + offset.z;
      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < outer*outer ){
	double sfe0, sff0;
	if( drr < inner*inner){
	  sfe0 = 1.0;
	  sff0 = 0.0;
	}
	else {
	  rc.Smooth( sqrt(drr), sfe0,sff0 );
	}
	ListItem li;
	li.d = d0;
	li.pair1 = i;
	li.pair2 = j;
	li.sfe = sfe0;
	li.sff = sff0;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}



int TruncPair::list_offsetpbc(
  const Mols& m1,
  const Mols& m2,
  const Vector3&     offset,
  const Box& box,
  const Truncation& rc
  )
{
  const int nmol1 = m1.Size();
  const int nmol2 = m2.Size();
  const double inner = rc.GetInner();
  const double outer = rc.GetOuter();
  for( int i=0; i<nmol1;i++ ){
    const Vector3& v1 = m1.Position(i);
    for( int j=0; j<nmol2; j++ ){
      const Vector3& v2 = m2.Position(j);
      Vector3 d0;
      d0.x = v1.x - v2.x + offset.x;
      d0.y = v1.y - v2.y + offset.y;
      d0.z = v1.z - v2.z + offset.z;
      //cerr << d0.x << " " << d0.y << " " << d0.z << endl;
      box.RelativePosition( d0 );
      //cerr << d0.x << " " << d0.y << " " << d0.z << endl;

      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < outer*outer ){
	double sfe0, sff0;
	if( drr < inner*inner){
	  sfe0 = 1.0;
	  sff0 = 0.0;
	}
	else {
	  rc.Smooth( sqrt(drr), sfe0,sff0 );
	}
	ListItem li;
	li.d = d0;
	li.pair1 = i;
	li.pair2 = j;
	li.sfe = sfe0;
	li.sff = sff0;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}



int TruncPair::list_pbc(
  const Mols& m1,
  const Mols& m2,
  const Box&     box,
  const Truncation& rc
  )
{
  const int nmol1 = m1.Size();
  const int nmol2 = m2.Size();
  const double inner = rc.GetInner();
  const double outer = rc.GetOuter();
  for( int i=0; i<nmol1;i++ ){
    const Vector3& v1 = m1.Position(i);
    for( int j=0; j<nmol2; j++ ){
      const Vector3& v2 = m2.Position(j);
      Vector3 d0;
      d0.x = v1.x - v2.x;
      d0.y = v1.y - v2.y;
      d0.z = v1.z - v2.z;
      box.RelativePosition( d0 );

      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < outer*outer ){
	double sfe0, sff0;
	if( drr < inner*inner){
	  sfe0 = 1.0;
	  sff0 = 0.0;
	}
	else {
	  rc.Smooth( sqrt(drr), sfe0,sff0 );
	}
	ListItem li;
	li.d = d0;
	li.pair1 = i;
	li.pair2 = j;
	li.sfe = sfe0;
	li.sff = sff0;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}



int TruncPair::list_pbc(
  const Mols& m,
  const Box&     box,
  const Truncation& rc
  )
{
  const int nmol = m.Size();
  const double inner = rc.GetInner();
  const double outer = rc.GetOuter();
  for( int i=0; i<nmol;i++ ){
    const Vector3& v1 = m.Position(i);
    for( int j=i+1; j<nmol; j++ ){
      const Vector3& v2 = m.Position(j);
      Vector3 d0;
      d0.x = v1.x - v2.x;
      d0.y = v1.y - v2.y;
      d0.z = v1.z - v2.z;              
      box.RelativePosition( d0 );

      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < outer*outer ){
	double sfe0, sff0;
	if( drr < inner*inner){
	  sfe0 = 1.0;
	  sff0 = 0.0;
	}
	else {
	  rc.Smooth( sqrt(drr), sfe0,sff0 );
	}
	ListItem li;
	li.d = d0;
	li.pair1 = i;
	li.pair2 = j;
	li.sfe = sfe0;
	li.sff = sff0;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}



int TruncPair::list_simple(
  const Mols& m1,
  const Mols& m2,
  const Truncation& rc
  )
{
  const int nmol1 = m1.Size();
  const int nmol2 = m2.Size();
  const double inner = rc.GetInner();
  const double outer = rc.GetOuter();
  for( int i=0; i<nmol1;i++ ){
    const Vector3& v1 = m1.Position( i );
    for( int j=0; j<nmol2; j++ ){
      const Vector3& v2 = m2.Position( j );
      Vector3 d0;
      //d0.x = r1.com.center.coord.x - r2.com.center.coord.x;
      //d0.y = r1.com.center.coord.y - r2.com.center.coord.y;
      //d0.z = r1.com.center.coord.z - r2.com.center.coord.z;
      d0.x = v1.x - v2.x;
      d0.y = v1.y - v2.y;
      d0.z = v1.z - v2.z;
      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < outer*outer ){
	double sfe0, sff0;
	if( drr < inner*inner){
	  sfe0 = 1.0;
	  sff0 = 0.0;
	}
	else {
	  rc.Smooth( sqrt(drr), sfe0,sff0 );
	}
	ListItem li;
	li.d = d0;
	li.pair1 = i;
	li.pair2 = j;
	li.sfe = sfe0;
	li.sff = sff0;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}



int TruncPair::list_simple(
  const Mols& m,
  const Truncation& rc
  )
{
  const int nmol = m.Size();
  const double inner = rc.GetInner();
  const double outer = rc.GetOuter();
  for( int i=0; i<nmol;i++ ){
    const Vector3& v1 = m.Position( i );
    for( int j=i+1; j<nmol; j++ ){
      const Vector3& v2 = m.Position( j );
      Vector3 d0;
      d0.x = v1.x - v2.x;
      d0.y = v1.y - v2.y;
      d0.z = v1.z - v2.z;
      double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
      if( drr < outer*outer ){
	double sfe0, sff0;
	if( drr < inner*inner){
	  sfe0 = 1.0;
	  sff0 = 0.0;
	}
	else {
	  rc.Smooth( sqrt(drr), sfe0,sff0 );
	}
	ListItem li;
	li.d = d0;
	li.pair1 = i;
	li.pair2 = j;
	li.sfe = sfe0;
	li.sff = sff0;
	listitems.push_back(li);
      }
    }
  }
  return listitems.size();
}




//使われていない?たぶんintel CPUだと圧縮しても速度が変わらない。
TruncPair*
SimplePair::Compress(
		     const Mols& m1,
		     const Mols& m2,
		     const Box& box,
		     const Truncation& rc
		     )
{
  TruncPair* truncpair = new TruncPair;
  const int npair = listitems.size();
  const double inner = rc.GetInner();
  const double outer = rc.GetOuter();
  for(int k=0; k<npair; k++){
    const int i = listitems[k].pair1;
    const int j = listitems[k].pair2;
    const Vector3& v1 = m1.Position(i);
    const Vector3& v2 = m2.Position(j);
    Vector3 d0;
    d0.x = v1.x - v2.x;
    d0.y = v1.y - v2.y;
    d0.z = v1.z - v2.z;
    box.RelativePosition( d0 );

    double drr = (d0.x*d0.x + d0.y*d0.y + d0.z*d0.z);
    if( drr < outer*outer ){
      double sfe0, sff0;
      if( drr < inner*inner){
	sfe0 = 1.0;
	sff0 = 0.0;
      }
      else {
	rc.Smooth( sqrt(drr), sfe0,sff0 );
      }
      ListItem li;
      li.d = d0;
      li.pair1 = i;
      li.pair2 = j;
      li.sfe = sfe0;
      li.sff = sff0;
      truncpair->push_back(li);
    }
  }
  return truncpair;
}
