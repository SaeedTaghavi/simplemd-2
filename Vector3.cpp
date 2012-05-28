#include <string>
#include <cassert>
#include <iostream>
#include <cstdio>
#include "Vector3.hpp"
#include "debug.hpp"

using namespace SI;

/*Class of 3D vector*/
Vector3::Vector3( double _x, double _y, double _z )
{
  Set( _x, _y, _z );
}



Vector3::Vector3()
{
  //Set(0,0,0);
}



void Vector3::Set( double _x, double _y, double _z ){
  x = _x;
  y = _y;
  z = _z;
}



void Vector3::Scale( double ratio )
{
  x *= ratio;
  y *= ratio;
  z *= ratio;
}



void Vector3::Scale( const Vector3& ratio )
{
  x *= ratio.x;
  y *= ratio.y;
  z *= ratio.z;
}



Box::Box() : Vector3( 1.0L, 1.0L, 1.0L ), volume( 1.0L ){}



Box::Box( double xx, double yy, double zz ) : Vector3( xx, yy, zz ), volume( xx*yy*zz ){}



void
Box::Set( double xx, double yy, double zz )
{
  Vector3::Set( xx,yy,zz );
  volume = xx*yy*zz;
}



void
Box::Set( double newvol )
{
  double ratio = pow( newvol / volume, 1.0/3.0 );
  x *= ratio;
  y *= ratio;
  z *= ratio;
  volume = newvol;
}



void
Box::RelativePosition( Vector3 &coord ) const
{
  coord.x = applyPBC0( coord.x, x / 2 );
  coord.y = applyPBC0( coord.y, y / 2 );
  coord.z = applyPBC0( coord.z, z / 2 );
}






bool Box::Contains( const Vector3 &coord ) const
{
  return ( 0       <= coord.x &&
	   coord.x <  x       &&
	   0       <= coord.y &&
	   coord.y <  y       &&
	   0       <= coord.z &&
	   coord.z <  z       );
}






void
Box::Scale( double ratio )
{
  Vector3::Scale( ratio );
  volume = x*y*z;
}



void
Box::Scale( const Vector3& ratio )
{
  Vector3::Scale( ratio );
  volume = x*y*z;
}



SemiperiodicBox::SemiperiodicBox( const int isP[3] )
{
  for( int i=0; i<3; i++ ){
    isPeriodic[i] = isP[i];
  }
}



void 
SemiperiodicBox::BoxCoordinate( Vector3 &coord ) const
{
  /*
    replaced by slow BoxCoordinate3 2006-5-10!!!
  coord.x = BoxCoordinate2( coord.x, x );
  coord.y = BoxCoordinate2( coord.y, y );
  coord.z = BoxCoordinate2( coord.z, z );
  */
  if ( isPeriodic[0] )
    coord.x = BoxCoordinate3( coord.x, x );
  if ( isPeriodic[1] )
    coord.y = BoxCoordinate3( coord.y, y );
  if ( isPeriodic[2] )
    coord.z = BoxCoordinate3( coord.z, z );
}



void
SemiperiodicBox::RelativePosition( Vector3 &coord ) const
{
  if ( isPeriodic[0] )
    coord.x = applyPBC0( coord.x, x / 2 );
  if ( isPeriodic[1] )
    coord.y = applyPBC0( coord.y, y / 2 );
  if ( isPeriodic[2] )
    coord.z = applyPBC0( coord.z, z / 2 );
}



bool
SemiperiodicBox::Contains( const Vector3& coord ) const
{
  bool isOutside = 0;
  //cerr << "SemiperiodicBox::Contains ";
  if ( ! isPeriodic[0] ){
    isOutside = isOutside || ( coord.x < 0 || x <= coord.x );
  }
  if ( ! isPeriodic[1] ){
    isOutside = isOutside || ( coord.y < 0 || y <= coord.y );
  }
  if ( ! isPeriodic[2] ){
    isOutside = isOutside || ( coord.z < 0 || z <= coord.z );
  }
  //cerr << coord.x << " " << coord.y << " " << coord.z;
  if ( isOutside ){
    //cerr << " outsize";
  }
  //cerr << endl;
  return ! isOutside;
}



void
Box::WriteBOX3( const Unit& u, ostream& to ) const
{
  to << "@BOX3\n"
     << x / ( Angstro * u.m )
     << " " << y / ( Angstro * u.m )
     << " " << z / ( Angstro * u.m ) << '\n';
}



//read initial values from file
void
Box::ReadBOX3( const Unit& u, FILE* file )
{
#define BUFSIZE 1024
  char buf[BUFSIZE];
  double x,y,z;
  
  fgets(buf,sizeof(buf),file);
  sscanf( buf, "%lf %lf %lf", &x, &y, &z );
  x *= Angstro * u.m;
  y *= Angstro * u.m;
  z *= Angstro * u.m;
  Set( x,y,z );
}



//read initial values from file
void
Box::ReadBXLA( const Unit& u, FILE* file )
{
#define BUFSIZE 1024
  char buf[BUFSIZE];
  double x,y,z;
  
  fgets(buf,sizeof(buf),file);
  sscanf( buf, "%lf", &x );
  x *= Angstro * u.m;
  y = z = x;
  Set( x,y,z );
}



/*nph.sh 1.18sec user; inline*/
double applyPBC3(double x,double bxlh)
{
  if ( x > bxlh ){
    x -= bxlh*2;
  }
  else if ( x < -bxlh )
    x += bxlh*2;
  return x;
}

/*nph.sh 1.26sec user; */
double applyPBC0(double x,double bxlh)
{
  while ( x > bxlh )
    x -= bxlh*2;
  while ( x < -bxlh )
    x += bxlh*2;
  return x;
}

/*nph.sh 17.5sec user; */
double applyPBC1(double x,double bxlh)
{
  return x - rint( x / (bxlh*2) ) * bxlh*2;
}


/*nph.sh 2.98sec user; */
double applyPBC2(double x,double bxlh)
{
  return x - (int)( x / (bxlh*2)  + 0.5 ) * bxlh*2;
}

/*nph.sh 2.98sec user; */
/*it is not correct when x is negative */
double BoxCoordinate2(double x,double bxl)
{
  return x - (int)( x / bxl ) * bxl;
}


double BoxCoordinate3(double x,double bxl)
{
  while( x < 0 ){
    x += bxl;
  }
  while( bxl <= x ){
    x -= bxl;
  }
  return x;
}



/*nph.sh 2.98sec user; */
double ImageOffset(double x,double bxl)
{
  return (int)( x / bxl ) * bxl;
}




Matrix33::Matrix33()
{
  Clear();
}



void Matrix33::Clear()
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      v[i][j] = 0;
    }
  }
}



void Matrix33::Set( double qa, double qb, double qc, double qd )
{
  v[0][0] = qa*qa+qb*qb-(qc*qc+qd*qd);
  v[0][1] =-2.0*(qa*qd+qb*qc);
  v[0][2] = 2.0*(qb*qd-qa*qc);
  v[1][0] = 2.0*(qa*qd-qb*qc);
  v[1][1] = qa*qa+qc*qc-(qb*qb+qd*qd);
  v[1][2] =-2.0*(qa*qb+qc*qd);
  v[2][0] = 2.0*(qa*qc+qb*qd);
  v[2][1] = 2.0*(qa*qb-qc*qd);
  v[2][2] = qa*qa+qd*qd-(qb*qb+qc*qc);
}



void
Matrix33::Quat( Quaternion& q )
{
  int m=0;
  //最も精度が稼げる射影を探す。2006-6-23
  double aa = (1.0+v[0][0]+v[1][1]+v[2][2]);
  double bb = (1.0+v[0][0]-v[1][1]-v[2][2]);
  double cc = (1.0-v[0][0]+v[1][1]-v[2][2]);
  double dd = (1.0-v[0][0]-v[1][1]+v[2][2]);
  double max=aa;
  if ( max < bb ){
    m=1;
    max = bb;
  }
  if ( max < cc ){
    m=2;
    max = cc;
  }
  if ( max < dd ){
    m=3;
    max = dd;
  }
  if ( m == 0 ){
    q.a=sqrt(aa/4);
    double ab = (v[2][1]-v[1][2])/4;
    double ac = (v[2][0]-v[0][2])/4;
    double ad = (v[1][0]-v[0][1])/4;
    q.b = ab/q.a;
    q.c = ac/q.a;
    q.d = ad/q.a;
    //cerr << q.a << " " << q.b << " " << q.c << " " << q.d << " " << endl;
  }
  else if ( m == 3 ){
    q.d=sqrt(dd/4);
    double ad = (v[1][0]-v[0][1])/4;
    double bd = (v[2][0]+v[0][2])/4;
    double cd =-(v[2][1]+v[1][2])/4;
    q.a = ad / q.d;
    q.b = bd / q.d;
    q.c = cd / q.d;
    //cerr << q.a << " " << q.b << " " << q.c << " " << q.d << " " << endl;
  }
  else if ( m == 1 ){
    q.b=sqrt(bb/4);
    double ab = (v[2][1]-v[1][2])/4;
    double bc =-(v[1][0]+v[0][1])/4;
    double bd = (v[2][0]+v[0][2])/4;
    q.a = ab/q.b;
    q.c = bc/q.b;
    q.d = bd/q.b;
    //cerr << q.a << " " << q.b << " " << q.c << " " << q.d << " " << endl;
  }
  else { // ( m == 2 ){
    q.c=sqrt(cc/4);
    double ac = (v[2][0]-v[0][2])/4;
    double bc =-(v[1][0]+v[0][1])/4;
    double cd =-(v[2][1]+v[1][2])/4;
    q.a = ac / q.c;
    q.b = bc / q.c;
    q.d = cd / q.c;
    //cerr << q.a << " " << q.b << " " << q.c << " " << q.d << " " << endl;
  }
  //cerr << q.a << " " << q.b << " " << q.c << " " << q.d << " " <<
  //  q.a*q.a+q.b*q.b+q.c*q.c+q.d*q.d << endl;
}
