#include <cstring>
#include <cstdio>
#include <cassert>
#include "SingleMol/RigidBody.hpp"
#include "debug.hpp"
#include "Interaction/Truncation.hpp"
#include "Interaction/PotVir.hpp"
#include "Interaction/Combination.hpp"
#include "Interaction/ListVector.hpp"


RigidBody::~RigidBody()
{
  mesg("~RigidBody\n");
}



//Copy constructor
/*
RigidBody::RigidBody( const RigidBody& m )
{
  *this = m;
}
*/


RigidBody::RigidBody( int nsite, const MolPropertyHandle& p )
  : PolyatomicMol( nsite, vector<double>(nsite) ), com( MonatomicMol( p->GetMass() ) )
{
  prop = p;
  order = -1;
  Rigid* r = dynamic_cast<Rigid*> ( p.get() );
  assert( r != 0 );
  in = r->GetInertia();
  //intra.resize( nsite );
  w.Set(0,0,0);
  dq.Set(0,0,0,0);
  torque.Set(0,0,0);
}



const Vector3&
RigidBody::Position() const
{
  return com.center.coord;
}



void RigidBody::Translate( const Vector3& offset )
{
    com.Translate( offset );
}



void RigidBody::resetsiteforce()
{
  com.center.force.Set(0,0,0);
  torque.Set(0,0,0);
  int nsite = atom.size();
  for( int m=0; m<nsite; m++ )
    atom[m].center.force.Set(0,0,0);
}



void RigidBody::preparerotationmatrix()
{
  double a = q.a;
  double b = q.b;
  double c = q.c;
  double d = q.d;
  rot.v[0][0] = (a*a+b*b-(c*c+d*d));
  rot.v[0][1] = -2.0*(a*d+b*c);
  rot.v[0][2] = 2.0*(b*d-a*c);
  rot.v[1][0] = 2.0*(a*d-b*c);
  rot.v[1][1] = a*a+c*c-(b*b+d*d);
  rot.v[1][2] = -2.0*(a*b+c*d);
  rot.v[2][0] = 2.0*(a*c+b*d);
  rot.v[2][1] = 2.0*(a*b-c*d);
  rot.v[2][2] = a*a+d*d-(b*b+c*c);
  //cerr << w.x << " " << rot.v[1][2] << " (1)\n";
}



void RigidBody::prepareintra()
{
  Rigid* p = dynamic_cast<Rigid*> ( prop.get() );
  assert( p != 0 );
  int nsite = atom.size();
  for( int s=0; s<nsite; s++ ){
    atom[s].center.coord.x = rot.v[0][0]*p->site[s].x + rot.v[0][1]*p->site[s].y + rot.v[0][2]*p->site[s].z;
  }
  for( int s=0; s<nsite; s++ ){
    atom[s].center.coord.y = rot.v[1][0]*p->site[s].x + rot.v[1][1]*p->site[s].y + rot.v[1][2]*p->site[s].z;
  }
  for( int s=0; s<nsite; s++ ){
    atom[s].center.coord.z = rot.v[2][0]*p->site[s].x + rot.v[2][1]*p->site[s].y + rot.v[2][2]*p->site[s].z;
  }
}



void RigidBody::Preforce()
{
  resetsiteforce();
  preparerotationmatrix();
  prepareintra();
  //cerr << order << " " << com.center.coord.x << " com.x\n";
}



void RigidBody::CollectForce()
{
  int nsite = atom.size();
  //collect force and calculate torque
  Vector3 torq,forc;
  forc.x = 0;
  forc.y = 0;
  forc.z = 0;
  for( int s=0; s<nsite; s++ ){
    //double x = atom[s].center.force.x;
    //double y = atom[s].center.force.y;
    //double z = atom[s].center.force.z;
    //printf("%f %f %f\n",x,y,z);
    forc.x += atom[s].center.force.x;
    forc.y += atom[s].center.force.y;
    forc.z += atom[s].center.force.z;
  }

  com.center.force.x += forc.x;
  com.center.force.y += forc.y;
  com.center.force.z += forc.z;
  
  torq.x = 0;
  torq.y = 0;
  torq.z = 0;
  for( int s=0; s<nsite; s++ ){
    SiteOfAction& soa = atom[s].center;
    torq.x += soa.coord.z*soa.force.y - soa.coord.y*soa.force.z;
    torq.y += soa.coord.x*soa.force.z - soa.coord.z*soa.force.x;
    torq.z += soa.coord.y*soa.force.x - soa.coord.x*soa.force.y;
  }
  torque.x = rot.v[0][0]*torq.x + rot.v[1][0]*torq.y + rot.v[2][0]*torq.z;
  torque.y = rot.v[0][1]*torq.x + rot.v[1][1]*torq.y + rot.v[2][1]*torq.z;
  torque.z = rot.v[0][2]*torq.x + rot.v[1][2]*torq.y + rot.v[2][2]*torq.z;
}



void RigidBody::Postforce( double mass )
{
  CollectForce();
}



void RigidBody::ScaleVelocity( const Vector3& r )
{
    com.ScaleVelocity( r );
}



void RigidBody::ScaleVelocity( double r )
{
  com.ScaleVelocity( r );
  w.x *= r;
  w.y *= r;
  w.z *= r;
}



void RigidBody::ScalePosition( const Vector3& r )
{
    com.ScalePosition( r );
}



void RigidBody::AddVelocity( const Vector3& v )
{
    com.AddVelocity( v );
}



void
RigidBody::BoxCoordinate( const Box& box )
{
  box.BoxCoordinate( com.center.coord );
}



void
RigidBody::ProgressPosition( double dt )
{
  com.ProgressPosition( dt );
  q.a += dq.a * dt;
  q.b += dq.b * dt;
  q.c += dq.c * dt;
  q.d += dq.d * dt;
  q.normalize();
}



void
RigidBody::ProgressMomentum( double dt )
{
  Vector3   dw;
  com.ProgressMomentum( dt );
  //6.36  Goldstein P.268 (really? !!!)
  dw.x = ( -torque.x + ( in.y - in.z ) * w.y * w.z ) / in.x;
  dw.y = ( -torque.y + ( in.z - in.x ) * w.z * w.x ) / in.y;
  dw.z = ( -torque.z + ( in.x - in.y ) * w.x * w.y ) / in.z;
  w.x += dw.x * dt;
  w.y += dw.y * dt;
  w.z += dw.z * dt;  
  dq.c = (-w.x*q.d - w.y*q.a + w.z*q.b) * 0.5;
  dq.b = ( w.x*q.a - w.y*q.d - w.z*q.c) * 0.5;
  dq.d = ( w.x*q.c + w.y*q.b + w.z*q.a) * 0.5;
  dq.a = (-w.x*q.b + w.y*q.c - w.z*q.d) * 0.5;
  //cerr << in.x << " " << torque.x << " " << dt << " " << dw.x << "\n";
}




double
Potential(
	  const RigidBody& r1,
	  const RigidBody& r2,
	  const Intersite& im,
	  const ListItem& li
	  )
{
  double epsum = 0;
  int npair = im.intr.size();
  for( int k=0; k<npair; k++ ){
    int s1 = im.site1[k];
    int s2 = im.site2[k];
    const SiteOfAction& soa1 = r1.atom[s1].center;
    const SiteOfAction& soa2 = r2.atom[s2].center;
    Vector3 delta;
    double pot;
    delta.x  = li.d.x + soa1.coord.x - soa2.coord.x;
    delta.y  = li.d.y + soa1.coord.y - soa2.coord.y;
    delta.z  = li.d.z + soa1.coord.z - soa2.coord.z;
    double fo;//dummy
    im.intr[k]->Force( delta, fo, pot );
    epsum += pot;
  }
  return li.sfe * epsum;
}



double
Potential(
	  const RigidBody& r1,
	  const MonatomicMol& r2,
	  const Intersite& im,
	  const ListItem& li
	  )
{
  double epsum = 0;
  int npair = im.intr.size();
  for( int k=0; k<npair; k++ ){
    int s1 = im.site1[k];
    //int s2 = im.site2[k];
    const SiteOfAction& soa1 = r1.atom[s1].center;
    //SiteOfAction& soa2 = r2.center;
    Vector3 delta;
    double pot;
    delta.x  = li.d.x + soa1.coord.x;
    delta.y  = li.d.y + soa1.coord.y;
    delta.z  = li.d.z + soa1.coord.z;
    double fo;//dummy
    im.intr[k]->Force( delta, fo, pot );
    epsum += pot;
  }
  return li.sfe * epsum;
}



double
Potential(
	  const SingleMolEntity& r1,
	  const SingleMolEntity& r2,
	  const Intersite& im,
	  const ListItem& li
	  )
{
  const RigidBody& rr1 = dynamic_cast<const RigidBody&> ( r1 );
  assert( &rr1 != 0 );
  const RigidBody& rr2 = dynamic_cast<const RigidBody&> ( r2 );
  if ( &rr2 != 0 ){
    return Potential( rr1, rr2, im, li );
  }
  else {
    const MonatomicMol& rr2 = dynamic_cast<const MonatomicMol&> ( r2 );
    assert ( &rr2 != 0 );
    return Potential( rr1, rr2, im, li );
  }
}

    


/*
bool
force_core_core(
		RigidBody& r1,
		RigidBody& r2,
		const Intersite& im,
		const Vector3& d0,
		PotVir &pv,
		double sfe,
		double sff
		)
{
  double ep=0;
  vector<Vector3> force1(im.nsite1);
  vector<Vector3> force2(im.nsite2);
  int npair = im.intr.size();
  for( int k=0; k<npair; k++ ){
    Vector3 delta;
    double  potential, fo;
    int s1 = im.site1[k];
    int s2 = im.site2[k];
    delta.x  = d0.x + r1.atom[s1].center.coord.x - r2.atom[s2].center.coord.x;
    delta.y  = d0.y + r1.atom[s1].center.coord.y - r2.atom[s2].center.coord.y;
    delta.z  = d0.z + r1.atom[s1].center.coord.z - r2.atom[s2].center.coord.z;
    im.intr[k]->Force( delta, fo, potential );
    ep += potential;
    double ffx = delta.x*fo*sfe;
    double ffy = delta.y*fo*sfe;
    double ffz = delta.z*fo*sfe;
    force1[s1].x += ffx;
    force1[s1].y += ffy;
    force1[s1].z += ffz;
    force2[s2].x -= ffx;
    force2[s2].y -= ffy;
    force2[s2].z -= ffz;
#ifdef DEBUG
    //printf("# %d %d %f\n", order[i], order[j], r);
    ::count++;
#endif
  }
  Vector3 comf;
  comf.x = -sff * ep * d0.x;
  comf.y = -sff * ep * d0.y;
  comf.z = -sff * ep * d0.z;
  r1.com.center.force.x += comf.x;
  r1.com.center.force.y += comf.y;
  r1.com.center.force.z += comf.z;
  r2.com.center.force.x -= comf.x;
  r2.com.center.force.y -= comf.y;
  r2.com.center.force.z -= comf.z;
  double epsum = ep * sfe;
  //printf("%f pv.ep\n", pv.ep);
  for(int s1=0;s1<im.nsite1;s1++){
    r1.atom[s1].center.force.x += force1[s1].x;
    r1.atom[s1].center.force.y += force1[s1].y;
    r1.atom[s1].center.force.z += force1[s1].z;
    comf.x += force1[s1].x;
    comf.y += force1[s1].y;
    comf.z += force1[s1].z;
  }
  for(int s2=0;s2<im.nsite2;s2++){
    r2.atom[s2].center.force.x += force2[s2].x;
    r2.atom[s2].center.force.y += force2[s2].y;
    r2.atom[s2].center.force.z += force2[s2].z;
  }
  Matrix33 vrsum;
  vrsum.v[0][0] = d0.x * comf.x;
  vrsum.v[0][1] = d0.x * comf.y;
  vrsum.v[0][2] = d0.x * comf.z;
  vrsum.v[1][0] = d0.y * comf.x;
  vrsum.v[1][1] = d0.y * comf.y;
  vrsum.v[1][2] = d0.y * comf.z;
  vrsum.v[2][0] = d0.z * comf.x;
  vrsum.v[2][1] = d0.z * comf.y;
  vrsum.v[2][2] = d0.z * comf.z;
  pv.Set( epsum, vrsum );
  return 1;
}
*/





void
RigidBody::WriteWTG5( const Unit& u, ostream& to )
{
  const Vector3& velocity = com.GetVelocity();
  //運動量ではなく速度を出力したいぞ。
  to << com.center.coord.x / ( Angstro * u.m )
     << " " <<com.center.coord.y / ( Angstro * u.m )
     << " " << com.center.coord.z / ( Angstro * u.m )
     << " " << q.a
     << " " << q.b
     << " " << q.c
     << " " << q.d
    //convert from internal unit to A/ps
     << " " << velocity.x / ( Angstro * u.m / (SI::pico * u.sec) )
     << " " << velocity.y / ( Angstro * u.m / (SI::pico * u.sec) )
     << " " << velocity.z / ( Angstro * u.m / (SI::pico * u.sec) )
    //   unit is /ps
     << " " << dq.a / ( 1 / (SI::pico * u.sec) )
     << " " << dq.b / ( 1 / (SI::pico * u.sec) )
     << " " << dq.c / ( 1 / (SI::pico * u.sec) )
     << " " << dq.d / ( 1 / (SI::pico * u.sec) )
    //  unit is rad/ps
     << " " << w.x / ( u.radian / (SI::pico * u.sec) )
     << " " << w.y / ( u.radian / (SI::pico * u.sec) )
     << " " << w.z / ( u.radian / (SI::pico * u.sec) )
    //  unit is kJ/mol/A
     << " " << com.center.force.x / ( SI::kilo * u.J / (mol * Angstro * u.m ) )
     << " " << com.center.force.y / ( SI::kilo * u.J / (mol * Angstro * u.m ) )
     << " " << com.center.force.z / ( SI::kilo * u.J / (mol * Angstro * u.m ) )
    //  unit is kJ/mol/A * A == kJ/mol(!)
     << " " << torque.x / ( SI::kilo * u.J / mol )
     << " " << torque.y / ( SI::kilo * u.J / mol )
     << " " << torque.z / ( SI::kilo * u.J / mol )
     << " " << order << '\n';
}



void
RigidBody::WriteNX4A( const Unit& u, ostream& to )
{
  to << com.center.coord.x / ( Angstro * u.m )
     << " " <<com.center.coord.y / ( Angstro * u.m )
     << " " << com.center.coord.z / ( Angstro * u.m );
  //Quaternion q;
  //rot.Quat( q );
  to << " " << q.a
     << " " << q.b
     << " " << q.c
     << " " << q.d << "\n";
}





void
RigidBody::ReadWTG5( int o, const Unit& u, FILE* file )
{
  char buf[1000];
  order = o;
  
  fgets(buf,sizeof(buf),file);
  //split
  com.center.coord.x = atof(strtok(buf," \t")) * Angstro * u.m;
  com.center.coord.y = atof(strtok(NULL," \t")) * Angstro * u.m;
  com.center.coord.z = atof(strtok(NULL," \t")) * Angstro * u.m;
  q.a    = atof(strtok(NULL," \t"));
  q.b    = atof(strtok(NULL," \t"));
  q.c    = atof(strtok(NULL," \t"));
  q.d    = atof(strtok(NULL," \t"));
  Vector3 velocity;
  velocity.x = atof(strtok(NULL," \t")) * ( Angstro * u.m / (SI::pico * u.sec) );
  velocity.y = atof(strtok(NULL," \t")) * ( Angstro * u.m / (SI::pico * u.sec) );
  velocity.z = atof(strtok(NULL," \t")) * ( Angstro * u.m / (SI::pico * u.sec) );
  com.SetVelocity( velocity );
  dq.a    = atof(strtok(NULL," \t")) * ( 1 / (SI::pico * u.sec) );
  dq.b    = atof(strtok(NULL," \t")) * ( 1 / (SI::pico * u.sec) );
  dq.c    = atof(strtok(NULL," \t")) * ( 1 / (SI::pico * u.sec) );
  dq.d    = atof(strtok(NULL," \t")) * ( 1 / (SI::pico * u.sec) );
  w.x    = atof(strtok(NULL," \t")) * ( u.radian / (SI::pico * u.sec) );
  w.y    = atof(strtok(NULL," \t")) * ( u.radian / (SI::pico * u.sec) );
  w.z    = atof(strtok(NULL," \t")) * ( u.radian / (SI::pico * u.sec) );
  com.center.force.x    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / (mol * Angstro * u.m ) );
  com.center.force.y    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / (mol * Angstro * u.m ) );
  com.center.force.z    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / (mol * Angstro * u.m ) );
  torque.x    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / mol );
  torque.y    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / mol );
  torque.z    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / mol );
  rot.Set( q.a, q.b, q.c, q.d );
}



void
RigidBody::ReadNX4A( int o, const Unit& u, FILE* file )
{
  char buf[1000];
  order = o;
  
  fgets(buf,sizeof(buf),file);
  //split
  com.center.coord.x = atof(strtok(buf," \t")) * Angstro * u.m;;
  com.center.coord.y = atof(strtok(NULL," \t")) * Angstro * u.m;;
  com.center.coord.z = atof(strtok(NULL," \t")) * Angstro * u.m;;
  q.a    = atof(strtok(NULL," \t"));
  q.b    = atof(strtok(NULL," \t"));
  q.c    = atof(strtok(NULL," \t"));
  q.d    = atof(strtok(NULL," \t"));
  rot.Set( q.a, q.b, q.c, q.d );
}



double*
RigidBody::unserialize( double* const p )
{
  com.unserialize( p );
  q.Set( p[3], p[4], p[5], p[6] );
  q.normalize();
  return p+7;
}



double*
RigidBody::serialize( double* p ) const
{
  com.serialize( p );
  p[3] = q.a;
  p[4] = q.b;
  p[5] = q.c;
  p[6] = q.d;
  return p+7;
}



double*
RigidBody::serializeforce( double* xi ) const
{
  com.serializeforce( xi );
  /*
  xi[0] = 0;
  xi[1] = 0;
  xi[2] = 0;
  */
  //quaternion is set
  //get torque
  const Vector3& in = GetMomentOfInertia();
  double dwx = torque.x / in.x;
  double dwy = torque.y / in.y;
  double dwz = torque.z / in.z;
  xi[3] = -(-dwx*q.b + dwy*q.c - dwz*q.d) * 0.5;
  xi[4] = -( dwx*q.a - dwy*q.d - dwz*q.c) * 0.5;
  xi[5] = -(-dwx*q.d - dwy*q.a + dwz*q.b) * 0.5;
  xi[6] = -( dwx*q.c + dwy*q.b + dwz*q.a) * 0.5;
  /*
  xi[3] = 0;
  xi[4] = 0;
  xi[5] = 0;
  xi[6] = 0;
  cout << q.a 
       << "/" << q.b
       << "/" << q.c
       << "/" << q.d << endl;
  */
  return xi+7;
}



double RigidBody::GetEk() const
{
  double ek = com.GetEk();
  double vx,vy,vz;
  vx = w.x * w.x;
  vy = w.y * w.y;
  vz = w.z * w.z;
  ek += 0.5 * (vx*in.x + vy*in.y + vz*in.z);
  return ek;
}
