#include <cstdio>
#include "SingleMol/RigidBody2.hpp"
#include "debug.hpp"

RigidBody2::~RigidBody2()
{
  mesg("~RigidBody2\n");
}



RigidBody2::RigidBody2( int nsite, const MolPropertyHandle& p )
  : RigidBody( nsite, p )
{
}






void RigidBody2::Preforce()
{
  resetsiteforce();
  prepareintra();
}










void
RigidBody2::ProgressMomentum( double dt )
{
  Vector3   dw;
  com.ProgressMomentum( dt );
  dw.x = torque.x / in.x;
  dw.y = torque.y / in.y;
  dw.z = torque.z / in.z;
  w.x += dw.x * dt;
  w.y += dw.y * dt;
  w.z += dw.z * dt;
  //cerr << in.x << " " << torque.x << " " << dt << " " << dw.x << "\n";
}



void
RigidBody2::ProgressPosition( double dt )
{
  com.ProgressPosition( dt );
  RotateX( dt / 2 );
  RotateY( dt / 2 );
  RotateZ( dt );
  RotateY( dt / 2 );
  RotateX( dt / 2 );
  //cerr << w.x << " " << rot.v[1][2] << " (2)\n";
}



void
RigidBody2::RotateX( double dt )
{
  double theta = dt * w.x;
  double costh = cos(theta);
  double sinth = sin(theta);
  //rotate angular velocity
  double w0 = costh*w.y - sinth*w.z;
  double w1 = sinth*w.y + costh*w.z;
  w.y = w0;
  w.z = w1;
  //rotate matrix
  w0 = costh*rot.v[0][1] - sinth*rot.v[0][2];
  w1 = sinth*rot.v[0][1] + costh*rot.v[0][2];
  rot.v[0][1] = w0;
  rot.v[0][2] = w1;
  w0 = costh*rot.v[1][1] - sinth*rot.v[1][2];
  w1 = sinth*rot.v[1][1] + costh*rot.v[1][2];
  rot.v[1][1] = w0;
  rot.v[1][2] = w1;
  w0 = costh*rot.v[2][1] - sinth*rot.v[2][2];
  w1 = sinth*rot.v[2][1] + costh*rot.v[2][2];
  rot.v[2][1] = w0;
  rot.v[2][2] = w1;
}



void
RigidBody2::RotateY( double dt )
{
  double theta = dt * w.y;
  double costh = cos(theta);
  double sinth = sin(theta);
  //rotate angular velocity
  double w0 = costh*w.z - sinth*w.x;
  double w1 = sinth*w.z + costh*w.x;
  w.z = w0;
  w.x = w1;
  //rotate matrix
  w0 = costh*rot.v[0][2] - sinth*rot.v[0][0];
  w1 = sinth*rot.v[0][2] + costh*rot.v[0][0];
  rot.v[0][2] = w0;
  rot.v[0][0] = w1;
  w0 = costh*rot.v[1][2] - sinth*rot.v[1][0];
  w1 = sinth*rot.v[1][2] + costh*rot.v[1][0];
  rot.v[1][2] = w0;
  rot.v[1][0] = w1;
  w0 = costh*rot.v[2][2] - sinth*rot.v[2][0];
  w1 = sinth*rot.v[2][2] + costh*rot.v[2][0];
  rot.v[2][2] = w0;
  rot.v[2][0] = w1;
}



void
RigidBody2::RotateZ( double dt )
{
  double theta = dt * w.z;
  double costh = cos(theta);
  double sinth = sin(theta);
  //rotate angular velocity
  double w0 = costh*w.x - sinth*w.y;
  double w1 = sinth*w.x + costh*w.y;
  w.x = w0;
  w.y = w1;
  //rotate matrix
  w0 = costh*rot.v[0][0] - sinth*rot.v[0][1];
  w1 = sinth*rot.v[0][0] + costh*rot.v[0][1];
  rot.v[0][0] = w0;
  rot.v[0][1] = w1;
  w0 = costh*rot.v[1][0] - sinth*rot.v[1][1];
  w1 = sinth*rot.v[1][0] + costh*rot.v[1][1];
  rot.v[1][0] = w0;
  rot.v[1][1] = w1;
  w0 = costh*rot.v[2][0] - sinth*rot.v[2][1];
  w1 = sinth*rot.v[2][0] + costh*rot.v[2][1];
  rot.v[2][0] = w0;
  rot.v[2][1] = w1;
}



void
RigidBody2::WriteWTG6( const Unit& u, ostream& to )
{
  const Vector3& velocity = com.GetVelocity();
  to << com.center.coord.x / ( Angstro * u.m )
     << " " <<com.center.coord.y / ( Angstro * u.m )
     << " " << com.center.coord.z / ( Angstro * u.m )
     << " " << rot.v[0][0]
     << " " << rot.v[0][1]
     << " " << rot.v[0][2]
     << " " << rot.v[1][0]
     << " " << rot.v[1][1]
     << " " << rot.v[1][2]
     << " " << rot.v[2][0]
     << " " << rot.v[2][1]
     << " " << rot.v[2][2]
    //convert from internal unit to A/ps
     << " " << velocity.x / ( Angstro * u.m / (SI::pico * u.sec) )
     << " " << velocity.y / ( Angstro * u.m / (SI::pico * u.sec) )
     << " " << velocity.z / ( Angstro * u.m / (SI::pico * u.sec) )
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
RigidBody2::WriteNX4A( const Unit& u, ostream& to )
{
  to << com.center.coord.x / ( Angstro * u.m )
     << " " <<com.center.coord.y / ( Angstro * u.m )
     << " " << com.center.coord.z / ( Angstro * u.m );
  Quaternion q;
  rot.Quat( q );
  to << " " << q.a
     << " " << q.b
     << " " << q.c
     << " " << q.d << "\n";
}



void
RigidBody2::ReadWTG6( int o, const Unit& u, FILE* file )
{
  char buf[1000];
  order = o;
  
  fgets(buf,sizeof(buf),file);
  //split
  com.center.coord.x = atof(strtok(buf," \t")) * Angstro * u.m;
  com.center.coord.y = atof(strtok(NULL," \t")) * Angstro * u.m;
  com.center.coord.z = atof(strtok(NULL," \t")) * Angstro * u.m;
  rot.v[0][0] = atof(strtok(NULL," \t"));
  rot.v[0][1] = atof(strtok(NULL," \t"));
  rot.v[0][2] = atof(strtok(NULL," \t"));
  rot.v[1][0] = atof(strtok(NULL," \t"));
  rot.v[1][1] = atof(strtok(NULL," \t"));
  rot.v[1][2] = atof(strtok(NULL," \t"));
  rot.v[2][0] = atof(strtok(NULL," \t"));
  rot.v[2][1] = atof(strtok(NULL," \t"));
  rot.v[2][2] = atof(strtok(NULL," \t"));
  Vector3 velocity;
  velocity.x = atof(strtok(NULL," \t")) * ( Angstro * u.m / (SI::pico * u.sec) );
  velocity.y = atof(strtok(NULL," \t")) * ( Angstro * u.m / (SI::pico * u.sec) );
  velocity.z = atof(strtok(NULL," \t")) * ( Angstro * u.m / (SI::pico * u.sec) );
  com.SetVelocity( velocity );
  w.x    = atof(strtok(NULL," \t")) * ( u.radian / (SI::pico * u.sec) );
  w.y    = atof(strtok(NULL," \t")) * ( u.radian / (SI::pico * u.sec) );
  w.z    = atof(strtok(NULL," \t")) * ( u.radian / (SI::pico * u.sec) );
  com.center.force.x    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / (mol * Angstro * u.m ) );
  com.center.force.y    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / (mol * Angstro * u.m ) );
  com.center.force.z    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / (mol * Angstro * u.m ) );
  torque.x    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / mol );
  torque.y    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / mol );
  torque.z    = atof(strtok(NULL," \t")) * ( SI::kilo * u.J / mol );
}



void
RigidBody2::ReadNX4A( int o, const Unit& u, FILE* file )
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
RigidBody2::unserialize( double* const p )
{
  com.unserialize( p );
  q.Set( p[3], p[4], p[5], p[6] );
  q.normalize();
  return p+7;
}



double*
RigidBody2::serialize( double* p ) const
{
  com.serialize( p );
  p[3] = q.a;
  p[4] = q.b;
  p[5] = q.c;
  p[6] = q.d;
  return p+7;
}



double*
RigidBody2::serializeforce( double* xi ) const
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
