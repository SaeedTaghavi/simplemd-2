#include "SingleMol/PolyatomicMol.hpp"
#include "debug.hpp"

PolyatomicMol::~PolyatomicMol()
{
  mesg("~PolyatomicMol\n");
}



PolyatomicMol::PolyatomicMol( int natom, vector<double> mass )
{
  numAtom = natom;
  atom    = vector<MonatomicMol> ( mass.begin(), mass.end() );
  order = -1;
}



void
PolyatomicMol::SetCom()
{
  com.Set(0,0,0);
  for(int i=0; i< numAtom; i++ ){
    com.x += atom[i].center.coord.x;
    com.y += atom[i].center.coord.y;
    com.z += atom[i].center.coord.z;
  }
  com.x /= numAtom;
  com.y /= numAtom;
  com.z /= numAtom;
}



const Vector3&
PolyatomicMol::Position() const
{
  //must be calculated in SetCom()....
  return com;
}



void PolyatomicMol::Translate( const Vector3& offset )
{
  for(int i=0; i< numAtom; i++ ){
    atom[i].Translate( offset );
  }
}



void PolyatomicMol::ScaleVelocity( const Vector3& r )
{
  for( int i=0; i<numAtom; i++ ){
    atom[i].ScaleVelocity( r );
  }
}



void PolyatomicMol::ScaleVelocity( double r )
{
  for( int i=0; i<numAtom; i++ ){
    atom[i].ScaleVelocity( r );
  }
}



void PolyatomicMol::ScalePosition( const Vector3& r )
{
  assert(0);
}



void PolyatomicMol::AddVelocity( const Vector3& v )
{
  for( int i=0; i<numAtom; i++ ){
    atom[i].AddVelocity( v );
  }
}



void PolyatomicMol::Preforce()
{
  for( int i=0; i<numAtom; i++ ){
    atom[i].Preforce();
  }
  //set the center-of-mass
  SetCom();
}



void PolyatomicMol::Postforce( //double vesseld1, double volume, 
                               vector<double> mass )
{
  for( int i=0; i<numAtom; i++ ){
      atom[i].Postforce( //vesseld1, volume, 
                         mass[i] );
  }
}



void PolyatomicMol::ProgressPosition( double dt )
{
  for( int i=0; i<numAtom; i++ ){
    atom[i].ProgressPosition( dt );
  }
}



void PolyatomicMol::ProgressMomentum2( double dt, vector<double> massi )
{
  for( int i=0; i<numAtom; i++ ){
    atom[i].ProgressMomentum( dt );
  }
}



void
PolyatomicMol::BoxCoordinate( const Box& box )
{
  assert(0); //だめ。boxの種類によって処理を変えなければいけないが・・・
  Vector3 pos = Position();
  pos.x = - ImageOffset( pos.x, box.x );
  pos.y = - ImageOffset( pos.y, box.y );
  pos.z = - ImageOffset( pos.z, box.z );
  Translate( pos );
}
  


double PolyatomicMol::GetEk() const
{
  double ek = 0;
  for( int i=0; i<numAtom; i++ ){
    ek += atom[i].GetEk();
  }
  return ek;
}
