#ifndef RIGIDBODY2_HPP
#define RIGIDBODY2_HPP
#include "SingleMol/RigidBody.hpp"

/*
 *RigidBodyではquaternionを使っていたが、RigidBody2では回転行列を直接時間発展させる。
 *
 *DLM97の方法。
 */



/*Position and force of a polyatomic molecule*/
class RigidBody2 : public RigidBody
{
public:
  RigidBody2( int nsite, const MolPropertyHandle& p );
  ~RigidBody2();
  void Preforce();
  void ProgressPosition( double dt );

  //local extension
  void ProgressMomentum( double dt );

  //I/O
  void WriteWTG6( const Unit& unit, ostream& to );
  void WriteNX4A( const Unit& unit, ostream& to );
  void ReadWTG6( int o, const Unit& u, FILE* file );
  void ReadNX4A( int o, const Unit& u, FILE* file );

  //for sorting
  int operator<( const RigidBody2& that ) const { return (order < that.order);}

  //Serialize/Unserialize configuration.
  double* unserialize( double* const p );
  double* serialize( double* const p ) const;
  double* serializeforce( double* const xi ) const;

private:
  void RotateX( double dt );
  void RotateY( double dt );
  void RotateZ( double dt );

};


#endif
