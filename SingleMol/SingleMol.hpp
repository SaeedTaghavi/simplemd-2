#ifndef SINGLEMOL_HPP
#define SINGLEMOL_HPP

#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "Vector3.hpp"
#include <cstring>

using namespace std;
using namespace boost;

/*abstract class*/
class SingleMolEntity
{
public:
  SingleMolEntity(){};
  virtual ~SingleMolEntity();
  virtual const Vector3& Position() const = 0;
  virtual void Preforce() = 0;

  //"Get" functons
  virtual int  GetOrder() const = 0;
  virtual double GetEk() const = 0;

  virtual void BoxCoordinate( const Box& box ) = 0;
  virtual void ProgressPosition( double ) = 0;

  //features for plugin
  virtual void AddVelocity( const Vector3& velocity ) = 0;
  virtual void Translate( const Vector3& ) = 0;
  virtual void ScaleVelocity( const Vector3& r ) = 0;
  virtual void ScaleVelocity( double r ) = 0;
  virtual void ScalePosition( const Vector3& r ) = 0;
  
  //Serialize/Unserialize configuration.
  virtual double* unserialize( double* const p ) = 0;
  virtual double* serialize( double* const p ) const = 0;
  virtual double* serializeforce( double* const xi ) const = 0;
private:
  //SingleMolEntity( const SingleMolEntity & );
  //SingleMolEntity& operator=( const SingleMolEntity & );
};



typedef shared_ptr<SingleMolEntity> SingleMolHandle;
typedef vector<SingleMolHandle>     SingleMolArray;


#endif
