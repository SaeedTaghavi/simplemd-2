#ifndef POLYATOMICMOL_HPP
#define POLYATOMICMOL_HPP
#include "SingleMol/MonatomicMol.hpp"

//多原子分子の一般的な型。この型は実体化させない。
class PolyatomicMol : public SingleMolEntity
{
private:
  //center-of-mass, used only in Position() and SetCom();
  Vector3 com;
public:
  //in PolyatomicMol, com is in absolute coordinate
  //in RigidBody, com is in coordinate relative to center of mass of the mol.
  vector<MonatomicMol> atom;
  int                   order;
  int                   numAtom;


  PolyatomicMol( int natom, vector<double> mass );

  const Vector3& Position() const;
  //void    Report();

  //"Get" functions
  double GetEk( ) const;
  int GetOrder() const {return order;}

  void Predict();
  void Preforce();
  void Postforce( //double vesseld1, double volume, 
                  vector<double> mass );
  void Correct();
  void BoxCoordinate( const Box& box );
  void ProgressPosition( double dt );
  void ProgressMomentum2( double dt, vector<double> massi );


  //features for plugin
  void ScaleVelocity( const Vector3& r );
  void ScaleVelocity( double r );
  void ScalePosition( const Vector3& r );
  void Translate( const Vector3& );
  void AddVelocity( const Vector3& velocity );

  //Serialize/Unserialize configuration.
  double* unserialize( double* const p ){ return 0; };
  double* serialize( double* const p ) const { return 0; };
  double* serializeforce( double* const xi ) const { return 0; };

private:
  void SetCom();
  //PolyatomicMol( const PolyatomicMol & );
  //PolyatomicMol& operator=( const PolyatomicMol & );

protected:
  ~PolyatomicMol();
};



//柔軟な分子の型。
//class FlexibleMol: public PolyatomicMol
//{
//public:
//protected:
//};


#endif
