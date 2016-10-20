#ifndef VECTOR3_HPP
#define VECTOR3_HPP
#include <cmath>
#include <vector>
#include <iostream>
//#include <boost/shared_ptr.hpp>
#include "Unit.hpp"

using namespace std;
//using namespace boost;

/*
 *単純な3次元ベクトル。いずれはvalarrayで記述することになるだろう。
 *いまのところ演算子を定義するつもりはない。
 */



/*Class of 3D vector*/
class Vector3
{
public:
  double x,y,z;
  Vector3( double _x, double _y, double _z );
  Vector3();
  double square() const { return x*x + y*y + z*z; }
  virtual void Set( double _x, double _y, double _z );
  virtual ~Vector3(){};
  virtual void Scale( double ratio );
  virtual void Scale( const Vector3& ratio );
  virtual Vector3* clone() const { return new Vector3(*this); }
  virtual string print() const { return to_string(x) + string(",") + to_string(y) + string(",") + to_string(z); }
};



class Quaternion
{
public:
  double a,b,c,d;
  void normalize()
  {
    double coeff = 1.0 / sqrt( a*a + b*b + c*c + d*d );
    a *= coeff;
    b *= coeff;
    c *= coeff;
    d *= coeff;
  };
  void Set( double aa, double bb, double cc, double dd )
  {
    a = aa;
    b = bb;
    c = cc;
    d = dd;
  }
};



struct SiteOfAction
{
  Vector3 coord;
  Vector3 force;  //also used as acceleration
};



class Box : public Vector3
{
private:
  double volume;
public:
  Box();
  Box(double x, double y, double z);
  void Set( double volume );
  double Volume(){ return volume; }
  void Set( double _x, double _y, double _z );
  virtual void RelativePosition( Vector3 &coord ) const;
  virtual void BoxCoordinate( Vector3 &coord ) const 
    {
      //cerr << "nop\n";
    }; // do nothing, because Box is not a periodic box.
  virtual bool Contains( const Vector3 &coord ) const;
  void Scale( double ratio );
  void Scale( const Vector3& ratio );
  virtual Box* clone() const { return new Box(*this); }


  //I/O
  void WriteBOX3( const Unit& unit, ostream& to ) const;
  void ReadBOX3( const Unit& u, FILE* file );
  void ReadBXLA( const Unit& u, FILE* file );
  //Coding Standard #50: Destructor should be public virtual or protected non-virtual.
  virtual ~Box(){
    //cerr << "~Box" << endl;
  }
};


//typedef std::shared_ptr<Box> BoxHandle;



class SemiperiodicBox : public Box
{
private:
  int isPeriodic[3]; // flags for each axis
public:
  SemiperiodicBox( const int isP[3] );
  void BoxCoordinate( Vector3 &coord ) const;
  bool Contains( const Vector3 &coord ) const;
  void RelativePosition( Vector3 &coord ) const;
  bool IsPeriodic() const{
    return isPeriodic[0] || isPeriodic[1] || isPeriodic[2];
  }
  virtual SemiperiodicBox* clone() const { return new SemiperiodicBox(*this); }
};


typedef std::shared_ptr<Box> BoxHandle;




class Matrix33
{
public:
  double v[3][3];
  Matrix33();
  void Clear();
  void Set( double qa, double qb, double qc, double qd );
  void Quat( Quaternion& q );
};



double applyPBC3(double x,double bxlh);
double applyPBC0(double x,double bxlh);
double applyPBC1(double x,double bxlh);
double applyPBC2(double x,double bxlh);
double BoxCoordinate2(double x,double bxl);
double BoxCoordinate3(double x,double bxl);
double ImageOffset(double x,double bxl);

#endif
