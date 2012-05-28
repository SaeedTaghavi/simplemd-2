#ifndef MOLPROPERTY_H
#define MOLPROPERTY_H
#include <vector>
#include <string>
#include "Interaction/Interaction.hpp"
#include "Vector3.hpp"

using namespace std;
using namespace boost;

/*Common molecular properties*/
class MolProperty
{
public:
  string id08;

  MolProperty(string id): id08(id) {}
  virtual ~MolProperty(){};

  //"Get" functions
  virtual int GetNumSite() const { return 1; };
  virtual int GetDoF() = 0;
  virtual double GetMass() const = 0;

protected:
private:
  //disable copy constructors
  MolProperty( const MolProperty & );
  MolProperty& operator=( const MolProperty & );
};


/*
 *単原子分子もこれを流用する。
 */
class Flexible : public MolProperty
{
public:
  Flexible( unsigned int nsite, const IntrParamsArray&, int DoF );
  Flexible( unsigned int nsite, const IntrParamsArray& i, vector<double>& mass_, vector<string>& name, string id, int DoF );
  Flexible();
  ~Flexible();
  //int  total_size() const;
  //void SetSize( int nmol ){ totalNum = nmol; }
  int GetNumSite() const { return numSite; }
  int GetDoF(){ return dof; }
  const IntrParamsArray& GetIntrParams() const;
  const vector<double>&   GetMassArray() const;
  double GetMass() const;
  const vector<double>&   GetMassI() const;
  string GetAtomName( int i ) const { return atomName[i]; }
  void ReadMonatomic( FILE* file );
  virtual void Write( const Unit& unit, ostream& to ) const;
protected:
  //number of sites in a molecules
  int numSite;
  IntrParamsArray  intr;  // IntrParams( nsite )
  vector<double>    mass;          // Mass ( nsite )
  vector<double>    massi;         // Mass ( nsite )
  void WriteDEFP(  const Unit& unit, ostream& to ) const;
  vector<string> atomName;
  //name of the molecule in 8 letters
  int dof;
private:
  //disable copy constructors
  Flexible( const Flexible & );
  Flexible& operator=( const Flexible & );
};



class Rigid : public Flexible
{
public:
  vector<Vector3> site;

  Rigid( int nsite, const IntrParamsArray&, vector<double>& mass, vector<Vector3>& site, vector<string>& atomName, string id, int DoF );
  Rigid();
  ~Rigid();
  //"Get" functions
  double GetMass() const;
  const IntrParamsArray& GetIntrParams() const;
  const Vector3& GetInertia() const;
  const Vector3& GetSite( int i ) const { return site[i]; }

  //"Set" functions
  void SetInertia( const Vector3& v );

  void ReadMonatomic( FILE* file );
  //coordinates of the sites in intramolecular coordinate.
  virtual void Write( const Unit& unit, ostream& to ) const;
private:
  void WriteDEFR( const Unit& unit, ostream& to ) const;
  //disable copy constructors
  Rigid( const Rigid & );
  Rigid& operator=( const Rigid & );
  //inertia
  Vector3 inertia;
  double total_mass;
};



typedef shared_ptr<MolProperty> MolPropertyHandle;
typedef shared_ptr<Flexible> FlexibleHandle;
typedef shared_ptr<Rigid> RigidHandle;
typedef vector<MolPropertyHandle> MolPropertyArray;
typedef vector<FlexibleHandle> FlexibleArray;

#endif
