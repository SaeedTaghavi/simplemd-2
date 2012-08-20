#ifndef SIMPLEMD_GUARD
#define SIMPLEMD_GUARD

#include <map>

#include "System/System.hpp"
#include "Unit.hpp"
#include "MolCollection.hpp"
#include "Vector3.hpp"
#include "Interaction/Truncation.hpp"
#include "Interaction/PotVir.hpp"
#include "Plugin/Plugin.hpp"

/*
 *Integratorの種類にあわせ、MDの進め方が違ってくる。その違いを吸収する
 *ためのクラス。
 *
 *MDのメインプログラムからは、ここのクラスのひとつを選択し、
 *初期化フェーズではinitialize()を、
 *実行フェーズではrunnungを、
 *終了フェーズではterminateを呼びだす。
 *
 *MDの過程で、他の処理をやらせたければ、これらのクラスを継承するか、包
 *含した新クラスを作ればよい。
 */



class SimpleMD : public VirtualMD
{
  //friend class ProcessPlugin;
  friend double func( double* p );
private:
protected:
  double dt;
  int step;
  double absoluteTime;
  map<string, FlexibleHandle> moldict;
  // Molecule container ( separated from the vessel )
  MolCollection* molecules;
  Truncation rc;
  Unit* unit;
  int dof;
  // for file I/O
  FILE*  input;
  ProcessPluginArray plugins;


  //Write final data for continuation
  virtual void write();

  //Calculate pressure ( in internal pressure unit )
  virtual double Pressure( double idealpv );

  //Simplest integrator
  virtual double Symplectic2nd();

  ostream* fout;
  ostream* ferr;
  PotVir pv;
public:

  explicit SimpleMD(
                     double absoluteTime0,
                     double interval,
                     MolCollection* mol,
                     Unit* const unit0,
                     Truncation rc0,
                     map<string, FlexibleHandle>& dict,
                     const ProcessPluginArray& pl,
		     ostream* fout,
		     ostream* ferr
    );
  virtual ~SimpleMD();
  virtual void initialize();
  virtual double OneStep();
  virtual bool running();
  virtual void terminate();
  virtual double GetEk();
  virtual void GetEkt( Vector3& );
  // Calculate force
  virtual void Force();

  virtual int GetStep() const
    {
      return step;
    }
  virtual const Unit& GetUnit() const
    {
      return *unit;
    }
  virtual const PotVir& GetPv()
    {
      return pv;
    }
  virtual MolCollection* GetMolCollection() const
    {
      return molecules;
    }
  //Write snapshot (coord only)
  void SnapShot( ostream& out=cerr, int mode=0 );

  //One-line log for monitor.
  virtual string log();
  
  //Calculate pressure ( in internal energy unit )
  virtual double Temperature( double ek );

  //obtain protected variables
  virtual ostream* GetFout()
  {
    return fout;
  }
  virtual ostream* GetFerr()
  {
    return ferr;
  }

};



typedef shared_ptr<SimpleMD> SimpleMDHandle;



class Simple4th : public SimpleMD
{
private:
  double Symplectic4th();
public:
  Simple4th(
             double absoluteTime0,
             double interval,
             MolCollection* mol,
             Unit* const unit0,
             Truncation rc0,
             map<string, FlexibleHandle>& dict,
             const ProcessPluginArray& pl,
	     ostream* fout,
	     ostream* ferr
    );
  double OneStep();
};












#endif
