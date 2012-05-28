#ifndef NOSEPOINCAREANDERSENMD_HPP
#define NOSEPOINCAREANDERSENMD_HPP

#include "System/SimpleMD.hpp"
class NosePoincareAndersen;


class NosePoincareAndersenMD : public SimpleMD
{
protected:
  NosePoincareAndersen* nosepoincareandersen;
  virtual double Symplectic2nd();
  void write();
public:
  NosePoincareAndersenMD(
		  double absoluteTime0,
		  double interval,
		  MolCollection* mol,
		  Unit* const unit0,
		  Truncation rc0,
		  map<string, FlexibleHandle>& dict,
                          const ProcessPluginArray& pl,
		  NosePoincareAndersen* nosepoincareandersen0,
			  bool forceloaded,
			  ostream* fout,
			  ostream* ferr
    );
  virtual ~NosePoincareAndersenMD();
  double GetEk();
  void GetEkt( Vector3& );
  string log();
  void terminate();
};



/*expands only in z-axis.*/
class NosePoincareAndersenZMD : public NosePoincareAndersenMD
{
  double Symplectic2nd();
  void write();
public:
  NosePoincareAndersenZMD(
			  double absoluteTime0,
			  double interval,
			  MolCollection* mol,
			  Unit* const unit0,
			  Truncation rc0,
			  map<string, FlexibleHandle>& dict,
			  const ProcessPluginArray& pl,
			  NosePoincareAndersen* nosepoincareandersen0,
			  bool forceloaded,
			  ostream* fout,
			  ostream* ferr
    );
};



#endif
