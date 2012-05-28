#ifndef NOSEPOINCAREMD_HPP
#define NOSEPOINCAREMD_HPP


#include "System/SimpleMD.hpp"
class NosePoincare;


class NosePoincareMD : public SimpleMD
{
protected:
  NosePoincare* nosepoincare;
private:
  double Symplectic2nd();
  void write();
public:
  NosePoincareMD(
		  double absoluteTime0,
		  double interval,
		  MolCollection* mol,
		  Unit* const unit0,
		  Truncation rc0,
		  map<string, FlexibleHandle>& dict,
                  const ProcessPluginArray& pl,
		  NosePoincare* nosepoincare0,
		  bool forceloaded,
		  ostream* fout,
		  ostream* ferr
    );

  virtual ~NosePoincareMD();
  double GetEk();
  void GetEkt( Vector3& );
  string log( );
  void terminate();
};



#endif
