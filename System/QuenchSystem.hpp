#ifndef QUENCHSYSTEM_HPP
#define QUENCHSYSTEM_HPP

//using GSL

#include <gsl/gsl_multimin.h>
#include "System/SimpleMD.hpp"

class SimpleQuench : public SimpleMD 
{
private:
  MolCollection2* mol2;
  double coeff;
  gsl_multimin_fdfminimizer *s;
  gsl_vector *x;
  const gsl_multimin_fdfminimizer_type *T;
  int iter;
  gsl_multimin_function_fdf my_func;
  double par[2];

public:
  explicit SimpleQuench(
			MolCollection2* mol,
			Unit* const unit0,
			Truncation rc0,
			map<string, FlexibleHandle>& dict,
			const ProcessPluginArray& pl,
                        double coeff,
			ostream* fout0,
			ostream* ferr0
			);
  MolCollection2* GetMolCollection2(){ return mol2; }
  double GetCoeff(){ return coeff; }
  void initialize();
  bool running();
  void terminate();
};



class SimpleQuench2 : public SimpleMD 
{
private:
  MolCollection2* mol2;
  double coeff;
  int iter;
public:
  explicit SimpleQuench2(
			MolCollection2* mol,
			Unit* const unit0,
			Truncation rc0,
			map<string, FlexibleHandle>& dict,
			const ProcessPluginArray& pl,
                        double coeff,
                        int nstep,
			ostream* fout0,
			ostream* ferr0
			);
  MolCollection2* GetMolCollection2(){ return mol2; }
  double GetCoeff(){ return coeff; }
  void initialize();
  bool running();
  void terminate();
};


#endif
