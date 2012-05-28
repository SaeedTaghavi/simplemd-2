#ifndef EXCHANGE_HPP
#define EXCHANGE_HPP
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include "Random.hpp"


using namespace boost;

class ExchangeDriver;
typedef shared_ptr<ExchangeDriver> ExchangeHandle;

class ExchangeDriver
{
};


class Exchange : ExchangeDriver
{
protected:
  RandomD01Handle random;
  //[0]:kT, [1]:Energy(Enthalpy), [2..]:aux
  vector<double> consts;
  vector<double> vars;

public:
  Exchange( const RandomD01Handle& r );
  bool Trial( Exchange& target );
  vector<double>& GetVars();
  const vector<double>& GetConsts() const;
  void Set( vector<double>& c, vector<double>& v );
  void exchange( Exchange& target );
};


typedef shared_ptr<ExchangeDriver> ExchangeHandle;




#endif
