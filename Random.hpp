#ifndef _RANDOM_HPP
#define _RANDOM_HPP
#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace boost;

class RandomD01
{
public:
  RandomD01( uint32_t seed )
    : rand( lagged_fibonacci607(seed), uniform_real<>(0,1) ){}
  double get(){ return rand(); }
private:
  variate_generator< lagged_fibonacci607, uniform_real<> > rand;
};

typedef shared_ptr<RandomD01> RandomD01Handle;

#endif
