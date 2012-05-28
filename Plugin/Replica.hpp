#ifndef REPLICA_HPP
#define REPLICA_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"
#include "Exchange.hpp"
#include "Random.hpp"
#include "Integrator/NosePoincare.hpp"
extern "C" {
#include "FakeMPI.h"
}

class SimpleMD;

class Replica : public ProcessPlugin
{
  SimpleMD* system;

  sFakeMPI* fmpi;
  ostream* fout;
  ostream* ferr;
  int ntrial;
  bool alter;
  NosePoincare* nosepoincare;
  vector<int> node_in_charge;
  vector<Exchange> replica;
  RandomD01Handle random;
  Unit* unit;

public:
  explicit Replica( sFakeMPI* f, int seed, int ntrial,  Unit* const u, NosePoincare* np );
  void initialize( SimpleMD* sys );
  bool running();
  ~Replica()
  {
    cerr << "~Plugin::Replica" << endl;
  }
};


#endif
