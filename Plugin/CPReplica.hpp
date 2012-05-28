#ifndef CPREPLICA_HPP
#define CPREPLICA_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"
#include "ExchangePT.hpp"
#include "Random.hpp"
#include "Integrator/NosePoincareAndersen.hpp"
extern "C" {
#include "FakeMPI.h"
}

class SimpleMD;

class CPReplica : public ProcessPlugin
{
  SimpleMD* system;

  sFakeMPI* fmpi;
  ostream* fout;
  ostream* ferr;
  int ntrial;
  bool alter;
  NosePoincareAndersen* nosepoincareandersen;
  vector<int> node_in_charge;
  vector<ExchangePT> replica;
  RandomD01Handle random;
  Unit* unit;

public:
  explicit CPReplica( sFakeMPI* f, int seed, int ntrial, Unit* const u, NosePoincareAndersen* npa );
  void initialize( SimpleMD* sys );
  bool running();
  ~CPReplica()
  {
    cerr << "~Plugin::CPReplica" << endl;
  }
};


#endif
