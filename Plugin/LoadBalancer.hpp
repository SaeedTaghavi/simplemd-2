#ifndef LOADBALANCER_HPP
#define LOADBALANCER_HPP

#include <iostream>
#include <string>
#include <boost/utility.hpp>
#include <boost/timer.hpp>
extern "C" {
#include "FakeMPI.h"
}
#include "System/SimpleMD.hpp"

using namespace std;
using namespace boost;


class LoadBalancer : public ProcessPlugin
{
protected:
  sFakeMPI* fmpi;
  int innerloop;
  int count;
  ProcessPlugin* plugin;

  int nextloop;
  int last;
  timer t;
  int step;
public:
  explicit LoadBalancer( sFakeMPI* f, int i, int c, ProcessPlugin* p=0 ) : fmpi(f), innerloop(i), count(c), plugin( p ) {};
  void set( ProcessPlugin* p )
  {
    plugin = p;
  }
  void initialize( SimpleMD* sys );
  bool running();
  void terminate();
  ~LoadBalancer()
  {
    delete plugin;
    cerr << "~Plugin::LoadBalancer" << endl;
  }
};

#endif
