#ifndef LOADBALANCER2_HPP
#define LOADBALANCER2_HPP

#include <iostream>
#include <string>
#include <ctime>
#include <boost/utility.hpp>
#include <boost/timer.hpp>
extern "C" {
#include "FakeMPI.h"
}
#include "System/SimpleMD.hpp"
#include "Plugin/LoadBalancer.hpp"

using namespace std;
using namespace boost;


class LoadBalancer2 : public LoadBalancer
{
public:
  explicit LoadBalancer2( sFakeMPI* f, int i, int c, ProcessPlugin* p ) : LoadBalancer( f, i, c, p ), interval(20) {};
  explicit LoadBalancer2( sFakeMPI* f, int i, int c ) : LoadBalancer( f, i, c ), interval(20) {};
  void initialize( SimpleMD* sys ){
    LoadBalancer::initialize( sys );
    interval = 20;
    starttime = time( NULL );
  }
  bool running();
  ~LoadBalancer2()
  {
    cerr << "~Plugin::LoadBalancer2" << endl;
  }
protected:
  int interval;
  int starttime; //absolute unix time from epoch
};

#endif
