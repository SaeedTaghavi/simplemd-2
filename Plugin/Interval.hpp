#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;
//a plugin to execute a plugin periodically.

class Interval : public ProcessPlugin
{
  SimpleMD* system;
  int interval;
  ProcessPlugin* plugin;
public:
  explicit Interval( int intv, ProcessPlugin* p );
  void initialize( SimpleMD* sys );
  bool running();
  void terminate();
  ~Interval()
  {
    delete plugin;
    cerr << "~Plugin::Interval" << endl;
  }
};


#endif
