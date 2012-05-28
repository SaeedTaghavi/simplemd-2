#ifndef STOPDRIFT_HPP
#define STOPDRIFT_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;

class StopDrift : public ProcessPlugin
{
  SimpleMD* system;
public:
  explicit StopDrift();
  void initialize( SimpleMD* sys );
  bool running();
  ~StopDrift()
  {
    cerr << "~Plugin::StopDrift" << endl;
  }
};


#endif
