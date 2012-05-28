#ifndef TRESET_HPP
#define TRESET_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;

class TReset : public ProcessPlugin
{
  double bathTemp;
  SimpleMD* system;
public:
  TReset( double temp );
  void initialize( SimpleMD* sys );
  bool running();
  void terminate();
  ~TReset()
  {
    cerr << "~Plugin::TReset" << endl;
  }
};


#endif
