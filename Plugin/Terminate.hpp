#ifndef TERMINATE_HPP
#define TERMINATE_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class Terminate : public ProcessPlugin
{
  int steps;
  SimpleMD* system;
public:
  explicit Terminate( int s );
  void initialize( SimpleMD* sys );
  bool running();
  ~Terminate()
    {
      cerr << "~Plugin::Terminate" << endl;
    }
};


#endif
