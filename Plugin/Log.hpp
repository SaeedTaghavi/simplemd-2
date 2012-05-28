#ifndef LOG_HPP
#define LOG_HPP

#include <iostream>
#include "Plugin/Plugin.hpp"

class SimpleMD;

class Log : public ProcessPlugin
{
  SimpleMD* system;

  ostream* fout;
  ostream* ferr;
  string   tag;

public:
  explicit Log( ostream* fout, ostream* ferr );
  explicit Log();
  void initialize( SimpleMD* sys );
  bool running();
  void terminate();
  ~Log()
  {
    cerr << "~Plugin::Log" << endl;
  }
};


#endif
