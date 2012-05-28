#ifndef PLUGIN_HPP
#define PLUGIN_HPP
/*
 * Sample implementation of the "Plug-in"
 *
 * Plug-ins are the derived class of Plugin == Main
 */
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
//#include "System.hpp"

using namespace std;
using namespace boost;

class SimpleMD;

class ProcessPlugin : noncopyable
{
public:
  virtual ~ProcessPlugin(){}

  // plugin may be constructed before the System is initialized,
  // so the pointer to the System is given here.
  // This pointer MUST NOT be freed.
  virtual void initialize( SimpleMD* sys ) = 0;

  //running() is called just after the main process of the System.
  //its order follows the order in the array.
  virtual bool running()
  {
    return true;
  }

  // ForceHook() is called inside SimpleMD::Force() 
  virtual void ForceHook()
  {
  }

  virtual void terminate()
  {
  }

  virtual void setpriority( int p )
  {
    priority = p;
  }
  int operator<( const ProcessPlugin& that ) const
  { 
    return (priority < that.priority);
  }
private:
  int priority;
};

typedef shared_ptr<ProcessPlugin> ProcessPluginHandle;
typedef vector<ProcessPluginHandle> ProcessPluginArray;

#endif
