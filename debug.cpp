#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <cstdio>
#include <iostream>

using namespace std;

void report_memory()
{
  struct rusage usage;
  int result = getrusage( RUSAGE_SELF, &usage );
  /*
  printf("%d result\n%ld shm\n%ld data\n%ld stack\n",
	 result,
	 usage.ru_ixrss,
	 usage.ru_idrss,
	 usage.ru_isrss );
  */
  cout << result          << "result" << endl
       << usage.ru_maxrss << "max"    << endl     
       << usage.ru_ixrss  << "shm"    << endl     
       << usage.ru_idrss  << "data"   << endl     
       << usage.ru_isrss  << "stack"  << endl;
}
