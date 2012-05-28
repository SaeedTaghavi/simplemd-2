#include <iostream>
#include "Plugin/LoadBalancer.hpp"

void
LoadBalancer::initialize( SimpleMD* sys )
{
  plugin->initialize( sys );

  nextloop = innerloop;
  last = 0;
  t.restart();
  step= 0;
}



bool
LoadBalancer::running()
{
  if ( step++ == last + nextloop ){
    last += nextloop;
    float dt = t.elapsed();
    cerr << dt << " elapsed\n";
    t.restart();
    int processloop[fmpi->nprocs];
    float processtime[fmpi->nprocs];
    fakempi_gather( fmpi, 4, &nextloop, (char*)processloop );
    fakempi_gather( fmpi, 4, &dt,       (char*)processtime );
    if ( fmpi->myrank == 0 ){
      cerr << "Loadbalancing...";
      float avgtime = 0;
      for( int i=0; i<fmpi->nprocs; i++ ){
	cerr << i 
	     << ":" << processloop[i]
	     << "," << processtime[i] << endl;
	avgtime += innerloop * processtime[i] / processloop[i];
      }
      avgtime /= fmpi->nprocs;
      cerr << avgtime << " avgtime\n";
      for( int i=0; i<fmpi->nprocs; i++ ){
	dt = processtime[i];
	if ( 100000 < dt ) dt = 100000;
	processloop[i] = avgtime * processloop[i] / processtime[i];
      }
      cerr << "Done.\n";
    }
    fakempi_scatter( fmpi, 4, &nextloop, (char*)processloop );
    
    dt = t.elapsed();
    cerr << dt << " seconds for load balancing.\n";
    cerr << nextloop << " new interval.\n";

    //execute the plugin
    bool result = plugin->running();
    
    //restart the interval timer
    t.restart();

    //decrement the loop counter and terminate when zero
    count --;
    result = (0 < count) && result;
    return result;
  }
  return true;
}



void
LoadBalancer::terminate()
{
  plugin->terminate();
}
