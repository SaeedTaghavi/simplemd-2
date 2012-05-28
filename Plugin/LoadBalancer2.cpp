#include <iostream>
#include "Plugin/LoadBalancer2.hpp"

bool
LoadBalancer2::running()
{
  step ++;

  int  t1 = starttime + interval;
  int  t2 = time(NULL);
  //bool finished = ( starttime + interval < time( NULL ) );
  bool finished = ( t1 < t2 );
  if ( last == 0 ){
    //初回のみ、指定された回数ループする。
    finished = ( nextloop < step );
  }
  if ( finished ){
    nextloop = step - last;
    last = step;
    //float dt = t.elapsed();
    int tt1 = t2;
    int tt2 = starttime;
    int intr = tt1 - tt2;
    fprintf(stderr, "%d %d %d intr\n", tt1, tt2, intr );
    cerr << t2 << "-" << starttime << ":" << intr << " elapsed " << nextloop << " loops\n";
    t.restart();
    int processloop[fmpi->nprocs];
    int processtime[fmpi->nprocs];
    fakempi_gather( fmpi, 4, &intr,       (char*)processtime );
    fakempi_gather( fmpi, 4, &nextloop, (char*)processloop );
    if ( fmpi->myrank == 0 ){
      cerr << "Loadbalancing...";
      //1ステップに必要な平均時間を求める。
      int totalloop = 0;
      int totaltime = 0;
      for( int i=0; i<fmpi->nprocs; i++ ){
	cerr << i 
	     << ":" << processloop[i]
	     << "," << processtime[i] << endl;
        totalloop += processloop[i];
        totaltime += processtime[i];
      }
      //innerloop回の計算に必要な、平均時間を求める。
      int avgtime = innerloop * totaltime / totalloop;

      //あまり極端に変化するのも困るので、feedbackで徐々に追従させる。
      int delta = (totaltime / fmpi->nprocs - avgtime) / 10;
      if ( delta < -60 ) delta = -60;
      if ( 60 < delta ) delta = 60;
      interval = totaltime / fmpi->nprocs - delta;

      cerr << avgtime << " avgtime " << interval << " next interval\n";

      for( int i=0; i<fmpi->nprocs; i++ ){
	processtime[i] = interval;
      }
      cerr << "Done.\n";
    }
    fakempi_scatter( fmpi, 4, &interval, (char*)processtime );
    if ( fmpi->myrank == 2 ){
      cerr << fmpi->myrank << ":" << interval << " next interval\n";
    }
    
    nextloop = 0;

    //execute the plugin
    bool result = plugin->running();
    
    //restart the interval timer
    t.restart();
    //どうもelapsedは実時間じゃないみたい。精度が悪くてもいいからtime()を使う。
    starttime = time( NULL );

    //decrement the loop counter and terminate when zero
    //outerloop counter
    count --;
    result = (0 < count) && result;
    return result;
  }
  return true;
}
