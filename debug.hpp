#ifndef DEBUG_H
#define DEBUG_H

#define mesg(x)    //fprintf( stdout, x)
#define mesg2(x,y) //fprintf( stdout, x, y )
#ifdef DEBUG
extern int count;
#endif

void report_memory();

#endif
