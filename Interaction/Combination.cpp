#include "Interaction/Combination.hpp"




void
Intersite::Set( int s1, const IntrParams& i1,
			int s2, const IntrParams& i2 )

{
  LJC lb = LJC( i1, i2 );
  int failed = 0;
  if ( lb.coulomb.sqrch == 0 ){
    if ( lb.lj.eps == 0 ){
      failed = 1;
    }
    else{
      intr.push_back( InteractionHandle( new LJ( lb.lj.eps, lb.lj.sig ) ) );
      //fprintf( stderr, "LJ\n" );
    }
  }
  else{
    if ( lb.lj.eps == 0 ){
      intr.push_back( InteractionHandle( new Coulomb( lb.coulomb.sqrch ) ) );
      //fprintf( stderr, "C\n" );
    }
    else{
      intr.push_back( InteractionHandle( new LJC( lb.lj.eps, lb.lj.sig, lb.coulomb.sqrch ) ) );
      //fprintf( stderr, "LJC\n" );
    }
  }
  if ( ! failed ){
    site1.push_back( s1 );
    site2.push_back( s2 );          
  }
  else{
    //fprintf( stderr, "None\n" );
  }
}


void
Intersite::Set( const FlexibleHandle& p1,
			const FlexibleHandle& p2 )
{
  const IntrParamsArray& ia1 = p1->GetIntrParams();
  const IntrParamsArray& ia2 = p2->GetIntrParams();
  nsite1 = p1->GetNumSite();
  nsite2 = p2->GetNumSite();
  for( int s1=0; s1<nsite1; s1++ ){
    for( int s2=0; s2<nsite2; s2++ ){
      //fprintf( stderr, "%d-%d:", s1, s2 );
      Set( s1, ia1[s1], s2, ia2[s2] );
    }
  }
}



int
Combination::push_back( const FlexibleHandle& ph )
{
  prop.push_back( ph );
  int ncompo = prop.size();
  int c1 = ncompo - 1;
  for( int c2=0; c2<ncompo; c2++ ){
    Intersite* im = new Intersite();
    im->Set( ph, prop[c2] );
    /*
     *Component c1 and c2 interact with parameters im.
     */
    compo1.push_back( c1 );
    compo2.push_back( c2 );
    ima.push_back( IntersiteHandle( im ) );
  }

  return ncompo;
}
    
