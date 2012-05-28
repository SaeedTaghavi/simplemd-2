#include "Interaction/PotVir.hpp"

PotVir::PotVir()
{
  //ep = 0;
  //vr.Clear();
}



PotVir::PotVir( double ep0, const Matrix33& vr0 ) : ep(ep0), vr(vr0)
{
}




void
PotVir::Add( const PotVir& pv )
{
  ep += pv.ep;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      vr.v[i][j] += pv.vr.v[i][j];
    }
  }
}



void PotVir::Clear()
{
  ep = 0;
  vr.Clear();
}



void PotVir::Set( double ep0, const Matrix33& vr0 )
{
  ep = ep0;
  vr = vr0;
}

  

double
PotVir::Tr() const
{
  return vr.v[0][0] + vr.v[1][1] + vr.v[2][2];
}
