#include "QuenchSystem.hpp"
#include "Plugin/CollectorPlugin.hpp"

//with Numerical Recipe.

SimpleQuench2::SimpleQuench2(
			   MolCollection2* mol,
			   Unit* const unit0,
			   Truncation rc0,
			   map<string, FlexibleHandle>& dict,
			   const ProcessPluginArray& pl,
                           double c,
                           int nstep,
			   ostream* fout0,
			   ostream* ferr0
  ) : SimpleMD( 0, 0, mol, unit0, rc0, dict, pl, fout0, ferr0 ), mol2( mol ), coeff( c ), iter(nstep)
{
}







static SimpleQuench2* sq;


double
my_f (double* vm)
{
  double* v = vm + 1;
  MolCollection2* mol2 = sq->GetMolCollection2();
  mol2->unserialize( v );
  //unserializeすると規格化されるので、それを呼び戻すことでvを規格化する。
  mol2->serialize( v );
  //size_t d = mol2->qdof();
  QDoFPlugin qdof;
  qdof.Execute0( mol2 );
  size_t d = qdof.Result();
  //for(int i=0;i<d;i++)
  //  printf( "%d %f v\n", i, v->data[i] );
  //
  sq->Force();
  //
  const PotVir& pv = sq->GetPv();
  double ep = pv.GetEp();
  //cerr << ep << " ep" << endl;
  return ep;
}



/* The gradient of f, df = (df/dx, df/dy). */
void
my_df (double* vm, double* dfm)
{
  double* v = vm + 1;
  double* df = dfm + 1;
  //you must set coord and calc force 
  MolCollection2* mol2 = sq->GetMolCollection2();

  //Recipeを使う場合は、dfuncは必ずfuncの直後に呼ばれるので、力を再計算する必要はない。
  //mol2->unserialize( v );
  //
  //sq->Force();
  mol2->serializeforce( df );
  //size_t d = mol2->qdof();
  QDoFPlugin qdof;
  qdof.Execute0( mol2 );
  size_t d = qdof.Result();
  //gradientの方向はforceの逆
  double coeff = sq->GetCoeff();
  for(int i=0;i<d;i++)
    df[i] *= - coeff;
  //for(int i=1;i<=d;i+=i)
  //  printf( "%d %f df %f v\n", i-1, dfm[i], vm[i] );
}






extern "C" {
  void frprmn(double p[], int n, double ftol, int *iter, double *fret,
              double (*func)(double []), void (*dfunc)(double [], double []));
}




void
SimpleQuench2::initialize()
{
  sq = this;

  MolCollection2* mol2 = sq->GetMolCollection2();
  //size_t d = 2; //mol2->qdof();
  //size_t d = mol2->qdof();
  QDoFPlugin qdof;
  qdof.Execute0( mol2 );
  size_t d = qdof.Result();
  cerr << d << " dof\n";

  //set coordinate to x
  double x[d];
  mol2->serialize( x );
  double* xm = x;
  xm --;
  double ftol = 1e-8;
  //int iter = 1000;
  double fret;
  frprmn( xm, d, ftol, &iter, &fret, &my_f, &my_df );
}



bool
SimpleQuench2::running()
{
  return false;
}



void
SimpleQuench2::terminate()
{
  *fout << "@ETOT" << endl 
//        << sq->pv.GetEp() / ( kilo * unit->J / mol ) << endl;
        << sq->pv.GetEp() / ( unit->J / mol / Kelvin ) << endl;
  SnapShot( *fout );
}
