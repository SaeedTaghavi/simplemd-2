#include "QuenchSystem.hpp"
#include "Plugin/CollectorPlugin.hpp"

// with GSL

SimpleQuench::SimpleQuench(
			   MolCollection2* mol,
			   Unit* const unit0,
			   Truncation rc0,
			   map<string, FlexibleHandle>& dict,
			   const ProcessPluginArray& pl,
                           double c,
			   ostream* fout0,
			   ostream* ferr0
  ) : SimpleMD( 0, 0, mol, unit0, rc0, dict, pl, fout0, ferr0 ), mol2( mol ), coeff(c)
{
}







static SimpleQuench* sq;


double
my_f (const gsl_vector *v, void *params)
{
  MolCollection2* mol2 = sq->GetMolCollection2();
  mol2->unserialize( v->data );
  //MolCollection* m = mol2;
  QDoFPlugin qdof;
  qdof.HookL3( mol2 );
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
my_df (const gsl_vector *v, void *params,
       gsl_vector *df)
{
  //you must set coord and calc force 
  MolCollection2* mol2 = sq->GetMolCollection2();
  mol2->unserialize( v->data );
  //
  sq->Force();
  mol2->serializeforce( df->data );
  //size_t d = mol2->qdof();
  QDoFPlugin qdof;
  qdof.HookL3( mol2 );
  size_t d = qdof.Result();
  double coeff =sq->GetCoeff();
  //gradientの方向はforceの逆
  for(int i=0;i<d;i++)
    df->data[i] *= - coeff;
  for(int i=1;i<=d;i+=i)
    printf( "%d %f df %f v\n", i-1, df->data[i-1], v->data[i-1] );
}



void
my_fdf (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
  *f = my_f(x, params);
  MolCollection2* mol2 = sq->GetMolCollection2();
  mol2->serializeforce( df->data );
  //size_t d = mol2->qdof();
  QDoFPlugin qdof;
  qdof.HookL3( mol2 );
  size_t d = qdof.Result();
  double coeff =sq->GetCoeff();
  //gradientの方向はforceの逆
  for(int i=0;i<d;i++)
    df->data[i] *= - coeff;
  for(int i=1;i<=d;i+=i)
    printf( "%d %f df %f x\n", i-1, df->data[i-1], x->data[i-1] );
  //cout << *f << " ep" << endl;
}



extern "C" {
  void frprmn(double p[], int n, double ftol, int *iter, double *fret,
              double (*func)(double []), void (*dfunc)(double [], double []));
}





double
sample_f (const gsl_vector *v, void *params)
{
  double x, y;
  double *dp = (double *)params;
  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
  
  return 10.0 * (x - dp[0]) * (x - dp[0]) +
    20.0 * (y - dp[1]) * (y - dp[1]) + 30.0;
}
     
/* The gradient of f, df = (df/dx, df/dy). */
void
sample_df (const gsl_vector *v, void *params,
       gsl_vector *df)
{
  double x, y;
  double *dp = (double *)params;
     
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
     
  gsl_vector_set(df, 0, 20.0 * (x - dp[0]));
  gsl_vector_set(df, 1, 40.0 * (y - dp[1]));
}
     
/* Compute both f and df together. */
void
sample_fdf (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
  *f = sample_f(x, params);
  sample_df(x, params, df);
  cout << *f << " sample_fdf\n";
}




void
SimpleQuench::initialize()
{
  sq = this;

  MolCollection2* mol2 = sq->GetMolCollection2();
  //size_t d = 2; //mol2->qdof();
  //size_t d = mol2->qdof();
  QDoFPlugin qdof;
  qdof.HookL3( mol2 );
  size_t d = qdof.Result();
  cout << d << " dof\n";

  par[0] = 1;
  par[1] = 2;
  //my_func.f = &my_f;
  //my_func.df = &my_df;
  //my_func.f = &sample_f;
  //my_func.df = &sample_df;
  //my_func.fdf = &sample_fdf;
  my_func.f = &my_f;
  my_func.df = &my_df;
  my_func.fdf = &my_fdf;
  my_func.n = d;
  my_func.params = &par;

  //set coordinate to x
  x = gsl_vector_alloc (d);
  mol2->serialize( x->data );
  //gsl_vector_set (x, 0, 5.0);
  //gsl_vector_set (x, 1, 7.0);

  T = gsl_multimin_fdfminimizer_conjugate_pr;
  //T = gsl_multimin_fdfminimizer_conjugate_fr;
  //T = gsl_multimin_fdfminimizer_steepest_descent;
  //T = gsl_multimin_fdfminimizer_vector_bfgs;
  s = gsl_multimin_fdfminimizer_alloc (T, d);

  double step_size = 1;
  double tol = 1e-8;
  gsl_multimin_fdfminimizer_set (s, &my_func, x, step_size, tol );

  iter = 0;
}



bool
SimpleQuench::running()
{
  iter ++;
  int status = gsl_multimin_fdfminimizer_iterate (s);
  //cout << status << " status1\n";

  /*
  if ( iter==300 ){
    //途中で切り替えてもうまくいかない。
    size_t d = mol2->qdof();
    MolCollection2* mol2 = sq->GetMolCollection2();
    mol2->unserialize( s->x->data );
    mol2->serialize( x->data );
    T = gsl_multimin_fdfminimizer_conjugate_fr;
    gsl_multimin_fdfminimizer_free (s);
    s = gsl_multimin_fdfminimizer_alloc (T, d);
    gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);
    status = 0;
  }
  */
  
  if (status)
    return false;

  status = gsl_multimin_test_gradient (s->gradient, 1e-8);
  //cout << status << " status2\n";

  if (status == GSL_SUCCESS){
    fprintf ( stderr, "Minimum found at:\n");
  }
  fprintf ( stderr, "%5d %.5f %.5f %10.5f\n", iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1), 
          s->f);

  return (status == GSL_CONTINUE);
}



void
SimpleQuench::terminate()
{
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
  SnapShot( *fout );
}
