#include "KappaSolver.h"

extern double MAX_KAPPA;

KappaSolver::KappaSolver(
  double N, double constant, double rbar
) : N(N), constant(constant), rbar(rbar)
{}

double KappaSolver::minimize()
{
  int num_params = 1;

  Vector x(num_params);  
  nlopt::opt opt(nlopt::LN_COBYLA,num_params);

  Vector lb(num_params,TOLERANCE);
  //Vector ub(num_params,MAX_KAPPA);

  opt.set_lower_bounds(lb);
  //opt.set_upper_bounds(ub);
  opt.set_xtol_rel(1e-4);

  x[0] = banerjee_approx(rbar);
  if (x[0] < TOLERANCE) x[0] = TOLERANCE;
  //cout << "banerjee: " << x[0] << endl;

  double minf;

  ObjectiveFunction function(N,constant);
  opt.set_min_objective(ObjectiveFunction::wrap, &function);
  nlopt::result result = opt.optimize(x, minf);

  //cout << "kappa: (" << x[0] << ")\n";
  //cout << "minf: " << minf << endl;

  return x[0];
}

