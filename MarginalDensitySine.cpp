#include "MarginalDensitySine.h"

extern int CONSTRAIN_KAPPA;
extern double MAX_KAPPA;

MarginalDensitySine::MarginalDensitySine(
                        double mu1, double kappa1, double kappa2, double lambda
                     ) : mu1(mu1), kappa1(kappa1), kappa2(kappa2), lambda(lambda)
{}

// unimodal marginal density
Vector MarginalDensitySine::minimize_unimodal_objective()
{
  Vector x(2);  // x[0]: theta, x[1]: kappa
  nlopt::opt opt(nlopt::LN_COBYLA,2);

  Vector lb(2,TOLERANCE);
  Vector ub(2,0);
  ub[0] = 2 * PI; ub[1] = MAX_KAPPA;

  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.set_xtol_rel(1e-4);

  x[0] = mu1;
  x[1] = kappa1;

  double minf;

  UnimodalSineObjectiveFunction obj_function(mu1,kappa1,kappa2,lambda);
  opt.set_min_objective(UnimodalSineObjectiveFunction::wrap, &obj_function);
  nlopt::result result = opt.optimize(x, minf);

  cout << "solution: (" << x[0] << ", " << x[1] << ")\n";
  cout << "minf: " << minf << endl;

  Vector y = x;
  y.push_back(minf);

  return y;
}

// bimodal marginal density
double MarginalDensitySine::solve_custom_function()
{
  CustomFunctionSine function(mu1,kappa1,kappa2,lambda);
}

