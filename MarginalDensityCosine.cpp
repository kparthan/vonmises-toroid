#include "MarginalDensityCosine.h"

extern double MAX_KAPPA;

MarginalDensityCosine::MarginalDensityCosine(
                        double mu1, double kappa1, double kappa2, double kappa3
                       ) : mu1(mu1), kappa1(kappa1), kappa2(kappa2), kappa3(kappa3)
{}

// unimodal marginal density
Vector MarginalDensityCosine::minimize_unimodal_objective()
{
  int num_params = 1;

  Vector x(num_params);  // x[0]: theta, x[1]: kappa
  nlopt::opt opt(nlopt::LN_COBYLA,num_params);

  Vector lb(num_params,TOLERANCE);
  Vector ub(num_params,0);
  ub[0] = 2 * PI; //ub[1] = 10000;

  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.set_xtol_rel(1e-4);

  x[0] = mu1;
  //x[1] = kappa1;

  if (x[0] < TOLERANCE) x[0] = TOLERANCE;

  double minf;

  UnimodalCosineObjectiveFunction function(mu1,kappa1,kappa2,kappa3);
  opt.set_min_objective(UnimodalCosineObjectiveFunction::wrap, &function);
  nlopt::result result = opt.optimize(x, minf);

  //cout << "solution: (" << x[0]*180/PI << ", " << x[1] << ")\n";
  //cout << "solution: (" << x[0]*180/PI << ")\n";
  //cout << "minf: " << minf << endl;

  Vector y = x;
  y.push_back(minf);

  return y;
}

// bimodal marginal density
double MarginalDensityCosine::solve_custom_function()
{
  //double left = mu1 - PI/2;
  //double right = mu1 + PI/2;
  double left = 0;
  double right = PI;

  CustomFunctionCosine function(mu1,kappa1,kappa2,kappa3);
  std::pair<double, double> result 
      = boost::math::tools::bisect(
          boost::bind(&CustomFunctionCosine::solve,&function,_1),
          left,right, TerminationCondition()
        );
  double root = (result.first + result.second) / 2;  
  return root;
}

// bimodal marginal density
Vector MarginalDensityCosine::minimize_bimodal_objective(double mode1, double mode2)
{
  int num_params = 1;

  Vector x(num_params);  // x[0]: theta, x[1]: kappa
  nlopt::opt opt(nlopt::LN_COBYLA,num_params);

  Vector lb(num_params,TOLERANCE);
  Vector ub(num_params,0);
  ub[0] = 2 * PI; //ub[1] = MAX_KAPPA;

  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.set_xtol_rel(1e-4);

  x[0] = mu1;
  //x[1] = kappa1;

  if (x[0] < TOLERANCE) x[0] = TOLERANCE;

  double minf;

  BimodalCosineObjectiveFunction function(mu1,kappa1,kappa2,kappa3,mode1,mode2);
  opt.set_min_objective(BimodalCosineObjectiveFunction::wrap, &function);
  nlopt::result result = opt.optimize(x, minf);

  //cout << "solution: (" << x[0] << ", " << x[1] << ")\n";
  cout << "solution: (" << x[0]*180/PI << ")\n";
  cout << "minf: " << minf << endl;

  Vector y = x;
  y.push_back(minf);

  return y;
}

