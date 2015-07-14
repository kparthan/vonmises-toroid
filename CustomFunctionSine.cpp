#include "CustomFunctionSine.h"
#include "Support.h"

CustomFunctionSine::CustomFunctionSine(
                        double mu1, double kappa1, double kappa2, double lambda
                    ) : mu1(mu1), kappa1(kappa1), kappa2(kappa2), lambda(lambda)
{}

double CustomFunctionSine::solve()
{

/*
  UnimodalSineObjectiveFunction obj_function(mu1,kappa1,kappa2,lambda);
  opt.set_min_objective(UnimodalSineObjectiveFunction::wrap, &obj_function);
  nlopt::result result = opt.optimize(x, minf);
*/
}

