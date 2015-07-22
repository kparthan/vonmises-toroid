#include "OptimizeSineSine.h"

extern int CONSTRAIN_KAPPA;
extern double MAX_KAPPA;
extern int ESTIMATION;

OptimizeSine::OptimizeSine(string type)
{
  if (type.compare("MLE") == 0) {
    ESTIMATION = MLE;
  } else if (type.compare("MAP") == 0) {
    ESTIMATION = MAP;
  } else if (type.compare("MML") == 0) {
    ESTIMATION = MML;
  } 
}

void OptimizeSine::initialize(
  double m1, double m2, double k1, double k2, double lam
) {
  mu1 = m1; mu2 = m2; kappa1 = k1; kappa2 = k2; lambda = lam;
  
  /*if (CONSTRAIN_KAPPA == SET) {
    if (kappa >= MAX_KAPPA) {
      double e = 2 * beta / kappa;
      kappa = MAX_KAPPA - TOLERANCE;
      beta = kappa * e / 2;
    }
  }*/
}

Vector OptimizeSine::minimize(struct SufficientStatisticsSine &suff_stats)
{
  int num_params = 5;

  Vector x(num_params);
  nlopt::opt opt(nlopt::LN_COBYLA, num_params);

  Vector lb(num_params,TOLERANCE);
  lb[4] = -HUGE_VAL; // lambda

  Vector ub(num_params,0);
  ub[0] = 2*PI;     // mu1
  ub[1] = 2*PI;     // mu2
  ub[2] = HUGE_VAL; // kappa1
  ub[3] = HUGE_VAL; // kappa2
  ub[4] = HUGE_VAL; // lambda

  double LIMIT = 1e-4;
  double minf;

  switch(ESTIMATION) {
    case MLE:
    {
      opt.set_lower_bounds(lb);
      opt.set_upper_bounds(ub);

      ML_Sine mle(suff_stats);
      opt.set_min_objective(ML_Sine::wrap, &mle);
      opt.set_xtol_rel(LIMIT);

      x[0] = mu1; x[1] = mu2; x[2] = kappa1; x[3] = kappa2; x[4] = lambda;
      for (int i=0; i<4; i++) {
        if (x[i] < TOLERANCE) x[i] = TOLERANCE;
      }
      nlopt::result result = opt.optimize(x,minf);
      assert(!boost::math::isnan(minf));
      break;
    }

    default:
      break;
  }
  //cout << "solution: "; print(cout,x,3); 
  //cout << "solution: (" << x[0] << ", " << x[1] << ")\n";
  //cout << "minf: " << minf << endl;
  return x;
}

