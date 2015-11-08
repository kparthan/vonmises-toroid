#include "OptimizeInd.h"

extern int CONSTRAIN_KAPPA;
extern double MAX_KAPPA;
extern int ESTIMATION;

OptimizeInd::OptimizeInd(string type)
{
  if (type.compare("MLE") == 0) {
    ESTIMATION = MLE;
  } else if (type.compare("MAP") == 0) {
    ESTIMATION = MAP;
  } else if (type.compare("MML") == 0) {
    ESTIMATION = MML;
  } 
}

void OptimizeInd::initialize(
  double m1, double m2, double k1, double k2
) {
  mu1 = m1; mu2 = m2; kappa1 = k1; kappa2 = k2;
  
  /*if (CONSTRAIN_KAPPA == SET) {
    if (kappa >= MAX_KAPPA) {
      double e = 2 * beta / kappa;
      kappa = MAX_KAPPA - TOLERANCE;
      beta = kappa * e / 2;
    }
  }*/
}

struct EstimatesInd OptimizeInd::minimize(
  struct SufficientStatisticsInd &suff_stats
) {
  struct EstimatesInd estimates;
  int num_params = 4;
  double LIMIT = 1e-4;
  double minf;

  Vector x(num_params);
  nlopt::opt opt(nlopt::LN_COBYLA, num_params);

  Vector lb(num_params,TOLERANCE);
  Vector ub(num_params,0);
  ub[0] = 2*PI;     // mu1
  ub[1] = 2*PI;     // mu2
  ub[2] = HUGE_VAL; // kappa1
  ub[3] = HUGE_VAL; // kappa2

  x[0] = mu1; x[1] = mu2; x[2] = kappa1; x[3] = kappa2; 
  for (int i=0; i<4; i++) {
    if (x[i] < TOLERANCE) x[i] = TOLERANCE;
  }

  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.set_xtol_rel(LIMIT);
  opt.add_inequality_constraint(ConstraintInd, NULL, TOLERANCE);

  switch(ESTIMATION) {
    case MLE:
    {
      ML_Ind mle(suff_stats);
      opt.set_min_objective(ML_Ind::wrap, &mle);
      nlopt::result result = opt.optimize(x,minf);
      break;
    }

    case MAP:
    {
      MAP_Ind map(suff_stats);
      opt.set_min_objective(MAP_Ind::wrap, &map);
      nlopt::result result = opt.optimize(x,minf);
      break;
    }

    case MML:
    {
      MML_Ind mml(suff_stats);
      opt.set_min_objective(MML_Ind::wrap, &mml);
      nlopt::result result = opt.optimize(x,minf);
      break;
    }

    default:
      break;
  }
  assert(!boost::math::isnan(minf));

  //cout << "solution: "; print(cout,x,3); 
  //cout << "solution: (" << x[0] << ", " << x[1] << ")\n";
  //cout << "minf: " << minf << endl;
  estimates.mu1 = x[0];
  estimates.mu2 = x[1];
  estimates.kappa1 = x[2];
  estimates.kappa2 = x[3];

  return estimates;
}

