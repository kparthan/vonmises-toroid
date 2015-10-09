#include "OptimizeCosine.h"

extern int CONSTRAIN_KAPPA;
extern double MAX_KAPPA;
extern int ESTIMATION;

OptimizeCosine::OptimizeCosine(string type)
{
  if (type.compare("MLE") == 0) {
    ESTIMATION = MLE;
  } else if (type.compare("MAP") == 0) {
    ESTIMATION = MAP;
  } else if (type.compare("MAP_TRANSFORM") == 0) {
    ESTIMATION = MAP_TRANSFORM;
  } else if (type.compare("MML") == 0) {
    ESTIMATION = MML;
  } 
}

void OptimizeCosine::initialize(
  double m1, double m2, double k1, double k2, double k3 
) {
  mu1 = m1; mu2 = m2; kappa1 = k1; kappa2 = k2; kappa3 = k3;
  
  /*if (CONSTRAIN_KAPPA == SET) {
    if (kappa >= MAX_KAPPA) {
      double e = 2 * beta / kappa;
      kappa = MAX_KAPPA - TOLERANCE;
      beta = kappa * e / 2;
    }
  }*/
}
/*
struct EstimatesCosine OptimizeCosine::minimize(
  std::vector<Vector> &data
) {
  int num_params = 5;

  Vector x(num_params);
  nlopt::opt opt(nlopt::LN_COBYLA, num_params);

  Vector lb(num_params,TOLERANCE);
  lb[4] = -HUGE_VAL; // kappa3

  Vector ub(num_params,0);
  ub[0] = 2*PI;     // mu1
  ub[1] = 2*PI;     // mu2
  ub[2] = HUGE_VAL; // kappa1
  ub[3] = HUGE_VAL; // kappa2
  ub[4] = HUGE_VAL; // kappa3

  double LIMIT = 1e-4;
  double minf;

  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.set_xtol_rel(LIMIT);
  x[0] = mu1; x[1] = mu2; x[2] = kappa1; x[3] = kappa2; x[4] = kappa3;
  for (int i=0; i<4; i++) {
    if (x[i] < TOLERANCE) x[i] = TOLERANCE;
  }
  opt.add_inequality_constraint(ConstraintCosine, NULL, TOLERANCE);

  // ESTIMATION = PMLE
  PML_Cosine pmle(data);
  opt.set_min_objective(PML_Cosine::wrap, &pmle);
  nlopt::result result = opt.optimize(x,minf);

  assert(!boost::math::isnan(minf));

  struct EstimatesCosine estimates;
  estimates.mu1 = x[0];
  estimates.mu2 = x[1];
  estimates.kappa1 = x[2];
  estimates.kappa2 = x[3];
  estimates.kappa3 = x[4];
  return estimates;
}*/

struct EstimatesCosine OptimizeCosine::minimize(
  struct SufficientStatisticsCosine &suff_stats
) {
  struct EstimatesCosine estimates;
  int num_params = 5;
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

  if (ESTIMATION == MLE || ESTIMATION == MAP || ESTIMATION == MML) {
    ub[4] = HUGE_VAL; // kappa3
    x[4] = kappa3;
    opt.add_inequality_constraint(ConstraintCosine, NULL, TOLERANCE);
  } else if (ESTIMATION == MAP_TRANSFORM) {
    ub[4] = 1; // rho 
    x[4] = kappa3 * (kappa1 + kappa2)/ (kappa1 * kappa2);
  }

  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.set_xtol_rel(LIMIT);

  switch(ESTIMATION) {
    case MLE:
    {
      ML_Cosine mle(suff_stats);
      opt.set_min_objective(ML_Cosine::wrap, &mle);
      nlopt::result result = opt.optimize(x,minf);
      break;
    }

    case MAP:
    {
      MAP_Cosine map(suff_stats);
      opt.set_min_objective(MAP_Cosine::wrap, &map);
      nlopt::result result = opt.optimize(x,minf);
      break;
    }

    case MAP_TRANSFORM:
    {
      MAP_Transform_Cosine map2(suff_stats);
      opt.set_min_objective(MAP_Transform_Cosine::wrap, &map2);
      nlopt::result result = opt.optimize(x,minf);
      break;
    }

    case MML:
    {
      MML_Cosine mml(suff_stats);
      opt.set_min_objective(MML_Cosine::wrap, &mml);
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

  if (ESTIMATION != MAP_TRANSFORM) {
    estimates.kappa3 = x[4];
    estimates.rho = estimates.kappa3 * (estimates.kappa1 + estimates.kappa2) / (estimates.kappa1 * estimates.kappa2);
  } else if (ESTIMATION == MAP_TRANSFORM) {
    estimates.rho = x[4];
    estimates.kappa3 = estimates.rho * (estimates.kappa1 * estimates.kappa2) / (estimates.kappa1 + estimates.kappa2);
  }

  return estimates;
}

