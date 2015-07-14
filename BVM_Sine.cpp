#include "BVM_Sine.h"
#include "vMC.h"
#include "MarginalDensitySine.h"
//#include "CustomFunctionSine.h"

BVM_Sine::BVM_Sine()
{
  mu1 = 0; mu2 = 0;
  kappa1 = 1; kappa2 = 1;
  lambda = 0; // independent
  computed = UNSET;
}

BVM_Sine::BVM_Sine(double kappa1, double kappa2, double lambda) :
                   kappa1(kappa1), kappa2(kappa2), lambda(lambda)
{
  computed = UNSET;
}

BVM_Sine::BVM_Sine(double mu1, double mu2, double kappa1, double kappa2, double lambda) :
          mu1(mu1), mu2(mu2), kappa1(kappa1), kappa2(kappa2), lambda(lambda)
{
  computed = UNSET;
}

BVM_Sine BVM_Sine::operator=(const BVM_Sine &source)
{
  if (this != &source) {
    mu1 = source.mu1;
    mu2 = source.mu2;
    kappa1 = source.kappa1;
    kappa2 = source.kappa2;
    lambda = source.lambda;
    constants = source.constants;
    computed = source.computed;
  }
  return *this;
}

std::vector<Vector> BVM_Sine::generate(int sample_size)
{
  std::vector<Vector> angle_pairs(sample_size);

  if (lambda != 0) {
    /* first, generate from marginal density */
    MarginalDensitySine marginal(mu1,kappa1,kappa2,lambda);
    
    /* check if marginal density is unimodal or bimodal */
    double ratio_bessel = computeRatioBessel(kappa2);
    double rhs = (kappa1 * kappa2) / (lambda * lambda);

    Vector thetas(sample_size,0);
    int unimodal;
    if (ratio_bessel <= rhs) unimodal = SET;
    else unimodal = UNSET;

    if (unimodal == SET) { /* marginal density is unimodal */
      Vector solution = marginal.minimize_unimodal_objective();
      double optimal_kappa = solution[1];
      vMC proposal(mu1,optimal_kappa);
      double log_max = -solution[2];
      if (log_max < 0) log_max = 0;

      int rem_sample_size = sample_size;
      int count = 0;
      int N;
      Vector proposal_thetas;
      repeat1:
      N = (int) (1.5 * rem_sample_size);
      proposal_thetas = proposal.generate(N);
      for (int i=0; i<proposal_thetas.size(); i++) {
        if (count == sample_size) break;
        double u = uniform_random();
        double log_fg = accept_reject_fval_unimodal_marginal_sine(
                          proposal_thetas[i],optimal_kappa,mu1,kappa1,kappa2,lambda
                        );
        if (log_fg > log(u) + log_max) {  /* accept */
          thetas[count++] = proposal_thetas[i];
        }
      } // for()
      if (count < sample_size) {
        rem_sample_size = sample_size - count;
        goto repeat1;
      }
    } else {  /* marginal density is bimodal */
      /* compute theta* to determine the two modes */
      //CustomFunctionSine function(mu1,kappa1,kappa2,lambda);
      //double theta_shift = function.solve();
      
      double theta_shift = marginal.solve_custom_function();
    }

    /* now, generate from conditional density which is a vMC */
    Vector pair(2,0);
    for (int i=0; i<thetas.size(); i++) {
      pair[0] = thetas[i];
      double lambda_sine = lambda * sin(thetas[i] - mu1);
      double beta = atan2(lambda_sine,kappa2);
      if (beta < 0) beta += (2 * PI);
      double m = mu2 + beta;

      double asq = kappa2 * kappa2 + lambda_sine * lambda_sine;
      double k = sqrt(asq);
      vMC vmc(m,k);
      Vector x1 = vmc.generate(1);
      pair[1] = x1[0];

      angle_pairs[i] = pair;
    } // for()
  } else {
    /* independently generate from two constituent vMC distributions */
    vMC vmc1(mu1,kappa1);
    Vector angles1 = vmc1.generate(sample_size);
    vMC vmc2(mu2,kappa2);
    Vector angles2 = vmc2.generate(sample_size);

    Vector pair(2,0);
    for (int i=0; i<sample_size; i++) {
      pair[0] = angles1[i];
      pair[1] = angles2[i];
      angle_pairs[i] = pair;
    } // for ()
  } // if (lambda != 0) {

  return angle_pairs;
}

std::vector<Vector> BVM_Sine::generate_cartesian(int sample_size)
{
  std::vector<Vector> angle_pairs = generate(sample_size);

  Vector cartesian(3,0);
  std::vector<Vector> random_sample(sample_size);
  double r1 = 3;
  double r2 = 1;
  for (int i=0; i<sample_size; i++) {
    double theta1 = angle_pairs[i][0];
    double theta2 = angle_pairs[i][1];
    cartesian[0] = (r1 + r2 * cos(theta2)) * cos(theta1); // x
    cartesian[1] = (r1 + r2 * cos(theta2)) * sin(theta1); // y
    cartesian[2] = r2 * sin(theta2);  // z
    random_sample[i] = cartesian;
  } // for()

  return random_sample;
}

BVM_Sine::Constants BVM_Sine::getConstants()
{
  if (computed != SET) {
    computeExpectation();
  }
  return constants;
}

void BVM_Sine::computeExpectation()
{
  computeConstants();

  computed = SET;
}

void BVM_Sine::computeConstants()
{
  constants.log_c = computeLogNormalizationConstant();
}

/*!
 *  logarithm of normalization constant
 */
double BVM_Sine::computeLogNormalizationConstant()
{
}

