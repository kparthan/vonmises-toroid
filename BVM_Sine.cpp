#include "BVM_Sine.h"
#include "vMC.h"
#include "Mixture_vMC.h"
#include "MarginalDensitySine.h"
#include "KappaSolver.h"
#include "OptimizeSine.h"

extern int ESTIMATION;

BVM_Sine::BVM_Sine()
{
  mu1 = 0; mu2 = 0;
  kappa1 = 1; kappa2 = 1;
  lambda = 0; // independent
  computed = UNSET;
  computed_lognorm = UNSET;
}

BVM_Sine::BVM_Sine(double kappa1, double kappa2, double lambda) :
                   kappa1(kappa1), kappa2(kappa2), lambda(lambda)
{
  mu1 = 0; mu2 = 0;
  computed = UNSET;
  computed_lognorm = UNSET;
}

BVM_Sine::BVM_Sine(
  double mu1, double mu2, double kappa1, double kappa2, double lambda
) : mu1(mu1), mu2(mu2), kappa1(kappa1), kappa2(kappa2), lambda(lambda)
{
  computed = UNSET;
  computed_lognorm = UNSET;
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
    computed_lognorm = source.computed_lognorm;
  }
  return *this;
}

double BVM_Sine::Mean1()
{
  return mu1;
}

double BVM_Sine::Mean2()
{
  return mu2;
}

double BVM_Sine::Kappa1()
{
  return kappa1;
}

double BVM_Sine::Kappa2()
{
  return kappa2;
}

double BVM_Sine::Lambda()
{
  return lambda;
}

// generates <theta1,theta2> pairs such that theta1,theta2 \in [0,2pi)
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
      //cout << "unimodal: \n";
      Vector solution = marginal.minimize_unimodal_objective();
      //double optimal_kappa = solution[1];
      //cout << "optimal_kappa: " << optimal_kappa << endl;
      double optimal_kappa = 1;

      vMC proposal(mu1,optimal_kappa);
      double log_max = -solution[1];
      if (log_max < 0) log_max = 0;
      //cout << "log_max: " << log_max << endl;

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
      cout << "bimodal: \n";
      //double theta_shift = marginal.solve_custom_function();
      double theta_shift = 0;
      cout << "mu1: " << mu1 * 180/PI << " degrees ...\n";
      cout << "theta_shift: " << theta_shift * 180/PI << " degrees ...\n";
      double m1 = mu1 - theta_shift;
      if (m1 < 0) m1 += (2 * PI);
      double m2 = mu2 + theta_shift;
      if (m2 < 0) m2 += (2 * PI);

      Vector solution = marginal.minimize_bimodal_objective(m1,m2);
      //double optimal_kappa = solution[1];
      //cout << "optimal_kappa: " << optimal_kappa << endl;
      double optimal_kappa = 1;

      vMC vmc1(m1,optimal_kappa);
      vMC vmc2(m2,optimal_kappa);
      Vector weights(2,0.5);
      std::vector<vMC> components;
      components.push_back(vmc1); components.push_back(vmc2);
      Mixture_vMC proposal(2,weights,components);
      double log_max = -solution[1];
      if (log_max < 0) log_max = 0;
      cout << "log_max: " << log_max << endl;

      int rem_sample_size = sample_size;
      int count = 0;
      int N;
      Vector proposal_thetas;
      repeat2:
      N = (int) (1.5 * rem_sample_size);
      proposal_thetas = proposal.generate(N);
      for (int i=0; i<proposal_thetas.size(); i++) {
        if (count == sample_size) break;
        double u = uniform_random();
        double log_fg = accept_reject_fval_bimodal_marginal_sine(
                          proposal_thetas[i],optimal_kappa,mu1,kappa1,kappa2,lambda,m1,m2
                        );
        if (log_fg > log(u) + log_max) {  /* accept */
          thetas[count++] = proposal_thetas[i];
        }
      } // for()
      if (count < sample_size) {
        rem_sample_size = sample_size - count;
        goto repeat2;
      }
    } // if (unimodal == SET)

    /* now, generate from conditional density which is a vMC */
    Vector pair(2,0);
    for (int i=0; i<thetas.size(); i++) {
      pair[0] = thetas[i];
      /*double lambda_sine = lambda * sin(thetas[i] - mu1);
      double beta = atan2(lambda_sine,kappa2);
      double m = mu2 + beta;
      if (m < 0) m += (2 * PI);

      double asq = kappa2 * kappa2 + lambda_sine * lambda_sine;
      double k = sqrt(asq);
      vMC vmc(m,k);*/
      vMC vmc = getConditionalDensitySine(thetas[i],mu1,mu2,kappa2,lambda);
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

  std::vector<Vector> random_sample = generate_cartesian(angle_pairs);

  return random_sample;
}

// convert theta,phi --> cartesian coordinates
std::vector<Vector> BVM_Sine::generate_cartesian(
  std::vector<Vector> &angle_pairs
) {
  int sample_size = angle_pairs.size();
  Vector cartesian(3,0);
  std::vector<Vector> random_sample(sample_size);
  for (int i=0; i<sample_size; i++) {
    toroid2cartesian(angle_pairs[i],cartesian);
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

  double log_ans;

  // E[ cos(t1-mu1) ] = (1/c) (dc/dk1)      > 0
  log_ans = -constants.log_c + constants.log_dc_dk1;
  constants.ck1_c = exp(log_ans);

  // E[ cos(t2-mu2) ] = (1/c) (dc/dk2)      > 0
  log_ans = -constants.log_c + constants.log_dc_dk2;
  constants.ck2_c = exp(log_ans);

  // E[ cos(t1-mu1) cos(t2-mu2) ]            > 0
  // = (1/c) (d^2 c/dk1 dk2)
  log_ans = -constants.log_c + constants.log_d2c_dk1dk2;
  constants.ck1k2_c = exp(log_ans);

  // E[ cos(t1-mu1)^2 ] = (1/c) (d^2 c/dk1^2)      > 0
  log_ans =  -constants.log_c + constants.log_d2c_dk1dk1;
  constants.ck1k1_c = exp(log_ans);

  // E[ cos(t2-mu2)^2 ] = (1/c) (d^2 c/dk2^2)      > 0
  log_ans =  -constants.log_c + constants.log_d2c_dk2dk2;
  constants.ck2k2_c = exp(log_ans);

  // E[ sin(t1-mu1) sin(t2-mu2) ]   *check sign (of lambda)*
  // = (1/c) (dc/dl)
  log_ans = -constants.log_c + constants.log_dc_dl;
  constants.cl_c = sign(lambda) * exp(log_ans);
  //cout << "cl_c: " << constants.cl_c << endl;

  // E[ sin(t1-mu1)^2 sin(t2-mu2)^2 ]  = (1/c) (d^2 c/dl^2) > 0
  log_ans = -constants.log_c + constants.log_d2c_dldl;
  constants.cll_c = exp(log_ans);

  // E[ cos(t1-mu1) sin(t1-mu1) sin(t2-mu2) ]   *check sign (of lambda)*
  // = (1/c) (d^2 c/dk1 dl)
  log_ans = -constants.log_c + constants.log_d2c_dk1dl;
  constants.ck1l_c = sign(lambda) * exp(log_ans);

  // E[ cos(t2-mu2) sin(t1-mu1) sin(t2-mu2) ]   *check sign (of lambda)*
  // = (1/c) (d^2 c/dk2 dl)
  log_ans = -constants.log_c + constants.log_d2c_dk2dl;
  constants.ck2l_c = sign(lambda) * exp(log_ans);

  computed = SET;
}

void BVM_Sine::computeConstants()
{
  constants.log_c = computeLogNormalizationConstant();  // > 0
  constants.log_dc_dk1 = compute_series_A(1,0);         // > 0
  constants.log_dc_dk2 = compute_series_A(0,1);         // > 0
  constants.log_d2c_dk1dk2 = compute_series_A(1,1);     // > 0

  double log_sum = compute_series_A(2,0);
  double log_diff = log_sum - constants.log_dc_dk1;
  double diff = exp(log_diff);
  constants.log_d2c_dk1dk1 = constants.log_dc_dk1 + log(diff+1/kappa1);

  log_sum = compute_series_A(0,2);
  log_diff = log_sum - constants.log_dc_dk2;
  diff = exp(log_diff);
  constants.log_d2c_dk2dk2 = constants.log_dc_dk2 + log(diff+1/kappa2);

  constants.log_dc_dl = compute_series_B(0,0);          // > 0
  constants.log_d2c_dk1dl = compute_series_B(1,0);      // > 0
  constants.log_d2c_dk2dl = compute_series_B(0,1);      // > 0
  constants.log_d2c_dldl = compute_series_C();          // > 0
}

double BVM_Sine::getLogNormalizationConstant()
{
  if (computed_lognorm == UNSET) {
    constants.log_c = computeLogNormalizationConstant();  
  }
  return constants.log_c;
}

/*!
 *  logarithm of normalization constant
 */
double BVM_Sine::computeLogNormalizationConstant()
{
  double log_series_sum = compute_series_A(0,0);
  //cout << "log_series_sum: " << log_series_sum << endl;
  computed_lognorm = SET;
  return log_series_sum;
}

// used to compute logarithm of: 
//  dc/dk1, dc/dk2, d2c/dk1dk2
//  d^2c/dk1^2, d^2c/dk2^2
double BVM_Sine::compute_series_A(double a, double b)
{
  // const = (lambda^2) / (2 * k1 * k2)
  double log_const = 2 * log(fabs(lambda)) - log(2) - log(kappa1) - log(kappa2);

  double log_bessel1_prev = computeLogModifiedBesselFirstKind(a,kappa1);
  double log_bessel2_prev = computeLogModifiedBesselFirstKind(b,kappa2);
  double log_f0 = log_bessel1_prev + log_bessel2_prev;
  double log_fj_prev = log_f0;
  double j = 0;
  double series_sum = 1;

  while (1) {
    double log_bessel1_current = computeLogModifiedBesselFirstKind(j+1+a,kappa1);
    double log_bessel2_current = computeLogModifiedBesselFirstKind(j+1+b,kappa2);
    double log_fj_current = log_const + log(2*j+1) - log(j+1)
                            + log_bessel1_current + log_bessel2_current
                            - log_bessel1_prev - log_bessel2_prev
                            + log_fj_prev;
    double log_diff = log_fj_current - log_f0;
    double current = exp(log_diff); // tj
    series_sum += current;
    if (current/series_sum <= 1e-6) {
      /*cout << "current/series_sum: " << current/series_sum << endl;
      cout << "j: " << j << "; series_sum: " << series_sum << endl;
      cout << "log_diff: " << log_diff << endl;
      cout << "log_const: " << log_const << endl;
      cout << "log_fj_prev: " << log_fj_prev << endl;
      cout << "log_fj_current: " << log_fj_current << endl;
      cout << "log_bessel1_prev: " << log_bessel1_prev << endl;
      cout << "log_bessel2_prev: " << log_bessel2_prev << endl;
      cout << "log_bessel1_current: " << log_bessel1_current << endl;
      cout << "log_bessel2_current: " << log_bessel2_current << endl;
      cout << "current: " << current << endl << endl;*/
      break;
    }
    j++;
    log_bessel1_prev = log_bessel1_current; 
    log_bessel2_prev = log_bessel2_current; 
    log_fj_prev = log_fj_current;
  } // while(1)
  //cout << "j: " << j << endl;
  double ans = 2*log(2*PI) + log_f0 + log(series_sum);
  return ans;
}

double BVM_Sine::compute_series_B(double a, double b)
{
  // const = (lambda^2) / (2 * k1 * k2)
  double log_const = 2 * log(fabs(lambda)) - log(2) - log(kappa1) - log(kappa2);

  double log_bessel1_prev = computeLogModifiedBesselFirstKind(1+a,kappa1);
  double log_bessel2_prev = computeLogModifiedBesselFirstKind(1+b,kappa2);
  double log_f1 = log_const + log_bessel1_prev + log_bessel2_prev;
  double log_fj_prev = log_f1;
  double j = 1;
  double series_sum = 1;
  while (1) {
    double log_bessel1_current = computeLogModifiedBesselFirstKind(j+1+a,kappa1);
    double log_bessel2_current = computeLogModifiedBesselFirstKind(j+1+b,kappa2);
    double log_fj_current = log_const + log(2*j+1) - log(j)
                            + log_bessel1_current + log_bessel2_current
                            - log_bessel1_prev - log_bessel2_prev
                            + log_fj_prev;
    double log_diff = log_fj_current - log_f1;
    double current = exp(log_diff); // tj
    series_sum += current;
    if (current/series_sum <= 1e-6) {
      break;
    }
    j++;
    log_bessel1_prev = log_bessel1_current; 
    log_bessel2_prev = log_bessel2_current; 
    log_fj_prev = log_fj_current;
  } // while(1)
  //cout << "j: " << j << endl;
  double ans = 2*log(2*PI) + log(2) + log_f1 + log(series_sum) - log(fabs(lambda));
  return ans;
}

double BVM_Sine::compute_series_C()
{
  // const = (lambda^2) / (2 * k1 * k2)
  double log_const = 2 * log(fabs(lambda)) - log(2) - log(kappa1) - log(kappa2);

  double log_bessel1_prev = computeLogModifiedBesselFirstKind(1,kappa1);
  double log_bessel2_prev = computeLogModifiedBesselFirstKind(1,kappa2);
  double log_f1 = log_const + log_bessel1_prev + log_bessel2_prev;
  double log_fj_prev = log_f1;
  double j = 1;
  double series_sum = 1;
  while (1) {
    double log_bessel1_current = computeLogModifiedBesselFirstKind(j+1,kappa1);
    double log_bessel2_current = computeLogModifiedBesselFirstKind(j+1,kappa2);
    double log_fj_current = log_const + 2*log(2*j+1) - log(j) - log(2*j-1)
                            + log_bessel1_current + log_bessel2_current
                            - log_bessel1_prev - log_bessel2_prev
                            + log_fj_prev;
    double log_diff = log_fj_current - log_f1;
    double current = exp(log_diff); // tj
    series_sum += current;
    if (current/series_sum <= 1e-6) {
      break;
    }
    j++;
    log_bessel1_prev = log_bessel1_current; 
    log_bessel2_prev = log_bessel2_current; 
    log_fj_prev = log_fj_current;
  } // while(1)
  //cout << "j: " << j << endl;
  double ans = 2*log(2*PI) + log(2) + log_f1 + log(series_sum) - 2*log(fabs(lambda));
  return ans;
}

double BVM_Sine::computeLogParametersProbability(double Neff)
{
  double log_prior_density = computeLogParametersPriorDensity();
  double log_expected_fisher = computeLogFisherInformation(Neff);
  double logp = -log_prior_density + 0.5 * log_expected_fisher;
  return logp;
}

// prior density of parameters ... h(kappa1,kappa2,lambda)
double BVM_Sine::computeLogParametersPriorDensity()
{
  // prior on means (uniform)
  double log_prior_means = -2 * log(2*PI);

  // joint prior on kappas and lambda (vMC)
  double log_scale = 0.5 * log(kappa1) - 1.5 * log(1+kappa1*kappa1);
  log_scale += (0.5 * log(kappa2) - 1.5 * log(1+kappa2*kappa2));
  log_scale -= log(2);

  double log_prior = log_prior_means + log_scale;
  return log_prior;
}

// prior density of parameters ... h(kappa1,kappa2,rho)
double BVM_Sine::computeLogParametersPriorDensityTransform()
{
  // prior on means (uniform)
  double log_prior_means = -2 * log(2*PI);

  // joint prior on kappas and rho (correlation) 
  double log_scale = log(kappa1) - 1.5 * log(1+kappa1*kappa1);
  log_scale += (log(kappa2) - 1.5 * log(1+kappa2*kappa2));
  log_scale -= log(2);

  double log_prior = log_prior_means + log_scale;
  return log_prior;
}

double BVM_Sine::computeLogFisherInformation(double N)
{
  double log_fisher = computeLogFisherInformation_Single(); 
  log_fisher += (5 * log(N));
  assert(!boost::math::isnan(log_fisher));
  return log_fisher;
}

double BVM_Sine::computeLogFisherInformation_Single()
{
  if (computed == UNSET) {
    computeExpectation();
  }
  double log_det_axes = computeLogFisherAxes();
  double log_det_scale = computeLogFisherScale();
  return log_det_axes + log_det_scale; 
}

double BVM_Sine::computeLogFisherAxes()
{
  // L: -ve log-likelihood

  // E [d^2 L / d mu1^2]
  double t1 = kappa1 * constants.ck1_c + lambda * constants.cl_c;

  // E [d^2 L / d mu2^2]
  double t2 = kappa2 * constants.ck2_c + lambda * constants.cl_c;

  // E [d^2 L / dmu1 dmu2]
  double t3 = -lambda * constants.ck1k2_c;

  double det = fabs(t1 * t2 - t3 * t3);
  //assert(det > 0);
  //cout << "det: " << det << endl;
  return log(det);
}

double BVM_Sine::computeLogFisherScale()
{
  Matrix fisher_scale = ZeroMatrix(3,3);

  // L: -ve log-likelihood

  // E [d^2 L / d k1^2]
  fisher_scale(0,0) = constants.ck1k1_c - (constants.ck1_c * constants.ck1_c);

  // E [d^2 L / dk1 dk2]
  fisher_scale(0,1) = constants.ck1k2_c - (constants.ck1_c * constants.ck2_c);

  // E [d^2 L / dk1 dl]
  fisher_scale(0,2) = constants.ck1l_c - (constants.ck1_c * constants.cl_c);

  // E [d^2 L / d k2^2]
  fisher_scale(1,1) = constants.ck2k2_c - (constants.ck2_c * constants.ck2_c);

  // E [d^2 L / dk2 dl]
  fisher_scale(1,2) = constants.ck2l_c - (constants.ck2_c * constants.cl_c);

  // E [d^2 L / dl^2]
  fisher_scale(2,2) = constants.cll_c - (constants.cl_c * constants.cl_c);

  fisher_scale(1,0) = fisher_scale(0,1);
  fisher_scale(2,0) = fisher_scale(0,2);
  fisher_scale(2,1) = fisher_scale(1,2);

  double det = determinant_3d(fisher_scale);
  //cout << "det: " << det << endl; //exit(1);
  return log(fabs(det));
}

double BVM_Sine::log_density(Vector &angle_pair)
{
  assert(angle_pair.size() == 2);
  return log_density(angle_pair[0],angle_pair[1]);
}

// log(pdf)
double BVM_Sine::log_density(double &theta1, double &theta2)
{
  if (computed != SET) {
    computeExpectation();
  }
  double ans = 0;
  ans -= constants.log_c;
  ans += (kappa1 * cos(theta1-mu1));
  ans += (kappa2 * cos(theta2-mu2));
  ans += (lambda * sin(theta1-mu1) * sin(theta2-mu2));
  return ans;
}

// data = angle_pairs
double BVM_Sine::computeNegativeLogLikelihood(std::vector<Vector> &data)
{
  struct SufficientStatisticsSine suff_stats;
  computeSufficientStatisticsSine(data,suff_stats);
  return computeNegativeLogLikelihood(suff_stats);
}

double BVM_Sine::computeNegativeLogLikelihood(
  struct SufficientStatisticsSine &suff_stats
) {
  double cosm1 = cos(mu1);
  double sinm1 = sin(mu1);
  double cosm2 = cos(mu2);
  double sinm2 = sin(mu2);

  // compute log-likelihood
  double ans = -suff_stats.N * getLogNormalizationConstant();
  ans += (kappa1 * cosm1 * suff_stats.cost1);
  ans += (kappa1 * sinm1 * suff_stats.sint1);
  ans += (kappa2 * cosm2 * suff_stats.cost2);
  ans += (kappa2 * sinm2 * suff_stats.sint2);

  ans += (lambda * cosm1 * cosm2 * suff_stats.sint1sint2);
  ans -= (lambda * cosm1 * sinm2 * suff_stats.sint1cost2);
  ans -= (lambda * sinm1 * cosm2 * suff_stats.cost1sint2);
  ans += (lambda * sinm1 * sinm2 * suff_stats.cost1cost2);

  // return -ve log-likelihood
  return -ans;
}

double BVM_Sine::computeNegativeLogLikelihood(
  struct EstimatesSine &estimates, struct SufficientStatisticsSine &suff_stats
) {
  BVM_Sine bvm_sine(
    estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.lambda
  );
  return bvm_sine.computeNegativeLogLikelihood(suff_stats);
}

// data = angle_pairs
double BVM_Sine::computeMessageLength(std::vector<Vector> &data)
{
  struct SufficientStatisticsSine suff_stats;
  computeSufficientStatisticsSine(data,suff_stats);
  return computeMessageLength(suff_stats);
}

double BVM_Sine::computeMessageLength(
  struct SufficientStatisticsSine &suff_stats
) {
  double log_prior = computeLogParametersPriorDensity();
  double log_fisher = computeLogFisherInformation(suff_stats.N);
  double part1 = -6.455 - log_prior + 0.5 * log_fisher;
  double part2 = computeNegativeLogLikelihood(suff_stats) + 2.5
                 - 2 * suff_stats.N * log(AOM);
  double msglen = part1 + part2;
  return msglen/log(2);
}

double BVM_Sine::computeMessageLength(
  struct EstimatesSine &estimates,
  struct SufficientStatisticsSine &suff_stats
) {
  BVM_Sine bvm_sine(
    estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.lambda
  );
  return bvm_sine.computeMessageLength(suff_stats);
}

// data = angle_pairs
void BVM_Sine::computeAllEstimators(
  std::vector<Vector> &data, 
  std::vector<struct EstimatesSine> &all_estimates,
  int verbose,
  int compute_kldiv
) {
  struct SufficientStatisticsSine suff_stats;
  computeSufficientStatisticsSine(data,suff_stats);

  computeAllEstimators(
    data,suff_stats,all_estimates,verbose,compute_kldiv
  );
}

void BVM_Sine::computeAllEstimators(
  std::vector<Vector> &data, 
  struct SufficientStatisticsSine &suff_stats,
  std::vector<struct EstimatesSine> &all_estimates,
  int verbose,
  int compute_kldiv
) {
  double msglen,negloglike,kldiv,min_msg;
  int min_index;

  all_estimates.clear();
  all_estimates = std::vector<struct EstimatesSine>(NUM_METHODS);

  string type = "initial";
  struct EstimatesSine initial_est = computeInitialEstimates(suff_stats);
  print(type,initial_est);

  /*type = "PMLE";
  struct EstimatesSine pml_est = initial_est;
  OptimizeSine opt_pmle(type);
  opt_pmle.initialize(
    pml_est.mu1,pml_est.mu2,pml_est.kappa1,pml_est.kappa2,pml_est.lambda
  );
  pml_est = opt_pmle.minimize(data);
  pml_est.msglen = computeMessageLength(pml_est,suff_stats);
  pml_est.negloglike = computeNegativeLogLikelihood(pml_est,suff_stats);
  if (compute_kldiv) {
    pml_est.kldiv = computeKLDivergence(pml_est);
  }
  if (verbose) {
    print(type,pml_est);
    cout << fixed << "msglen: " << pml_est.msglen << endl;
    cout << "negloglike: " << pml_est.negloglike << endl;
    cout << "KL-divergence: " << pml_est.kldiv << endl << endl;
  }
  all_estimates[PMLE] = pml_est;
  //if (pml_est.msglen < min_msg) {
    min_index = PMLE;
    min_msg = pml_est.msglen;
  //}
  */

  type = "MLE";
  struct EstimatesSine ml_est = initial_est;
  OptimizeSine opt_mle(type);
  opt_mle.initialize(
    ml_est.mu1,ml_est.mu2,ml_est.kappa1,ml_est.kappa2,ml_est.lambda
  );
  ml_est = opt_mle.minimize(suff_stats);
  ml_est.msglen = computeMessageLength(ml_est,suff_stats);
  ml_est.negloglike = computeNegativeLogLikelihood(ml_est,suff_stats);
  if (compute_kldiv) {
    ml_est.kldiv = computeKLDivergence(ml_est);
  }
  if (verbose) {
    print(type,ml_est);
    cout << fixed << "msglen: " << ml_est.msglen << endl;
    cout << "negloglike: " << ml_est.negloglike << endl;
    cout << "KL-divergence: " << ml_est.kldiv << endl << endl;
  }
  all_estimates[MLE] = ml_est;
  //if (ml_est.msglen < min_msg) {
    min_index = MLE;
    min_msg = ml_est.msglen;
  //}

  type = "MAP";
  struct EstimatesSine map_est = initial_est;
  OptimizeSine opt_map(type);
  opt_map.initialize(
    map_est.mu1,map_est.mu2,map_est.kappa1,map_est.kappa2,map_est.lambda
  );
  map_est = opt_map.minimize(suff_stats);
  map_est.msglen = computeMessageLength(map_est,suff_stats);
  map_est.negloglike = computeNegativeLogLikelihood(map_est,suff_stats);
  if (compute_kldiv) {
    map_est.kldiv = computeKLDivergence(map_est);
  }
  if (verbose) {
    print(type,map_est);
    cout << fixed << "msglen: " << map_est.msglen << endl;
    cout << "negloglike: " << map_est.negloglike << endl;
    cout << "KL-divergence: " << map_est.kldiv << endl << endl;
  }
  all_estimates[MAP] = map_est;
  if (map_est.msglen < min_msg) {
    min_index = MAP;
    min_msg = map_est.msglen;
  }

  /***************************** MAP Variant *********************************/

  type = "MAP_TRANSFORM";
  struct EstimatesSine map2_est = initial_est;
  OptimizeSine opt_map2(type);
  opt_map2.initialize(
    map2_est.mu1,map2_est.mu2,map2_est.kappa1,map2_est.kappa2,map2_est.lambda
  );
  map2_est = opt_map2.minimize(suff_stats);
  map2_est.msglen = computeMessageLength(map2_est,suff_stats);
  map2_est.negloglike = computeNegativeLogLikelihood(map2_est,suff_stats);
  if (compute_kldiv) {
    map2_est.kldiv = computeKLDivergence(map2_est);
  }
  if (verbose) {
    print(type,map2_est);
    cout << fixed << "msglen: " << map2_est.msglen << endl;
    cout << "negloglike: " << map2_est.negloglike << endl;
    cout << "KL-divergence: " << map2_est.kldiv << endl << endl;
  }
  all_estimates[MAP_TRANSFORM] = map2_est;
  if (map2_est.msglen < min_msg) {
    min_index = MAP_TRANSFORM;
    min_msg = map2_est.msglen;
  }

  /***************************** MAP Variant *********************************/

  type = "MML";
  //struct EstimatesSine mml_est = initial_est;
  struct EstimatesSine mml_est = all_estimates[min_index];
  OptimizeSine opt_mml(type);
  opt_mml.initialize(
    mml_est.mu1,mml_est.mu2,mml_est.kappa1,mml_est.kappa2,mml_est.lambda
  );
  mml_est = opt_mml.minimize(suff_stats);
  mml_est.msglen = computeMessageLength(mml_est,suff_stats);
  mml_est.negloglike = computeNegativeLogLikelihood(mml_est,suff_stats);
  if (compute_kldiv) {
    mml_est.kldiv = computeKLDivergence(mml_est);
  }
  if (verbose) {
    print(type,mml_est);
    cout << fixed << "msglen: " << mml_est.msglen << endl;
    cout << "negloglike: " << mml_est.negloglike << endl;
    cout << "KL-divergence: " << mml_est.kldiv << endl << endl;
  }
  all_estimates[MML] = mml_est;
  if (mml_est.msglen < min_msg) {
    min_index = MML;
    min_msg = mml_est.msglen;
  }

  /* MAP TYPE 3 : ROSENBLATT TRANSFORM */
  // same as ML estimates ...
  // print the transformation z1,...,z5
  double z1 = all_estimates[MLE].mu1 / (2*PI);
  double z2 = all_estimates[MLE].mu2 / (2*PI);
  double z3 = 1 - cos(atan(all_estimates[MLE].kappa1));
  double z4 = 1 - cos(atan(all_estimates[MLE].kappa2));
  double z5 = (all_estimates[MLE].rho + 1) * 0.5;
  cout << "ROSENBLATT TRANSFORM ...\n";
  cout << "z1: " << z1 << "; z2: " << z2 << "; z3: " << z3 
       << "; z4: " << z4 << "; z5: " << z5 << endl;
}

struct EstimatesSine BVM_Sine::computeInitialEstimates(
  struct SufficientStatisticsSine &suff_stats
) {
  struct EstimatesSine estimates;
  
  double mu1_est = atan2(suff_stats.sint1,suff_stats.cost1); 
  if (mu1_est < 0) mu1_est += (2*PI);
  estimates.mu1 = mu1_est;
  double mu2_est = atan2(suff_stats.sint2,suff_stats.cost2); 
  if (mu2_est < 0) mu2_est += (2*PI);
  estimates.mu2 = mu2_est;

  double cos_sum = suff_stats.cost1;
  double sin_sum = suff_stats.sint1;
  double rbar = (sqrt(cos_sum * cos_sum + sin_sum * sin_sum))/suff_stats.N;
  //cout << "rbar: " << rbar << endl;
  double constant = cos(estimates.mu1) * suff_stats.cost1 
                    + sin(estimates.mu1) * suff_stats.sint1;
  KappaSolver solver1(suff_stats.N,constant,rbar);
  estimates.kappa1 = solver1.minimize();
  //cout << "kappa1_init: " << estimates.kappa1 << endl;

  cos_sum = suff_stats.cost2;
  sin_sum = suff_stats.sint2;
  rbar = (sqrt(cos_sum * cos_sum + sin_sum * sin_sum))/suff_stats.N;
  //cout << "rbar: " << rbar << endl;
  constant = cos(estimates.mu2) * suff_stats.cost2 
             + sin(estimates.mu2) * suff_stats.sint2;
  KappaSolver solver2(suff_stats.N,constant,rbar);
  estimates.kappa2 = solver2.minimize();
  //cout << "kappa2_init: " << estimates.kappa2 << endl;

  estimates.rho = 0.5;
  estimates.lambda = estimates.rho * sqrt(estimates.kappa1 * estimates.kappa2);
  //estimates.lambda = 0.5 * (estimates.kappa1 + estimates.kappa2);
  //estimates.lambda = 0.5 * sqrt(estimates.kappa1 * estimates.kappa2);
  //estimates.lambda = uniform_random() * sqrt(estimates.kappa1 * estimates.kappa2);

  return estimates;
}

void BVM_Sine::estimateParameters(std::vector<Vector> &data, Vector &weights)
{
  struct SufficientStatisticsSine suff_stats;
  computeSufficientStatisticsSine(data,suff_stats,weights);

  struct EstimatesSine estimates;
  struct EstimatesSine initial_est = computeInitialEstimates(suff_stats);

  string type;
  switch(ESTIMATION) {
    case MLE:
    {
      type = "MLE";
      struct EstimatesSine ml_est = initial_est;
      OptimizeSine opt_mle(type);
      opt_mle.initialize(
        ml_est.mu1,ml_est.mu2,ml_est.kappa1,ml_est.kappa2,ml_est.lambda
      );
      ml_est = opt_mle.minimize(suff_stats);
      estimates = ml_est;
      break;
    }

    case MAP:
    {
      type = "MAP";
      struct EstimatesSine map_est = initial_est;
      OptimizeSine opt_map(type);
      opt_map.initialize(
        map_est.mu1,map_est.mu2,map_est.kappa1,map_est.kappa2,map_est.lambda
      );
      map_est = opt_map.minimize(suff_stats);
      estimates = map_est;
      break;
    }

    case MML:
    {
      type = "MML";
      struct EstimatesSine mml_est = initial_est;
      OptimizeSine opt_mml(type);
      opt_mml.initialize(
        mml_est.mu1,mml_est.mu2,mml_est.kappa1,mml_est.kappa2,mml_est.lambda
      );
      mml_est = opt_mml.minimize(suff_stats);
      estimates = mml_est;
      break;
    }
  } // switch()
  
  updateParameters(estimates);
}

void BVM_Sine::updateParameters(struct EstimatesSine &estimates)
{
  mu1 = estimates.mu1;
  mu2 = estimates.mu2;
  kappa1 = estimates.kappa1;
  kappa2 = estimates.kappa2;
  lambda = estimates.lambda;
  assert(!boost::math::isnan(mu1));
  assert(!boost::math::isnan(mu2));
  assert(!boost::math::isnan(kappa1));
  assert(!boost::math::isnan(kappa2));
  assert(!boost::math::isnan(lambda));
  computeExpectation();
}

void BVM_Sine::printParameters(ostream &os)
{
  os << "[mus]: " << "(" << mu1*180/PI << ", " << mu2*180/PI << ")";
  os << "\t[kappas]: " << fixed << setprecision(3) 
     << "(" << kappa1 << ", " << kappa2 << ")";
  os << "\t[lambda]: " << fixed << setprecision(3) << lambda << endl;
}

double BVM_Sine::computeKLDivergence(BVM_Sine &other)
{
  if (computed != SET) {
    computeExpectation();
  }

  double log_norm_b = other.getLogNormalizationConstant();
  double ans = log_norm_b - constants.log_c;

  double mu1b = other.Mean1();
  double mu1_diff = mu1 - mu1b;
  double mu2b = other.Mean2();
  double mu2_diff = mu2 - mu2b;

  double kappa1b = other.Kappa1();
  double kappa2b = other.Kappa2();

  double kappa1_diff = kappa1 - kappa1b * cos(mu1_diff);
  ans += (kappa1_diff * constants.ck1_c);
  double kappa2_diff = kappa2 - kappa2b * cos(mu2_diff);
  ans += (kappa2_diff * constants.ck2_c);

  double lambdab = other.Lambda();
  double lambda_diff = lambda - (lambdab * cos(mu1_diff) * cos(mu2_diff));
  ans += (lambda_diff * constants.cl_c);

  double lambda_term = lambdab * sin(mu1_diff) * sin(mu2_diff);
  //cout << "lambda_term: " << lambda_term << endl;
  ans -= (lambda_term * constants.ck1k2_c);

  assert(ans >= 0);
  return ans/log(2);  // KL divergence (in bits)
}

/*double BVM_Sine::computeKLDivergence2(BVM_Sine &other)
{
  if (computed != SET) {
    computeExpectation();
  }

  double log_norm_b = other.getLogNormalizationConstant();
  double ans = log_norm_b - constants.log_c;

  double mu1b = other.Mean1();
  double mu1_diff = mu1 - mu1b;
  double mu2b = other.Mean2();
  double mu2_diff = mu2 - mu2b;

  double kappa1b = other.Kappa1();
  double kappa2b = other.Kappa2();

  double kappa1_diff = kappa1 - kappa1b * cos(mu1_diff);
  ans += (kappa1_diff * constants.ck1_c);
  double kappa2_diff = kappa2 - kappa2b * cos(mu2_diff);
  ans += (kappa2_diff * constants.ck2_c);

  double lambdab = other.Lambda();
  double lambda_diff = lambda - (lambdab * cos(mu1_diff) * cos(mu2_diff));
  ans += (lambda_diff * constants.cl_c);

  double lambda_term = lambdab * sin(mu1_diff) * sin(mu2_diff);
  //cout << "lambda_term: " << lambda_term << endl;
  ans += (lambda_term * constants.ck1k2_c);

  assert(ans >= 0);
  return ans/log(2);  // KL divergence (in bits)
}*/

double BVM_Sine::computeKLDivergence(struct EstimatesSine &estimates)
{
  BVM_Sine bvm_sine(
    estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.lambda
  );
  return computeKLDivergence(bvm_sine);
}

double BVM_Sine::computeKLDivergence(BVM_Sine &other, std::vector<Vector> &angle_pairs)
{
  double ans = 0;
  for (int i=0; i<angle_pairs.size(); i++) {
    double log_f = log_density(angle_pairs[i]);
    double log_g = other.log_density(angle_pairs[i]);
    ans += (log_f - log_g);
  } // for()
  //assert(ans >= 0);
  ans /= angle_pairs.size();
  return ans/log(2);  // KL divergence (in bits)
}

