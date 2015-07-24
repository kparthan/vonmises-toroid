#include "BVM_Sine.h"
#include "vMC.h"
#include "Mixture_vMC.h"
#include "MarginalDensitySine.h"
#include "KappaSolver.h"
#include "OptimizeSine.h"

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
  mu1 = 0; mu2 = 0;
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
      double lambda_sine = lambda * sin(thetas[i] - mu1);
      double beta = atan2(lambda_sine,kappa2);
      double m = mu2 + beta;
      if (m < 0) m += (2 * PI);

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
  double r1 = 2;
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
  // const = (lambda^2) / (2 * k1 * k2)
  double log_const = 2 * log(fabs(lambda)) - log(2) - log(kappa1) - log(kappa2);

  double log_bessel1_prev = computeLogModifiedBesselFirstKind(0,kappa1);
  double log_bessel2_prev = computeLogModifiedBesselFirstKind(0,kappa2);
  double log_f0 = log_bessel1_prev + log_bessel2_prev;
  double log_fj_prev = log_f0;
  double j = 0;
  double series_sum = 1;
  while (1) {
    double log_bessel1_current = computeLogModifiedBesselFirstKind(j+1,kappa1);
    double log_bessel2_current = computeLogModifiedBesselFirstKind(j+1,kappa2);
    double log_fj_current = log_const + log(2*j+1) - log(j+1)
                            + log_bessel1_current + log_bessel2_current
                            - log_bessel1_prev - log_bessel2_prev;
                            + log_fj_prev;
    double log_diff = log_fj_current - log_f0;
    double current = exp(log_diff); // tj
    series_sum += current;
    if (current/series_sum <= 1e-6) {
      /*cout << "j: " << j << "; series_sum: " << series_sum << endl;
      cout << "log_const: " << log_const << endl;
      cout << "log_fj_prev: " << log_fj_prev << endl;
      cout << "log_fj_current: " << log_fj_current << endl;
      cout << "log_bessel1_prev: " << log_bessel1_prev << endl;
      cout << "log_bessel2_prev: " << log_bessel2_prev << endl;
      cout << "log_bessel1_current: " << log_bessel1_current << endl;
      cout << "log_bessel2_current: " << log_bessel2_current << endl;
      cout << "current: " << current << endl;*/
      break;
    }
    j++;
    log_bessel1_prev = log_bessel1_current; 
    log_bessel2_prev = log_bessel2_current; 
    log_fj_prev = log_fj_current;
  } // while(1)
  //cout << "j: " << j << endl;
  double ans = 2*log(2*PI) + log_f0 + log(series_sum);
  //ans += computeLogModifiedBesselFirstKind(0,kappa1);
  //ans += computeLogModifiedBesselFirstKind(0,kappa2);
  return ans;
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
  if (computed != SET) {
    computeExpectation();
  }

  double cosm1 = cos(mu1);
  double sinm1 = sin(mu1);
  double cosm2 = cos(mu2);
  double sinm2 = sin(mu2);

  // compute log-likelihood
  struct SufficientStatisticsSine suff_stats;
  computeSufficientStatisticsSine(data,suff_stats);
  double ans = -data.size() * constants.log_c;
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
  if (computed != SET) {
    computeExpectation();
  }

  double cosm1 = cos(estimates.mu1);
  double sinm1 = sin(estimates.mu1);
  double cosm2 = cos(estimates.mu2);
  double sinm2 = sin(estimates.mu2);

  // compute log-likelihood
  double ans = -suff_stats.N * constants.log_c;
  ans += (estimates.kappa1 * cosm1 * suff_stats.cost1);
  ans += (estimates.kappa1 * sinm1 * suff_stats.sint1);
  ans += (estimates.kappa2 * cosm2 * suff_stats.cost2);
  ans += (estimates.kappa2 * sinm2 * suff_stats.sint2);

  ans += (estimates.lambda * cosm1 * cosm2 * suff_stats.sint1sint2);
  ans -= (estimates.lambda * cosm1 * sinm2 * suff_stats.sint1cost2);
  ans -= (estimates.lambda * sinm1 * cosm2 * suff_stats.cost1sint2);
  ans += (estimates.lambda * sinm1 * sinm2 * suff_stats.cost1cost2);

  // return -ve log-likelihood
  return -ans;
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
    suff_stats,all_estimates,verbose,compute_kldiv
  );
}

void BVM_Sine::computeAllEstimators(
  struct SufficientStatisticsSine &suff_stats,
  std::vector<struct EstimatesSine> &all_estimates,
  int verbose,
  int compute_kldiv
) {
  double msglen,negloglike,kldiv,min_msg;
  int min_index;

  all_estimates.clear();

  string type = "initial";
  struct EstimatesSine initial_est = computeInitialEstimates(suff_stats);
  //print(type,initial_est);

  type = "MLE";
  struct EstimatesSine ml_est = initial_est;
  OptimizeSine opt1(type);
  opt1.initialize(
    ml_est.mu1,ml_est.mu2,ml_est.kappa1,ml_est.kappa2,ml_est.lambda
  );
  ml_est = opt1.minimize(suff_stats);
  //ml_est.msglen = computeMessageLength(ml_est,suff_stats);
  ml_est.negloglike = computeNegativeLogLikelihood(ml_est,suff_stats);
  /*if (compute_kldiv) {
    ml_est.kldiv = computeKLDivergence(ml_est);
  }*/
  if (verbose) {
    print(type,ml_est);
    //cout << fixed << "msglen: " << ml_est.msglen << endl;
    //cout << "negloglike: " << ml_est.negloglike << endl;
    //cout << "KL-divergence: " << ml_est.kldiv << endl << endl;
  }
  all_estimates.push_back(ml_est);
  /*if (ml_est.msglen < min_msg) {
    min_index = MLE;
    min_msg = ml_est.msglen;
  }*/

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
  cout << "rbar: " << rbar << endl;
  double constant = cos(estimates.mu1) * suff_stats.cost1 
                    + sin(estimates.mu1) * suff_stats.sint1;
  KappaSolver solver1(suff_stats.N,constant,rbar);
  estimates.kappa1 = solver1.minimize();
  cout << "kappa1_init: " << estimates.kappa1 << endl;

  cos_sum = suff_stats.cost2;
  sin_sum = suff_stats.sint2;
  rbar = (sqrt(cos_sum * cos_sum + sin_sum * sin_sum))/suff_stats.N;
  cout << "rbar: " << rbar << endl;
  constant = cos(estimates.mu2) * suff_stats.cost2 
             + sin(estimates.mu2) * suff_stats.sint2;
  KappaSolver solver2(suff_stats.N,constant,rbar);
  estimates.kappa2 = solver2.minimize();
  cout << "kappa2_init: " << estimates.kappa2 << endl;

  //estimates.lambda = 0.5 * (estimates.kappa1 + estimates.kappa2);
  estimates.lambda = uniform_random() * sqrt(estimates.kappa1 * estimates.kappa2);

  return estimates;
}

