#include "BVM_Cosine.h"
#include "vMC.h"
#include "Mixture_vMC.h"
#include "MarginalDensityCosine.h"
#include "KappaSolver.h"
#include "OptimizeCosine.h"

extern int ESTIMATION;

BVM_Cosine::BVM_Cosine()
{
  mu1 = 0; mu2 = 0;
  kappa1 = 1; kappa2 = 1;
  kappa3 = 0; // independent
  computed = UNSET;
  computed_lognorm = UNSET;
}

BVM_Cosine::BVM_Cosine(double kappa1, double kappa2, double kappa3) :
                       kappa1(kappa1), kappa2(kappa2), kappa3(kappa3)
{
  mu1 = 0; mu2 = 0;
  computed = UNSET;
  computed_lognorm = UNSET;
}

BVM_Cosine::BVM_Cosine(
  double mu1, double mu2, double kappa1, double kappa2, double kappa3
) : mu1(mu1), mu2(mu2), kappa1(kappa1), kappa2(kappa2), kappa3(kappa3)
{
  computed = UNSET;
  computed_lognorm = UNSET;
}

BVM_Cosine BVM_Cosine::operator=(const BVM_Cosine &source)
{
  if (this != &source) {
    mu1 = source.mu1;
    mu2 = source.mu2;
    kappa1 = source.kappa1;
    kappa2 = source.kappa2;
    kappa3 = source.kappa3;
    constants = source.constants;
    computed = source.computed;
    computed_lognorm = source.computed_lognorm;
  }
  return *this;
}

double BVM_Cosine::Mean1()
{
  return mu1;
}

double BVM_Cosine::Mean2()
{
  return mu2;
}

double BVM_Cosine::Kappa1()
{
  return kappa1;
}

double BVM_Cosine::Kappa2()
{
  return kappa2;
}

double BVM_Cosine::Kappa3()
{
  return kappa3;
}

// generates <theta1,theta2> pairs such that theta1,theta2 \in [0,2pi)
std::vector<Vector> BVM_Cosine::generate(int sample_size)
{
  std::vector<Vector> angle_pairs(sample_size);

  if (kappa3 != 0) {
    /* first, generate from marginal density */
    MarginalDensityCosine marginal(mu1,kappa1,kappa2,kappa3);
    
    /* check if marginal density is unimodal or bimodal */
    double kappa_diff = fabs(kappa2 - kappa3);
    double ratio_bessel = computeRatioBessel(kappa_diff);
    double rhs = (kappa_diff * kappa1) / (kappa2 * kappa3);

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
        double log_fg = accept_reject_fval_unimodal_marginal_cosine(
                          proposal_thetas[i],optimal_kappa,mu1,kappa1,kappa2,kappa3
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
        double log_fg = accept_reject_fval_bimodal_marginal_cosine(
                          proposal_thetas[i],optimal_kappa,mu1,kappa1,kappa2,kappa3,m1,m2
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
      double diff = thetas[i] - mu1;
      double num = -kappa3 * sin(diff);
      double denom = kappa2 - kappa3 * cos(diff);
      double beta = atan2(num,denom);
      double m = mu2 + beta;
      if (m < 0) m += (2 * PI);

      double k23sq = kappa2 * kappa2 + kappa3 * kappa3 
                     - (2 * kappa2 * kappa3 * cos(diff));
      double k23 = sqrt(k23sq);

      vMC vmc(m,k23);
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
  } // if (kappa3 != 0) {

  return angle_pairs;
}

std::vector<Vector> BVM_Cosine::generate_cartesian(int sample_size)
{
  std::vector<Vector> angle_pairs = generate(sample_size);

  std::vector<Vector> random_sample = generate_cartesian(angle_pairs);

  return random_sample;
}

// convert theta,phi --> cartesian coordinates
std::vector<Vector> BVM_Cosine::generate_cartesian(
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

BVM_Cosine::Constants BVM_Cosine::getConstants()
{
  if (computed != SET) {
    computeExpectation();
  }
  return constants;
}

void BVM_Cosine::computeExpectation()
{
  computeConstants();

  computed = SET;
}

void BVM_Cosine::computeConstants()
{
  constants.log_c = computeLogNormalizationConstant();
}

double BVM_Cosine::getLogNormalizationConstant()
{
  if (computed_lognorm == UNSET) {
    constants.log_c = computeLogNormalizationConstant();  
  }
  return constants.log_c;
}

/*!
 *  logarithm of normalization constant
 */
/*
double BVM_Cosine::computeLogNormalizationConstant()
{
  double log_bessel_sum = computeLogModifiedBesselFirstKind(0,kappa1)
                          + computeLogModifiedBesselFirstKind(0,kappa2)
                          + computeLogModifiedBesselFirstKind(0,kappa3);
  int j = 1;
  double log_f1 = computeLogModifiedBesselFirstKind(j,kappa1)
                  + computeLogModifiedBesselFirstKind(j,kappa2)
                  + computeLogModifiedBesselFirstKind(j,kappa3);
  double series_sum = 1;

  while (1) {
    j++;
    double log_fj = computeLogModifiedBesselFirstKind(j,kappa1)
                    + computeLogModifiedBesselFirstKind(j,kappa2)
                    + computeLogModifiedBesselFirstKind(j,kappa3);
    double log_diff = log_fj - log_f1;
    double t_j = exp(log_diff);
    series_sum += t_j;
    if (t_j/series_sum <= 1e-6) {
      cout << "log_bessel_sum: " << log_bessel_sum << endl;
      cout << "j: " << j << "; series_sum: " << series_sum << endl;
      cout << "log_diff: " << log_diff << endl;
      cout << "t_j: " << t_j << endl;
      cout << "t_j/series_sum: " << t_j/series_sum << endl << endl;
      break;
    } // if()
  } // while()
  double log_sum_from_f1 = log_f1 + log(series_sum);
  double log_tmp = log_sum_from_f1 + log(2) - log_bessel_sum;
  double tmp2 = 1 + exp(log_tmp);
  double log_const = 2 * log(2*PI) + log_bessel_sum + log(tmp2);
  computed_lognorm = SET;
  return log_const;
}*/

double BVM_Cosine::computeLogNormalizationConstant()
{
  double log_bessel_sum = computeLogModifiedBesselFirstKind(0,kappa1)
                          + computeLogModifiedBesselFirstKind(0,kappa2)
                          + computeLogModifiedBesselFirstKind(0,kappa3);

  double log_sum_odd = compute_series_partial(1);
  double log_sum_even = compute_series_partial(2);

  assert(log_sum_odd > log_sum_even);
  double ratio = exp(log_sum_even - log_sum_odd);
  double log_tmp = log(2) - log_bessel_sum + log_sum_odd + log(1-ratio);
  double tmp2 = 1 - exp(log_tmp);

  double log_const = 2 * log(2*PI) + log_bessel_sum + log(tmp2);
  computed_lognorm = SET;
  return log_const;
}

double BVM_Cosine::compute_series_partial(int begin)
{
  int j = begin;
  double log_f1 = computeLogModifiedBesselFirstKind(j,kappa1)
                  + computeLogModifiedBesselFirstKind(j,kappa2)
                  + computeLogModifiedBesselFirstKind(j,kappa3);
  double series_sum = 1;

  while (1) {
    j += 2;
    double log_fj = computeLogModifiedBesselFirstKind(j,kappa1)
                    + computeLogModifiedBesselFirstKind(j,kappa2)
                    + computeLogModifiedBesselFirstKind(j,kappa3);
    double log_diff = log_fj - log_f1;
    double t_j = exp(log_diff);
    series_sum += t_j;
    if (t_j/series_sum <= 1e-15) {
      /*cout << "log_bessel_sum: " << log_bessel_sum << endl;
      cout << "j: " << j << "; series_sum: " << series_sum << endl;
      cout << "log_diff: " << log_diff << endl;
      cout << "t_j: " << t_j << endl;
      cout << "t_j/series_sum: " << t_j/series_sum << endl << endl;*/
      break;
    } // if()
  } // while()
  double log_sum_from_f1 = log_f1 + log(series_sum);
  return log_sum_from_f1;
}

double BVM_Cosine::computeLogParametersProbability(double Neff)
{
  double log_prior_density = computeLogParametersPriorDensity();
  double log_expected_fisher = computeLogFisherInformation(Neff);
  double logp = -log_prior_density + 0.5 * log_expected_fisher;
  return logp;
}

// prior density of parameters ... h(kappa1,kappa2,kappa3)
double BVM_Cosine::computeLogParametersPriorDensity()
{
  // prior on means (uniform)
  double log_prior_means = -2 * log(2*PI);

  // joint prior on kappas (vMC)
  double log_scale = log(kappa1 + kappa2);
  log_scale -= (1.5 * log(1+kappa1*kappa1));
  log_scale -= (1.5 * log(1+kappa2*kappa2));

  double log_prior = log_prior_means + log_scale;
  return log_prior;
}

// prior density of parameters ... h(kappa1,kappa2,rho)
double BVM_Cosine::computeLogParametersPriorDensityTransform()
{
  // prior on means (uniform)
  double log_prior_means = -2 * log(2*PI);

  // joint prior on kappas and rho (correlation) 
  double log_scale = log(kappa1) - 1.5 * log(1+kappa1*kappa1);
  log_scale += (log(kappa2) - 1.5 * log(1+kappa2*kappa2));

  double log_prior = log_prior_means + log_scale;
  return log_prior;
}

double BVM_Cosine::computeLogFisherInformation(double N)
{
  double log_fisher = computeLogFisherInformation_Single(); 
  log_fisher += (5 * log(N));
  assert(!boost::math::isnan(log_fisher));
  return log_fisher;
}

double BVM_Cosine::computeLogFisherInformation_Single()
{
  if (computed == UNSET) {
    computeExpectation();
  }
  double log_det_axes = computeLogFisherAxes();
  double log_det_scale = computeLogFisherScale();
  return log_det_axes + log_det_scale; 
}

double BVM_Cosine::computeLogFisherAxes()
{
  // L: -ve log-likelihood
/*
  // E [d^2 L / d mu1^2]
  double t1 = kappa1 * constants.ck1_c + lambda * constants.cl_c;

  // E [d^2 L / d mu2^2]
  double t2 = kappa2 * constants.ck2_c + lambda * constants.cl_c;

  // E [d^2 L / dmu1 dmu2]
  double t3 = -lambda * constants.ck1k2_c;

  double det = fabs(t1 * t2 - t3 * t3);
  //assert(det > 0);
  //cout << "det: " << det << endl;
  return log(det);*/
  return 0;
}

double BVM_Cosine::computeLogFisherScale()
{
  Matrix fisher_scale = ZeroMatrix(3,3);
/*
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
*/
  double det = determinant_3d(fisher_scale);
  //cout << "det: " << det << endl; //exit(1);
  return log(fabs(det));
}

double BVM_Cosine::log_density(Vector &angle_pair)
{
  assert(angle_pair.size() == 2);
  return log_density(angle_pair[0],angle_pair[1]);
}

// log(pdf)
double BVM_Cosine::log_density(double &theta1, double &theta2)
{
  if (computed != SET) {
    computeExpectation();
  }
  double ans = 0;
  ans -= constants.log_c;
  ans += (kappa1 * cos(theta1-mu1));
  ans += (kappa2 * cos(theta2-mu2));
  ans -= (kappa3 * cos(theta1-theta2-mu1+mu2));
  return ans;
}

// data = angle_pairs
double BVM_Cosine::computeNegativeLogLikelihood(std::vector<Vector> &data)
{
  struct SufficientStatisticsCosine suff_stats;
  computeSufficientStatisticsCosine(data,suff_stats);
  return computeNegativeLogLikelihood(suff_stats);
}

double BVM_Cosine::computeNegativeLogLikelihood(
  struct SufficientStatisticsCosine &suff_stats
) {
  double cosm1 = cos(mu1);
  double sinm1 = sin(mu1);
  double cosm2 = cos(mu2);
  double sinm2 = sin(mu2);
  double cosm1_m2 = cos(mu1-mu2);
  double sinm1_m2 = sin(mu1-mu2);

  // compute log-likelihood
  double ans = -suff_stats.N * getLogNormalizationConstant();
  ans += (kappa1 * cosm1 * suff_stats.cost1);
  ans += (kappa1 * sinm1 * suff_stats.sint1);
  ans += (kappa2 * cosm2 * suff_stats.cost2);
  ans += (kappa2 * sinm2 * suff_stats.sint2);

  ans -= (kappa3 * cosm1_m2 * suff_stats.cost1_t2);
  ans -= (kappa3 * sinm1_m2 * suff_stats.sint1_t2);

  // return -ve log-likelihood
  return -ans;
}

double BVM_Cosine::computeNegativeLogLikelihood(
  struct EstimatesCosine &estimates, struct SufficientStatisticsCosine &suff_stats
) {
  BVM_Cosine bvm_cosine(
    estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.kappa3
  );
  return bvm_cosine.computeNegativeLogLikelihood(suff_stats);
}

// data = angle_pairs
double BVM_Cosine::computeMessageLength(std::vector<Vector> &data)
{
  struct SufficientStatisticsCosine suff_stats;
  computeSufficientStatisticsCosine(data,suff_stats);
  return computeMessageLength(suff_stats);
}

double BVM_Cosine::computeMessageLength(
  struct SufficientStatisticsCosine &suff_stats
) {
  double log_prior = computeLogParametersPriorDensity();
  double log_fisher = computeLogFisherInformation(suff_stats.N);
  double part1 = -6.455 - log_prior + 0.5 * log_fisher;
  double part2 = computeNegativeLogLikelihood(suff_stats) + 2.5
                 - 2 * suff_stats.N * log(AOM);
  double msglen = part1 + part2;
  return msglen/log(2);
}

double BVM_Cosine::computeMessageLength(
  struct EstimatesCosine &estimates,
  struct SufficientStatisticsCosine &suff_stats
) {
  BVM_Cosine bvm_cosine(
    estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.kappa3
  );
  return bvm_cosine.computeMessageLength(suff_stats);
}

// data = angle_pairs
void BVM_Cosine::computeAllEstimators(
  std::vector<Vector> &data, 
  std::vector<struct EstimatesCosine> &all_estimates,
  int verbose,
  int compute_kldiv
) {
  struct SufficientStatisticsCosine suff_stats;
  computeSufficientStatisticsCosine(data,suff_stats);

  computeAllEstimators(
    data,suff_stats,all_estimates,verbose,compute_kldiv
  );
}

void BVM_Cosine::computeAllEstimators(
  std::vector<Vector> &data, 
  struct SufficientStatisticsCosine &suff_stats,
  std::vector<struct EstimatesCosine> &all_estimates,
  int verbose,
  int compute_kldiv
) {
  double msglen,negloglike,kldiv,min_msg;
  int min_index;

  all_estimates.clear();
  all_estimates = std::vector<struct EstimatesCosine>(NUM_METHODS);

  string type = "initial";
  struct EstimatesCosine initial_est = computeInitialEstimates(suff_stats);
  print(type,initial_est);

  type = "MLE";
  struct EstimatesCosine ml_est = initial_est;
  OptimizeCosine opt_mle(type);
  opt_mle.initialize(
    ml_est.mu1,ml_est.mu2,ml_est.kappa1,ml_est.kappa2,ml_est.kappa3
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
  struct EstimatesCosine map_est = initial_est;
  OptimizeCosine opt_map(type);
  opt_map.initialize(
    map_est.mu1,map_est.mu2,map_est.kappa1,map_est.kappa2,map_est.kappa3
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
  struct EstimatesCosine map2_est = initial_est;
  OptimizeCosine opt_map2(type);
  opt_map2.initialize(
    map2_est.mu1,map2_est.mu2,map2_est.kappa1,map2_est.kappa2,map2_est.kappa3
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
  //struct EstimatesCosine mml_est = initial_est;
  struct EstimatesCosine mml_est = all_estimates[min_index];
  OptimizeCosine opt_mml(type);
  opt_mml.initialize(
    mml_est.mu1,mml_est.mu2,mml_est.kappa1,mml_est.kappa2,mml_est.kappa3
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
}

struct EstimatesCosine BVM_Cosine::computeInitialEstimates(
  struct SufficientStatisticsCosine &suff_stats
) {
  struct EstimatesCosine estimates;
  
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
  estimates.kappa3 = estimates.rho * (estimates.kappa1 * estimates.kappa2) / (estimates.kappa1 + estimates.kappa2);

  return estimates;
}

void BVM_Cosine::estimateParameters(std::vector<Vector> &data, Vector &weights)
{
  struct SufficientStatisticsCosine suff_stats;
  computeSufficientStatisticsCosine(data,suff_stats,weights);

  struct EstimatesCosine estimates;
  struct EstimatesCosine initial_est = computeInitialEstimates(suff_stats);

  string type;
  switch(ESTIMATION) {
    case MLE:
    {
      type = "MLE";
      struct EstimatesCosine ml_est = initial_est;
      OptimizeCosine opt_mle(type);
      opt_mle.initialize(
        ml_est.mu1,ml_est.mu2,ml_est.kappa1,ml_est.kappa2,ml_est.kappa3
      );
      ml_est = opt_mle.minimize(suff_stats);
      estimates = ml_est;
      break;
    }

    case MAP:
    {
      type = "MAP";
      struct EstimatesCosine map_est = initial_est;
      OptimizeCosine opt_map(type);
      opt_map.initialize(
        map_est.mu1,map_est.mu2,map_est.kappa1,map_est.kappa2,map_est.kappa3
      );
      map_est = opt_map.minimize(suff_stats);
      estimates = map_est;
      break;
    }

    case MML:
    {
      type = "MML";
      struct EstimatesCosine mml_est = initial_est;
      OptimizeCosine opt_mml(type);
      opt_mml.initialize(
        mml_est.mu1,mml_est.mu2,mml_est.kappa1,mml_est.kappa2,mml_est.kappa3
      );
      mml_est = opt_mml.minimize(suff_stats);
      estimates = mml_est;
      break;
    }
  } // switch()
  
  updateParameters(estimates);
}

void BVM_Cosine::updateParameters(struct EstimatesCosine &estimates)
{
  mu1 = estimates.mu1;
  mu2 = estimates.mu2;
  kappa1 = estimates.kappa1;
  kappa2 = estimates.kappa2;
  kappa3 = estimates.kappa3;
  assert(!boost::math::isnan(mu1));
  assert(!boost::math::isnan(mu2));
  assert(!boost::math::isnan(kappa1));
  assert(!boost::math::isnan(kappa2));
  assert(!boost::math::isnan(kappa3));
  computeExpectation();
}

void BVM_Cosine::printParameters(ostream &os)
{
  os << "[mus]: " << "(" << mu1*180/PI << ", " << mu2*180/PI << ")";
  os << "\t[kappas]: " << fixed << setprecision(3) 
     << "(" << kappa1 << ", " << kappa2 << ")";
  os << "\t[kappa3]: " << fixed << setprecision(3) << kappa3 << endl;
}

double BVM_Cosine::computeKLDivergence(BVM_Cosine &other)
{
  if (computed != SET) {
    computeExpectation();
  }
/*
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
  ans += (lambda_term * constants.ck1k2_c);

  assert(ans >= 0);
  return ans/log(2);  // KL divergence (in bits)
*/
}

double BVM_Cosine::computeKLDivergence(struct EstimatesCosine &estimates)
{
  BVM_Cosine bvm_cosine(
    estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.kappa3
  );
  return computeKLDivergence(bvm_cosine);
}

