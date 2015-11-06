#include "BVM_Ind.h"
#include "vMC.h"
#include "KappaSolver.h"
#include "OptimizeInd.h"

extern int ESTIMATION;

BVM_Ind::BVM_Ind()
{
  mu1 = 0; mu2 = 0;
  kappa1 = 1; kappa2 = 1;
  computed = UNSET;
}

BVM_Ind::BVM_Ind(double kappa1, double kappa2) : kappa1(kappa1), kappa2(kappa2)
{
  mu1 = 0; mu2 = 0;
  computed = UNSET;
}

BVM_Ind::BVM_Ind(
  double mu1, double mu2, double kappa1, double kappa2
) : mu1(mu1), mu2(mu2), kappa1(kappa1), kappa2(kappa2)
{
  computed = UNSET;
}

BVM_Ind BVM_Ind::operator=(const BVM_Ind &source)
{
  if (this != &source) {
    mu1 = source.mu1;
    mu2 = source.mu2;
    kappa1 = source.kappa1;
    kappa2 = source.kappa2;
    log_c = source.log_c;
    computed = source.computed;
  }
  return *this;
}

double BVM_Ind::Mean1()
{
  return mu1;
}

double BVM_Ind::Mean2()
{
  return mu2;
}

double BVM_Ind::Kappa1()
{
  return kappa1;
}

double BVM_Ind::Kappa2()
{
  return kappa2;
}

// generates <theta1,theta2> pairs such that theta1,theta2 \in [0,2pi)
std::vector<Vector> BVM_Ind::generate(int sample_size)
{
  std::vector<Vector> angle_pairs(sample_size);

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

  return angle_pairs;
}

std::vector<Vector> BVM_Ind::generate_cartesian(int sample_size)
{
  std::vector<Vector> angle_pairs = generate(sample_size);

  std::vector<Vector> random_sample = generate_cartesian(angle_pairs);

  return random_sample;
}

// convert theta,phi --> cartesian coordinates
std::vector<Vector> BVM_Ind::generate_cartesian(
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

double BVM_Ind::getLogNormalizationConstant()
{
  if (computed == UNSET) {
    log_c = computeLogNormalizationConstant();  
  }
  return log_c;
}

/*!
 *  logarithm of normalization constant
 */
double BVM_Ind::computeLogNormalizationConstant()
{
  double log_bessel1 = computeLogModifiedBesselFirstKind(0,kappa1);
  double log_bessel2 = computeLogModifiedBesselFirstKind(0,kappa2);
  log_c = 2 * log(2*PI) + log_bessel1 + log_bessel2;
  computed = SET;
  return log_c;
}

double BVM_Ind::computeLogParametersProbability(double Neff)
{
  double log_prior_density = computeLogParametersPriorDensity();
  double log_expected_fisher = computeLogFisherInformation(Neff);
  double logp = -log_prior_density + 0.5 * log_expected_fisher;
  return logp;
}

// prior density of parameters ... h(kappa1,kappa2)
double BVM_Ind::computeLogParametersPriorDensity()
{
  // prior on means (uniform)
  double log_prior_means = -2 * log(2*PI);

  // joint prior on kappas (vMC)
  double log_scale = log(kappa1) - 1.5 * log(1+kappa1*kappa1);
  log_scale += (log(kappa2) - 1.5 * log(1+kappa2*kappa2));

  double log_prior = log_prior_means + log_scale;
  return log_prior;
}

double BVM_Ind::computeLogFisherInformation(double N)
{
  double log_fisher = computeLogFisherInformation_Single(); 
  log_fisher += (4 * log(N));
  assert(!boost::math::isnan(log_fisher));
  return log_fisher;
}

double BVM_Ind::computeLogFisherInformation_Single()
{
  double log_a2k = computeLogRatioBessel(2,kappa1);
  double a2k = exp(log_a2k);
  double a2k_der = computeDerivativeOfRatioBessel(kappa1,a2k);
  double log_f1 = log(kappa1) + log_a2k + log(a2k_der)

  log_a2k = computeLogRatioBessel(2,kappa2);
  a2k = exp(log_a2k);
  a2k_der = computeDerivativeOfRatioBessel(kappa2,a2k);
  double log_f2 = log(kappa2) + log_a2k + log(a2k_der)
 
  double log_f = log_f1 + log_f2;
  return log_f;
}

double BVM_Ind::log_density(Vector &angle_pair)
{
  assert(angle_pair.size() == 2);
  return log_density(angle_pair[0],angle_pair[1]);
}

// log(pdf)
double BVM_Ind::log_density(double &theta1, double &theta2)
{
  if (computed != SET) {
    computeLogNormalizationConstant();
  }
  double ans = 0;
  ans -= constants.log_c;
  ans += (kappa1 * cos(theta1-mu1));
  ans += (kappa2 * cos(theta2-mu2));
  return ans;
}

// data = angle_pairs
double BVM_Ind::computeNegativeLogLikelihood(std::vector<Vector> &data)
{
  struct SufficientStatisticsInd suff_stats;
  computeSufficientStatisticsInd(data,suff_stats);
  return computeNegativeLogLikelihood(suff_stats);
}

double BVM_Ind::computeNegativeLogLikelihood(
  struct SufficientStatisticsInd &suff_stats
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

  // return -ve log-likelihood
  return -ans;
}

double BVM_Ind::computeNegativeLogLikelihood(
  struct EstimatesInd &estimates, struct SufficientStatisticsInd &suff_stats
) {
  BVM_Ind bvm_ind(
    estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2
  );
  return bvm_ind.computeNegativeLogLikelihood(suff_stats);
}

// data = angle_pairs
double BVM_Ind::computeMessageLength(std::vector<Vector> &data)
{
  struct SufficientStatisticsInd suff_stats;
  computeSufficientStatisticsInd(data,suff_stats);
  return computeMessageLength(suff_stats);
}

double BVM_Ind::computeMessageLength(
  struct SufficientStatisticsInd &suff_stats
) {
  double log_prior = computeLogParametersPriorDensity();
  double log_fisher = computeLogFisherInformation(suff_stats.N);
  double part1 = -5.138 - log_prior + 0.5 * log_fisher;
  double part2 = computeNegativeLogLikelihood(suff_stats) + 2
                 - 2 * suff_stats.N * log(AOM);
  double msglen = part1 + part2;
  return msglen/log(2);
}

double BVM_Ind::computeMessageLength(
  struct EstimatesInd &estimates,
  struct SufficientStatisticsInd &suff_stats
) {
  BVM_Ind bvm_ind(
    estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2
  );
  return bvm_ind.computeMessageLength(suff_stats);
}

// data = angle_pairs
void BVM_Ind::computeAllEstimators(
  std::vector<Vector> &data, 
  std::vector<struct EstimatesInd> &all_estimates,
  int verbose,
  int compute_kldiv
) {
  struct SufficientStatisticsInd suff_stats;
  computeSufficientStatisticsInd(data,suff_stats);

  computeAllEstimators(
    data,suff_stats,all_estimates,verbose,compute_kldiv
  );
}

void BVM_Ind::computeAllEstimators(
  std::vector<Vector> &data, 
  struct SufficientStatisticsInd &suff_stats,
  std::vector<struct EstimatesInd> &all_estimates,
  int verbose,
  int compute_kldiv
) {
  double msglen,negloglike,kldiv,min_msg;
  int min_index;

  all_estimates.clear();
  all_estimates = std::vector<struct EstimatesInd>(NUM_METHODS);

  string type = "initial";
  struct EstimatesInd initial_est = computeInitialEstimates(suff_stats);
  print(type,initial_est);

  type = "MLE";
  struct EstimatesInd ml_est = initial_est;
  OptimizeInd opt_mle(type);
  opt_mle.initialize(
    ml_est.mu1,ml_est.mu2,ml_est.kappa1,ml_est.kappa2
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
  struct EstimatesInd map_est = initial_est;
  OptimizeInd opt_map(type);
  opt_map.initialize(
    map_est.mu1,map_est.mu2,map_est.kappa1,map_est.kappa2
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

  type = "MML";
  //struct EstimatesInd mml_est = initial_est;
  struct EstimatesInd mml_est = all_estimates[min_index];
  OptimizeInd opt_mml(type);
  opt_mml.initialize(
    mml_est.mu1,mml_est.mu2,mml_est.kappa1,mml_est.kappa2
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

struct EstimatesInd BVM_Ind::computeInitialEstimates(
  struct SufficientStatisticsInd &suff_stats
) {
  struct EstimatesInd estimates;
  
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

  return estimates;
}

void BVM_Ind::estimateParameters(std::vector<Vector> &data, Vector &weights)
{
  struct SufficientStatisticsInd suff_stats;
  computeSufficientStatisticsInd(data,suff_stats,weights);

  struct EstimatesInd estimates;
  struct EstimatesInd initial_est = computeInitialEstimates(suff_stats);

  string type;
  switch(ESTIMATION) {
    case MLE:
    {
      type = "MLE";
      struct EstimatesInd ml_est = initial_est;
      OptimizeInd opt_mle(type);
      opt_mle.initialize(
        ml_est.mu1,ml_est.mu2,ml_est.kappa1,ml_est.kappa2
      );
      ml_est = opt_mle.minimize(suff_stats);
      estimates = ml_est;
      break;
    }

    case MAP:
    {
      type = "MAP";
      struct EstimatesInd map_est = initial_est;
      OptimizeInd opt_map(type);
      opt_map.initialize(
        map_est.mu1,map_est.mu2,map_est.kappa1,map_est.kappa2
      );
      map_est = opt_map.minimize(suff_stats);
      estimates = map_est;
      break;
    }

    case MML:
    {
      type = "MML";
      struct EstimatesInd mml_est = initial_est;
      OptimizeInd opt_mml(type);
      opt_mml.initialize(
        mml_est.mu1,mml_est.mu2,mml_est.kappa1,mml_est.kappa2
      );
      mml_est = opt_mml.minimize(suff_stats);
      estimates = mml_est;
      break;
    }
  } // switch()
  
  updateParameters(estimates);
}

void BVM_Ind::updateParameters(struct EstimatesInd &estimates)
{
  mu1 = estimates.mu1;
  mu2 = estimates.mu2;
  kappa1 = estimates.kappa1;
  kappa2 = estimates.kappa2;
  assert(!boost::math::isnan(mu1));
  assert(!boost::math::isnan(mu2));
  assert(!boost::math::isnan(kappa1));
  assert(!boost::math::isnan(kappa2));
  computeLogNormalizationConstant();
}

void BVM_Ind::printParameters(ostream &os)
{
  os << "[mus]: " << "(" << mu1*180/PI << ", " << mu2*180/PI << ")";
  os << "\t[kappas]: " << fixed << setprecision(3) 
     << "(" << kappa1 << ", " << kappa2 << ")\n";
}

double BVM_Ind::computeKLDivergence(BVM_Ind &other)
{
  if (computed != SET) {
    computeLogNormalizationConstant();
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

  assert(ans >= 0);
  return ans/log(2);  // KL divergence (in bits)
}

double BVM_Ind::computeKLDivergence(struct EstimatesInd &estimates)
{
  BVM_Ind bvm_ind(
    estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2
  );
  return computeKLDivergence(bvm_ind);
}

double BVM_Ind::computeKLDivergence(BVM_Ind &other, std::vector<Vector> &angle_pairs)
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

