#include "Test.h"
#include "Support.h"
#include "Normal.h"
#include "vMC.h"
#include "Mixture_vMC.h"
#include "BVM_Sine.h"
#include "BVM_Cosine.h"

extern int ENABLE_DATA_PARALLELISM;
extern int NUM_THREADS;
extern int ESTIMATION,CRITERION;

void Test::bessel()
{
  double d,k,log_bessel;
  
  d = 3;
  k = 730;

  log_bessel = computeLogModifiedBesselFirstKind(d,k);
  cout << "log I(" << d << "," << k << "): " << log_bessel << endl;

  log_bessel = log(cyl_bessel_i(d,k));
  cout << "Boost log I(" << d << "," << k << "): " << log_bessel << endl;
}

void Test::generate_vmc()
{
  int N = 100;
  double mu = 210; // in degrees
  double kappa = 100;

  mu *= PI / 180;
  vMC vmc(mu,kappa);
  std::vector<Vector> random_sample = vmc.generate_cartesian(N);

  writeToFile("vmc.dat",random_sample);
}

void Test::generate_mix_vmc()
{
  std::vector<vMC> components;
  int N = 100;
  double mu = 210; // in degrees
  double kappa = 100;

  mu *= PI / 180;
  vMC vmc1(mu,kappa);
  components.push_back(vmc1);

  mu = 290;
  mu *= PI / 180;
  vMC vmc2(mu,kappa);
  components.push_back(vmc2);

  mu = 90; kappa = 1000;
  mu *= PI / 180;
  vMC vmc3(mu,kappa);
  components.push_back(vmc3);

  Vector weights(3,1.0/3);
  
  Mixture_vMC mix(3,weights,components);
  std::vector<Vector> random_sample = mix.generate_cartesian(N);

  writeToFile("mix_vmc.dat",random_sample);
}

/* sine model related */
void Test::generate_bvm_sine()
{
  double mu1,mu2,kappa1,kappa2,lambda;
  int N = 5000;

  mu1 = 88.1434; mu2 = 82.768; kappa1 = 64.9647; kappa2 = 7.57064; 
  lambda = 0.000195011;

  mu1 *= PI/180; mu2 *= PI/180;
  BVM_Sine bvm_sine(mu1,mu2,kappa1,kappa2,lambda);

  std::vector<Vector> angle_pairs = bvm_sine.generate(N);
  writeToFile("angle_pairs.dat",angle_pairs);

  std::vector<Vector> random_sample = bvm_sine.generate_cartesian(angle_pairs);
  writeToFile("bvm_sine.dat",random_sample);
}

void Test::bvm_sine_normalization_constant()
{
  double mu1,mu2,kappa1,kappa2,lambda;
  int N = 1000;

  mu1 = 90; mu2 = 90; kappa1 = 100; kappa2 = 20; lambda = -30;

  mu1 *= PI/180; mu2 *= PI/180;
  BVM_Sine bvm_sine(mu1,mu2,kappa1,kappa2,lambda);
  double log_norm = bvm_sine.computeLogNormalizationConstant();
  cout << "log_norm: " << log_norm << endl;
}

void Test::bvm_sine_constants()
{
  double mu1,mu2,kappa1,kappa2,lambda;
  int N = 1000;

  mu1 = 90; mu2 = 90; kappa1 = 100; kappa2 = 20; lambda = -30;

  mu1 *= PI/180; mu2 *= PI/180;
  BVM_Sine bvm_sine(mu1,mu2,kappa1,kappa2,lambda);

  BVM_Sine::Constants constants = bvm_sine.getConstants();
  cout << "log_norm: " << constants.log_c << endl;
  cout << "log_dc_dk1: " << constants.log_dc_dk1 << endl;
  cout << "log_dc_dk2: " << constants.log_dc_dk2 << endl;
  cout << "log_d2c_dk1dk2: " << constants.log_d2c_dk1dk2 << endl;
  cout << "log_d2c_dk1dk1: " << constants.log_d2c_dk1dk1 << endl;
  cout << "log_d2c_dk2dk2: " << constants.log_d2c_dk2dk2 << endl << endl;

  cout << "log_dc_dl: " << constants.log_dc_dl << endl;
  cout << "log_d2c_dk1dl: " << constants.log_d2c_dk1dl << endl;
  cout << "log_d2c_dk2dl: " << constants.log_d2c_dk2dl << endl;
  cout << "log_d2c_dldl: " << constants.log_d2c_dldl << endl;
}

void Test::sanity_check()
{
  double kappa1 = 100, kappa2 = 20, lambda = -30;

  double log_const = 2 * log(fabs(lambda)) - log(4) - log(kappa1) - log(kappa2);

  double log_bessel1 = computeLogModifiedBesselFirstKind(0,kappa1);
  double log_bessel2 = computeLogModifiedBesselFirstKind(0,kappa2);
  double log_f = log_bessel1 + log_bessel2; 
  cout << "log_f0: " << log_f << endl;

  log_bessel1 = computeLogModifiedBesselFirstKind(1,kappa1);
  log_bessel2 = computeLogModifiedBesselFirstKind(1,kappa2);
  log_f = log(2) + log_const + log_bessel1 + log_bessel2; 
  cout << "log_f1: " << log_f << endl;

  log_bessel1 = computeLogModifiedBesselFirstKind(2,kappa1);
  log_bessel2 = computeLogModifiedBesselFirstKind(2,kappa2);
  log_f = log(6) + 2 * log_const + log_bessel1 + log_bessel2; 
  cout << "log_f2: " << log_f << endl;
}

void Test::check_sufficient_stats_sine()
{
  double mu1,mu2,kappa1,kappa2,lambda;
  int N = 1000;

  mu1 = 90; mu2 = 90; kappa1 = 100; kappa2 = 10; lambda = -3;

  mu1 *= PI/180; mu2 *= PI/180;
  BVM_Sine bvm_sine(mu1,mu2,kappa1,kappa2,lambda);

  std::vector<Vector> angle_pairs = bvm_sine.generate(N);
  writeToFile("angle_pairs.dat",angle_pairs);

  cout << "Normal:\n";
  struct SufficientStatisticsSine suff_stats;
  computeSufficientStatisticsSine(angle_pairs,suff_stats);

  cout << "\nParallel:\n";
  struct SufficientStatisticsSine suff_stats_parallel;
  computeSufficientStatisticsSine(angle_pairs,suff_stats_parallel);
}

void Test::bvm_sine_kldiv()
{
  double mu1,mu2,kappa1,kappa2,lambda;
  int N = 10;

  mu1 = 90; mu2 = 90; kappa1 = 1; kappa2 = 1; lambda = 0.9;
  mu1 *= PI/180; mu2 *= PI/180;

  BVM_Sine bvm_sine1(mu1,mu2,kappa1,kappa2,lambda);

  mu1 = 100; mu2 = 100; kappa1 = 1; kappa2 = 1; lambda = 0.9;
  mu1 *= PI/180; mu2 *= PI/180;
  BVM_Sine bvm_sine2(mu1,mu2,kappa1,kappa2,lambda);

  double kldiv = bvm_sine1.computeKLDivergence(bvm_sine2);
  cout << "kldiv: " << kldiv << endl;
}

void Test::bvm_sine_kldiv2()
{
  double mu1,mu2,kappa1,kappa2,lambda;
  int N = 100000;

  mu1 = 60; mu2 = 90; kappa1 = 100; kappa2 = 1; lambda = 0.5;
  mu1 *= PI/180; mu2 *= PI/180;

  BVM_Sine bvm_sine(mu1,mu2,kappa1,kappa2,lambda);
  std::vector<Vector> angle_pairs = bvm_sine.generate(N);
  //writeToFile("angle_pairs.dat",angle_pairs);

  mu1 = 50; mu2 = 40; kappa1 = 1; kappa2 = 10; lambda = 0.4;
  mu1 *= PI/180; mu2 *= PI/180;

  BVM_Sine bvm_sine2(mu1,mu2,kappa1,kappa2,lambda);

  /*std::vector<struct EstimatesSine> estimates;
  bvm_sine.computeAllEstimators(angle_pairs,estimates,1,1);

  struct EstimatesSine mle = estimates[MLE];
  BVM_Sine bvm_sine2(mle.mu1,mle.mu2,mle.kappa1,mle.kappa2,mle.lambda);*/

  double kldiv = bvm_sine.computeKLDivergence(bvm_sine2);
  cout << "\nkldiv (actual -): " << kldiv << endl;
  /*kldiv = bvm_sine.computeKLDivergence2(bvm_sine2);
  cout << "kldiv (actual +): " << kldiv << endl;*/

  kldiv = bvm_sine.computeKLDivergence(bvm_sine2,angle_pairs);
  cout << "kldiv (empirical): " << kldiv << endl;
}

void Test::bvm_sine_ml_estimation()
{
  double mu1,mu2,kappa1,kappa2,lambda;
  int N = 1000;

  mu1 = 90; mu2 = 90; kappa1 = 100; kappa2 = 10; lambda = -3;

  mu1 *= PI/180; mu2 *= PI/180;
  BVM_Sine bvm_sine(mu1,mu2,kappa1,kappa2,lambda);
  //cout << "log_norm: " << bvm_sine.computeLogNormalizationConstant() << endl;

  std::vector<Vector> angle_pairs = bvm_sine.generate(N);
  std::vector<struct EstimatesSine> estimates;
  bvm_sine.computeAllEstimators(angle_pairs,estimates,1,1);

  std::vector<Vector> random_sample = bvm_sine.generate_cartesian(angle_pairs);
  writeToFile("bvm_sine.dat",random_sample);
  double negloglike_true = bvm_sine.computeNegativeLogLikelihood(angle_pairs);
  //cout << "\nloglike true: " << -negloglike_true << endl;

  struct EstimatesSine mle = estimates[MLE];
  BVM_Sine bvm_mle(mle.mu1,mle.mu2,mle.kappa1,mle.kappa2,mle.lambda);
  std::vector<Vector> random_sample2 = bvm_mle.generate_cartesian(N);
  writeToFile("bvm_sine2.dat",random_sample2);
  double negloglike_est = bvm_mle.computeNegativeLogLikelihood(angle_pairs);
  //cout << "loglike est: " << -negloglike_est << endl;
}

void Test::bvm_sine_all_estimation()
{
  double mu1,mu2,kappa1,kappa2,lambda;
  int N = 10;

  mu1 = 90; mu2 = 90; kappa1 = 10; kappa2 = 10; lambda = 9;

  mu1 *= PI/180; mu2 *= PI/180;
  BVM_Sine bvm_sine(mu1,mu2,kappa1,kappa2,lambda);
  //cout << "log_norm: " << bvm_sine.computeLogNormalizationConstant() << endl;

  std::vector<Vector> angle_pairs = bvm_sine.generate(N);
  writeToFile("angle_pairs.dat",angle_pairs);
  std::vector<struct EstimatesSine> all_estimates;
  bvm_sine.computeAllEstimators(angle_pairs,all_estimates,1,1);

  Vector statistics,pvalues;
  chisquare_hypothesis_testing(all_estimates,statistics,pvalues);

  std::vector<Vector> random_sample = bvm_sine.generate_cartesian(angle_pairs);
  writeToFile("bvm_sine.dat",random_sample);
}

void Test::testing_sample_empirical_distribution()
{
  int N = 10000;
  double res = 1;
  std::vector<std::vector<int> > true_bins;
  std::vector<Vector> sampled_angle_pairs 
  = sample_empirical_distribution(N,res,true_bins);
  
  /*std::vector<std::vector<int> > 
  sampled_bins = updateBins(sampled_angle_pairs,res);
  outputBins(sampled_bins,res);*/

  string sampled_angle_pairs_density = "./visualize/sampled_data/sampled_density.dat";
  ofstream out(sampled_angle_pairs_density.c_str());
  int row,col;
  for (int i=0; i<N; i++) {
    Vector angle_pair = sampled_angle_pairs[i];
    double theta = angle_pair[0] * 180 / PI;  // convert to degrees
    if (fabs(theta) <= ZERO) {
      row = 0;
    } else {
      row = (int)(ceil(theta/res) - 1);
    }
    double phi = angle_pair[1] * 180 / PI;    // convert to degrees
    if (fabs(phi) <= ZERO) {
      col = 0;
    } else {
      col = (int)(ceil(phi/res) - 1);
    }
    out << fixed << scientific << angle_pair[0] << "\t"
        << angle_pair[1] << "\t"
        << true_bins[row][col] << endl;
  }
  out.close();
}

/* cosine model related */
void Test::generate_bvm_cosine()
{
  double mu1,mu2,kappa1,kappa2,kappa3;
  int N = 1000;

  //mu1 = 90; mu2 = 90; kappa1 = 200; kappa2 = 200; kappa3 = 100;
  //mu1 = 360-66.2; mu2 = 149.6; kappa1 = 65.4; kappa2 = 45.7; kappa3 = 17.3;
  mu1 = 67.4; mu2 = 96.2; kappa1 = 3.6; kappa2 = 1.9; kappa3 = 0.8;

  mu1 *= PI/180; mu2 *= PI/180;
  BVM_Cosine bvm_cosine(mu1,mu2,kappa1,kappa2,kappa3);

  std::vector<Vector> angle_pairs = bvm_cosine.generate(N);
  writeToFile("angle_pairs.dat",angle_pairs);

  std::vector<Vector> random_sample = bvm_cosine.generate_cartesian(angle_pairs);
  writeToFile("bvm_cosine.dat",random_sample);
}


void Test::bvm_cosine_normalization_constant()
{
  double mu1,mu2,kappa1,kappa2,kappa3;

  mu1 = 90; mu2 = 90; kappa1 = 100; kappa2 = 50; kappa3 = -10;
  mu1 *= PI/180; mu2 *= PI/180;

  BVM_Cosine bvm_cosine(mu1,mu2,kappa1,kappa2,kappa3);

  double log_norm = bvm_cosine.computeLogNormalizationConstant();
  cout << "log_norm: " << log_norm << endl;
}

void Test::bvm_cosine_constants()
{
  double mu1,mu2,kappa1,kappa2,kappa3;
  int N = 1000;

  mu1 = 0; mu2 = 0; kappa1 = 100; kappa2 = 50; kappa3 = -30;
  mu1 *= PI/180; mu2 *= PI/180;

  BVM_Cosine bvm_cosine(mu1,mu2,kappa1,kappa2,kappa3);
  cout << "correlation: " << cosine_correlation(kappa1,kappa2,kappa3) << endl;
  BVM_Cosine::Constants constants = bvm_cosine.getConstants();

  double norm2 = integration(mu1,mu2,kappa1,kappa2,kappa3);

 // cout << "log_norm: " << constants.log_c << endl;
 // cout << "log_dc_dk1: " << constants.log_dc_dk1 << endl;
 // cout << "log_dc_dk2: " << constants.log_dc_dk2 << endl;
 // cout << "log_dc_dk3: " << constants.log_dc_dk3 << endl;
  cout << "norm: " << exp(constants.log_c) << endl;
  cout << "log_norm: " << constants.log_c << endl;
//  cout << "dc_dk1: " << exp(constants.log_dc_dk1) << endl;
//  cout << "dc_dk2: " << exp(constants.log_dc_dk2) << endl;
//  cout << "dc_dk3: " << exp(constants.log_dc_dk3) << endl;
//  cout << "d2c_dk1dk1: " << exp(constants.log_d2c_dk1dk1) << endl;
//  cout << "d2c_dk2dk2: " << exp(constants.log_d2c_dk2dk2) << endl;
//  cout << "d2c_dk3dk3: " << exp(constants.log_d2c_dk3dk3) << endl;
//  cout << "d2c_dk1dk2: " << exp(constants.log_d2c_dk1dk2) << endl;
//  cout << "d2c_dk1dk3: " << exp(constants.log_d2c_dk1dk3) << endl;
//  cout << "d2c_dk2dk3: " << exp(constants.log_d2c_dk2dk3) << endl;
}

