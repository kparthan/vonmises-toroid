#include "Test.h"
#include "Support.h"
#include "Normal.h"
#include "vMC.h"
#include "Mixture_vMC.h"
#include "BVM_Sine.h"
#include "BVM_Cosine.h"

extern Vector XAXIS,YAXIS,ZAXIS;
extern int ENABLE_DATA_PARALLELISM;
extern int NUM_THREADS;
extern int ESTIMATION,CRITERION;

void Test::load_data()
{
  string data_file = "random_sample.dat";
  std::vector<Vector> data = load_data_table(data_file);
  print(cout,data[10],3); cout << endl;
  string output = "copy.dat";
  writeToFile(output,data);
}

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

  std::vector<Vector> random_sample = bvm_sine.generate_cartesian(N);
  writeToFile("bvm_sine2.dat",random_sample);
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
  cout << "log_dc_dl: " << constants.log_dc_dl << endl;
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

void Test::bvm_sine_ml_estimation()
{
  double mu1,mu2,kappa1,kappa2,lambda;
  int N = 1000;

  mu1 = 90; mu2 = 90; kappa1 = 100; kappa2 = 10; lambda = -30;

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

  struct EstimatesSine mle = estimates[0];
  BVM_Sine bvm_mle(mle.mu1,mle.mu2,mle.kappa1,mle.kappa2,mle.lambda);
  std::vector<Vector> random_sample2 = bvm_mle.generate_cartesian(N);
  writeToFile("bvm_sine2.dat",random_sample2);
  double negloglike_est = bvm_mle.computeNegativeLogLikelihood(angle_pairs);
  //cout << "loglike est: " << -negloglike_est << endl;
}

/* cosine model related */
void Test::generate_bvm_cosine()
{
  double mu1,mu2,kappa1,kappa2,kappa3;
  int N = 1000;

  mu1 = 90; mu2 = 90; kappa1 = 200; kappa2 = 200; kappa3 = 100;

  mu1 *= PI/180; mu2 *= PI/180;
  BVM_Cosine bvm_cosine(mu1,mu2,kappa1,kappa2,kappa3);

  std::vector<Vector> random_sample = bvm_cosine.generate_cartesian(N);
  writeToFile("bvm_cosine.dat",random_sample);
}


