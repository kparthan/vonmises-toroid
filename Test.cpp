#include "Test.h"
#include "Support.h"
#include "Normal.h"
#include "vMC.h"
#include "Mixture_vMC.h"
#include "BVM_Sine.h"

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

void Test::generate_bvm_sine()
{
  double mu1,mu2,kappa1,kappa2,lambda;
  int N = 1000;

  mu1 = 30; mu2 = 60; kappa1 = 25; kappa2 = 100; lambda = 70;

  mu1 *= PI/180; mu2 *= PI/180;
  BVM_Sine bvm_sine(mu1,mu2,kappa1,kappa2,lambda);

  std::vector<Vector> random_sample = bvm_sine.generate_cartesian(N);
  writeToFile("bvm_sine.dat",random_sample);
}
