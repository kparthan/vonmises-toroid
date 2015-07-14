#include "Test.h"
#include "Support.h"
#include "Normal.h"
#include "vMC.h"

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

  writeToFile("random_sample.dat",random_sample);
}

