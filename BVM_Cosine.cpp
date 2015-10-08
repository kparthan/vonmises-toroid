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
double BVM_Cosine::computeLogNormalizationConstant()
{
}

