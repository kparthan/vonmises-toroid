#ifndef OPTIMIZE_SINE
#define OPTIMIZE_SINE

#include <nlopt.hpp>
#include "BVM_Sine.h"
#include "Support.h"

//extern int MOMENT_FAIL,MLE_FAIL,MAP_FAIL,MML_FAIL;
//extern bool FAIL_STATUS;

// MLE 
class ML_Sine
{
  private:
    struct SufficientStatisticsSine suff_stats_sine;

  public:
    ML_Sine(
      struct SufficientStatisticsSine &suff_stats_sine
    ) : suff_stats_sine(suff_stats_sine)
    {}

    static double wrap(
      const Vector &x, 
      Vector &grad, 
      void *data) {
        /*for (int i=0; i<x.size(); i++) {
          if(boost::math::isnan(x[i])) {  // return this sub-optimal state 
            if (FAIL_STATUS == 0) {
              MLE_FAIL++;
              FAIL_STATUS = 1;
            }
            return -HUGE_VAL;
          } // if() ends ...
        }*/ // for() ends ...
        return (*reinterpret_cast<ML_Sine*>(data))(x, grad); 
    }

    /*!
     *  minimize function: N * log c(k,b) - k  (m0' x) - b (mj' xx' mj -  mi xx' mi)
     *  \sum_x: sample mean
     *  \sum_xx' : dispersion matrix (S)
     *  k,b,m0,mj,mi are parameters
     */
    double operator() (const Vector &x, Vector &grad) {
      struct EstimatesSine estimates;
      estimates.mu1 = x[0];
      estimates.mu2 = x[1];
      estimates.kappa1 = x[2];
      estimates.kappa2 = x[3];
      estimates.lambda = x[4];

      BVM_Sine bvm(
        estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.lambda
      );
      double fval = bvm.computeNegativeLogLikelihood(estimates,suff_stats_sine)
                    - 2 * suff_stats_sine.N * log(AOM);
      return fval;
    }
};

class OptimizeSine
{
  private:
    double mu1,mu2,kappa1,kappa2,lambda;

  public:
    OptimizeSine(string);

    void initialize(double, double, double, double, double);

    Vector minimize(struct SufficientStatisticsSine &);
};

#endif

