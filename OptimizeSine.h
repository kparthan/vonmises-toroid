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
    struct SufficientStatisticsSine suff_stats;

  public:
    ML_Sine(
      struct SufficientStatisticsSine &suff_stats
    ) : suff_stats(suff_stats)
    {}

    static double wrap(
      const Vector &x, 
      Vector &grad, 
      void *data) {
        for (int i=0; i<x.size(); i++) {
          if(boost::math::isnan(x[i])) {  // return this sub-optimal state 
            /*if (FAIL_STATUS == 0) {
              MLE_FAIL++;
              FAIL_STATUS = 1;
            }*/
            return -HUGE_VAL;
          } // if() ends ...
        } // for() ends ...
        return (*reinterpret_cast<ML_Sine*>(data))(x, grad); 
    }

    /*!
     *  minimize function: 
     *   N * log c(k1,k2,lambda) - k1 * cos(t1-mu1) - k2 * cos(t2-mu2) 
     *                           - lambda * sin(t1-mu1) * sin(t2-mu2)
     */
    double operator() (const Vector &x, Vector &grad) {
      struct EstimatesSine estimates;
      estimates.mu1 = x[0];
      estimates.mu2 = x[1];
      //if (estimates.mu1 < 0) estimates.mu1 += 2*PI;
      //if (estimates.mu2 < 0) estimates.mu2 += 2*PI;
      estimates.kappa1 = x[2];
      estimates.kappa2 = x[3];
      estimates.lambda = x[4];

      BVM_Sine bvm(
        estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.lambda
      );
      double fval = bvm.computeNegativeLogLikelihood(estimates,suff_stats);
                    //- 2 * suff_stats.N * log(AOM);
      //cout << "-fval: " << -fval << "\t"; print(cout,x,3); cout << endl;
      return fval;
    }
};

// MAP 
class MAP_Sine
{
  private:
    struct SufficientStatisticsSine suff_stats;

  public:
    MAP_Sine(
      struct SufficientStatisticsSine &suff_stats
    ) : suff_stats(suff_stats)
    {}

    static double wrap(
      const Vector &x, 
      Vector &grad, 
      void *data) {
        for (int i=0; i<x.size(); i++) {
          if(boost::math::isnan(x[i])) {  // return this sub-optimal state 
            /*if (FAIL_STATUS == 0) {
              MLE_FAIL++;
              FAIL_STATUS = 1;
            }*/
            return -HUGE_VAL;
          } // if() ends ...
        } // for() ends ...
        return (*reinterpret_cast<MAP_Sine*>(data))(x, grad); 
    }

    /*!
     *  minimize function: 
     *   N * log c(k1,k2,lambda) - k1 * cos(t1-mu1) - k2 * cos(t2-mu2) 
     *                           - lambda * sin(t1-mu1) * sin(t2-mu2)
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
      double log_prior = bvm.computeLogParametersPriorDensity();
      double fval = -log_prior + bvm.computeNegativeLogLikelihood(estimates,suff_stats);
                    //- 2 * suff_stats.N * log(AOM);
      //cout << "-fval: " << -fval << "\t"; print(cout,x,3); cout << endl;
      return fval;
    }
};

// MML
class MML_Sine
{
  private:
    struct SufficientStatisticsSine suff_stats;

    double const_lattk;

  public:
    MML_Sine(
      struct SufficientStatisticsSine &suff_stats
    ) : suff_stats(suff_stats)
    {
      const_lattk = -6.455;
    }

    static double wrap(
      const std::vector<double> &x, 
      std::vector<double> &grad, 
      void *data) {
        for (int i=0; i<x.size(); i++) {
          //assert(!boost::math::isnan(x[i]));
          if(boost::math::isnan(x[i])) {  // return this sub-optimal state 
            /*if (FAIL_STATUS == 0) {
              FAIL_STATUS = 1;
            }*/
            return 0;
          } // if() ends ...
        } // for() ends ...
        return (*reinterpret_cast<MML_Sine*>(data))(x, grad); 
    }

    /*!
     *
     */
    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      struct EstimatesSine estimates;
      estimates.mu1 = x[0];
      estimates.mu2 = x[1];
      estimates.kappa1 = x[2];
      estimates.kappa2 = x[3];
      estimates.lambda = x[4];

      BVM_Sine bvm(
        estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.lambda
      );
      double log_prior = bvm.computeLogParametersPriorDensity();
      double log_fisher = bvm.computeLogFisherInformation(suff_stats.N);
      double part1 = const_lattk - log_prior + 0.5 * log_fisher;
      if (part1 < 0) {
        cout << "part1: " << part1;
        part1 = 0; 
        //part1 = fabs(part1);
        cout << " :Part 1 is negative ...";
        //MML_FAIL = 1;
      }
      double part2 = bvm.computeNegativeLogLikelihood(estimates,suff_stats) + 2.5
                     - 2 * suff_stats.N * log(AOM);
      double fval = part1 + part2;
      //cout << "fval: " << fval << "\t"; print(cout,x,3); cout << endl;
      //assert(!boost::math::isnan(fval));
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

    struct EstimatesSine minimize(struct SufficientStatisticsSine &);
};

#endif

