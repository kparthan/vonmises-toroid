#ifndef OPTIMIZE_COSINE
#define OPTIMIZE_COSINE

#include <nlopt.hpp>
#include "vMC.h"
#include "BVM_Cosine.h"
#include "Support.h"

// MLE 
class ML_Cosine
{
  private:
    struct SufficientStatisticsCosine suff_stats;

  public:
    ML_Cosine(
      struct SufficientStatisticsCosine &suff_stats
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
        return (*reinterpret_cast<ML_Cosine*>(data))(x, grad); 
    }

    /*!
     *  minimize function: 
     *   N * log c(k1,k2,kappa3) - k1 * cos(t1-mu1) - k2 * cos(t2-mu2) 
     *                           + kappa3 * cos((t1-mu1)-(t2-mu2))
     */
    double operator() (const Vector &x, Vector &grad) {
      struct EstimatesCosine estimates;
      estimates.mu1 = x[0];
      estimates.mu2 = x[1];
      //if (estimates.mu1 < 0) estimates.mu1 += 2*PI;
      //if (estimates.mu2 < 0) estimates.mu2 += 2*PI;
      estimates.kappa1 = x[2];
      estimates.kappa2 = x[3];
      estimates.kappa3 = x[4];

      BVM_Cosine bvm(
        estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.kappa3
      );
      double fval = bvm.computeNegativeLogLikelihood(estimates,suff_stats);
                    //- 2 * suff_stats.N * log(AOM);
      //cout << "-fval: " << -fval << "\t"; print(cout,x,3); cout << endl;
      return fval;
    }
};

// MAP 
class MAP_Cosine
{
  private:
    struct SufficientStatisticsCosine suff_stats;

  public:
    MAP_Cosine(
      struct SufficientStatisticsCosine &suff_stats
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
        return (*reinterpret_cast<MAP_Cosine*>(data))(x, grad); 
    }

    /*!
     *  minimize function: 
     *   N * log c(k1,k2,kappa3) - k1 * cos(t1-mu1) - k2 * cos(t2-mu2) 
     *                           + kappa3 * cos((t1-mu1)-(t2-mu2))
     */
    double operator() (const Vector &x, Vector &grad) {
      struct EstimatesCosine estimates;
      estimates.mu1 = x[0];
      estimates.mu2 = x[1];
      estimates.kappa1 = x[2];
      estimates.kappa2 = x[3];
      estimates.kappa3 = x[4];

      BVM_Cosine bvm(
        estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.kappa3
      );
      double log_prior = bvm.computeLogParametersPriorDensity();
      double fval = -log_prior + bvm.computeNegativeLogLikelihood(estimates,suff_stats);
                    //- 2 * suff_stats.N * log(AOM);
      //cout << "-fval: " << -fval << "\t"; print(cout,x,3); cout << endl;
      return fval;
    }
};

// MAP Transform 
class MAP_Transform_Cosine
{
  private:
    struct SufficientStatisticsCosine suff_stats;

  public:
    MAP_Transform_Cosine(
      struct SufficientStatisticsCosine &suff_stats
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
        return (*reinterpret_cast<MAP_Transform_Cosine*>(data))(x, grad); 
    }

    /*!
     *  minimize function: 
     *   N * log c(k1,k2,kappa3) - k1 * cos(t1-mu1) - k2 * cos(t2-mu2) 
     *                           + kappa3 * cos((t1-mu1)-(t2-mu2))
     */
    double operator() (const Vector &x, Vector &grad) {
      struct EstimatesCosine estimates;
      estimates.mu1 = x[0];
      estimates.mu2 = x[1];
      estimates.kappa1 = x[2];
      estimates.kappa2 = x[3];
      double rho = x[4];
      estimates.kappa3 = rho * (estimates.kappa1 * estimates.kappa2) / (estimates.kappa1 + estimates.kappa2);

      BVM_Cosine bvm(
        estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.kappa3
      );
      double log_prior = bvm.computeLogParametersPriorDensityTransform();
      double fval = -log_prior + bvm.computeNegativeLogLikelihood(estimates,suff_stats);
                    //- 2 * suff_stats.N * log(AOM);
      //cout << "-fval: " << -fval << "\t"; print(cout,x,3); cout << endl;
      return fval;
    }
};

// MML
class MML_Cosine
{
  private:
    struct SufficientStatisticsCosine suff_stats;

    double const_lattk;

  public:
    MML_Cosine(
      struct SufficientStatisticsCosine &suff_stats
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
        return (*reinterpret_cast<MML_Cosine*>(data))(x, grad); 
    }

    /*!
     *
     */
    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      struct EstimatesCosine estimates;
      estimates.mu1 = x[0];
      estimates.mu2 = x[1];
      estimates.kappa1 = x[2];
      estimates.kappa2 = x[3];
      estimates.kappa3 = x[4];

      BVM_Cosine bvm(
        estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2,estimates.kappa3
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

class OptimizeCosine
{
  private:
    double mu1,mu2,kappa1,kappa2,kappa3;

  public:
    OptimizeCosine(string);

    void initialize(double, double, double, double, double);

    struct EstimatesCosine minimize(struct SufficientStatisticsCosine &);
};

#endif

