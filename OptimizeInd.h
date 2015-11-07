#ifndef OPTIMIZE_IND
#define OPTIMIZE_IND

#include <nlopt.hpp>
#include "vMC.h"
#include "BVM_Ind.h"
#include "Support.h"

// MLE 
class ML_Ind
{
  private:
    struct SufficientStatisticsInd suff_stats;

  public:
    ML_Ind(
      struct SufficientStatisticsInd &suff_stats
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
        return (*reinterpret_cast<ML_Ind*>(data))(x, grad); 
    }

    /*!
     *  minimize function: 
     *   N * log c(k1,k2) - k1 * \sum cos(t1-mu1) - k2 * \sum cos(t2-mu2) 
     */
    double operator() (const Vector &x, Vector &grad) {
      struct EstimatesInd estimates;
      estimates.mu1 = x[0];
      estimates.mu2 = x[1];
      estimates.kappa1 = x[2];
      estimates.kappa2 = x[3];

      BVM_Ind bvm(
        estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2
      );
      double fval = bvm.computeNegativeLogLikelihood(estimates,suff_stats);
                    //- 2 * suff_stats.N * log(AOM);
      //cout << "-fval: " << -fval << "\t"; print(cout,x,3); cout << endl;
      return fval;
    }
};

// MAP 
class MAP_Ind
{
  private:
    struct SufficientStatisticsInd suff_stats;

  public:
    MAP_Ind(
      struct SufficientStatisticsInd &suff_stats
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
        return (*reinterpret_cast<MAP_Ind*>(data))(x, grad); 
    }

    /*!
     *  minimize function: 
     *   N * log c(k1,k2,lambda) - k1 * cos(t1-mu1) - k2 * cos(t2-mu2) 
     */
    double operator() (const Vector &x, Vector &grad) {
      struct EstimatesInd estimates;
      estimates.mu1 = x[0];
      estimates.mu2 = x[1];
      estimates.kappa1 = x[2];
      estimates.kappa2 = x[3];

      BVM_Ind bvm(
        estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2
      );
      double log_prior = bvm.computeLogParametersPriorDensity();
      double fval = -log_prior + bvm.computeNegativeLogLikelihood(estimates,suff_stats);
      //cout << "-fval: " << -fval << "\t"; print(cout,x,3); cout << endl;
      return fval;
    }
};

// MML
class MML_Ind
{
  private:
    struct SufficientStatisticsInd suff_stats;

    double const_lattk;

  public:
    MML_Ind(
      struct SufficientStatisticsInd &suff_stats
    ) : suff_stats(suff_stats)
    {
      const_lattk = -5.138;
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
        return (*reinterpret_cast<MML_Ind*>(data))(x, grad); 
    }

    /*!
     *
     */
    double operator() (const std::vector<double> &x, std::vector<double> &grad) {
      struct EstimatesInd estimates;
      estimates.mu1 = x[0];
      estimates.mu2 = x[1];
      estimates.kappa1 = x[2];
      estimates.kappa2 = x[3];

      BVM_Ind bvm(
        estimates.mu1,estimates.mu2,estimates.kappa1,estimates.kappa2
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
      double part2 = bvm.computeNegativeLogLikelihood(estimates,suff_stats) + 2
                     - 2 * suff_stats.N * log(AOM);
      double fval = part1 + part2;
      //cout << "fval: " << fval << "\t"; print(cout,x,3); cout << endl;
      //assert(!boost::math::isnan(fval));
      return fval;
    }
};

class OptimizeInd
{
  private:
    double mu1,mu2,kappa1,kappa2;

  public:
    OptimizeInd(string);

    void initialize(double, double, double, double);

    struct EstimatesInd minimize(std::vector<Vector> &);
    struct EstimatesInd minimize(struct SufficientStatisticsInd &);
};

#endif

