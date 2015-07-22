#ifndef BVM_SINE_H
#define BVM_SINE_H

#include "Header.h"

class BVM_Sine
{
  private:
    double mu1,mu2;

    double kappa1,kappa2,lambda;  // l^2 < k1 * k2 (unimodal)

    struct Constants {
      double log_c,log_ck1,log_ck2,log_cl;
      double ck1_c,ck2_c,cl_c;
    } constants;

    int computed;

  public:
		//! Constructor
		BVM_Sine();

		//! Constructor that sets value of parameters
		BVM_Sine(double, double, double);

		//! Constructor that sets value of parameters
		BVM_Sine(double, double, double, double, double);

    //! Assignment of an existing BVM_Sine distribution
    BVM_Sine operator=(const BVM_Sine &);

    //! Generate a random sample of (theta,phi) pairs
    std::vector<Vector> generate(int);

    //! Generate a random sample of 3D coordinates on the 2D-torus 
    std::vector<Vector> generate_cartesian(int);

    struct Constants getConstants();

    void computeExpectation();

    void computeConstants();

    double computeLogNormalizationConstant();

    double log_density(double &, double &);

    double computeNegativeLogLikelihood(std::vector<Vector> &);
    double computeNegativeLogLikelihood(
      struct EstimatesSine &, struct SufficientStatisticsSine &
    );

    void computeAllEstimators(
      std::vector<Vector> &, 
      std::vector<struct EstimatesSine> &,
      int, int
    );

    void computeAllEstimators(
      struct SufficientStatisticsSine &,
      std::vector<struct EstimatesSine> &,
      int, int
    );

    struct EstimatesSine computeInitialEstimates(
      struct SufficientStatisticsSine &
    );
};

#endif

