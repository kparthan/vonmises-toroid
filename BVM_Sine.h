#ifndef BVM_SINE_H
#define BVM_SINE_H

#include "Header.h"

class BVM_Sine
{
  friend class Test;

  private:
    double mu1,mu2;

    double kappa1,kappa2,lambda;  // l^2 < k1 * k2 (unimodal)

    struct Constants {
      double log_c,log_dc_dk1,log_dc_dk2,log_d2c_dk1dk2;
      double log_dc_dl;
      double E_cost1;       // E[ cos(t1-mu1) ]
      double E_cost2;       // E[ cos(t2-mu2) ]
      double E_sint1sint2;  // E[ sin(t1-mu1) sin(t2-mu2) ]
      double E_cost1cost2;  // E[ cos(t1-mu1) cos(t2-mu2) ]
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

    double Mean1();
    double Mean2();
    double Kappa1();
    double Kappa2();
    double Lambda();

    //! Generate a random sample of 3D coordinates on the 2D-torus 
    std::vector<Vector> generate_cartesian(int);
    std::vector<Vector> generate_cartesian(std::vector<Vector> &);

    struct Constants getConstants();

    void computeExpectation();

    void computeConstants();

    double computeLogNormalizationConstant();

    double compute_series_A(double, double);
    double compute_series_B();

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

    double computeKLDivergence(BVM_Sine &);
    double computeKLDivergence(struct EstimatesSine &);
};

#endif

