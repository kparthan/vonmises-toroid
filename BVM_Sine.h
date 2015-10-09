#ifndef BVM_SINE_H
#define BVM_SINE_H

#include "Header.h"

class BVM_Sine
{
  friend class Test;

  private:
    double mu1,mu2;

    double kappa1,kappa2,lambda;  // l^2 < k1 * k2 (unimodal)
                                  // rho = lambda/sqrt(k1*k2)
    struct Constants {
      double log_c,log_dc_dk1,log_dc_dk2;
      double log_dc_dl;
      double log_d2c_dk1dk2,log_d2c_dk1dk1,log_d2c_dk2dk2;
      double log_d2c_dldl,log_d2c_dk1dl,log_d2c_dk2dl;

      double ck1_c;   // E[ cos(t1-mu1) ]
      double ck2_c;   // E[ cos(t2-mu2) ]
      double cl_c;    // E[ sin(t1-mu1) sin(t2-mu2) ]
      double ck1k2_c; // E[ cos(t1-mu1) cos(t2-mu2) ]

      double ck1k1_c; // E[ cos(t1-mu1)^2 ]
      double ck2k2_c; // E[ cos(t2-mu2)^2 ]
      double ck1l_c;  // E[ cos(t1-mu1) sin(t1-mu1) sin(t2-mu2) ]
      double ck2l_c;  // E[ cos(t2-mu2) sin(t1-mu1) sin(t2-mu2) ]
      double cll_c;   // E[ sin(t1-mu1)^2 sin(t2-mu2)^2 ]
    } constants;

    int computed,computed_lognorm;

  public:
		//! Constructor
		BVM_Sine();

		//! Constructor that sets value of parameters
		BVM_Sine(double, double, double);

		//! Constructor that sets value of parameters
		BVM_Sine(double, double, double, double, double);

    //! Assignment of an existing BVM_Sine distribution
    BVM_Sine operator=(const BVM_Sine &);

    double Mean1();
    double Mean2();
    double Kappa1();
    double Kappa2();
    double Lambda();

    //! Generate a random sample of (theta,phi) pairs
    std::vector<Vector> generate(int);

    //! Generate a random sample of 3D coordinates on the 2D-torus 
    std::vector<Vector> generate_cartesian(int);
    std::vector<Vector> generate_cartesian(std::vector<Vector> &);

    struct Constants getConstants();

    void computeExpectation();

    void computeConstants();

    double getLogNormalizationConstant();
    double computeLogNormalizationConstant();

    double compute_series_A(double, double);
    double compute_series_B(double, double);
    double compute_series_C();

    double computeLogParametersProbability(double);
    double computeLogParametersPriorDensity();
    double computeLogParametersPriorDensityTransform();

    double computeLogFisherInformation(double);
    double computeLogFisherInformation_Single();
    double computeLogFisherAxes();
    double computeLogFisherScale();

    double log_density(Vector &);
    double log_density(double &, double &);

    double computeNegativeLogLikelihood(std::vector<Vector> &);
    double computeNegativeLogLikelihood(struct SufficientStatisticsSine &);
    double computeNegativeLogLikelihood(
      struct EstimatesSine &, struct SufficientStatisticsSine &
    );

    double computeMessageLength(std::vector<Vector> &);
    double computeMessageLength(struct SufficientStatisticsSine &);
    double computeMessageLength(
      struct EstimatesSine &, struct SufficientStatisticsSine &
    );

    void computeAllEstimators(
      std::vector<Vector> &, 
      std::vector<struct EstimatesSine> &,
      int, int
    );

    void computeAllEstimators(
      std::vector<Vector> &, 
      struct SufficientStatisticsSine &,
      std::vector<struct EstimatesSine> &,
      int, int
    );

    struct EstimatesSine computeInitialEstimates(
      struct SufficientStatisticsSine &
    );

    void estimateParameters(std::vector<Vector> &, Vector &);
    void updateParameters(struct EstimatesSine &);

    void printParameters(ostream &);

    double computeKLDivergence(BVM_Sine &);
    double computeKLDivergence(struct EstimatesSine &);
};

#endif

