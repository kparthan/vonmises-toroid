#ifndef BVM_COSINE_H
#define BVM_COSINE_H

#include "Header.h"

class BVM_Cosine
{
  friend class Test;

  private:
    double mu1,mu2;

    double kappa1,kappa2,kappa3;  // k3 < (k1 * k2) / (k1 + k2) (unimodal)
                                  // rho = k3 / (sqrt(k1-k3) * sqrt(k2-k3))

    struct Constants {
      double log_c,log_dc_dk1,log_dc_dk2,log_dc_dk3;
      double log_d2c_dk1dk1,log_d2c_dk2dk2,log_d2c_dk3dk3;
      double log_d2c_dk1dk2,log_d2c_dk1dk3,log_d2c_dk2dk3;

      double ck1_c;   // E[ cos(t1-mu1) ]
      double ck2_c;   // E[ cos(t2-mu2) ]
      double ck3_c;   // - E[ cos(t1-t2-mu1+mu2) ]

      double ck1k1_c; // E[ cos(t1-mu1)^2 ]
      double ck2k2_c; // E[ cos(t2-mu2)^2 ]
      double ck3k3_c; // E[ cos(t1-t2-mu1+mu2)^2 ]

      double ck1k2_c; // E[ cos(t1-mu1) cos(t2-mu2) ]
      double ck1k3_c; // - E[ cos(t1-mu1) cos(t1-t2-mu1+mu2) ]
      double ck2k3_c; // - E[ cos(t2-mu2) cos(t1-t2-mu1+mu2) ]
    } constants;

    int computed,computed_lognorm;

  public:
		//! Constructor
		BVM_Cosine();

		//! Constructor that sets value of parameters
		BVM_Cosine(double, double, double);

		//! Constructor that sets value of parameters
		BVM_Cosine(double, double, double, double, double);

    //! Assignment of an existing BVM_Cosine distribution
    BVM_Cosine operator=(const BVM_Cosine &);

    double Mean1();
    double Mean2();
    double Kappa1();
    double Kappa2();
    double Kappa3();

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

    double compute_log_dc_dk1();
    double compute_log_dc_dk2();
    double compute_log_dc_dk3();

    double series_type1_kappa1_partial(int, int);
    double series_type2_kappa1_partial(int, int);

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
    double computeNegativeLogLikelihood(struct SufficientStatisticsCosine &);
    double computeNegativeLogLikelihood(
      struct EstimatesCosine &, struct SufficientStatisticsCosine &
    );

    double computeMessageLength(std::vector<Vector> &);
    double computeMessageLength(struct SufficientStatisticsCosine &);
    double computeMessageLength(
      struct EstimatesCosine &, struct SufficientStatisticsCosine &
    );

    void computeAllEstimators(
      std::vector<Vector> &, 
      std::vector<struct EstimatesCosine> &,
      int, int
    );

    void computeAllEstimators(
      std::vector<Vector> &, 
      struct SufficientStatisticsCosine &,
      std::vector<struct EstimatesCosine> &,
      int, int
    );

    struct EstimatesCosine computeInitialEstimates(
      struct SufficientStatisticsCosine &
    );

    void estimateParameters(std::vector<Vector> &, Vector &);
    void updateParameters(struct EstimatesCosine &);

    void printParameters(ostream &);

    double computeKLDivergence(BVM_Cosine &);
    double computeKLDivergence(struct EstimatesCosine &);
};

#endif

