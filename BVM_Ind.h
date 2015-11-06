#ifndef BVM_IND_H
#define BVM_IND_H

#include "Header.h"

// f1 = 1/C1 exp{kappa1 * cos(theta1 - mu1)}, where C1 = 2 * pi * I_0(kappa1)
// f2 = 1/C2 exp{kappa2 * cos(theta2 - mu2)}, where C2 = 2 * pi * I_0(kappa2)
// pdf = f1 * f2
class BVM_Ind
{
  friend class Test;

  private:
    double mu1,mu2;

    double kappa1,kappa2;
                                  
    double log_c;

    int computed;

  public:
		//! Constructor
		BVM_Ind();

		//! Constructor that sets value of parameters
		BVM_Ind(double, double);

		//! Constructor that sets value of parameters
		BVM_Ind(double, double, double, double, double);

    //! Assignment of an existing BVM_Ind distribution
    BVM_Ind operator=(const BVM_Ind &);

    double Mean1();
    double Mean2();
    double Kappa1();
    double Kappa2();

    //! Generate a random sample of (theta,phi) pairs
    std::vector<Vector> generate(int);

    //! Generate a random sample of 3D coordinates on the 2D-torus 
    std::vector<Vector> generate_cartesian(int);
    std::vector<Vector> generate_cartesian(std::vector<Vector> &);

    double getLogNormalizationConstant();
    double computeLogNormalizationConstant();

    double computeLogParametersProbability(double);
    double computeLogParametersPriorDensity();

    double computeLogFisherInformation(double);
    double computeLogFisherInformation_Single();

    double log_density(Vector &);
    double log_density(double &, double &);

    double computeNegativeLogLikelihood(std::vector<Vector> &);
    double computeNegativeLogLikelihood(struct SufficientStatisticsInd &);
    double computeNegativeLogLikelihood(
      struct EstimatesInd &, struct SufficientStatisticsInd &
    );

    double computeMessageLength(std::vector<Vector> &);
    double computeMessageLength(struct SufficientStatisticsInd &);
    double computeMessageLength(
      struct EstimatesInd &, struct SufficientStatisticsInd &
    );

    void computeAllEstimators(
      std::vector<Vector> &, 
      std::vector<struct EstimatesInd> &,
      int, int
    );

    void computeAllEstimators(
      std::vector<Vector> &, 
      struct SufficientStatisticsInd &,
      std::vector<struct EstimatesInd> &,
      int, int
    );

    struct EstimatesInd computeInitialEstimates(
      struct SufficientStatisticsInd &
    );

    void estimateParameters(std::vector<Vector> &, Vector &);
    void updateParameters(struct EstimatesInd &);

    void printParameters(ostream &);

    double computeKLDivergence(BVM_Ind &);
    //double computeKLDivergence2(BVM_Ind &);
    double computeKLDivergence(struct EstimatesInd &);
    double computeKLDivergence(BVM_Ind &, std::vector<Vector> &);
};

#endif

