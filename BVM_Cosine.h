#ifndef BVM_COSINE_H
#define BVM_COSINE_H

#include "Header.h"

class BVM_Cosine
{
  private:
    double mu1,mu2;

    double kappa1,kappa2,kappa3;  // k3 < (k1 * k2) / (k1 + k2) (unimodal)

    struct Constants {
      double log_c,log_ck1,log_ck2,log_ck3;
      double ck1_c,ck2_c,ck3_c;
    } constants;

    int computed;

  public:
		//! Constructor
		BVM_Cosine();

		//! Constructor that sets value of parameters
		BVM_Cosine(double, double, double);

		//! Constructor that sets value of parameters
		BVM_Cosine(double, double, double, double, double);

    //! Assignment of an existing BVM_Cosine distribution
    BVM_Cosine operator=(const BVM_Cosine &);

    //! Generate a random sample of (theta,phi) pairs
    std::vector<Vector> generate(int);

    //! Generate a random sample of 3D coordinates on the 2D-torus 
    std::vector<Vector> generate_cartesian(int);
    std::vector<Vector> generate_cartesian(std::vector<Vector> &);

    struct Constants getConstants();

    void computeExpectation();

    void computeConstants();

    double computeLogNormalizationConstant();

};

#endif

