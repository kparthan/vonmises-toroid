#ifndef BVM_COSINE_H
#define BVM_COSINE_H

#include "Header.h"

class BVM_Cosine
{
  private:
    double mu1,mu2;

    double kappa1,kappa2,kappa3;  // k3 < (k1 * k2) / (k1 + k2) (unimodal)

    struct Constants {
      double log_c,log_dc_dk1,log_dc_dk2,log_dc_dk3;
      double log_d2c_dk1dk2,log_d2c_dk1dk1,log_d2c_dk2dk2;
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

    struct Constants {
      double log_c,log_ck1,log_ck2,log_ck3;
      double ck1_c,ck2_c,ck3_c;
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

    double computeLogNormalizationConstant();

};

#endif

