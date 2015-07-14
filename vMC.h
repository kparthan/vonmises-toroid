#ifndef VMC_H
#define VMC_H

#include "Header.h"

// pdf = 1/C exp{kappa * cos(theta - mu)}, where C = 2 * pi * I_0(kappa)
// see log_c (private variable) ...
class vMC
{
  private:
    //! angular mean of the distribution (in radians) [0, 2 pi)
		double mu;

    //! Concentration parameter 
		double kappa;

    //! Logarithm of normalization constant
    double log_c;

    int computed;

  public:
		//! Constructor
		vMC();

		//! Constructor that sets value of parameters
		vMC(double, double);

    //! Assignment of an existing vMC distribution
    vMC operator=(const vMC &);

		//! Gets the mean 
		double Mean();

    //! Gets the Kappa 
    double Kappa(); 

    //! Generate a random sample of 1D angles
    Vector generate(int);

    //! Generate a random sample of (x,y) coordinates on the unit circle 
    std::vector<Vector> generate_cartesian(int);

    double computeLogNormalizationConstant();

    // log(density)
    double log_density(double); // angle in radians
};

#endif

