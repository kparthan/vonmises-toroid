#ifndef KAPPA_SOLVER
#define KAPPA_SOLVER

#include <nlopt.hpp>
#include "Support.h"

class ObjectiveFunction
{
  private:
    double rbar;

  public:
    ObjectiveFunction(double rbar) : rbar(rbar)
    {}

    static double wrap(
      const Vector &x, 
      Vector &grad, 
      void *data
    ) {
        return (*reinterpret_cast<ObjectiveFunction*>(data))(x, grad); 
    }

    // fval = A_2(k) - rbar
    double operator() (
      const Vector &x, Vector &grad
    ) {
      double kappa = x[0];

      double ratio_bessel = computeRatioBessel(kappa);
      double fval = ratio_bessel - rbar;
      return fval;
    }
};

class KappaSolver
{
  private:
    double rbar;

  public:
    double solve();
};

#endif

