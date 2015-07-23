#ifndef KAPPA_SOLVER
#define KAPPA_SOLVER

#include <nlopt.hpp>
#include "Support.h"

class ObjectiveFunction
{
  private:
    double N,constant;

  public:
    ObjectiveFunction(double N, double constant) : N(N), constant(constant)
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

      double log_bessel = computeLogModifiedBesselFirstKind(0,kappa);
      double fval = N * log_bessel - kappa * constant;
      return fval;
    }
};

class KappaSolver
{
  private:
    double N,constant,rbar;

  public:
    KappaSolver(double, double, double);

    double minimize();
};

#endif

