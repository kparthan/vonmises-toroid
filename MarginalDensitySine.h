#ifndef MARGINAL_DENSITY_SINE
#define MARGINAL_DENSITY_SINE

#include <nlopt.hpp>
#include "vMC.h"
#include "Support.h"

class UnimodalSineObjectiveFunction
{
  private:
    double mu1;

    double kappa1,kappa2,lambda;

  public:
    UnimodalSineObjectiveFunction(
      double mu1, double kappa1, double kappa2, double lambda
    ) : mu1(mu1), kappa1(kappa1), kappa2(kappa2), lambda(lambda)
    {}

    static double wrap(
      const Vector &x, 
      Vector &grad, 
      void *data
    ) {
        return (*reinterpret_cast<UnimodalSineObjectiveFunction*>(data))(x, grad); 
    }

    double operator() (
      const Vector &x, Vector &grad
    ) {
      double theta = x[0];
      double kappa = x[1];

      double fval = accept_reject_fval_unimodal_marginal_sine(
                      theta,kappa,mu1,kappa1,kappa2,lambda
                    );
      return -fval;
    }
};

struct FunctionToApproximate  {
  double operator() (double x)  {
    return x*x - 3*x + 1;  // Replace with your function
  }
};

struct TerminationCondition  {
  bool operator() (double min, double max)  {
    return abs(min - max) <= 0.000001;
  }
};

class CustomFunctionSine
{
  private:
    double mu1;

    double kappa1,kappa2,lambda,rhs;

  public:
    CustomFunctionSine(
      double mu1, double kappa1, double kappa2, double lambda
    ) : mu1(mu1), kappa1(kappa1), kappa2(kappa2), lambda(lambda)
    {
      rhs = kappa1 / (lambda * lambda);
    }

    double solve (double x) {
      double tmp = lambda * sin(x - mu1);
      double asq = kappa2 * kappa2 + tmp * tmp;
      double a = sqrt(asq);
      double ratio_bessel = computeRatioBessel(a);
      
      double fval = cos(x - mu1) * ratio_bessel / a;
      fval -= rhs;
      return fval;
    }
};

class MarginalDensitySine
{
  private:
    double mu1;

    double kappa1,kappa2,lambda;

  public:
    MarginalDensitySine(double, double, double, double);

    Vector minimize_unimodal_objective();

    double solve_custom_function();
};

#endif

