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
      //double kappa = x[1];
      //double kappa = kappa1;
      double kappa = 1;

      double fval = accept_reject_fval_unimodal_marginal_sine(
                      theta,kappa,mu1,kappa1,kappa2,lambda
                    );
      return -fval;
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

    double solve (double x) 
    {
      double tmp = lambda * sin(x - mu1);
      double asq = kappa2 * kappa2 + tmp * tmp;
      double a = sqrt(asq);
      double ratio_bessel = computeRatioBessel(a);
      
      double fval = cos(x - mu1) * ratio_bessel / a;
      fval -= rhs;
      return fval;
    }
};

class BimodalSineObjectiveFunction
{
  private:
    double mu1,mode1,mode2;

    double kappa1,kappa2,lambda;

  public:
    BimodalSineObjectiveFunction(
      double mu1, double kappa1, double kappa2, double lambda, double mode1, double mode2
    ) : mu1(mu1), kappa1(kappa1), kappa2(kappa2), lambda(lambda), mode1(mode1), mode2(mode2)
    {}

    static double wrap(
      const Vector &x, 
      Vector &grad, 
      void *data
    ) {
        return (*reinterpret_cast<BimodalSineObjectiveFunction*>(data))(x, grad); 
    }

    double operator() (
      const Vector &x, Vector &grad
    ) {
      double theta = x[0];
      //double kappa = x[1];
      double kappa = 1;

      double fval = accept_reject_fval_bimodal_marginal_sine(
                      theta,kappa,mu1,kappa1,kappa2,lambda,mode1,mode2
                    );
      return -fval;
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

    Vector minimize_bimodal_objective(double, double);
};

#endif

