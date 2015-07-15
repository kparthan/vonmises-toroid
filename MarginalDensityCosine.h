#ifndef MARGINAL_DENSITY_COSINE
#define MARGINAL_DENSITY_COSINE

#include <nlopt.hpp>
#include "vMC.h"
#include "Support.h"

class UnimodalCosineObjectiveFunction
{
  private:
    double mu1;

    double kappa1,kappa2,kappa3;

  public:
    UnimodalCosineObjectiveFunction(
      double mu1, double kappa1, double kappa2, double kappa3
    ) : mu1(mu1), kappa1(kappa1), kappa2(kappa2), kappa3(kappa3)
    {}

    static double wrap(
      const Vector &x, 
      Vector &grad, 
      void *data
    ) {
        return (*reinterpret_cast<UnimodalCosineObjectiveFunction*>(data))(x, grad); 
    }

    double operator() (
      const Vector &x, Vector &grad
    ) {
      double theta = x[0];
      //double kappa = x[1];
      //double kappa = kappa1;
      double kappa = 1;

      double fval = accept_reject_fval_unimodal_marginal_cosine(
                      theta,kappa,mu1,kappa1,kappa2,kappa3
                    );
      return -fval;
    }
};

class CustomFunctionCosine
{
  private:
    double mu1;

    double kappa1,kappa2,kappa3,rhs;

  public:
    CustomFunctionCosine(
      double mu1, double kappa1, double kappa2, double kappa3
    ) : mu1(mu1), kappa1(kappa1), kappa2(kappa2), kappa3(kappa3)
    {
      rhs = kappa1 / (kappa2 * kappa3);
    }

    double solve (double x) {
      double k23sq = kappa2 * kappa2 + kappa3 * kappa3 
                     - (2 * kappa2 * kappa3 * cos(x-mu1));
      double k23 = sqrt(k23sq);

      double ratio_bessel = computeRatioBessel(k23);
      
      double fval = ratio_bessel / k23;
      fval -= rhs;
      return fval;
    }
};

class BimodalCosineObjectiveFunction
{
  private:
    double mu1,mode1,mode2;

    double kappa1,kappa2,kappa3;

  public:
    BimodalCosineObjectiveFunction(
      double mu1, double kappa1, double kappa2, double kappa3, double mode1, double mode2
    ) : mu1(mu1), kappa1(kappa1), kappa2(kappa2), kappa3(kappa3), mode1(mode1), mode2(mode2)
    {}

    static double wrap(
      const Vector &x, 
      Vector &grad, 
      void *data
    ) {
        return (*reinterpret_cast<BimodalCosineObjectiveFunction*>(data))(x, grad); 
    }

    double operator() (
      const Vector &x, Vector &grad
    ) {
      double theta = x[0];
      //double kappa = x[1];
      double kappa = 1;

      double fval = accept_reject_fval_bimodal_marginal_cosine(
                      theta,kappa,mu1,kappa1,kappa2,kappa3,mode1,mode2
                    );
      return -fval;
    }
};

class MarginalDensityCosine
{
  private:
    double mu1;

    double kappa1,kappa2,kappa3;

  public:
    MarginalDensityCosine(double, double, double, double);

    Vector minimize_unimodal_objective();

    double solve_custom_function();

    Vector minimize_bimodal_objective(double, double);
};

#endif

