#ifndef CUSTOM_FUNCTION_SINE
#define CUSTOM_FUNCTION_SINE

#include "Header.h"

class CustomFunctionSine
{
  private:
    double mu1;

    double kappa1,kappa2,lambda;

  public:
    CustomFunctionSine(double, double, double, double);

    double solve();
};

#endif

