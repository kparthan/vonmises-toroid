#ifndef MIXTURE_VMC_H
#define MIXTURE_VMC_H

#include "vMC.h"

class Mixture_vMC
{
  private:
    int K;

    Vector weights;

    std::vector<vMC> components;

  public:
    Mixture_vMC(int, Vector &, std::vector<vMC> &);

    Mixture_vMC operator=(const Mixture_vMC &);

    int randomComponent();

    Vector generate(int);

    std::vector<Vector> generate_cartesian(int);

    double log_density(double);
};

#endif

