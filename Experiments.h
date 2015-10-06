#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "Header.h"

class Experiments
{
  private:
    int iterations;

  public:
    Experiments();

    void fisher_uncertainty();

    void simulate_sine(struct Parameters &);
};

#endif

