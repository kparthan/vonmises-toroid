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

    void simulate_sine(int iterations=1);
};

#endif

