#ifndef TEST_H
#define TEST_H

#include "Header.h"

class Test
{
  private:

  public:
    void load_data();

    void bessel();

    void generate_vmc();
    void generate_mix_vmc();

    /* sine model related */
    void generate_bvm_sine();

    void bvm_sine_normalization_constant();
    void bvm_sine_constants();

    void sanity_check();

    void check_sufficient_stats_sine();

    void bvm_sine_kldiv();
    void bvm_sine_kldiv2();

    void bvm_sine_ml_estimation();
    void bvm_sine_all_estimation();

    void testing_sample_empirical_distribution();

    /* cosine model related */
    void generate_bvm_cosine();

    void bvm_cosine_normalization_constant();
    void bvm_cosine_constants();
};

#endif

