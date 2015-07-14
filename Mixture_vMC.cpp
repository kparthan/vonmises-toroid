#include "Mixture_vMC.h"
#include "Support.h"

Mixture_vMC::Mixture_vMC(int K, Vector &weights, std::vector<vMC> &components) :
                         K(K), weights(weights), components(components)
{
  assert(components.size() == K);
}

Mixture_vMC Mixture_vMC::operator=(const Mixture_vMC &source)
{
  if (this != &source) {
    K = source.K;
    weights = source.weights;
    components = source.components; 
  }
  return *this;
}

int Mixture_vMC::randomComponent()
{
  double random = uniform_random();
  double previous = 0;
  for (int i=0; i<weights.size(); i++) {
    if (random <= weights[i] + previous) {
      return i;
    } // if()
    previous += weights[i];
  } // for()
}

Vector Mixture_vMC::generate(int num_samples)
{
  Vector sample_size(K,0);
  for (int i=0; i<num_samples; i++) {
    int k = randomComponent();
    sample_size[k]++;
  }

  Vector random_sample;
  for (int i=0; i<K; i++) {
    int comp_num_samples = (int)sample_size[i];
    Vector sample = components[i].generate(comp_num_samples);
    for (int j=0; j<sample.size(); j++) {
      random_sample.push_back(sample[j]);
    } // for(j) 
  } // for(i)

  return random_sample;
}

std::vector<Vector> Mixture_vMC::generate_cartesian(int sample_size)
{
  Vector thetas = generate(sample_size);

  Vector v(2,0);
  std::vector<Vector> cartesian_sample(sample_size);
  for (int i=0; i<sample_size; i++) {
    v[0] = cos(thetas[i]);
    v[1] = sin(thetas[i]);
    cartesian_sample[i] = v;
  }
  return cartesian_sample;
}

// theta in radians ...
double Mixture_vMC::log_density(double theta)
{
  Vector log_densities(K,0);
  for (int j=0; j<K; j++) {
    log_densities[j] = components[j].log_density(theta);
    assert(!boost::math::isnan(log_densities[j]));
  }
  int max_index = maximumIndex(log_densities);
  double max_log_density = log_densities[max_index];
  double density = 0;
  for (int j=0; j<K; j++) {
    log_densities[j] -= max_log_density;
    density += weights[j] * exp(log_densities[j]);
  }
  return max_log_density + log(density);
}

