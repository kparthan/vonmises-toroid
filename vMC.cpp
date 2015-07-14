#include "vMC.h"
#include "Support.h"

/*!
 *  \brief This is a constructor module
 */
vMC::vMC()
{
  mu = 0;
  kappa = 1;
  computed = UNSET;
}

/*!
 *  \brief constructor function which sets the value of mean and 
 *  kappa of the distribution
 *  \param mu a reference to a vector<double>
 *  \param kappa a double
 */
vMC::vMC(double mu, double kappa) : mu(mu), kappa(kappa)
{
  computed = UNSET;
}

/*!
 *  \brief This function assigns a source vMC distribution.
 *  \param source a reference to a vMC
 */
vMC vMC::operator=(const vMC &source)
{
  if (this != &source) {
    mu = source.mu;
    kappa = source.kappa;
    log_c = source.log_c;
    computed = source.computed;
  }
  return *this;
}

/*!
 *  \brief This function returns the mean of the distribution
 *  \return the mean of the distribution
 */
double vMC::Mean(void)
{
	return mu;
}

/*!
 *  \brief This function returns the kappa of the distribution
 *  \return the kappa of the distribution
 */
double vMC::Kappa(void)
{
	return kappa;
}

// returns a random sample of angles
Vector vMC::generate(int sample_size)
{
  double tmp = 1 + (4 * kappa * kappa);
  double a = 1 + sqrt(tmp);

  tmp = a - sqrt(2*a);
  double b = tmp / (2*kappa);

  tmp = 1 + (b*b);
  double r = tmp / (2*b);

  double u1,u2,u3,z,f,c,check1,check2,angle;
  Vector thetas(sample_size,0);
  for (int i=0; i<sample_size; i++) {
    repeat:
    u1 = uniform_random();
    z = cos(PI*u1);
    f = (1+r*z) / (r+z);
    c = kappa * (r-f);

    u2 = uniform_random();
    check1 = (c * (2-c)) - u2;
    if (check1 > 0) {
      goto accept;
    } else if (check1 <= 0) {
      check2 = log(c) - log(u2) + 1 - c;
      if (check2 < 0) goto repeat;
      else if (check2 >= 0) goto accept;
    }
    accept:
    u3 = uniform_random();
    angle = mu + (sign(u3-0.5)*acos(f));
    thetas[i] = angle;
  }
  //assert(thetas.size() == sample_size);
  return thetas;
}

std::vector<Vector> vMC::generate_cartesian(int sample_size)
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

double vMC::computeLogNormalizationConstant()
{
  log_c = log(2*PI);
  log_c += computeLogModifiedBesselFirstKind(0,kappa);
  computed = SET;
  return log_c;
}

// theta in radians ...
double vMC::log_density(double theta)
{
  if (computed == UNSET) {
    computeLogNormalizationConstant();
  }

  double log_pdf = -log_c + (kappa * cos(theta - mu));

  return log_pdf;
}

