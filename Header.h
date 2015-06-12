#ifndef HEADER_H
#define HEADER_H

#include <iostream>
#include <memory>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <ctime>
#include <cassert>
#include <omp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//#include <gsl/gsl_math.h>
//#include <gsl/gsl_integration.h>
//#include <gsl/gsl_monte.h>
//#include <gsl/gsl_monte_plain.h>
//#include <gsl/gsl_monte_miser.h>
//#include <gsl/gsl_monte_vegas.h>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/random.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/detail/erf_inv.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::program_options;
using namespace boost::filesystem;
using namespace boost::math;
using namespace boost::numeric::ublas;
using namespace boost::numeric::odeint;

typedef std::vector<double> Vector;
typedef boost::numeric::ublas::matrix<double> Matrix;
typedef boost::numeric::ublas::identity_matrix<double> IdentityMatrix;
typedef boost::numeric::ublas::zero_matrix<double> ZeroMatrix;

// numeric constants
#define AOM 0.001
#define LARGE_NUMBER 1000000000
#define PI boost::math::constants::pi<double>()
#define LOG_PI log(PI)
#define ZERO std::numeric_limits<double>::epsilon()
#define TOLERANCE 1e-6

#define SET 1 
#define UNSET 0

#endif

