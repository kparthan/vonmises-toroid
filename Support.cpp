#include "Support.h"
#include "Mixture_vMC.h"
#include "Test.h"
#include "Structure.h"
#include "BVM_Sine.h"
#include "UniformRandomNumberGenerator.h"
#include "Experiments.h"

int MIXTURE_ID = 1;
int MIXTURE_SIMULATION;
int INFER_COMPONENTS;
int ENABLE_DATA_PARALLELISM;
int NUM_THREADS;
int ESTIMATION,CRITERION;
double MAX_KAPPA;
double IMPROVEMENT_RATE;
int CONSTRAIN_KAPPA;
UniformRandomNumberGenerator *uniform_generator;
int VERBOSE,COMPUTE_KLDIV;
int IGNORE_SPLIT;
double MIN_N;
int SPLITTING = 0;
string EM_LOG_FOLDER;
int DISTRIBUTION;
int MSGLEN_FAIL;
double SERIES_TOLERANCE=1e-100;

struct stat st = {0};

////////////////////// GENERAL PURPOSE FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

/*!
 *  \brief This function checks to see if valid arguments are given to the 
 *  command line output.
 *  \param argc an integer
 *  \param argv an std::vector of strings
 *  \return the parameters of the model
 */
struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;
  string constrain,pdf,estimation_method,criterion;
  double improvement_rate;
  int prior;

  bool noargs = 1;

  cout << "Checking command-line input ..." << endl;
  options_description desc("Allowed options");
  desc.add_options()
       ("help","produce help component")
       ("test","run some test cases")
       ("experiments","run experiments")
       ("iter",value<int>(&parameters.iterations),"number of iterations")
       ("k1",value<double>(&parameters.k1),"kappa 1")
       ("k2",value<double>(&parameters.k2),"kappa 2")
       ("rho",value<double>(&parameters.rho),"rho")
       ("size",value<int>(&parameters.exp_sample_size),"rho")
       ("profile",value<string>(&parameters.profile_file),"path to the profile")
       ("profiles",value<string>(&parameters.profiles_dir),"path to all profiles")
       ("constrain",value<string>(&constrain),"to constrain kappa")
       ("max_kappa",value<double>(&parameters.max_kappa),"maximum value of kappa allowed")
       ("mixture","flag to do mixture modelling")
       ("pdf",value<string>(&pdf),"type of distribution")
       ("k",value<int>(&parameters.fit_num_components),"number of components")
       ("infer_components","flag to infer the number of components")
       ("begin",value<int>(&parameters.start_from),"# of components to begin inference from")
       ("max_k",value<int>(&parameters.max_components),"max components to infer")
       ("log",value<string>(&parameters.infer_log),"log file")
       ("continue","flag to continue inference from some state")
       ("simulate","to simulate a mixture model")
       ("load",value<string>(&parameters.mixture_file),"mixture file")
       ("components",value<int>(&parameters.simulated_components),"# of simulated components")
       ("samples",value<int>(&parameters.sample_size),"sample size generated")
       ("bins","parameter to generate heat maps")
       ("res",value<double>(&parameters.res),"resolution used in heat map images")
       ("mt",value<int>(&parameters.num_threads),"flag to enable multithreading")
       ("estimation",value<string>(&estimation_method),"type of estimation")
       ("criterion",value<string>(&criterion),"type of criterion")
       ("estimate_all","flag to estimate using all methods")
       ("responsibility","flag to compute responsibility matrix")
       ("improvement",value<double>(&improvement_rate),"improvement rate")
  ;
  variables_map vm;
  store(command_line_parser(argc,argv).options(desc).run(),vm);
  notify(vm);

  if (vm.count("help")) {
    Usage(argv[0],desc);
  }

  if (vm.count("test")) {
    parameters.test = SET;
  } else {
    parameters.test = UNSET;
  }

  if (vm.count("experiments")) {
    parameters.experiments = SET;
    if (!vm.count("iter")) {
      parameters.iterations = 1;
    }
  } else {
    parameters.experiments = UNSET;
  }

  if (vm.count("bins")) {
    parameters.heat_map = SET;
    if (!vm.count("res")) {
      parameters.res = DEFAULT_RESOLUTION;
    }
  } else {
    parameters.heat_map = UNSET;
  }

  if (vm.count("profiles") || vm.count("profile")) {
    parameters.read_profiles = SET;
  } else {
    parameters.read_profiles = UNSET;
  }

  if (vm.count("constrain")) {
    CONSTRAIN_KAPPA = SET;
    if (!vm.count("max_kappa")) {
      MAX_KAPPA = DEFAULT_MAX_KAPPA;
    } else {
      MAX_KAPPA = parameters.max_kappa;
    }
  } else {
    CONSTRAIN_KAPPA = UNSET;
    MAX_KAPPA = DEFAULT_MAX_KAPPA;
  }

  if (vm.count("pdf")) {
    if (pdf.compare("sine") == 0) {
      DISTRIBUTION = SINE;
    } else if (pdf.compare("cosine") == 0) {
      DISTRIBUTION = COSINE;
    } else {
      cout << "Unsupported distribution ...\n";
      Usage(argv[0],desc);
    }
  } else {
    DISTRIBUTION = SINE;  // default distribution
  }

  if (vm.count("mixture")) {
    parameters.mixture_model = SET;
    if (!vm.count("k")) {
      parameters.fit_num_components = DEFAULT_FIT_COMPONENTS;
    }
    if (vm.count("infer_components")) {
      parameters.infer_num_components = SET;
      INFER_COMPONENTS = SET;
      if (!vm.count("begin")) {
        parameters.start_from = 1;
      }
      if (!vm.count("max_k")) {
        parameters.max_components = -1;
        if (vm.count("continue")) {
          parameters.continue_inference = SET;
        } else {
          parameters.continue_inference = UNSET;
        }
      }
    } else {
      parameters.infer_num_components = UNSET;
      INFER_COMPONENTS = UNSET;
    }
  } else {
    parameters.mixture_model = UNSET;
  }

  if (vm.count("simulate")) {
    parameters.simulation = SET;
    MIXTURE_SIMULATION = SET;
    if (!vm.count("samples")) {
      parameters.sample_size = DEFAULT_SAMPLE_SIZE;
    }
    if (vm.count("load")) {
      parameters.load_mixture = SET;
    } else {
      parameters.load_mixture = UNSET;
      if (!vm.count("components")) {
        parameters.simulated_components = DEFAULT_SIMULATE_COMPONENTS;
      }
    }
  } else {
    parameters.simulation = UNSET;
    MIXTURE_SIMULATION = UNSET;
  }

  if (vm.count("mt")) {
    NUM_THREADS = parameters.num_threads;
    ENABLE_DATA_PARALLELISM = SET;
  } else {
    ENABLE_DATA_PARALLELISM = UNSET;
    NUM_THREADS = 1;
  }

  if (vm.count("improvement")) {
    IMPROVEMENT_RATE = improvement_rate;
  } else {
    IMPROVEMENT_RATE = 1e-5; // 0.001 % default
  }

  if (vm.count("estimate_all")) {
    parameters.estimate_all = SET;
  } else {
    parameters.estimate_all = UNSET;
  }

  if (vm.count("responsibility")) {
    parameters.compute_responsibility_matrix = SET;
  } else {
    parameters.compute_responsibility_matrix = UNSET;
  }

  if (vm.count("estimation")) {
    if (estimation_method.compare("mle") == 0) {
      ESTIMATION = MLE;
    } else if (estimation_method.compare("map") == 0) {
      ESTIMATION = MAP;
    } else if (estimation_method.compare("map_transform") == 0) {
      ESTIMATION = MAP_TRANSFORM;
    } else if (estimation_method.compare("mml") == 0) {
      ESTIMATION = MML;
    } else {
      cout << "Invalid estimation method ...\n";
      Usage(argv[0],desc);
    }
  } else {  // default is MML estimation ...
    estimation_method = "mml";
    ESTIMATION = MML;
  }

  if (vm.count("criterion")) {
    if (criterion.compare("aic") == 0) {
      CRITERION = AIC;
    } else if (criterion.compare("bic") == 0) {
      CRITERION = BIC;
    } else if (criterion.compare("icl") == 0) {
      CRITERION = ICL;
    } else if (criterion.compare("mml") == 0) {
      CRITERION = MMLC;
    }
  } else {
    criterion = "mmlc";
    CRITERION = MMLC;
  }

  return parameters;
}

/*!
 *  \brief This module prints the acceptable input format to the program
 *  \param exe a reference to a const char
 *  \param desc a reference to a options_description object
 */
void Usage(const char *exe, options_description &desc)
{
  cout << "Usage: " << exe << " [options]" << endl;
  cout << desc << endl;
  exit(1);
}

/*!
 *  \brief This module checks whether the input file exists or not.
 *  \param file_name a reference to a string
 *  \return true or false depending on whether the file exists or not.
 */
bool checkFile(string &file_name)
{
  ifstream file(file_name.c_str());
  return file;
}

/*!
 *  \brief This module prints the elements of a std::vector<Vector<> > to a file
 *  \param v a reference to std::vector<Vector>
 *  \param file_name a pointer to a const char
 */
void writeToFile(const char *file_name, std::vector<Vector> &v, int precision)
{
  ofstream file(file_name);
  for (int i=0; i<v.size(); i++) {
    for (int j=0; j<v[i].size(); j++) {
      file << fixed << setw(10) << setprecision(precision) << v[i][j];
    }
    file << endl;
  }
  file.close(); 
}

void writeToFile(const char *file_name, std::vector<Vector> &v)
{
  ofstream file(file_name);
  for (int i=0; i<v.size(); i++) {
    for (int j=0; j<v[i].size(); j++) {
      file << setw(15) << scientific << v[i][j];
    }
    file << endl;
  }
  file.close(); 
}

void writeToFile(string &file_name, std::vector<Vector> &v)
{
  ofstream file(file_name.c_str());
  for (int i=0; i<v.size(); i++) {
    for (int j=0; j<v[i].size(); j++) {
      file << setw(15) << scientific << v[i][j];
    }
    file << endl;
  }
  file.close(); 
}

void writeToFile(string &file_name, std::vector<std::vector<int> > &v)
{
  ofstream file(file_name.c_str());
  for (int i=0; i<v.size(); i++) {
    for (int j=0; j<v[i].size(); j++) {
      file << fixed << setw(15) << v[i][j];
    }
    file << endl;
  }
  file.close(); 
}

/*!
 *  \brief This module extracts the file name from the path
 *  \param file a reference to a string
 *  \return the extracted portion of the file name
 */
string extractName(string &file)
{
  unsigned pos1 = file.find_last_of("/");
  unsigned pos2 = file.find(".");
  int length = pos2 - pos1 - 1;
  string sub = file.substr(pos1+1,length);
  return sub;
}

/*!
 *  \brief This function prints the elements of an Vector.
 *  \param os a reference to a ostream
 *  \param v a reference to a Vector
 */
void print(ostream &os, const Vector &v, int precision)
{
  if (precision == 0) {
    if (v.size() == 1) {
      os << scientific << "(" << v[0] << ")";
    } else if (v.size() > 1) {
      os << scientific << "(" << v[0] << ", ";
      for (int i=1; i<v.size()-1; i++) {
        os << scientific << v[i] << ", ";
      }
      os << scientific << v[v.size()-1] << ")\t";
    } else {
      os << "No elements in v ...";
    }
  } else if (precision != 0) { // scientific notation
    if (v.size() == 1) {
      os << fixed << setprecision(precision) << "(" << v[0] << ")";
    } else if (v.size() > 1) {
      os << fixed << setprecision(precision) << "(" << v[0] << ", ";
      for (int i=1; i<v.size()-1; i++) {
        os << fixed << setprecision(precision) << v[i] << ", ";
      }
      os << fixed << setprecision(precision) << v[v.size()-1] << ")\t";
    } else {
      os << "No elements in v ...";
    }
  }
}

void print(string &type, struct EstimatesSine &estimates)
{
  cout << "\nTYPE: " << type << endl;
  cout << "m1_est: " << estimates.mu1 * 180/PI << "; ";
  cout << "m2_est: " << estimates.mu2 * 180/PI << "; ";
  cout << "k1_est: " << estimates.kappa1 << "; ";
  cout << "k2_est: " << estimates.kappa2 << "; ";
  cout << "lambda_est: " << estimates.lambda << "; ";
  cout << "rho_est: " << estimates.rho << endl;
}

void print(struct SufficientStatisticsSine &suff_stats)
{
  cout << "sufficient stats sine:\n";
  cout << "N: " << suff_stats.N << endl;
  cout << "cost1: " << suff_stats.cost1 << endl;
  cout << "sint1: " << suff_stats.sint1 << endl;
  cout << "cost2: " << suff_stats.cost2 << endl;
  cout << "sint2: " << suff_stats.sint2 << endl;
  cout << "sint1 sint2: " << suff_stats.sint1sint2 << endl;
  cout << "sint1 cost2: " << suff_stats.sint1cost2 << endl;
  cout << "cost1 sint2: " << suff_stats.cost1sint2 << endl;
  cout << "cost1 cost2: " << suff_stats.cost1cost2 << endl;
}

void print(string &type, struct EstimatesCosine &estimates)
{
  cout << "\nTYPE: " << type << endl;
  cout << "m1_est: " << estimates.mu1 * 180/PI << "; ";
  cout << "m2_est: " << estimates.mu2 * 180/PI << "; ";
  cout << "k1_est: " << estimates.kappa1 << "; ";
  cout << "k2_est: " << estimates.kappa2 << "; ";
  cout << "k3_est: " << estimates.kappa3 << "; ";
  cout << "rho_est: " << estimates.rho << endl;
}

void print(struct SufficientStatisticsCosine &suff_stats)
{
  cout << "sufficient stats cosine:\n";
  cout << "N: " << suff_stats.N << endl;
  cout << "cost1: " << suff_stats.cost1 << endl;
  cout << "sint1: " << suff_stats.sint1 << endl;
  cout << "cost2: " << suff_stats.cost2 << endl;
  cout << "sint2: " << suff_stats.sint2 << endl;
  cout << "cos(t1-t2): " << suff_stats.cost1_t2 << endl;
  cout << "sin(t1-t2): " << suff_stats.sint1_t2 << endl;
}

void check_and_create_directory(string &directory)
{
  if (stat(directory.c_str(), &st) == -1) {
    mkdir(directory.c_str(), 0700);
  }
}

////////////////////// MATH FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

double scale_to_aom(double &x)
{
  int aom_inv = 1.0/AOM;

  int tmp = x * aom_inv;

  double y = tmp * AOM;

  return y;
}

std::vector<Vector> scale_to_aom(std::vector<Vector> &sample)
{
  int N = sample.size();
  int D = sample[0].size();

  Vector emptyvec(D,0);
  std::vector<Vector> scaled_sample(N,emptyvec);
  for (int i=0; i<N; i++) {
    for (int j=0; j<D; j++) {
      scaled_sample[i][j] = scale_to_aom(sample[i][j]);
    }
  }
  return scaled_sample;
}

/*!
 *  \brief This module returns the sign of a number.
 *  \param number a double
 *  \return the sign
 */
int sign(double number)
{
  if (fabs(number) <= ZERO) {
    return 0;
  } else if (number > 0) {
    return 1;
  } else {
    return -1;
  }
}

/*!
 *  \brief Normalizes a Vector
 *  \param x a reference to a Vector
 *  \param unit a reference to a Vector
 *  \return the norm of the Vector
 */
double normalize(Vector &x, Vector &unit)
{
  assert(x.size() == unit.size());
  double l2norm = norm(x);
  for (int i=0; i<x.size(); i++) {
    unit[i] = x[i] / l2norm;
  }
  return l2norm;
}

double norm(Vector &v)
{
  double normsq = 0;
  for (int i=0; i<v.size(); i++) {
    normsq += v[i] * v[i];
  }
  return sqrt(normsq);
}

/*!
 *  \brief This function converts the cartesian coordinates into spherical.
 *  (theta with +X and phi with +Y)
 *  \param cartesian a reference to a Vector 
 *  \param spherical a reference to a Vector 
 */
void cartesian2spherical(Vector &cartesian, Vector &spherical)
{
  Vector unit(3,0);
  double r = normalize(cartesian,unit);

  double x = unit[0];
  double y = unit[1];
  double z = unit[2];

  // theta \in [0,PI]: angle with X-axis
  double theta = acos(x);

  // phi \in[0,2 PI]: angle with positive Y-axis
  double ratio = y/sin(theta);
  if (ratio > 1) {
    ratio = 1;
  } else if (ratio < -1) {
    ratio = -1;
  }
  double angle = acos(ratio);
  double phi = 0;
  if (y == 0 && z == 0) {
    phi = 0;
  } else if (y == 0) {
    if (z > 0) {
      phi = angle;
    } else {
      phi = 2 * PI - angle;
    }
  } else if (z >= 0) {
    phi = angle;
  } else if (z < 0) {
    phi = 2 * PI - angle;
  }

  spherical[0] = r;
  spherical[1] = theta;
  spherical[2] = phi;
}

void spherical2cartesian(Vector &spherical, Vector &cartesian)
{
  cartesian[0] = spherical[0] * cos(spherical[1]);
  cartesian[1] = spherical[0] * sin(spherical[1]) * cos(spherical[2]);
  cartesian[2] = spherical[0] * sin(spherical[1]) * sin(spherical[2]);
}

// angle_pair[0] = theta1, angle_pair[1] = theta2 \in [0,2 \pi)
void toroid2cartesian(Vector &angle_pair, Vector &cartesian)
{
  double r1 = 2;
  double r2 = 1;

  double theta1 = angle_pair[0];
  double theta2 = angle_pair[1];

  cartesian[0] = (r1 + r2 * cos(theta2)) * cos(theta1); // x
  cartesian[1] = (r1 + r2 * cos(theta2)) * sin(theta1); // y
  cartesian[2] = r2 * sin(theta2);  // z
}

/*!
 *  \brief This function computes the dot product between two Vectors.
 *  \param v1 a reference to a Vector
 *  \param v2 a reference to a Vector
 *  \return the dot product
 */
double computeDotProduct(Vector &v1, Vector &v2) 
{
  assert(v1.size() == v2.size());
  double dot_product = 0;
  for (int i=0; i<v1.size(); i++) {
    dot_product += v1[i] * v2[i];
  }
  return dot_product;
}

Vector crossProduct(Vector &v1, Vector &v2) 
{
  Vector ans(3,0);
  ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
  ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
  ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return ans;
}

double computeSum(Vector &data)
{
  double sum = 0;
  for (int i=0; i<data.size(); i++) {
    sum += data[i];
  }
  return sum;
}

/*!
 *  \brief This function computes the surface area of nd-sphere
 *  Surface area = Gamma(d/2+1)/d \pi^(d/2)
 *  \param d an integer
 *  \return the surface area
 */
double computeLogSurfaceAreaSphere(int d)
{
  double log_num = log(d) + ((d/2.0) * log(PI));
  double log_denom = boost::math::lgamma<double>(d/2.0+1);
  return (log_num - log_denom);
}

double uniform_random()
{
  return (*uniform_generator)();
  //return rand()/(double)RAND_MAX;
}

double computeLogModifiedBesselFirstKind(double alpha, double z)
{

  if (!(alpha >= 0 && fabs(z) >= 0)) {
    cout << "Error logModifiedBesselFirstKind: (alpha,x) = (" << alpha << "," << z << ")\n";
    exit(1);
  }
  if (fabs(z) <= TOLERANCE) {
    return -LARGE_NUMBER;
  } 

  double x = fabs(z);

  // constant term log(x^2/4)
  double log_x2_4 = 2.0 * log(x/2.0);
  double four_x2 = 4.0 / (x * x);

  double m = 1;
  double log_sm_prev = -boost::math::lgamma<double>(alpha+1); // log(t0)
  double log_sm_current;
  double log_tm_prev = log_sm_prev; // log(t0)
  double log_tm_current;
  double cm_prev = 0,cm_current; 
  double ratio = (alpha+1) * four_x2;  // t0/t1
  while(ratio < 1) {
    cm_current = (cm_prev + 1) * ratio;
    log_tm_current = log_tm_prev - log(ratio); 
    log_sm_prev = log_tm_current + log(cm_current + 1);
    m++;
    ratio = m * (m+alpha) * four_x2;
    log_tm_prev = log_tm_current;
    cm_prev = cm_current;
  } // while() ends ...
  double k = m;
  log_tm_current = log_tm_prev - log(ratio);  // log(tk)
  double c = log_tm_current - log_sm_prev;
  double tk_sk_1 = exp(c);
  double y = log_sm_prev;
  double zm = 1;
  double log_rm_prev = 0,log_rm_current,rm;
  while(1) {
    log_sm_current = y + log(1 + tk_sk_1 * zm);
    m++;
    log_rm_current = (log_x2_4 - log(m) - log(m+alpha)) + log_rm_prev;
    rm = exp(log_rm_current);
    zm += rm;
    if (rm/zm < 1e-16)  break;
    log_sm_prev = log_sm_current;
    log_rm_prev = log_rm_current;
  } // while() ends ...
  log_sm_current = y + log(1 + tk_sk_1 * zm);
  return (log_sm_current + (alpha * log(x/2.0)));

//  return log(cyl_bessel_i(alpha,z));
}

// A_d(k) = I_{d/2}(k) / I_{d/2-1}(k)
double computeLogRatioBessel(double d, double kappa)
{
  double index = d / 2.0;
  double log_num = computeLogModifiedBesselFirstKind(index,kappa);
  index -= 1;
  double log_denom = computeLogModifiedBesselFirstKind(index,kappa);
  double log_Ad = log_num - log_denom;
  return log_Ad;
}

// A_2(k) = I1 / I0
double computeRatioBessel(double kappa)
{
  /*double log_I0 = computeLogModifiedBesselFirstKind(0,kappa);
  double log_I1 = computeLogModifiedBesselFirstKind(1,kappa);
  double log_Ad = log_I1 - log_I0;*/
  double log_Ad = computeLogRatioBessel(2,kappa);
  return exp(log_Ad);
}

////////////////////// GEOMETRY FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

std::vector<Vector> load_data_table(string &file_name, int dim)
{
  std::vector<Vector> sample;
  ifstream file(file_name.c_str());
  string line;
  Vector numbers(dim,0),unit_vector(dim,0);
  int i;
  while (getline(file,line)) {
    boost::char_separator<char> sep(" \t");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    i = 0;
    BOOST_FOREACH(const string &t, tokens) {
      istringstream iss(t);
      double x;
      iss >> x;
      numbers[i++] = x;
    }
    normalize(numbers,unit_vector);
    sample.push_back(unit_vector);
  }
  file.close();
  return sample;
}

/*!
 *  v1 and v2 are considered to be a column std::vectors
 *  output: v1 * v2' (the outer product matrix)
 */
Matrix outer_prod(Vector &v1, Vector &v2)
{
  assert(v1.size() == v2.size());
  int m = v1.size();
  Matrix ans(m,m);
  for (int i=0; i<m; i++) {
    for (int j=0; j<m; j++) {
      ans(i,j) = v1[i] * v2[j];
    }
  }
  return ans;
}

/*!
 *  v is considered to be a column std::vector
 *  output: m * v (a row std::vector)
 */
Vector prod(Matrix &m, Vector &v)
{
  assert(m.size2() == v.size());
  Vector ans(m.size1(),0);
  for (int i=0; i<m.size1(); i++) {
    for (int j=0; j<m.size2(); j++) {
      ans[i] += m(i,j) * v[j];
    }
  }
  return ans;
}

/*!
 *  v is considered to be a column std::vector
 *  output: v' * m (a row std::vector)
 */
Vector prod(Vector &v, Matrix &m)
{
  assert(m.size1() == v.size());
  Vector ans(m.size2(),0);
  for (int i=0; i<m.size2(); i++) {
    for (int j=0; j<m.size1(); j++) {
      ans[i] += v[j] * m(j,i);
    }
  }
  return ans;
}

/*!
 *  v is considered to be a column std::vector
 *  output: v' M v
 */
double prod_vMv(Vector &v, Matrix &M)
{
  Vector vM = prod(v,M);
  return computeDotProduct(vM,v);
}

/*!
 *  x,y are considered to be a column std::vectors
 *  output: x' M y
 */
double prod_xMy(Vector &x, Matrix &M, Vector &y)
{
  Vector xM = prod(x,M);
  return computeDotProduct(xM,y);
}

/*!
 *  determinant of 3 X 3 matrix
 */
double determinant_3d(Matrix &m)
{
  double det = 0,subdet;
  subdet = m(1,1) * m(2,2) - m(1,2) * m(2,1);
  det += m(0,0) * subdet;

  subdet = m(1,0) * m(2,2) - m(1,2) * m(2,0);
  det -= m(0,1) * subdet;

  subdet = m(1,0) * m(2,1) - m(1,1) * m(2,0);
  det += m(0,2) * subdet;

  return det;
}

/*!
 *  Computes \sum x (x is a vector)
 */
Vector computeVectorSum(std::vector<Vector> &sample) 
{
  int d = sample[0].size();
  Vector resultant(d,0);  // resultant direction

  std::vector<Vector> _resultants;
  for (int i=0; i<NUM_THREADS; i++) {
    _resultants.push_back(resultant);
  }
  int tid;
  #pragma omp parallel if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) private(tid) 
  {
    tid = omp_get_thread_num();
    #pragma omp for
    for (int i=0; i<sample.size(); i++) {
      for (int j=0; j<d; j++) {
        _resultants[tid][j] += sample[i][j];
      }
    } // i loop ends ...
  }

  for (int i=0; i<NUM_THREADS; i++) {
    for (int j=0; j<d; j++) {
      resultant[j] += _resultants[i][j];
    }
  }
  return resultant;
}

// assigns/overrides value to Neff = \sum weights[i]
Vector computeVectorSum(std::vector<Vector> &sample, Vector &weights, double &Neff) 
{
  int d = sample[0].size();
  Vector resultant(d,0);  // resultant direction

  std::vector<Vector> _resultants;
  for (int i=0; i<NUM_THREADS; i++) {
    _resultants.push_back(resultant);
  }
  int tid;
  double sum_neff = 0;
  #pragma omp parallel if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) private(tid) reduction(+:sum_neff) 
  {
    tid = omp_get_thread_num();
    #pragma omp for
    for (int i=0; i<sample.size(); i++) {
      for (int j=0; j<d; j++) {
        _resultants[tid][j] += sample[i][j] * weights[i];
      }
      sum_neff += weights[i];
    } // i loop ends ...
  }
  Neff = sum_neff;

  for (int i=0; i<NUM_THREADS; i++) {
    for (int j=0; j<d; j++) {
      resultant[j] += _resultants[i][j];
    }
  }
  return resultant;
}

/*!
 *  Computes \sum x / N (x is a vector)
 *  *** not a unit vector ... mere division ***
 */
Vector computeNormalizedVectorSum(std::vector<Vector> &sample) 
{
  Vector sum = computeVectorSum(sample);
  for (int j=0; j<sum.size(); j++) {
    sum[j] /= sample.size();
  }
  return sum;
}

/*!
 *  Computes \sum x * x' (x is a vector)
 */
Matrix computeDispersionMatrix(std::vector<Vector> &sample)
{
  int d = sample[0].size();
  Matrix dispersion = ZeroMatrix(d,d);
  for (int i=0; i<sample.size(); i++) {
    dispersion += outer_prod(sample[i],sample[i]);
  }
  return dispersion;
}

Matrix computeDispersionMatrix(std::vector<Vector> &sample, Vector &weights)
{
  int d = sample[0].size();
  Matrix dispersion = ZeroMatrix(d,d);
  Matrix tmp;
  for (int i=0; i<sample.size(); i++) {
    tmp = outer_prod(sample[i],sample[i]);
    dispersion += (weights[i] * tmp);
  }
  return dispersion;
}

/*!
 *  Computes \sum x * x' / N (x is a vector)
 */
Matrix computeNormalizedDispersionMatrix(std::vector<Vector> &sample)
{
  Matrix dispersion = computeDispersionMatrix(sample);
  return dispersion/sample.size();
}

void computeMeanAndCovariance(
  std::vector<Vector> &data, 
  Vector &weights,
  Vector &mean, 
  Matrix &cov
) {
  int N = data.size();
  int D = data[0].size();

  double Neff;
  mean = computeVectorSum(data,weights,Neff);
  for (int i=0; i<D; i++) {
    mean[i] /= Neff;
  }

  std::vector<Vector> x_mu(N);
  Vector diff(D,0);
  for (int i=0; i<N; i++) {
    for (int j=0; j<D; j++) {
      diff[j] = data[i][j] - mean[j];
    }
    x_mu[i] = diff;
  }
  Matrix S = computeDispersionMatrix(x_mu,weights);
  if (Neff > 1) {
    cov = S / (Neff - 1);
  } else {
    cov = S / Neff;
  }
}

/*!
 *  Matrix inverse C++ Boost::ublas
 */
bool invertMatrix(const Matrix &input, Matrix &inverse)
{
  typedef permutation_matrix<std::size_t> pmatrix;

  // create a working copy of the input
  Matrix A(input);

  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  int res = lu_factorize(A, pm);
  if (res != 0)
    return false;

  // create identity matrix of "inverse"
  inverse.assign(IdentityMatrix (A.size1()));

  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);

  return true;
}

/*!
 *  Eigen decomposition
 *  Inputs:
 *                m -- a symmetric matrix
 *    eigen_vectors -- an identity matrix
 *  Outputs:
 *    eigen_vectors -- each column is a unit eigen vector
 */
void eigenDecomposition(
  Matrix m, 
  Vector &eigen_values,
  Matrix &eigen_vectors
) {
  // check if m is symmetric
  int num_rows = m.size1();
  int num_cols = m.size2();
  if (num_rows != num_cols) {
    cout << "Error: rows: " << num_rows << " != cols: " << num_cols << endl;
    exit(1);
  }
  for (int i=0; i<num_rows; i++) {
    for (int j=0; j<num_cols; j++) {
      if (fabs(m(i,j)-m(j,i)) >= TOLERANCE) {
        cout << "Error: Matrix is not symmetric ...\n";
        cout << "m: " << m << endl;
        cout << "m(" << i << "," << j << ") != m(" << j << "," << i << ")\n";
        exit(1);
      }
    }
  }

  // matrix is now validated ...
  int MAX_ITERATIONS = 100;
  for (int i=0; i < MAX_ITERATIONS; i++) {
    //find the largest off-diagonal 
    int max_row = 0, max_col = 1;
    int cur_row, cur_col;
    double max_val = m(max_row,max_col);
    for (cur_row = 0; cur_row < num_rows-1; ++cur_row) {
      for (cur_col = cur_row + 1; cur_col < num_cols; ++cur_col) {
        if (fabs(m(cur_row,cur_col)) > max_val) {
          max_row = cur_row;
          max_col = cur_col;
          max_val = fabs(m(cur_row,cur_col));
        }
      }
    }

    if (max_val <= ZERO) {
      break; //finished
    }

    jacobiRotateMatrix(m,eigen_vectors,max_row,max_col);
  }

  for (int i = 0; i < num_cols; i++) {
    eigen_values[i] = m(i,i);
  }

  //cout << "eigen_values: "; print(cout,eigen_values,0); cout << endl;
  //cout << "eigen_vectors: " << eigen_vectors << endl;
}

void jacobiRotateMatrix(
  Matrix &m,
  Matrix &eigen_vectors, 
  int max_row, 
  int max_col
) {
  double diff = m(max_col,max_col) - m(max_row,max_row);
  double phi, t, c, s, tau, temp;
  int i;
  
  phi = diff / (2.0 * m(max_row,max_col));
  t = 1.0 / (std::fabs(phi) + std::sqrt((phi*phi) + 1.0));
  if(phi < 0){ t = -t; }

  c = 1.0 / std::sqrt(t*t + 1.0);
  s = t*c;
  tau = s/(1.0 + c);

  temp = m(max_row,max_col);
  m(max_row,max_col) = 0;
  m(max_row,max_row) = m(max_row,max_row) - (t*temp);
  m(max_col,max_col) = m(max_col,max_col) + (t*temp);
  for(i = 0; i < max_row; i++){ // Case i < max_row
    temp = m(i,max_row);
    m(i,max_row) = temp - (s*(m(i,max_col) + (tau*temp)));
    m(i,max_col) = m(i,max_col) + (s*(temp - (tau*m(i,max_col))));
  }
  for(i = max_row + 1; i < max_col; i++){ // Case max_row < i < max_col
    temp = m(max_row,i);
    m(max_row,i) = temp - (s*(m(i,max_col) + (tau*m(max_row,i))));
    m(i,max_col) = m(i,max_col) + (s*(temp - (tau*m(i,max_col))));
  }
  for(i = max_col + 1; i < m.size2(); i++){ // Case i > max_col
    temp = m(max_row,i);
    m(max_row,i) = temp - (s*(m(max_col,i) + (tau*temp)));
    m(max_col,i) = m(max_col,i) + (s*(temp - (tau*m(max_col,i))));
  }

  for (i = 0; i < eigen_vectors.size1(); i++) { // update the transformation matrix
    temp = eigen_vectors(i,max_row);
    eigen_vectors(i,max_row) = temp
      - (s*(eigen_vectors(i,max_col) + (tau*temp)));
    eigen_vectors(i,max_col) = eigen_vectors(i,max_col)
      + (s*(temp - (tau*eigen_vectors(i,max_col))));
  }
  return;
}

double computeSquaredEuclideanDistance(Vector &p1, Vector &p2)
{
  int D = p1.size();
  double distsq = 0;
  for (int i=0; i<D; i++) {
    distsq += (p1[i] - p2[i]) * (p1[i] - p2[i]);
  }
  //return sqrt(distsq);
  return distsq;
}

double compute_aic(int k, int n, double neg_log_likelihood)
{
  /*double ans = 2 * k * n / (double) (n - k -1);
  ans += 2 * neg_log_likelihood;*/

  double ans = neg_log_likelihood + k;
  return ans / (2*log(2));
}

double compute_bic(int k, int n, double neg_log_likelihood)
{
  /*double ans = 2 * neg_log_likelihood;
  ans += k * (log(n) + log(2*PI));*/

  double ans = neg_log_likelihood + (0.5 * k * log(n));
  return ans / (log(2));
}

////////////////////// MIXTURE FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

/*!
 *  \brief This function computes the approximation of the constant term for
 *  the constant term in the message length expression (pg. 257 Wallace)
 *  \param d an integer
 *  \return the constant term
 */
double computeConstantTerm(int d)
{
  double ad = 0;
  ad -= 0.5 * d * log(2 * PI);
  ad += 0.5 * log(d * PI);
  return ad;
}

double logLatticeConstant(int d)
{
  double cd = computeConstantTerm(d);
  double ans = -1;
  ans += (2.0 * cd / d);
  return ans;
}

/*!
 *  \brief This function bins the sample data 
 *  \param res a double
 *  \param unit_coordinates a reference to a vector<vector<double> > 
 */
std::vector<std::vector<int> > updateBins(
  std::vector<Vector> &angle_pairs, double res
) {
  std::vector<std::vector<int> > bins;
  int num_rows = 360 / res;
  int num_cols = 360 / res;
  std::vector<int> tmp(num_cols,0);
  for (int i=0; i<num_rows; i++) {
    bins.push_back(tmp);
  }

  cout << "angle_pairs.size(): " << angle_pairs.size() << endl;
  int row,col;
  for (int i=0; i<angle_pairs.size(); i++) {
    //cout << "i: " << i << endl; 
    double theta = angle_pairs[i][0] * 180 / PI;  // convert to degrees
    if (fabs(theta) <= ZERO) {
      row = 0;
    } else {
      row = (int)(ceil(theta/res) - 1);
    }
    double phi = angle_pairs[i][1] * 180 / PI;    // convert to degrees
    if (fabs(phi) <= ZERO) {
      col = 0;
    } else {
      col = (int)(ceil(phi/res) - 1);
    }
    if (row >= bins.size() || col >= bins[0].size()) {
      cout << "outside bounds: " << row << " " << col << "\n";
      cout << "theta: " << theta << " phi: " << phi << endl;
      cout << "angle_pairs[i][0]: " << angle_pairs[i][0] 
           << " angle_pairs[i][1]: " << angle_pairs[i][1] << endl;
      fflush(stdout);
    }
    bins[row][col]++;
  } // for (i)
  return bins;
}

/*!
 *  \brief This function outputs the bin data.
 *  \param bins a reference to a std::vector<std::vector<int> >
 */
void outputBins(std::vector<std::vector<int> > &bins, double res)
{
  double theta=0,phi;
  string fbins2D_file,fbins3D_file;
  fbins2D_file = "./visualize/sampled_data/bins2D.dat";
  fbins3D_file = "./visualize/sampled_data/bins3D.dat";
  ofstream fbins2D(fbins2D_file.c_str());
  ofstream fbins3D(fbins3D_file.c_str());
  //Vector cartesian(3,0);
  Vector angle_pair(2,0);
  for (int i=0; i<bins.size(); i++) {
    phi = 0;
    angle_pair[0] = theta * PI / 180;
    for (int j=0; j<bins[i].size(); j++) {
      fbins2D << fixed << setw(10) << bins[i][j];
      phi += res;
      angle_pair[1] = phi * PI / 180;
      /*toroid2cartesian(angle_pair,cartesian);
      for (int k=0; k<3; k++) {
        fbins3D << fixed << setw(10) << setprecision(4) << cartesian[k];
      }*/ // for(k)
      fbins3D << fixed << scientific << setprecision(6) 
              << angle_pair[0] << "\t" << angle_pair[1] << "\t";
      fbins3D << fixed << setw(10) << bins[i][j] << endl;
    } // for(j)
    theta += res;
    fbins2D << endl;
  } // for(i)
  fbins2D.close();
  fbins3D.close();
}

std::vector<Vector> sample_empirical_distribution(
  int N, double res, std::vector<std::vector<int> > &true_bins
) {
  struct Parameters parameters;
  parameters.profile_file = "./data/dihedral_angles.dat";

  std::vector<Vector> data;
  gatherData(parameters,data);
  cout << "data.size(): " << data.size() << endl;
  //std::vector<std::vector<int> > true_bins = updateBins(data,res);
  true_bins = updateBins(data,res);
  //string true_bins_file = "true_bins.dat";  // integers
  //writeToFile(true_bins_file,true_bins);

  int num_rows = true_bins.size();
  int num_cols = true_bins[0].size();
  int num_bins = num_rows * num_cols;
  cout << "num_bins: " << num_bins << endl;

  Vector emptyvec(num_cols,0);
  std::vector<Vector> prob_bins(num_rows,emptyvec);
  Vector elements(num_bins,0);
  int count = 0;
  for (int i=0; i<num_rows; i++) {
    for (int j=0; j<num_cols; j++) {
      assert(!boost::math::isnan(true_bins[i][j]));
      prob_bins[i][j] = true_bins[i][j] / (double) data.size();
      elements[count++] = prob_bins[i][j];
    } // for (j)
  } // for (i)
  //string prob_bins_file = "prob_bins.dat";  // fractional values
  //writeToFile(prob_bins_file,prob_bins);

  std::vector<int> sorted_index;
  Vector sorted_elements = sort(elements,sorted_index);
  Vector cumsum(num_bins,0);
  cumsum[0] = sorted_elements[0];
  for (int i=1; i<num_bins; i++) {
    cumsum[i] = cumsum[i-1] + sorted_elements[i];
    //cout << sorted_index[i] << "\t\t" << cumsum[i] << endl;
  }

  std::vector<Vector> random_sample;
  Vector angle_pair(2,0);
  for (int i=0; i<N; i++) {
    double cdf = uniform_random();
    int bin;
    for (int j=0; j<num_bins; j++) {
      if (cdf <= cumsum[j]) {
        bin = sorted_index[j];
        break;
      } // if ()
    } // for (j)
    int row = bin / num_cols;
    double theta = (row + uniform_random()) * res;  // in degrees
    //double theta = row * res;
    int col = bin % num_cols;
    double phi = (col + uniform_random()) * res;   // in degrees`
    //double phi = col * res;
    angle_pair[0] = theta * PI/180;
    angle_pair[1] = phi * PI/180;
    random_sample.push_back(angle_pair);
  } // for (i)

  return random_sample;
}

/*!
 *  \brief This function is used to read the angular profiles and use this data
 *  to estimate parameters of a Von Mises distribution.
 *  \param parameters a reference to a struct Parameters
 */
void computeEstimators(struct Parameters &parameters)
{
  std::vector<Vector> angle_pairs;
  bool success = gatherData(parameters,angle_pairs);
  if (parameters.heat_map == SET) {
    std::vector<std::vector<int> > bins = updateBins(angle_pairs,parameters.res);
    outputBins(bins,parameters.res);
  }
  if (success && parameters.mixture_model == UNSET) {  // no mixture modelling
    modelOneComponent(parameters,angle_pairs);
  } else if (success && parameters.mixture_model == SET) { // mixture modelling
    modelMixture(parameters,angle_pairs);
  }
}

bool gatherData(struct Parameters &parameters, std::vector<Vector> &angle_pairs)
{
  if (parameters.profile_file.compare("") == 0) {
    path p(parameters.profiles_dir);
    cout << "path: " << p.string() << endl;
    if (exists(p)) { 
      if (is_directory(p)) { 
        std::vector<path> files; // store paths,
        copy(directory_iterator(p), directory_iterator(), back_inserter(files));
        cout << "# of profiles: " << files.size() << endl;
        int tid;
        std::vector<std::vector<Vector> > _angle_pairs(NUM_THREADS);
        #pragma omp parallel num_threads(NUM_THREADS) private(tid)
        {
          tid = omp_get_thread_num();
          if (tid == 0) {
            cout << "# of threads: " << omp_get_num_threads() << endl;
          }
          #pragma omp for 
          for (int i=0; i<files.size(); i++) {
            Structure structure;
            structure.load(files[i]);
            std::vector<Vector> coords = structure.getAnglePairs();
            for (int j=0; j<coords.size(); j++) {
              _angle_pairs[tid].push_back(coords[j]);
            }
          }
        }
        for (int i=0; i<NUM_THREADS; i++) {
          for (int j=0; j<_angle_pairs[i].size(); j++) {
            angle_pairs.push_back(_angle_pairs[i][j]);
          }
        }
        cout << "# of profiles read: " << files.size() << endl;
        return 1;
      } else {
        cout << p << " exists, but is neither a regular file nor a directory\n";
      }
    } else {
      cout << p << " does not exist\n";
    }
    return 0;
  } else if (parameters.profiles_dir.compare("") == 0) {
    if (checkFile(parameters.profile_file)) {
      // read a single profile
      Structure structure;
      structure.load(parameters.profile_file);
      angle_pairs = structure.getAnglePairs();
      return 1;
    } else {
      cout << "Profile " << parameters.profile_file << " does not exist ...\n";
      return 0;
    }
  }
}

void modelOneComponent(struct Parameters &parameters, std::vector<Vector> &angle_pairs)
{
  cout << "Sample size: " << angle_pairs.size() << endl;
  Vector weights(angle_pairs.size(),1);
  if (DISTRIBUTION == SINE) {
    BVM_Sine bvm_sine;
    std::vector<struct EstimatesSine> all_estimates;
    bvm_sine.computeAllEstimators(angle_pairs,all_estimates,1,0);
  } else if (DISTRIBUTION == COSINE) {
  }
}

void modelMixture(struct Parameters &parameters, std::vector<Vector> &data)
{
  Vector data_weights(data.size(),1);
  // if the optimal number of components need to be determined
  if (DISTRIBUTION == SINE) {
    if (parameters.infer_num_components == SET) {
      Mixture_Sine mixture;
      if (parameters.continue_inference == UNSET) {
        Mixture_Sine m(parameters.start_from,data,data_weights);
        mixture = m;
        mixture.estimateParameters();
      } else if (parameters.continue_inference == SET) {
        mixture.load(parameters.mixture_file,data,data_weights);
      } // continue_inference
      ofstream log(parameters.infer_log.c_str());
      Mixture_Sine stable = inferComponents(mixture,data.size(),log);
      log.close();
    } else if (parameters.infer_num_components == UNSET) {
      // for a given value of number of components
      // do the mixture modelling
      Mixture_Sine mixture(parameters.fit_num_components,data,data_weights);
      mixture.estimateParameters();
    }
  } else if (DISTRIBUTION == COSINE) {
  }
}

void simulateMixtureModel(struct Parameters &parameters)
{
  std::vector<Vector> data;
  if (DISTRIBUTION == SINE) {
    if (parameters.load_mixture == SET) {
      Mixture_Sine original;
      original.load(parameters.mixture_file);
      bool save = 1;
      if (parameters.read_profiles == SET) {
        bool success = gatherData(parameters,data);
        if (!success) {
          cout << "Error in reading data...\n";
          exit(1);
        }
      } else if (parameters.read_profiles == UNSET) {
        data = original.generate(parameters.sample_size,save);
      }
      double msglen = original.compress(data);
      if (parameters.heat_map == SET) {
        original.generateHeatmapData(parameters.res);
        std::vector<std::vector<int> > bins = updateBins(data,parameters.res);
        outputBins(bins,parameters.res);
        /* save mixture_density if not already done */
        if (parameters.read_profiles == SET) {
          string mix_density_file = "./visualize/sampled_data/mixture_density.dat";
          ofstream mix(mix_density_file.c_str());
          double mix_density;
          for (int i=0; i<data.size(); i++) {
            mix_density = exp(original.log_probability(data[i]));
            for (int j=0; j<data[i].size(); j++) {
              mix << scientific << setprecision(6) << data[i][j] << "\t\t";
            } // j
            mix << scientific << setprecision(6) << mix_density << endl;
          } // i
          mix.close();
        } // if (parameters.read_profiles == SET) {
      } // if (parameters.heat_map == SET)
    } else if (parameters.load_mixture == UNSET) {
      int k = parameters.simulated_components;
      //srand(time(NULL));
      Vector weights = generateFromSimplex(k);
      std::vector<BVM_Sine> components = generateRandomComponents(k);
      Mixture_Sine original(k,components,weights);
      bool save = 1;
      data = original.generate(parameters.sample_size,save);
      // save the simulated mixture
      ofstream file("./simulation/simulated_mixture");
      for (int i=0; i<k; i++) {
        file << fixed << setw(10) << setprecision(5) << weights[i];
        file << "\t";
        components[i].printParameters(file);
      }
      file.close();
    }
  } else if (DISTRIBUTION == COSINE) {
  }

  // model a mixture using the original data
  if (parameters.heat_map == UNSET) {
    if (parameters.mixture_model == UNSET) {
      modelOneComponent(parameters,data);
    } else if (parameters.mixture_model == SET) {
      modelMixture(parameters,data);
    }
  }
}

Vector generateFromSimplex(int K)
{
  Vector values(K,0);
  double random,sum = 0;
  for (int i=0; i<K; i++) {
    // generate a random value in (0,1)
    random = uniform_random();
    assert(random > 0 && random < 1);
    // sampling from an exponential distribution with \lambda = 1
    values[i] = -log(1-random);
    sum += values[i];
  }
  for (int i=0; i<K; i++) {
    values[i] /= sum;
  }
  return values;
}

std::vector<BVM_Sine> generateRandomComponents(int num_components)
{
  std::vector<BVM_Sine> components;
  for (int i=0; i<num_components; i++) {
    double mu1 = uniform_random() * 2 * PI;
    double mu2 = uniform_random() * 2 * PI;
    double kappa1 = uniform_random() * MAX_KAPPA;
    double kappa2 = uniform_random() * MAX_KAPPA;
    double pho = (uniform_random() * 2 - 1);
    double lambda = pho * sqrt(kappa1 * kappa2);
    BVM_Sine bvm_sine(mu1,mu2,kappa1,kappa2,lambda);
    components.push_back(bvm_sine);
  }
  return components;
}

Mixture_Sine inferComponents(std::vector<Vector> &data, string &log_file)
{
  Vector data_weights(data.size(),1.0);
  Mixture_Sine m(1,data,data_weights);
  Mixture_Sine mixture = m;
  mixture.estimateParameters();
  ofstream log(log_file.c_str());
  Mixture_Sine inferred = inferComponents(mixture,data.size(),log);
  log.close();
  return inferred;
}

Mixture_Sine inferComponents(Mixture_Sine &mixture, int N, ostream &log)
{
  int K,iter = 0;
  std::vector<BVM_Sine> components;
  Mixture_Sine modified,improved,parent;
  Vector sample_size;

  double null_msglen = mixture.computeNullModelMessageLength();
  log << "Null model encoding: " << null_msglen << " bits."
      << "\t(" << null_msglen/N << " bits/point)\n\n";

  improved = mixture;

  while (1) {
    parent = improved;
    iter++;
    log << "Iteration #" << iter << endl;
    log << "Parent:\n";
    parent.printParameters(log,1);
    parent.printIndividualMsgLengths(log);
    components = parent.getComponents();
    sample_size = parent.getSampleSize();
    K = components.size();
    for (int i=0; i<K; i++) { // split() ...
      if (sample_size[i] > MIN_N) {
        IGNORE_SPLIT = 0;
        cout << "Splitting component " << i+1 << " ...\n";
        modified = parent.split(i,log);
        if (IGNORE_SPLIT == 0) {
          updateInference(modified,improved,N,log,SPLIT);
        } else {
          log << "\t\tIGNORING SPLIT ... \n\n";
        }
      }
    } // split() ends ... 
    if (K >= 2) {  // kill() ...
      for (int i=0; i<K; i++) {
        cout << "Deleting component " << i+1 << " ...\n";
        modified = parent.kill(i,log);
        updateInference(modified,improved,N,log,KILL);
      } // killing each component
    } // if (K > 2) loop
    if (K > 1) {  // join() ...
      for (int i=0; i<K; i++) {
        int j = parent.getNearestComponent(i); // closest component
        cout << "Merging components " << i+1 << "and " << "j+1" << " ...\n";
        modified = parent.join(i,j,log);
        updateInference(modified,improved,N,log,JOIN);
      } // join() ing nearest components
    } // if (K > 1) loop
    if (improved == parent) goto finish;
  } // if (improved == parent || iter%2 == 0) loop

  finish:
  string inferred_mixture_file = "./simulation/inferred_mixture";
  parent.printParameters(inferred_mixture_file);
  return parent;
}

/*!
 *  \brief Updates the inference
 *  \param modified a reference to a Mixture_Sine
 *  \param current a reference to a Mixture_Sine
 *  \param log a reference to a ostream
 *  \param operation an integer
 */
void updateInference(
  Mixture_Sine &modified, 
  Mixture_Sine &current, 
  int N, 
  ostream &log, 
  int operation
) {
  double modified_value,current_value,improvement_rate;

  /*switch(CRITERION) {
    case AIC:
      modified_value = modified.getAIC();
      current_value = current.getAIC();
      break;

    case BIC:
      modified_value = modified.getBIC();
      current_value = current.getBIC();
      break;

    case ICL:
      modified_value = modified.getICL();
      current_value = current.getICL();
      break;

    case MMLC:
      modified_value = modified.getMinimumMessageLength();
      current_value = current.getMinimumMessageLength();
      break;
  }*/

/*
  if (current_value > modified_value) {
    improvement_rate = (current_value - modified_value) / fabs(current_value);
    //if (operation == KILL || operation == JOIN || operation == SPLIT) {
    if (operation == KILL || operation == JOIN) {
      log << "\t ... IMPROVEMENT ... (+ " << fixed << setprecision(3) 
          << 100 * improvement_rate << " %) ";
      log << "\t\t[ACCEPT]\n\n";
      current = modified;
    } // kill | join 
    if (operation == SPLIT) {
      if (improvement_rate >= IMPROVEMENT_RATE) {
        log << "\t ... IMPROVEMENT ... (+ " << fixed << setprecision(3) 
            << 100 * improvement_rate << " %) ";
        log << "\t\t[ACCEPT]\n\n";
        current = modified;
      } else {
        log << "\t ... IMPROVEMENT ... (+ " << fixed << setprecision(3) 
            << 100 * improvement_rate << " %) < (" 
            << 100 * IMPROVEMENT_RATE << " %)" ;
        log << "\t\t[REJECT]\n\n";
      } // if-else() 
    } // split
  } else {
    log << "\t ... NO IMPROVEMENT\t\t\t[REJECT]\n\n";
  } // if (current > modified)
*/

  modified_value = modified.getMinimumMessageLength();
  current_value = current.getMinimumMessageLength();

  if (operation == KILL || operation == JOIN || operation == SPLIT) {
    if (current_value > modified_value) {
      improvement_rate = (current_value - modified_value) / fabs(current_value);
      log << "\t ... IMPROVEMENT ... (+ " << fixed << setprecision(3) 
          << 100 * improvement_rate << " %) ";
      log << "\t\t[ACCEPT]\n\n";
      current = modified;
    } else {
      log << "\t ... NO IMPROVEMENT\t\t\t[REJECT]\n\n";
    } // if (current_value > modified_value)
  } // if (operation)
}

////////////////////// TESTING FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

void TestFunctions(void)
{
  Test test;

  //test.generate_vmc();

  //test.generate_mix_vmc();

  //test.generate_bvm_sine();

  //test.bvm_sine_normalization_constant();

  //test.bvm_sine_constants();

  //test.sanity_check();

  //test.check_sufficient_stats_sine();

  //test.bvm_sine_kldiv();

  //test.bvm_sine_kldiv2();

  //test.bvm_sine_ml_estimation();

  //test.bvm_sine_all_estimation();

  //test.testing_sample_empirical_distribution();

  //test.generate_bvm_cosine();

  //test.bvm_cosine_normalization_constant();

  test.bvm_cosine_constants();
}

////////////////////// EXPERIMENTS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

void RunExperiments(struct Parameters &parameters)
{
  Experiments experiments;

  //experiments.fisher_uncertainty();

  experiments.simulate_sine(parameters);
}

/*!
 *  \brief This function sorts the elements in the list
 *  \param list a reference to a vector<double>
 *  \return the sorted list (in increasing order)
 */
Vector sort(Vector &list)
{
  int num_samples = list.size();
	Vector sortedList(list);
  std::vector<int> index(num_samples,0);
	for(int i=0; i<num_samples; i++) {
			index[i] = i;
  }
	quicksort(sortedList,index,0,num_samples-1);
  return sortedList;
}

Vector sort(Vector &list, std::vector<int> &sorted_index)
{
  int num_samples = list.size();
	Vector sortedList(list);
  sorted_index = std::vector<int>(num_samples,0);
	for(int i=0; i<num_samples; i++) {
			sorted_index[i] = i;
  }
	quicksort(sortedList,sorted_index,0,num_samples-1);
  return sortedList;
}

/*!
 *  This is an implementation of the classic quicksort() algorithm to sort a
 *  list of data values. The module uses the overloading operator(<) to 
 *  compare two Point<T> objects. 
 *  Pivot is chosen as the right most element in the list(default)
 *  This function is called recursively.
 *  \param list a reference to a Vector
 *	\param index a reference to a std::vector<int>
 *  \param left an integer
 *  \param right an integer
 */
void quicksort(Vector &list, std::vector<int> &index, int left, int right)
{
	if(left < right)
	{
		int pivotNewIndex = partition(list,index,left,right);
		quicksort(list,index,left,pivotNewIndex-1);
		quicksort(list,index,pivotNewIndex+1,right);
	}
}

/*!
 *  This function is called from the quicksort() routine to compute the new
 *  pivot index.
 *  \param list a reference to a Vector
 *	\param index a reference to a std::vector<int>
 *  \param left an integer
 *  \param right an integer
 *  \return the new pivot index
 */
int partition(Vector &list, std::vector<int> &index, int left, int right)
{
	double temp,pivotPoint = list[right];
	int storeIndex = left,temp_i;
	for(int i=left; i<right; i++) {
		if(list[i] < pivotPoint) {
			temp = list[i];
			list[i] = list[storeIndex];
			list[storeIndex] = temp;
			temp_i = index[i];
			index[i] = index[storeIndex];
			index[storeIndex] = temp_i;
			storeIndex += 1;	
		}
	}
	temp = list[storeIndex];
	list[storeIndex] = list[right];
	list[right] = temp;
	temp_i = index[storeIndex];
	index[storeIndex] = index[right];
	index[right] = temp_i;
	return storeIndex;
}

std::vector<Vector> flip(std::vector<Vector> &table)
{
  int num_rows = table.size();
  Vector empty_vector(num_rows,0);
  int num_cols = table[0].size();
  std::vector<Vector> inverted_table(num_cols,empty_vector);
  for (int i=0; i<num_cols; i++) {
    for (int j=0; j<num_rows; j++) {
      inverted_table[i][j] = table[j][i];
    }
  }
  return inverted_table;
}

/*!
 *  \brief This module computes the median of a sorted set of samples
 *  \param list a reference to a std::vector<double>
 *  \return the median value
 */
double computeMedian(Vector &list)
{
  Vector sorted_list = sort(list);
  int n = sorted_list.size();
  if (n % 2 == 1) {
    return sorted_list[n/2];
  } else {
    return (sorted_list[n/2-1]+sorted_list[n/2])/2;
  }
}

Vector computeMedians(std::vector<Vector> &table)
{
  std::vector<Vector> inverted_table = flip(table);
  int num_cols = table[0].size();
  Vector medians(num_cols,0);
  for (int i=0; i<num_cols; i++) {
    medians[i] = computeMedian(inverted_table[i]);
  }
  return medians;
}

/*!
 *  \brief This module computes the mean of a set of samples
 *  \param list a reference to a std::vector<double>
 *  \return the mean value
 */
double computeMean(Vector &list)
{
  double sum = 0;
  for (int i=0; i<list.size(); i++) {
    sum += list[i];
  }
  return sum / (double)list.size();
}

Vector computeMeans(std::vector<Vector> &table)
{
  std::vector<Vector> inverted_table = flip(table);
  int num_cols = table[0].size();
  Vector means(num_cols,0);
  for (int i=0; i<num_cols; i++) {
    means[i] = computeMean(inverted_table[i]);
  }
  return means;
}

/*!
 *  \brief Computes the variance
 */
double computeVariance(Vector &list)
{
  double mean = computeMean(list);
  double sum = 0;
  for (int i=0; i<list.size(); i++) {
    sum += (list[i]-mean) * (list[i]-mean);
  }
  return sum / (double) (list.size()-1);
}

int minimumIndex(Vector &values)
{
  int min_index = 0;
  double min_val = values[0];
  for (int i=1; i<values.size(); i++) { 
    if (values[i] <= min_val) {
      min_index = i;
      min_val = values[i];
    }
  }
  return min_index;
}

int maximumIndex(Vector &values)
{
  int max_index = 0;
  double max_val = values[0];
  for (int i=1; i<values.size(); i++) { 
    if (values[i] > max_val) {
      max_index = i;
      max_val = values[i];
    }
  }
  return max_index;
}

void chisquare_hypothesis_testing(
  std::vector<struct EstimatesSine> &all_bvm_sine_estimates,
  Vector &statistics,
  Vector &pvalues
) {
  int num_methods = all_bvm_sine_estimates.size();
  statistics = Vector(num_methods,0);
  pvalues = Vector(num_methods,0);

  chi_squared chisq(5);

  BVM_Sine bvm_sine_ml(
    all_bvm_sine_estimates[MLE].mu1,
    all_bvm_sine_estimates[MLE].mu2,
    all_bvm_sine_estimates[MLE].kappa1,
    all_bvm_sine_estimates[MLE].kappa2,
    all_bvm_sine_estimates[MLE].lambda
  );
  double bvm_sine_ml_negloglike = all_bvm_sine_estimates[MLE].negloglike;
  //cout << "BVM_Sine (MLE) negloglike: " << bvm_sine_ml_negloglike << endl;

  double bvm_sine_negloglike,log_ratio_statistic,pvalue;
  for (int i=0; i<num_methods; i++) {
    //if (i != PMLE || i != MLE) {
    if (i != MLE) {
      // null: Kent(i)
      BVM_Sine bvm_sine_ml(
        all_bvm_sine_estimates[i].mu1,
        all_bvm_sine_estimates[i].mu2,
        all_bvm_sine_estimates[i].kappa1,
        all_bvm_sine_estimates[i].kappa2,
        all_bvm_sine_estimates[i].lambda
      );
      bvm_sine_negloglike = all_bvm_sine_estimates[i].negloglike;
      log_ratio_statistic = 2 * (bvm_sine_negloglike - bvm_sine_ml_negloglike);
      if (log_ratio_statistic < 0) {
        //assert(fabs(log_ratio_statistic) < 1e-4);
        log_ratio_statistic = fabs(log_ratio_statistic);
      }
      pvalue = 1 - cdf(chisq,log_ratio_statistic);
    } else { 
      log_ratio_statistic = 0;
      pvalue = 1;
    } // if()
    cout << "(null: BVM_Sine[" << i <<"]) negloglike: " << bvm_sine_negloglike
         << "; log_ratio_statistic: " << log_ratio_statistic
         << "; pvalue: " << pvalue << endl;
    statistics[i] = log_ratio_statistic;
    pvalues[i] = pvalue;
  } // for()
}

////////////////////// BVM SPECIFIC FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

double accept_reject_fval_unimodal_marginal_sine(
  double &theta, 
  double &kappa, 
  double &mu1, 
  double &kappa1, 
  double &kappa2, 
  double &lambda
) {
  double tmp = lambda * sin(theta - mu1);
  double asq = kappa2 * kappa2 + tmp * tmp;
  double a = sqrt(asq);

  double fval = 0;
  fval += computeLogModifiedBesselFirstKind(0,a);
  fval += computeLogModifiedBesselFirstKind(0,kappa);
  fval += ((kappa1 - kappa) * cos(theta - mu1));

  return fval;
}

double accept_reject_fval_bimodal_marginal_sine(
  double &theta, 
  double &kappa, 
  double &mu1, 
  double &kappa1, 
  double &kappa2, 
  double &lambda,
  double &mode1,
  double &mode2
) {
  double tmp = lambda * sin(theta - mu1);
  double asq = kappa2 * kappa2 + tmp * tmp;
  double a = sqrt(asq);

  double fval = 0;
  fval += computeLogModifiedBesselFirstKind(0,a);
  fval += (kappa1 * cos(theta - mu1));

  Vector weights(2,0.5);
  vMC vmc1(mode1,kappa);
  vMC vmc2(mode2,kappa);
  std::vector<vMC> components;
  components.push_back(vmc1); components.push_back(vmc2);
  Mixture_vMC mix(2,weights,components);
  
  fval -= mix.log_density(theta);

  return fval;
}

double accept_reject_fval_unimodal_marginal_cosine(
  double &theta, 
  double &kappa, 
  double &mu1, 
  double &kappa1, 
  double &kappa2, 
  double &kappa3
) {
  double k23sq = kappa2 * kappa2 + kappa3 * kappa3 
                 - (2 * kappa2 * kappa3 * cos(theta-mu1));
  double k23 = sqrt(k23sq);

  double fval = 0;
  fval += computeLogModifiedBesselFirstKind(0,k23);
  fval += computeLogModifiedBesselFirstKind(0,kappa);
  fval += ((kappa1 - kappa) * cos(theta - mu1));

  return fval;
}

double accept_reject_fval_bimodal_marginal_cosine(
  double &theta, 
  double &kappa, 
  double &mu1, 
  double &kappa1, 
  double &kappa2, 
  double &kappa3,
  double &mode1,
  double &mode2
) {
  double k23sq = kappa2 * kappa2 + kappa3 * kappa3 
                 - (2 * kappa2 * kappa3 * cos(theta-mu1));
  double k23 = sqrt(k23sq);

  double fval = 0;
  fval += computeLogModifiedBesselFirstKind(0,k23);
  fval += (kappa1 * cos(theta - mu1));

  Vector weights(2,0.5);
  vMC vmc1(mode1,kappa);
  vMC vmc2(mode2,kappa);
  std::vector<vMC> components;
  components.push_back(vmc1); components.push_back(vmc2);
  Mixture_vMC mix(2,weights,components);
  
  fval -= mix.log_density(theta);

  return fval;
}

double banerjee_approx(double &rbar)
{
  assert(rbar < 1);

  double rbarsq = rbar * rbar;
  double denom = 1 - rbarsq;
  double num = rbar * (denom + 1);
  return num/denom;
}

/* Sufficient statistics for the Sine model */
// data = angle_pairs
void computeSufficientStatisticsSineNotParallel(
  std::vector<Vector> &data,
  struct SufficientStatisticsSine &suff_stats
) {
  suff_stats.N = data.size();

  suff_stats.cost1 = 0; suff_stats.cost2 = 0;
  suff_stats.sint1 = 0; suff_stats.sint2 = 0;
  suff_stats.sint1sint2 = 0; suff_stats.sint1cost2 = 0;
  suff_stats.cost1sint2 = 0; suff_stats.cost1cost2 = 0;

  for (int i=0; i<data.size(); i++) {
    double t1 = data[i][0];
    double cost1 = cos(t1) ;
    double sint1 = sin(t1) ;
    double t2 = data[i][1];
    double cost2 = cos(t2) ;
    double sint2 = sin(t2) ;
 
    suff_stats.cost1 += cost1;
    suff_stats.cost2 += cost2;
    suff_stats.sint1 += sint1;
    suff_stats.sint2 += sint2;

    suff_stats.sint1sint2 += sint1 * sint2;
    suff_stats.sint1cost2 += sint1 * cost2;
    suff_stats.cost1sint2 += cost1 * sint2;
    suff_stats.cost1cost2 += cost1 * cost2;
  } // for()

  //print(suff_stats);
}

// data = angle_pairs
void computeSufficientStatisticsSine(
  std::vector<Vector> &data,
  struct SufficientStatisticsSine &suff_stats
) {

  if (NUM_THREADS == 1) {
    return computeSufficientStatisticsSineNotParallel(data,suff_stats);
  }

  // empty the existing sufficient stats
  suff_stats.N = 0;
  suff_stats.cost1 = 0; suff_stats.cost2 = 0;
  suff_stats.sint1 = 0; suff_stats.sint2 = 0;
  suff_stats.sint1sint2 = 0; suff_stats.sint1cost2 = 0;
  suff_stats.cost1sint2 = 0; suff_stats.cost1cost2 = 0;

  std::vector<struct SufficientStatisticsSine> suff_stats_vector(NUM_THREADS);
  for (int i=0; i<NUM_THREADS; i++) {
    suff_stats_vector[i] = suff_stats;
  } // for (i)

  int tid;
  #pragma omp parallel if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) private(tid) 
  {
    tid = omp_get_thread_num();
    #pragma omp for
    for (int i=0; i<data.size(); i++) {
      suff_stats_vector[tid].N += 1;

      double t1 = data[i][0];
      double cost1 = cos(t1) ;
      double sint1 = sin(t1) ;
      double t2 = data[i][1];
      double cost2 = cos(t2) ;
      double sint2 = sin(t2) ;

      suff_stats_vector[tid].cost1 += cost1;
      suff_stats_vector[tid].cost2 += cost2;
      suff_stats_vector[tid].sint1 += sint1;
      suff_stats_vector[tid].sint2 += sint2;

      suff_stats_vector[tid].sint1sint2 += sint1 * sint2;
      suff_stats_vector[tid].sint1cost2 += sint1 * cost2;
      suff_stats_vector[tid].cost1sint2 += cost1 * sint2;
      suff_stats_vector[tid].cost1cost2 += cost1 * cost2;
    } // for (i)
  } // pragma omp parallel

  for (int i=0; i<NUM_THREADS; i++) {
    suff_stats.N += suff_stats_vector[i].N;

    suff_stats.cost1 += suff_stats_vector[i].cost1;
    suff_stats.cost2 += suff_stats_vector[i].cost2;
    suff_stats.sint1 += suff_stats_vector[i].sint1;
    suff_stats.sint2 += suff_stats_vector[i].sint2;

    suff_stats.sint1sint2 += suff_stats_vector[i].sint1sint2;
    suff_stats.sint1cost2 += suff_stats_vector[i].sint1cost2;
    suff_stats.cost1sint2 += suff_stats_vector[i].cost1sint2;
    suff_stats.cost1cost2 += suff_stats_vector[i].cost1cost2;

    //cout << "suff_stats[" << i << "]:\n"; print(suff_stats_vector[i]);
  } // for (i)

  //print(suff_stats);
}

// data = angle_pairs
void computeSufficientStatisticsSineNotParallel(
  std::vector<Vector> &data,
  struct SufficientStatisticsSine &suff_stats,
  Vector &weights
) {
  suff_stats.N = 0;
  suff_stats.cost1 = 0; suff_stats.cost2 = 0;
  suff_stats.sint1 = 0; suff_stats.sint2 = 0;
  suff_stats.sint1sint2 = 0; suff_stats.sint1cost2 = 0;
  suff_stats.cost1sint2 = 0; suff_stats.cost1cost2 = 0;

  for (int i=0; i<data.size(); i++) {
    suff_stats.N += weights[i];

    double t1 = data[i][0];
    double cost1 = cos(t1) ;
    double sint1 = sin(t1) ;
    double t2 = data[i][1];
    double cost2 = cos(t2) ;
    double sint2 = sin(t2) ;
 
    suff_stats.cost1 += (weights[i] * cost1);
    suff_stats.cost2 += (weights[i] * cost2);
    suff_stats.sint1 += (weights[i] * sint1);
    suff_stats.sint2 += (weights[i] * sint2);

    suff_stats.sint1sint2 += (weights[i] * sint1 * sint2);
    suff_stats.sint1cost2 += (weights[i] * sint1 * cost2);
    suff_stats.cost1sint2 += (weights[i] * cost1 * sint2);
    suff_stats.cost1cost2 += (weights[i] * cost1 * cost2);
  } // for()

  //print(suff_stats);
}

// data = angle_pairs
void computeSufficientStatisticsSine(
  std::vector<Vector> &data,
  struct SufficientStatisticsSine &suff_stats,
  Vector &weights
) {

  if (NUM_THREADS == 1) {
    return computeSufficientStatisticsSineNotParallel(data,suff_stats,weights);
  }

  // empty the existing sufficient stats
  suff_stats.N = 0;
  suff_stats.cost1 = 0; suff_stats.cost2 = 0;
  suff_stats.sint1 = 0; suff_stats.sint2 = 0;
  suff_stats.sint1sint2 = 0; suff_stats.sint1cost2 = 0;
  suff_stats.cost1sint2 = 0; suff_stats.cost1cost2 = 0;

  std::vector<struct SufficientStatisticsSine> suff_stats_vector(NUM_THREADS);
  for (int i=0; i<NUM_THREADS; i++) {
    suff_stats_vector[i] = suff_stats;
  } // for (i)

  int tid;
  #pragma omp parallel if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) private(tid) 
  {
    tid = omp_get_thread_num();
    #pragma omp for
    for (int i=0; i<data.size(); i++) {
      suff_stats_vector[tid].N += weights[i];

      double t1 = data[i][0];
      double cost1 = cos(t1) ;
      double sint1 = sin(t1) ;
      double t2 = data[i][1];
      double cost2 = cos(t2) ;
      double sint2 = sin(t2) ;

      suff_stats_vector[tid].cost1 += (weights[i] * cost1);
      suff_stats_vector[tid].cost2 += (weights[i] * cost2);
      suff_stats_vector[tid].sint1 += (weights[i] * sint1);
      suff_stats_vector[tid].sint2 += (weights[i] * sint2);

      suff_stats_vector[tid].sint1sint2 += (weights[i] * sint1 * sint2);
      suff_stats_vector[tid].sint1cost2 += (weights[i] * sint1 * cost2);
      suff_stats_vector[tid].cost1sint2 += (weights[i] * cost1 * sint2);
      suff_stats_vector[tid].cost1cost2 += (weights[i] * cost1 * cost2);
    } // for (i)
  } // pragma omp parallel

  for (int i=0; i<NUM_THREADS; i++) {
    suff_stats.N += suff_stats_vector[i].N;

    suff_stats.cost1 += suff_stats_vector[i].cost1;
    suff_stats.cost2 += suff_stats_vector[i].cost2;
    suff_stats.sint1 += suff_stats_vector[i].sint1;
    suff_stats.sint2 += suff_stats_vector[i].sint2;

    suff_stats.sint1sint2 += suff_stats_vector[i].sint1sint2;
    suff_stats.sint1cost2 += suff_stats_vector[i].sint1cost2;
    suff_stats.cost1sint2 += suff_stats_vector[i].cost1sint2;
    suff_stats.cost1cost2 += suff_stats_vector[i].cost1cost2;

    //cout << "suff_stats[" << i << "]:\n"; print(suff_stats_vector[i]);
  } // for (i)

  //print(suff_stats);
}

double ConstraintSine(const Vector &x, std::vector<double> &grad, void *data)
{
    double k1 = x[2];
    double k2 = x[3];
    double lam = x[4];
    return (lam*lam - k1*k2);
    //return (2 * x[1] - x[0]);
}

// g(theta2|theta1)
vMC getConditionalDensitySine(
  double &theta1, double &mu1, double &mu2, double &kappa2, double &lambda
) {
  double lambda_sine = lambda * sin(theta1 - mu1);
  double beta = atan2(lambda_sine,kappa2);
  double m = mu2 + beta;
  if (m < 0) m += (2 * PI);

  double asq = kappa2 * kappa2 + lambda_sine * lambda_sine;
  double k = sqrt(asq);
  vMC vmc(m,k);
  return vmc;
}

/* Sufficient statistics for the Cosine model */

double cosine_correlation(double kappa1, double kappa2, double kappa3)
{
  double diff1 = (kappa1 - kappa3);
  double diff2 = (kappa2 - kappa3);
  double product = diff1 * diff2;
  return -(kappa3/sqrt(product));
}

// data = angle_pairs
void computeSufficientStatisticsCosineNotParallel(
  std::vector<Vector> &data,
  struct SufficientStatisticsCosine &suff_stats
) {
  suff_stats.N = data.size();

  suff_stats.cost1 = 0; suff_stats.cost2 = 0;
  suff_stats.sint1 = 0; suff_stats.sint2 = 0;
  suff_stats.cost1_t2 = 0; suff_stats.sint1_t2 = 0;

  for (int i=0; i<data.size(); i++) {
    double t1 = data[i][0];
    double cost1 = cos(t1) ;
    double sint1 = sin(t1) ;
    double t2 = data[i][1];
    double cost2 = cos(t2) ;
    double sint2 = sin(t2) ;
 
    suff_stats.cost1 += cost1;
    suff_stats.cost2 += cost2;
    suff_stats.sint1 += sint1;
    suff_stats.sint2 += sint2;

    suff_stats.cost1_t2 += cos(t1-t2);
    suff_stats.sint1_t2 += sin(t1-t2);
  } // for()

  //print(suff_stats);
}

// data = angle_pairs
void computeSufficientStatisticsCosine(
  std::vector<Vector> &data,
  struct SufficientStatisticsCosine &suff_stats
) {

  if (NUM_THREADS == 1) {
    return computeSufficientStatisticsCosineNotParallel(data,suff_stats);
  }

  // empty the existing sufficient stats
  suff_stats.N = 0;
  suff_stats.cost1 = 0; suff_stats.cost2 = 0;
  suff_stats.sint1 = 0; suff_stats.sint2 = 0;
  suff_stats.cost1_t2 = 0; suff_stats.sint1_t2 = 0;

  std::vector<struct SufficientStatisticsCosine> suff_stats_vector(NUM_THREADS);
  for (int i=0; i<NUM_THREADS; i++) {
    suff_stats_vector[i] = suff_stats;
  } // for (i)

  int tid;
  #pragma omp parallel if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) private(tid) 
  {
    tid = omp_get_thread_num();
    #pragma omp for
    for (int i=0; i<data.size(); i++) {
      suff_stats_vector[tid].N += 1;

      double t1 = data[i][0];
      double cost1 = cos(t1) ;
      double sint1 = sin(t1) ;
      double t2 = data[i][1];
      double cost2 = cos(t2) ;
      double sint2 = sin(t2) ;

      suff_stats_vector[tid].cost1 += cost1;
      suff_stats_vector[tid].cost2 += cost2;
      suff_stats_vector[tid].sint1 += sint1;
      suff_stats_vector[tid].sint2 += sint2;

      suff_stats_vector[tid].cost1_t2 += cos(t1-t2);
      suff_stats_vector[tid].sint1_t2 += sin(t1-t2);
    } // for (i)
  } // pragma omp parallel

  for (int i=0; i<NUM_THREADS; i++) {
    suff_stats.N += suff_stats_vector[i].N;

    suff_stats.cost1 += suff_stats_vector[i].cost1;
    suff_stats.cost2 += suff_stats_vector[i].cost2;
    suff_stats.sint1 += suff_stats_vector[i].sint1;
    suff_stats.sint2 += suff_stats_vector[i].sint2;

    suff_stats.cost1_t2 += suff_stats_vector[i].cost1_t2;
    suff_stats.sint1_t2 += suff_stats_vector[i].sint1_t2;

    //cout << "suff_stats[" << i << "]:\n"; print(suff_stats_vector[i]);
  } // for (i)

  //print(suff_stats);
}

// data = angle_pairs
void computeSufficientStatisticsCosineNotParallel(
  std::vector<Vector> &data,
  struct SufficientStatisticsCosine &suff_stats,
  Vector &weights
) {
  suff_stats.N = 0;
  suff_stats.cost1 = 0; suff_stats.cost2 = 0;
  suff_stats.sint1 = 0; suff_stats.sint2 = 0;
  suff_stats.cost1_t2 = 0; suff_stats.sint1_t2 = 0;

  for (int i=0; i<data.size(); i++) {
    suff_stats.N += weights[i];

    double t1 = data[i][0];
    double cost1 = cos(t1) ;
    double sint1 = sin(t1) ;
    double t2 = data[i][1];
    double cost2 = cos(t2) ;
    double sint2 = sin(t2) ;
 
    suff_stats.cost1 += (weights[i] * cost1);
    suff_stats.cost2 += (weights[i] * cost2);
    suff_stats.sint1 += (weights[i] * sint1);
    suff_stats.sint2 += (weights[i] * sint2);

    suff_stats.cost1_t2 += (weights[i] * cos(t1-t2));
    suff_stats.sint1_t2 += (weights[i] * sin(t1-t2));
  } // for()

  //print(suff_stats);
}

// data = angle_pairs
void computeSufficientStatisticsCosine(
  std::vector<Vector> &data,
  struct SufficientStatisticsCosine &suff_stats,
  Vector &weights
) {

  if (NUM_THREADS == 1) {
    return computeSufficientStatisticsCosineNotParallel(data,suff_stats,weights);
  }

  // empty the existing sufficient stats
  suff_stats.N = 0;
  suff_stats.cost1 = 0; suff_stats.cost2 = 0;
  suff_stats.sint1 = 0; suff_stats.sint2 = 0;
  suff_stats.cost1_t2 = 0; suff_stats.sint1_t2 = 0;

  std::vector<struct SufficientStatisticsCosine> suff_stats_vector(NUM_THREADS);
  for (int i=0; i<NUM_THREADS; i++) {
    suff_stats_vector[i] = suff_stats;
  } // for (i)

  int tid;
  #pragma omp parallel if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) private(tid) 
  {
    tid = omp_get_thread_num();
    #pragma omp for
    for (int i=0; i<data.size(); i++) {
      suff_stats_vector[tid].N += weights[i];

      double t1 = data[i][0];
      double cost1 = cos(t1) ;
      double sint1 = sin(t1) ;
      double t2 = data[i][1];
      double cost2 = cos(t2) ;
      double sint2 = sin(t2) ;

      suff_stats_vector[tid].cost1 += (weights[i] * cost1);
      suff_stats_vector[tid].cost2 += (weights[i] * cost2);
      suff_stats_vector[tid].sint1 += (weights[i] * sint1);
      suff_stats_vector[tid].sint2 += (weights[i] * sint2);

      suff_stats_vector[tid].cost1_t2 += (weights[i] * cos(t1-t2));
      suff_stats_vector[tid].sint1_t2 += (weights[i] * sin(t1-t2));
    } // for (i)
  } // pragma omp parallel

  for (int i=0; i<NUM_THREADS; i++) {
    suff_stats.N += suff_stats_vector[i].N;

    suff_stats.cost1 += suff_stats_vector[i].cost1;
    suff_stats.cost2 += suff_stats_vector[i].cost2;
    suff_stats.sint1 += suff_stats_vector[i].sint1;
    suff_stats.sint2 += suff_stats_vector[i].sint2;

    suff_stats.cost1_t2 += suff_stats_vector[i].cost1_t2;
    suff_stats.sint1_t2 += suff_stats_vector[i].sint1_t2;

    //cout << "suff_stats[" << i << "]:\n"; print(suff_stats_vector[i]);
  } // for (i)

  //print(suff_stats);
}

double ConstraintCosine(const Vector &x, std::vector<double> &grad, void *data)
{
    double k1 = x[2];
    double k2 = x[3];
    double k3 = x[4];
    return (k3*(k1+k2) - k1*k2);
}

double norm_constant_integral (double *k, size_t dim, void *params)
{
  double *params2 = (double *)params;
  /*cout << "params: ";
  for (int i=0; i<5; i++) {
    cout << params2[i] << "; ";
  } cout << endl;*/

  double mu1 = params2[0];
  double mu2 = params2[1];
  double kappa1 = params2[2];
  double kappa2 = params2[3];
  double kappa3 = params2[4];

  double x = k[0] - mu1;
  double y = k[1] - mu2;

  double t1 = kappa1 * cos(x);
  double t2 = kappa2 * cos(y);
  double t3 = kappa3 * cos(x-y);

  double exponent = t1 + t2 - t3;
  return exp(exponent);
}

double integration(double mu1, double mu2, double kappa1, double kappa2, double kappa3)
{
  double res, err;

  double params[5] = {mu1,mu2,kappa1,kappa2,kappa3};
  double xl[2] = { 0, 0};
  double xu[2] = { 2*PI, 2*PI};

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = { &norm_constant_integral, 2, params };

  size_t calls = 1000000;

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  {
    gsl_monte_plain_state *s = gsl_monte_plain_alloc (2);
    gsl_monte_plain_integrate (&G, xl, xu, 2, calls, r, s, 
                               &res, &err);
    gsl_monte_plain_free (s);

    cout << "integral value 1: " << res << endl;
  }

  {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (2);
    gsl_monte_miser_integrate (&G, xl, xu, 2, calls, r, s,
                               &res, &err);
    gsl_monte_miser_free (s);

    cout << "integral value 2: " << res << endl;
  }

  {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

    gsl_monte_vegas_integrate (&G, xl, xu, 2, 10000, r, s,
                               &res, &err);

    do
      {
        gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s,
                                   &res, &err);
        //printf ("result = % .6f sigma = % .6f "
        //        "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    cout << "integral value 3: " << res << endl;

    gsl_monte_vegas_free (s);
  }

  gsl_rng_free (r);
}

