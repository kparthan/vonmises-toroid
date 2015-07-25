#include "Support.h"
#include "Test.h"
#include "UniformRandomNumberGenerator.h"

Vector XAXIS,YAXIS,ZAXIS;
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
    if (estimation_method.compare("moment") == 0) {
      ESTIMATION = MOMENT;
    } else if (estimation_method.compare("mle") == 0) {
      ESTIMATION = MLE;
    } else if (estimation_method.compare("map") == 0) {
      ESTIMATION = MAP;
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
  cout << "m1_est: " << estimates.mu1 * 180/PI << endl;
  cout << "m2_est: " << estimates.mu2 * 180/PI << endl;
  cout << "k1_est: " << estimates.kappa1 << endl;
  cout << "k2_est: " << estimates.kappa2 << endl;
  cout << "lambda_est: " << estimates.lambda << endl;
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

double computeLogModifiedBesselFirstKind(double alpha, double x)
{
  if (!(alpha >= 0 && fabs(x) >= 0)) {
    cout << "Error logModifiedBesselFirstKind: (alpha,x) = (" << alpha << "," << x << ")\n";
    exit(1);
  }
  if (fabs(x) <= TOLERANCE) {
    return -LARGE_NUMBER;
  } 

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

// theta in radians ...
Matrix rotate_about_arbitrary_axis(Vector &axis, double theta)
{
  Matrix K = ZeroMatrix(3,3);
  K(0,1) = -axis[2];
  K(0,2) = axis[1];
  K(1,2) = -axis[0];
  K(1,0) = -K(0,1);
  K(2,0) = -K(0,2);
  K(2,1) = -K(1,2);

  Matrix Ksq = prod(K,K);
  Matrix I = IdentityMatrix(3,3);
  Matrix sinm = sin(theta) * K;
  Matrix cosm = (1-cos(theta)) * Ksq;
  Matrix tmp = sinm + cosm;
  Matrix R = I + tmp;
  return R;
}

// anti-clockwise rotation about +X
Matrix rotate_about_xaxis(double theta)
{
  Matrix r = IdentityMatrix(3,3);
  r(1,1) = cos(theta);
  r(1,2) = -sin(theta);
  r(2,1) = -r(1,2); // sin(theta)
  r(2,2) = r(1,1);  // cos(theta)
  return r;
}

// anti-clockwise rotation about +Y
Matrix rotate_about_yaxis(double theta)
{
  Matrix r = IdentityMatrix(3,3);
  r(0,0) = cos(theta);
  r(0,2) = sin(theta);
  r(2,0) = -r(0,2); // -sin(theta)
  r(2,2) = r(0,0);  // cos(theta)
  return r;
}

// anti-clockwise rotation about +Z
Matrix rotate_about_zaxis(double theta)
{
  Matrix r = IdentityMatrix(3,3);
  r(0,0) = cos(theta);
  r(0,1) = -sin(theta);
  r(1,0) = -r(0,1); // sin(theta)
  r(1,1) = r(0,0);  // cos(theta)
  return r;
}

/*
 *  \brief Transformation of x using T
 *  \param x a reference to a vector<vector<double> >
 *  \param T a reference to a Matrix<double>
 *  \return the transformed vector list
 */
std::vector<Vector> transform(
  std::vector<Vector> &x, 
  Matrix &T
) {
  std::vector<Vector> y(x.size());
  for (int i=0; i<x.size(); i++) {
    y[i] = prod(T,x[i]);
  }
  return y;
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

  test.bvm_sine_ml_estimation();

  //test.generate_bvm_cosine();
}

////////////////////// EXPERIMENTS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

void RunExperiments(int iterations)
{
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
  long double min_val = values[0];
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

// data = angle_pairs
void computeSufficientStatisticsSine(
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

/*
  cout << "sufficient stats sine:\n";
  cout << "cost1: " << suff_stats.cost1 << endl;
  cout << "sint1: " << suff_stats.sint1 << endl;
  cout << "cost2: " << suff_stats.cost2 << endl;
  cout << "sint2: " << suff_stats.sint2 << endl;
  cout << "sint1 sint2: " << suff_stats.sint1sint2 << endl;
  cout << "sint1 cost2: " << suff_stats.sint1cost2 << endl;
  cout << "cost1 sint2: " << suff_stats.cost1sint2 << endl;
  cout << "cost1 cost2: " << suff_stats.cost1cost2 << endl;
*/  
}

double ConstraintSine(const Vector &x, std::vector<double> &grad, void *data)
{
    double k1 = x[2];
    double k2 = x[3];
    double lam = x[4];
    return (lam*lam - k1*k2);
    //return (2 * x[1] - x[0]);
}

