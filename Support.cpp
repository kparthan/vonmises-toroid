#include "Support.h"
#include "Test.h"
#include "Experiments.h"
#include "Structure.h"
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
int PRIOR;
string tracking_file;

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
       ("prior",value<int>(&prior),"vMF Kappa prior (2D/3D)")
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
    if (pdf.compare("vmf") == 0) {
      DISTRIBUTION = VMF;
    } else if (pdf.compare("kent") == 0) {
      DISTRIBUTION = KENT;
    } else {
      cout << "Unsupported distribution ...\n";
      Usage(argv[0],desc);
    }
  } else {
    DISTRIBUTION = KENT;  // default distribution
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

  tracking_file = "kent_" + criterion + "_" + estimation_method;

  if (!vm.count("prior")) {
    PRIOR = 3;
  } else {
    PRIOR = prior;
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
void print(ostream &os, Vector &v, int precision)
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

void print(string &type, struct Estimates &estimates)
{
  //cout << "MIX_ID: " << MIXTURE_ID-1 << endl;
  cout << "TYPE: " << type << endl;
  Vector spherical(3,0);
  cartesian2spherical(estimates.mean,spherical);
  cout << "m0_est: "; print(cout,estimates.mean,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(estimates.major_axis,spherical);
  cout << "m1_est: "; print(cout,estimates.major_axis,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cartesian2spherical(estimates.minor_axis,spherical);
  cout << "m2_est: "; print(cout,estimates.minor_axis,3);
  cout << "\t(" << spherical[1]*180/PI << "," << spherical[2]*180/PI << ")\n";
  cout << "psi_est: " << estimates.psi * 180/PI << "; ";
  cout << "alpha_est: " << estimates.alpha * 180/PI << "; ";
  cout << "eta_est: " << estimates.eta * 180/PI << "; ";
  cout << "kappa_est: " << estimates.kappa << "; beta_est: " << estimates.beta 
       << "; ecc_est: " << 2 * estimates.beta / estimates.kappa << endl;
  // check:
  /*cout << fixed << scientific
       << "m0_est . m1_est = " << computeDotProduct(estimates.mean,estimates.major_axis)
       << "; m0_est . m2_est = " << computeDotProduct(estimates.mean,estimates.minor_axis)
       << "; m1_est . m2_est = " << computeDotProduct(estimates.major_axis,estimates.minor_axis)
       << endl;*/
}

void print(string &type, struct Estimates_vMF &estimates)
{
  cout << "TYPE: " << type << endl;
  Vector spherical(3,0);
  cout << "m0_est: "; print(cout,estimates.mean,3);
  cout << "\t(" << estimates.alpha*180/PI << "," << estimates.eta*180/PI << ")\n";
  cout << "kappa_est: " << estimates.kappa << endl;
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

void computeLambertProjection(Vector &cartesian, Vector &projection)
{
  Vector spherical(3,0);
  cartesian2spherical(cartesian,spherical);
  double theta = spherical[1];
  double phi = spherical[2];
  double rho = 2 * sin(theta/2);
  projection[0] = rho * cos(phi);
  projection[1] = rho * sin(phi);
}

void computeLambertProjection(std::vector<Vector> &sample)
{
  Vector proj(2,0);
  std::vector<Vector> projections(sample.size(),proj);
  for (int i=0; i<sample.size(); i++) {
    computeLambertProjection(sample[i],projections[i]);
  }
  string output = "./visualize/sampled_data/transformed.dat";
  writeToFile(output,projections);
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

/*!
 *  Find the roots of a quadratic: ax^2 + bx + c = 0
 */
void solveQuadratic(
  Vector &roots, 
  double a,
  double b, 
  double c
) {
  double D = sqrt(b*b-4*a*c);
  roots[0] = (-b + D) / (2*a);
  roots[1] = (-b - D) / (2*a);
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

////////////////////// GEOMETRY FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

std::vector<Vector> load_data_table(string &file_name)
{
  std::vector<Vector> sample;
  ifstream file(file_name.c_str());
  string line;
  Vector numbers(3,0),unit_vector(3,0);
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
double determinant(Matrix &m)
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
/*Vector computeVectorSum(std::vector<Vector> &sample) 
{
  int d = sample[0].size();
  Vector sum(d,0);
  for (int i=0; i<sample.size(); i++) {
    for (int j=0; j<d; j++) {
      sum[j] += sample[i][j];
    }
  }
  return sum;
}*/

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

/*Vector computeVectorSum(std::vector<Vector> &sample, Vector &weights, double &Neff) 
{
  int d = sample[0].size();
  Vector sum(d,0);
  Neff = 0;
  for (int i=0; i<sample.size(); i++) {
    for (int j=0; j<d; j++) {
      sum[j] += sample[i][j] * weights[i];
    }
    Neff += weights[i];
  }
  return sum;
}*/

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

/*!
 *  Matrix to rotate FROM the standard frame of reference.
 */
Matrix computeOrthogonalTransformation(
  Vector &mean,
  Vector &major_axis
) {
  Matrix r1 = align_vector_with_xaxis(mean);
  Matrix r1_inv = trans(r1);
  Vector mj_xy = prod(r1,major_axis);
  Vector spherical(3,0);
  cartesian2spherical(mj_xy,spherical);
  double psi = spherical[2];
  Matrix r2 = rotate_about_xaxis(psi);
  Matrix r = prod(r1_inv,r2);
  return r;
}

/*!
 *  Matrix to rotate FROM the standard frame of reference.
 */
Matrix computeOrthogonalTransformation(double psi, double alpha, double eta)
{
  Matrix r1 = rotate_about_xaxis(psi);
  Matrix r2 = rotate_about_zaxis(alpha);
  Matrix r3 = rotate_about_xaxis(eta);
  Matrix tmp = prod(r3,r2);
  Matrix r = prod(tmp,r1);
  return r;
}

void computeOrthogonalTransformation(
  Vector &mean, 
  Vector &major_axis,
  double &psi,
  double &alpha,
  double &eta
) {
  Vector spherical(3,0);
  cartesian2spherical(mean,spherical);
  alpha = spherical[1];
  eta = spherical[2];

  Matrix r1 = align_vector_with_xaxis(mean);
  Vector mj_xy = prod(r1,major_axis);
  cartesian2spherical(mj_xy,spherical);
  psi = spherical[2];
}

Matrix align_xaxis_with_vector(Vector &y)
{
  Matrix r = align_vector_with_xaxis(y);
  return trans(r);
}

Matrix align_vector_with_xaxis(Vector &y)
{
  Vector spherical(3,0);
  cartesian2spherical(y,spherical);
  double alpha = spherical[1];
  double eta = spherical[2];
  return align_vector_with_xaxis(alpha,eta);
}

// alpha = co-latitude
// eta = longitude
Matrix align_vector_with_xaxis(double alpha, double eta)
{
  Matrix r = ZeroMatrix(3,3);

  r(0,0) = cos(alpha);
  r(0,1) = sin(alpha) * cos(eta);
  r(0,2) = sin(alpha) * sin(eta);

  r(1,0) = -sin(alpha);
  r(1,1) = cos(alpha) * cos(eta);
  r(1,2) = cos(alpha) * sin(eta);

  r(2,0) = 0;
  r(2,1) = -sin(eta);
  r(2,2) = cos(eta);

  return r;
}

void generateRandomOrthogonalVectors(
  Vector &mean,
  Vector &major_axis,
  Vector &minor_axis
) {
  double psi = (2 * PI) * uniform_random();
  double alpha = PI * uniform_random();
  double eta = (2 * PI) * uniform_random();

  Matrix r = computeOrthogonalTransformation(psi,alpha,eta);

  mean = prod(r,XAXIS);
  major_axis = prod(r,YAXIS);
  minor_axis = prod(r,ZAXIS);
}

void generateRandomOrientations(
  double &psi, double &alpha, double &eta
) {
  psi = (2 * PI) * uniform_random();
  alpha = PI * uniform_random();
  eta = (2 * PI) * uniform_random();
}

Matrix generateRandomCovarianceMatrix(int D)
{
  // generate random orthogonal vectors
  Matrix V = ZeroMatrix(D,D);
  if (D == 2) {
    Vector x(2,0);
    double phi = (2 * PI) * uniform_random();
    x[0] = cos(phi); x[1] = sin(phi);
    V(0,0) = x[0]; V(0,1) = -x[1];
    V(1,0) = x[1]; V(1,1) = x[0];
  } else if (D == 3) {
    Vector m0,m1,m2;
    std::vector<Vector> axis(3);
    generateRandomOrthogonalVectors(axis[0],axis[1],axis[2]);
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        V(j,i) = axis[i][j];
      }
    }
  } else {
    cout << "D = " << D << " not supported ...\n";
    exit(1);
  }
  Matrix diag = ZeroMatrix(D,D);
  for (int i=0; i<D; i++) {
    diag(i,i) = uniform_random() * 10;
  }
  Matrix tmp = prod(V,diag);
  Matrix Vt = trans(V);
  Matrix cov = prod(tmp,Vt);
  return cov;
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

/*
 *  Dawson's integral required to compute FB4 normalization constant
 */
double computeDawsonsIntegral(double limit)
{
  //state_type x (1,0.0);
  std::vector<double> x (1,0.0);
  //integrate(rhs,x,0.0,limit,0.1,track);
  integrate(rhs,x,0.0,limit,0.1);
  //cout << "ans: " << x[0] << endl;
  return x[0];
}

void rhs(const std::vector<double> &x, std::vector<double> &dxdt, const double t)
{
    dxdt[0] = 1 - (2 * t * x[0]); 
}

void track(const std::vector<double> &x, const double t)
{
    cout << t << "\t" << x[0] << endl;
}

double Constraint2(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    //double k = x[0];
    //double b = x[1];
    return (2 * x[1] - x[0]);
}

double Constraint5(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    //double k = x[3];
    //double b = x[4];
    return (2 * x[4] - x[3]);
}

double Constraint5_2(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    //double k = x[3];
    //double b = x[4];
    return (TOLERANCE * x[3] - 2 * x[4]);
}

double test_function_integral(double x, void * params) 
{
  double alpha = *(double *) params;
  double f = log(alpha*x) / sqrt(x);
  return f;
}

double integrate_wrt_theta(double theta, void *params)
{
  //cout << "here in theta ...\n";
  cout << "";
  //struct IntegrationParams int_params = *(struct IntegrationParams *) params;
  struct IntegrationParams int_params = *reinterpret_cast<struct IntegrationParams *>(params);
  int_params.theta = theta;
  double fval = integrate_wrt_phi(int_params) * sin(theta);
  return fval;
}

double compute_bivariate_fval(double phi, void *params)
{
  //struct IntegrationParams int_params = *(struct IntegrationParams *) params;
  struct IntegrationParams int_params = *reinterpret_cast<struct IntegrationParams *>(params);

  Vector x(3,0);
  Vector spherical(3,1);
  spherical[1] = int_params.theta;
  spherical[2] = phi;
  spherical2cartesian(spherical,x);

  double c1 = computeDotProduct(int_params.mu,x);

  Matrix xx = outer_prod(x,x);
  double mj = prod_xMy(int_params.mj,xx,int_params.mj);
  double mi = prod_xMy(int_params.mi,xx,int_params.mi);
  double c2 = mj - mi;

  //double ans = -int_params.log_norm_const + int_params.kappa * c1 + int_params.beta * c2;
  double ans = int_params.kappa * c1 + int_params.beta * c2;
  return exp(ans);
}

double integrate_wrt_phi(struct IntegrationParams &int_params)
{
  //cout << "here in phi ...\n";
  cout << "";
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
  
  double result, error;

  gsl_function F;
  F.function = &compute_bivariate_fval;
  //F.params = &int_params;
  F.params = reinterpret_cast<void *>(&int_params);

  gsl_integration_qags (&F, int_params.phi_lower, int_params.phi_upper, 0, 1e-8, 1000, //2,
                        w, &result, &error); 

  gsl_integration_workspace_free (w);
}

void display_results(string &title, double result, double error)
{
  //double exact = 1.3932039296856768591842462603255;
  double exact = 1;
  cout << title << " ==================\n";
  cout << fixed << setprecision(6);
  cout << "result = " << result << endl;
  cout << "sigma  = " << error << endl; 
  cout << "exact  = " << exact << endl;
  cout << "error  = " << result - exact << " = " <<
          fabs (result - exact) / error << endl;
}

double test_function_integral2 (double *k, size_t dim, void *params)
{
  double A = 1.0 / (M_PI * M_PI * M_PI);
  return A / (1.0 - cos (k[0]) * cos (k[1]) * cos (k[2]));
}

double compute_bivariate_fval(double *k, size_t dim, void *params)
{
  struct IntegrationParams int_params = *reinterpret_cast<struct IntegrationParams *>(params);

  Vector x(3,0);
  Vector spherical(3,1);
  spherical[1] = k[0];
  spherical[2] = k[1];
  spherical2cartesian(spherical,x);

  double c1 = computeDotProduct(int_params.mu,x);

  Matrix xx = outer_prod(x,x);
  double mj = prod_xMy(int_params.mj,xx,int_params.mj);
  double mi = prod_xMy(int_params.mi,xx,int_params.mi);
  double c2 = mj - mi;

  double ans = -int_params.log_norm_const + int_params.kappa * c1 + int_params.beta * c2;
  return exp(ans) * sin(k[0]);
}


double computeTestStatistic( // vMF (null) against Kent (alternative)
  double kappa, 
  double eig_val_sq, 
  double rbar, 
  int N
) {
  double num = N * kappa * kappa * kappa * eig_val_sq;
  double denom = 4 * (kappa - 3* rbar);
  return num/denom;
}

// return the p-value of the test
double compute_pvalue(double t, chi_squared &chisq)
{
  //double critical_value = quantile(chisq,1-alpha);
  double pvalue = 1 - cdf(chisq,t);
  return pvalue;
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
  //return -d*log(12);
}

/*!
 *  \brief This function bins the sample data 
 *  \param res a double
 *  \param unit_coordinates a reference to a vector<vector<double> > 
 */
std::vector<std::vector<int> > updateBins(std::vector<Vector> &unit_coordinates, double res)
{
  std::vector<std::vector<int> > bins;
  int num_rows = 180 / res;
  int num_cols = 360 / res;
  std::vector<int> tmp(num_cols,0);
  for (int i=0; i<num_rows; i++) {
    bins.push_back(tmp);
  }

  cout << "unit_coordinates.size(): " << unit_coordinates.size() << endl;
  double theta,phi;
  int row,col;
  Vector spherical(3,0);
  //string output_file = "constrained_bins.dat";
  //ofstream out(output_file.c_str());
  for (int i=0; i<unit_coordinates.size(); i++) {
    //cout << "i: " << i << endl; 
    cartesian2spherical(unit_coordinates[i],spherical);
    theta = spherical[1] * 180 / PI;
    if (fabs(theta) <= ZERO) {
      row = 0;
    } else {
      row = (int)(ceil(theta/res) - 1);
    }
    phi = spherical[2] * 180 / PI;
    if (fabs(phi) <= ZERO) {
      col = 0;
    } else {
      col = (int)(ceil(phi/res) - 1);
    }
    if (row >= bins.size() || col >= bins[0].size()) {
      cout << "outside bounds: " << row << " " << col << "\n";
      cout << "theta: " << theta << " phi: " << phi << endl;
      cout << "spherical[1]: " << spherical[1] << " spherical[2]: " << spherical[2] << endl;
      cout << "unit_coordinates[i][1]: " << unit_coordinates[i][1] << " unit_coordinates[i][2]: " << unit_coordinates[i][2] << endl;
      fflush(stdout);
    }
    bins[row][col]++;
    /*if (bins[row][col] < 500) {
      bins[row][col]++;
      for (int j=0; j<3; j++) {
        out << fixed << setw(10) << setprecision(3) << unit_coordinates[i][j];
      } // for (j)
      out << endl;
    } // if ()*/
    //cout << "row,col: " << row << "," << col << endl;
  } // for (i)
  //out.close();
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
  Vector cartesian(3,0);
  Vector spherical(3,1);
  for (int i=0; i<bins.size(); i++) {
    phi = 0;
    spherical[1] = theta * PI / 180;
    for (int j=0; j<bins[i].size(); j++) {
      fbins2D << fixed << setw(10) << bins[i][j];
      phi += res;
      spherical[2] = phi * PI / 180;
      spherical2cartesian(spherical,cartesian);
      for (int k=0; k<3; k++) {
        fbins3D << fixed << setw(10) << setprecision(4) << cartesian[k];
      }
      fbins3D << fixed << setw(10) << bins[i][j] << endl;
    }
    theta += res;
    fbins2D << endl;
  }
  fbins2D.close();
  fbins3D.close();
}

/*!
 *  \brief This function is used to read the angular profiles and use this data
 *  to estimate parameters of a Von Mises distribution.
 *  \param parameters a reference to a struct Parameters
 */
void computeEstimators(struct Parameters &parameters)
{
  std::vector<Vector> unit_coordinates;
  bool success = gatherData(parameters,unit_coordinates);
  if (parameters.heat_map == SET) {
    std::vector<std::vector<int> > bins = updateBins(unit_coordinates,parameters.res);
    outputBins(bins,parameters.res);
  }
  if (success && parameters.mixture_model == UNSET) {  // no mixture modelling
    modelOneComponent(parameters,unit_coordinates);
  } else if (success && parameters.mixture_model == SET) { // mixture modelling
    modelMixture(parameters,unit_coordinates);
  }
}

void computeResponsibilityGivenMixture(struct Parameters &parameters)
{
  std::vector<Vector> unit_coordinates;
  bool success = gatherData(parameters,unit_coordinates);
  if (success) {
    Mixture_vMF mixture;
    mixture.load(parameters.mixture_file);
    mixture.computeResponsibilityMatrix(unit_coordinates,parameters.infer_log);
  } else {
    cout << "Something wrong in reading data ...\n";
    exit(1);
  }
}

/*!
 *  \brief This function reads through the profiles from a given directory
 *  and collects the data to do mixture modelling.
 *  \param parameters a reference to a struct Parameters
 *  \param unit_coordinates a reference to a std::vector<Vector>
 */
bool gatherData(struct Parameters &parameters, std::vector<Vector> &unit_coordinates)
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
        std::vector<std::vector<Vector> > _unit_coordinates(NUM_THREADS);
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
            std::vector<Vector> coords = structure.getUnitCoordinates();
            for (int j=0; j<coords.size(); j++) {
              _unit_coordinates[tid].push_back(coords[j]);
            }
          }
        }
        for (int i=0; i<NUM_THREADS; i++) {
          for (int j=0; j<_unit_coordinates[i].size(); j++) {
            unit_coordinates.push_back(_unit_coordinates[i][j]);
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
      unit_coordinates = structure.getUnitCoordinates();
      return 1;
    } else {
      cout << "Profile " << parameters.profile_file << " does not exist ...\n";
      return 0;
    }
  }
}

/*!
 *  \brief This function is used to simulate the mixture model.
 *  \param parameters a reference to a struct Parameters
 */
void simulateMixtureModel(struct Parameters &parameters)
{
  std::vector<Vector> data;
  if (DISTRIBUTION == KENT) {
    if (parameters.load_mixture == SET) {
      Mixture original;
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
          string mix_density_file = "./visualize/sampled_data/bins_kent/mixture_density.dat";
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
      std::vector<Kent> components = generateRandomComponents(k);
      Mixture original(k,components,weights);
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
  } else if (DISTRIBUTION == VMF) {
    if (parameters.load_mixture == SET) {
      Mixture_vMF original;
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
          string mix_density_file = "./visualize/sampled_data/bins_vmf/mixture_density.dat";
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
      std::vector<vMF> components = generateRandomComponents_vMF(k);
      Mixture_vMF original(k,components,weights);
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

/*!
 *  \brief This function is used to generate random weights such that
 *  0 < w_i < 1 and \sum_i w_i = 1
 *  \param K an integer
 *  \return the list of weights
 */
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

/*!
 *  \brief This function is used to generate random components.
 *  \param num_components an integer
 *  \param D an integer
 *  \return the list of components
 */
std::vector<Kent> generateRandomComponents(int num_components)
{
  Vector mean(3,0),major_axis(3,0),minor_axis(3,0);

  // generate random kappas
  Vector kappas = generateRandomKappas(num_components);

  // generate random betas 
  Vector betas = generateRandomBetas(kappas);

  std::vector<Kent> components;
  double psi,alpha,eta;
  for (int i=0; i<num_components; i++) {
    //generateRandomOrthogonalVectors(mean,major_axis,minor_axis);
    //Kent kent(mean,major_axis,minor_axis,kappas[i],betas[i]);

    generateRandomOrientations(psi,alpha,eta);
    Kent kent(psi,alpha,eta,kappas[i],betas[i]);
    components.push_back(kent);
  }
  return components;
}

std::vector<vMF> generateRandomComponents_vMF(int num_components)
{
  // generate random kappas
  Vector kappas = generateRandomKappas(num_components);

  Vector mean(3,0),spherical(3,1);
  double theta,phi;

  std::vector<vMF> components;
  for (int i=0; i<num_components; i++) {
    spherical[1] = PI * uniform_random();
    spherical[2] = (2 * PI) * uniform_random();
    spherical2cartesian(spherical,mean); 
    vMF vmf(mean,kappas[i]);
    components.push_back(vmf);
  }
  return components;
}

/*!
 *  \brief This function generates random kappas 
 *  \param num_components an integer
 *  \return the list of random kappas 
 */
Vector generateRandomKappas(int K)
{
  Vector random_kappas;
  double max_kappa;

  //max_kappa = MAX_KAPPA;
  max_kappa = 100;
  for (int i=0; i<K; i++) {
    double kappa = uniform_random() * max_kappa;
    random_kappas.push_back(kappa);
  }
  return random_kappas;
}

Vector generateRandomBetas(Vector &kappas)
{
  Vector random_betas;
  for (int i=0; i<kappas.size(); i++) {
    double beta = uniform_random() * 0.5 * kappas[i];
    random_betas.push_back(beta);
  }
  return random_betas;
}

/*!
 *  \brief This function models a single component.
 *  \param parameters a reference to a struct Parameters
 *  \param data a reference to a std::vector<Vector>
 */
void modelOneComponent(struct Parameters &parameters, std::vector<Vector> &data)
{
  cout << "Sample size: " << data.size() << endl;
  Vector weights(data.size(),1);
  if (DISTRIBUTION == KENT) {
    //Kent kent;
    //Kent kent(ZAXIS,XAXIS,YAXIS,100,45);
    double psi = 90; psi *= PI/180;
    double alpha = 90; alpha *= PI/180;
    double eta = 90; eta *= PI/180;
    double kappa = 10;
    double ecc = 0.5;
    double beta = 0.5 * kappa * ecc;
    Kent kent(psi,alpha,eta,kappa,beta);
    std::vector<struct Estimates> all_estimates;
    //kent.computeAllEstimators(data,all_estimates,1,1);
    kent.computeAllEstimators(data,all_estimates,1,0);
  } else if (DISTRIBUTION == VMF) {
    vMF vmf;
    vmf.estimateParameters(data,weights);
  }
}

/*!
 *  \brief This function models a mixture of several components.
 *  \param parameters a reference to a struct Parameters
 *  \param data a reference to a std::vector<std::vector<double,3> >
 */
void modelMixture(struct Parameters &parameters, std::vector<Vector> &data)
{
  Vector data_weights(data.size(),1);
  // if the optimal number of components need to be determined
  if (DISTRIBUTION == KENT) {
    if (parameters.infer_num_components == SET) {
      Mixture mixture;
      if (parameters.continue_inference == UNSET) {
        Mixture m(parameters.start_from,data,data_weights);
        mixture = m;
        mixture.estimateParameters();
      } else if (parameters.continue_inference == SET) {
        mixture.load(parameters.mixture_file,data,data_weights);
      } // continue_inference
      ofstream log(parameters.infer_log.c_str());
      Mixture stable = inferComponents(mixture,data.size(),log);
      log.close();
    } else if (parameters.infer_num_components == UNSET) {
      // for a given value of number of components
      // do the mixture modelling
      Mixture mixture(parameters.fit_num_components,data,data_weights);
      mixture.estimateParameters();
    }
  } else if (DISTRIBUTION == VMF) {
    if (parameters.infer_num_components == SET) {
      Mixture_vMF mixture;
      if (parameters.continue_inference == UNSET) {
        Mixture_vMF m(parameters.start_from,data,data_weights);
        mixture = m;
        mixture.estimateParameters();
      } else if (parameters.continue_inference == SET) {
        mixture.load(parameters.mixture_file,data,data_weights);
      } // continue_inference
      ofstream log(parameters.infer_log.c_str());
      Mixture_vMF stable = inferComponents_vMF(mixture,data.size(),log);
      /*double percent_remove = 0.2;
      Mixture_vMF jk2 = jackknife(data,stable,percent_remove,log);
      Mixture_vMF inferred2 = inferComponents_vMF(jk2,data.size(),log);*/
      log.close();
    } else if (parameters.infer_num_components == UNSET) {
      // for a given value of number of components
      // do the mixture modelling
      Mixture_vMF mixture(parameters.fit_num_components,data,data_weights);
      mixture.estimateParameters();
    }
  }
}

void update_tracking_file(int iter, Mixture &parent, ostream &out)
{
  out << fixed << setw(5) << iter << "\t";
  out << fixed << setw(5) << parent.getNumberOfComponents() << "\t";

  /*if (CRITERION == AIC) {
    out << fixed << scientific << setprecision(6) << parent.getAIC() << endl;
  } else if (CRITERION == BIC) {
    out << fixed << scientific << setprecision(6) << parent.getBIC() << endl;
  } else if (CRITERION == ICL) {
    out << fixed << scientific << setprecision(6) << parent.getICL() << endl;
  } else if (CRITERION == MMLC) {
    out << fixed << scientific << setprecision(6) << parent.first_part()
        << "\t\t" << parent.second_part()
        << "\t\t" << parent.getMinimumMessageLength() << endl;
  }*/

  double aic = parent.computeAIC();
  double bic = parent.computeBIC();
  double icl = parent.computeICL();
  double msglen = parent.computeMinimumMessageLength();

  out << fixed << scientific << setprecision(6) << aic << "\t";
  out << fixed << scientific << setprecision(6) << bic << "\t";
  out << fixed << scientific << setprecision(6) << icl << "\t";
  out << fixed << scientific << setprecision(6) << parent.first_part() << "\t";
  out << fixed << scientific << setprecision(6) << parent.second_part() << "\t";
  out << fixed << scientific << setprecision(6) << msglen << endl;
}

Mixture inferComponents(std::vector<Vector> &data, string &log_file)
{
  Vector data_weights(data.size(),1.0);
  Mixture m(1,data,data_weights);
  Mixture mixture = m;
  mixture.estimateParameters();
  ofstream log(log_file.c_str());
  Mixture inferred = inferComponents(mixture,data.size(),log);
  log.close();
  return inferred;
}

Mixture inferComponents(Mixture &mixture, int N, ostream &log)
{
  int K,iter = 0;
  std::vector<Kent> components;
  Mixture modified,improved,parent;
  Vector sample_size;

  double null_msglen = mixture.computeNullModelMessageLength();
  log << "Null model encoding: " << null_msglen << " bits."
      << "\t(" << null_msglen/N << " bits/point)\n\n";

  ofstream out(tracking_file.c_str());

  improved = mixture;

  while (1) {
    parent = improved;
    iter++;
    update_tracking_file(iter,parent,out);
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
        modified = parent.split(i,log);
        if (IGNORE_SPLIT == 0) {
          updateInference(modified,improved,N,log,SPLIT);
        } else {
          log << "\t\tIGNORING SPLIT ... \n\n";
        }
      }
    }
    if (K >= 2) {  // kill() ...
      for (int i=0; i<K; i++) {
        modified = parent.kill(i,log);
        updateInference(modified,improved,N,log,KILL);
      } // killing each component
    } // if (K > 2) loop
    if (K > 1) {  // join() ...
      for (int i=0; i<K; i++) {
        int j = parent.getNearestComponent(i); // closest component
        modified = parent.join(i,j,log);
        updateInference(modified,improved,N,log,JOIN);
      } // join() ing nearest components
    } // if (K > 1) loop
    if (improved == parent) goto finish;
  } // if (improved == parent || iter%2 == 0) loop

  finish:
  out.close();
  string inferred_mixture_file = "./simulation/inferred_mixture";
  parent.printParameters(inferred_mixture_file);
  return parent;
}

/*!
 *  \brief Updates the inference
 *  \param modified a reference to a Mixture
 *  \param current a reference to a Mixture
 *  \param log a reference to a ostream
 *  \param operation an integer
 */
void updateInference(Mixture &modified, Mixture &current, int N, ostream &log, int operation)
{
  double modified_value,current_value,improvement_rate;

  switch(CRITERION) {
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
  }

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

  /*if (operation == KILL || operation == JOIN || operation == SPLIT) {
    if (current_value > modified_value) {
      improvement_rate = (current_value - modified_value) / fabs(current_value);
      log << "\t ... IMPROVEMENT ... (+ " << fixed << setprecision(3) 
          << 100 * improvement_rate << " %) ";
      log << "\t\t[ACCEPT]\n\n";
      current = modified;
    } else {
      log << "\t ... NO IMPROVEMENT\t\t\t[REJECT]\n\n";
    }
  }*/
}

Mixture_vMF inferComponents_vMF(std::vector<Vector> &data, string &log_file)
{
  Vector data_weights(data.size(),1.0);
  Mixture_vMF m(1,data,data_weights);
  Mixture_vMF mixture = m;
  mixture.estimateParameters();
  ofstream log(log_file.c_str());
  Mixture_vMF inferred = inferComponents_vMF(mixture,data.size(),log);
  return inferred;

  /* jackknife */
  /*double percent_remove = 0.2;
  Mixture_vMF jk1 = jackknife(data,inferred,percent_remove,log);
  Mixture_vMF inferred1 = inferComponents_vMF(jk1,data.size(),log);
  Mixture_vMF jk2 = jackknife(data,inferred1,percent_remove,log);
  Mixture_vMF inferred2 = inferComponents_vMF(jk2,data.size(),log);
  log.close();
  return inferred2;*/
}

Mixture_vMF jackknife(
  std::vector<Vector> &data, Mixture_vMF &original, double percent_remove, ostream &log
) {
  log << " jackknife ...\n\n";
  int N = data.size();
  int K = original.getNumberOfComponents();

  if (K == 1) return original;

  int K_remove = ceil(percent_remove * K);
  log << "K: " << K << endl;
  log << "K_remove: " << K_remove << endl;
  int K_m = K - K_remove;
  log << "K_m: " << K_m << endl;

  std::vector<int> flags(K,0);  // = 1 => discard the component
  int count = 0;
  while (1) {
    int c = rand() % K;
    if (flags[c] == 0) {
      flags[c] = 1;
      count++;
    }
    if (count == K_remove) break;
  }

  // adjust weights
  Vector weights = original.getWeights();
  double residual_sum = 0;
  for (int i=0; i<K; i++) {
    if (flags[i] != 1) residual_sum += weights[i];
  }
  Vector weights_m(K_m,0);
  int index = 0;
  for (int i=0; i<K; i++) {
    if (flags[i] != 1) weights_m[index++] = weights[i] / residual_sum;
  }
  
  // adjust responsibility matrix
  std::vector<Vector> responsibility = original.getResponsibilityMatrix();
  Vector residual_sums(N,0);
  for (int i=0; i<N; i++) {
    for (int j=0; j<K; j++) {
      if (flags[j] != 1) {
        residual_sums[i] += responsibility[j][i];
      } // if()
    } // j
  } // i
  Vector resp(N,0);
  std::vector<Vector> responsibility_m(K_m,resp);
  index = 0;
  for (int i=0; i<K; i++) {
    if (flags[i] != 1) {
      #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) //private(residual_sum) 
      for (int j=0; j<N; j++) {
        if (residual_sums[j] <= 0.01) {
          responsibility_m[index][j] = 1.0 / K_m;
        } else {
          responsibility_m[index][j] = responsibility[i][j] / residual_sums[j];
        } // if()
        assert(responsibility_m[index][j] >= 0 && responsibility_m[index][j] <= 1);
      } // j
      index++;
    } // if ()
  } // i

  // adjust effective sample size
  Vector sample_size_m(K_m,0);
  for (int i=0; i<K_m; i++) {
    double sum = 0;
    #pragma omp parallel for if(ENABLE_DATA_PARALLELISM) num_threads(NUM_THREADS) reduction(+:sum) 
    for (int j=0; j<N; j++) {
      sum += responsibility_m[i][j];
    } // j
    sample_size_m[i] = sum;
  } // i

  // child components
  std::vector<vMF> components_m(K_m);
  std::vector<vMF> components = original.getComponents();
  index = 0;
  for (int i=0; i<K; i++) {
    if (flags[i] != 1) {
      components_m[index++] = components[i];
    } // if
  } // i

  log << "\t\tResidual:\n";
  Vector data_weights_m(N,1);
  Mixture_vMF modified(K_m,components_m,weights_m,sample_size_m,responsibility_m,data,data_weights_m);
  log << "\t\tBefore adjustment ...\n";
  modified.computeMinimumMessageLength();
  modified.printParameters(log,2);
  modified.EM();
  log << "\t\tAfter adjustment ...\n";
  modified.printParameters(log,2);
  modified.printIndividualMsgLengths(log);
  return modified;
}

Mixture_vMF inferComponents_vMF(Mixture_vMF &mixture, int N, ostream &log)
{
  int K,iter = 0;
  std::vector<vMF> components;
  Mixture_vMF modified,improved,parent;
  Vector sample_size;

  double null_msglen = mixture.computeNullModelMessageLength();
  log << "Null model encoding: " << null_msglen << " bits."
      << "\t(" << null_msglen/N << " bits/point)\n\n";

  //MIN_N = 1;

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
        modified = parent.split(i,log);
        if (IGNORE_SPLIT == 0) {
          updateInference_vMF(modified,improved,N,log,SPLIT);
        } else {
            log << "\t\tIGNORING SPLIT ... \n\n";
        }
      }
    }
    if (K >= 2) {  // kill() ...
      for (int i=0; i<K; i++) {
        modified = parent.kill(i,log);
        updateInference_vMF(modified,improved,N,log,KILL);
      } // killing each component
    } // if (K > 2) loop
    if (K > 1) {  // join() ...
      for (int i=0; i<K; i++) {
        int j = parent.getNearestComponent(i); // closest component
        modified = parent.join(i,j,log);
        updateInference_vMF(modified,improved,N,log,JOIN);
      } // join() ing nearest components
    } // if (K > 1) loop
    if (improved == parent) goto finish;
  } // if (improved == parent || iter%2 == 0) loop

  finish:
  //string inferred_mixture_file = "./simulation/inferred_mixture_vmf";
  string inferred_mixture_file = "inferred_mixture_vmf";
  parent.printParameters(inferred_mixture_file);
  return parent;
}

/*!
 *  \brief Updates the inference
 *  \param modified a reference to a Mixture_vMF
 *  \param current a reference to a Mixture_vMF
 *  \param log a reference to a ostream
 *  \param operation an integer
 */
void updateInference_vMF(Mixture_vMF &modified, Mixture_vMF &current, int N, ostream &log, int operation)
{
  double modified_value,current_value,improvement_rate;

  switch(CRITERION) {
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
  }

  if (current_value > modified_value) {
    improvement_rate = (current_value - modified_value) / fabs(current_value);
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
}

////////////////////// TESTING FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

void TestFunctions(void)
{
  Test test;

  //test.load_data();

  //test.bessel();

  //test.testing_cartesian2sphericalPoleXAxis();

  //test.parallel_sum_computation();

  //test.uniform_number_generation();

  //test.arbitrary_rotation();

  //test.matrixFunctions();

  //test.productMatrixVector();

  //test.dispersionMatrix();

  //test.numericalIntegration();

  //test.normalDistributionFunctions();

  //test.orthogonalTransformations();

  //test.orthogonalTransformations2();

  //test.randomSampleGeneration();

  //test.multivariate_normal();

  //test.acg();

  //test.bingham();

  //test.kent_bingham_generation();

  //test.normalization_constant();

  //test.optimization();

  //test.moment_estimation();

  //test.ml_estimation();

  //test.expectation();

  //test.kl_divergence();

  //test.fisher();

  //test.fisher2();

  //test.mml_estimation();

  //test.mml_estimation2();

  //test.plot_posterior_density();

  //test.vmf_all_estimation();

  //test.chi_square();

  //test.hypothesis_testing();

  //test.hypothesis_testing2();

  //test.confidence_region();

  //test.infer_mixture();

  //test.infer_mixture_vmf();

  test.contours();

  //test.gsl_numerical_integration();

  //test.gsl_numerical_kent_density_integration();

  //test.gsl_monte_carlo_integration();

  //test.gsl_monte_carlo_kent_density_integration();
}

////////////////////// EXPERIMENTS \\\\\\\\\\\\\\\\\\\\\\\\\\\\

void RunExperiments(int iterations)
{
  Experiments experiments(iterations);

  //experiments.fisher_uncertainty();

  //experiments.simulate();

  //experiments.exp1();

  //experiments.exp2();

  //experiments.exp3();

  //experiments.exp4();

  //experiments.exp5();

  //experiments.exp6();
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

void chisquare_hypothesis_testing(
  std::vector<Vector> &data,
  std::vector<struct Estimates> &all_kent_estimates,
  Vector &statistics,
  Vector &pvalues
) {
  int num_methods = all_kent_estimates.size();
  statistics = Vector(num_methods,0);
  pvalues = Vector(num_methods,0);

  chi_squared chisq(5);

  Kent kent_ml(
    all_kent_estimates[MLE].psi,
    all_kent_estimates[MLE].alpha,
    all_kent_estimates[MLE].eta,
    all_kent_estimates[MLE].kappa,
    all_kent_estimates[MLE].beta
  );
  double kent_ml_negloglike = kent_ml.computeNegativeLogLikelihood(data);
  cout << "Kent (MLE) negloglike: " << kent_ml_negloglike << endl;

  double kent_negloglike,log_ratio_statistic,pvalue;
  for (int i=0; i<num_methods; i++) {
    if (i != MLE) {
      // null: Kent(i)
      Kent kent(
        all_kent_estimates[i].psi,
        all_kent_estimates[i].alpha,
        all_kent_estimates[i].eta,
        all_kent_estimates[i].kappa,
        all_kent_estimates[i].beta
      );
      kent_negloglike = kent.computeNegativeLogLikelihood(data);
      log_ratio_statistic = 2 * (kent_negloglike - kent_ml_negloglike);
      if (log_ratio_statistic < 0) {
        //assert(fabs(log_ratio_statistic) < 1e-4);
        log_ratio_statistic = fabs(log_ratio_statistic);
      }
      pvalue = 1 - cdf(chisq,log_ratio_statistic);
    } else { 
      log_ratio_statistic = 0;
      pvalue = 1;
    } // if()
    cout << "(null: Kent[i]) negloglike: " << kent_negloglike
         << "; log_ratio_statistic: " << log_ratio_statistic
         << "; pvalue: " << pvalue << endl;
    statistics[i] = log_ratio_statistic;
    pvalues[i] = pvalue;
  } // for()
}

