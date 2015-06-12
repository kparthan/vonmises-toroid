#ifndef SUPPORT_H
#define SUPPORT_H

struct Parameters
{
  int test;                 // flag to test some modules
  int experiments;          // flag to run some experiments 
  int iterations;           // number of iterations
  string profile_file;      // path to a single profile
  string profiles_dir;      // path to the directory containing the profiles
  int heat_map;             // flag to generate heat map images
  double res;          // resolution used in heat map images
  int read_profiles;        // flag to read profile(s)
  int mixture_model;        // flag to model a mixture
  int fit_num_components;   // # of components in the mixture model
  int infer_num_components; // flag to infer # of components
  int max_components;       // max components to infer
  string infer_log;         // log file
  int continue_inference;   // flag to continue inference from some state
  int simulation;           // flag to run mixture model simulation
  int load_mixture;         // flag to read mixture from a file
  int simulated_components; // # of components to be simulated
  string mixture_file;      // file containing the mixture information
  int sample_size;          // sample size to be generated from the simulated mixture
  int num_threads;          // flag to enable multithreading
  double max_kappa;    // max value of kappa allowed
  int start_from;           // starting value of number of components
                            // during inference
  int estimate_all;         // estimate using all methods
  int compute_responsibility_matrix;  // flag
};

// general functions
//void Setup();
struct Parameters parseCommandLineInput (int, char **); 
void Usage (const char *, options_description &);
bool checkFile(string &);
void writeToFile(const char *, std::vector<Vector> &, int);
void writeToFile(const char *, std::vector<Vector> &);
void writeToFile(string &file_name, std::vector<Vector> &);
void writeToFile(string &file_name, std::vector<std::vector<int> > &);
string extractName(string &);
void print(ostream &, Vector &, int);
void print(string &, struct Estimates &);
void print(string &, struct Estimates_vMF &);
void check_and_create_directory(string &);

double scale_to_aom(double &);
std::vector<Vector> scale_to_aom(std::vector<Vector> &);
int sign(double);
double normalize(Vector &, Vector &);
double norm(Vector &);
//void cartesian2spherical2(Vector &, Vector &);
void cartesian2spherical(Vector &, Vector &);
//void spherical2cartesian2(Vector &, Vector &);
void spherical2cartesian(Vector &, Vector &);
void computeLambertProjection(Vector &, Vector &);
void computeLambertProjection(std::vector<Vector> &);
double computeDotProduct(Vector &, Vector &);
Vector crossProduct(Vector &, Vector &); 
double computeSum(Vector &);
double computeLogSurfaceAreaSphere(int);
void solveQuadratic(Vector &, double, double, double);
double uniform_random();
double computeLogModifiedBesselFirstKind(double, double);

std::vector<Vector> load_data_table(string &);
Matrix outer_prod(Vector &, Vector &);
Vector prod(Matrix &, Vector &);
Vector prod(Vector &, Matrix &);
double prod_vMv(Vector &, Matrix &);
double prod_xMy(Vector &, Matrix &, Vector &);
double determinant(Matrix &);
Vector computeVectorSum(std::vector<Vector> &);
Vector computeVectorSum(std::vector<Vector> &, Vector &, double &);
Vector computeNormalizedVectorSum(std::vector<Vector> &);
Matrix computeDispersionMatrix(std::vector<Vector> &);
Matrix computeDispersionMatrix(std::vector<Vector> &, Vector &);
Matrix computeNormalizedDispersionMatrix(std::vector<Vector> &);
Matrix rotate_about_arbitrary_axis(Vector &, double);
Matrix rotate_about_xaxis(double);
Matrix rotate_about_yaxis(double);
Matrix rotate_about_zaxis(double);
Matrix computeOrthogonalTransformation(Vector &, Vector &);
Matrix computeOrthogonalTransformation(double, double, double);
void computeOrthogonalTransformation(Vector &, Vector &, double &, double &, double &);
Matrix align_xaxis_with_vector(Vector &);
Matrix align_vector_with_xaxis(Vector &);
Matrix align_vector_with_xaxis(double, double);
void generateRandomOrthogonalVectors(Vector &, Vector &, Vector &);
void generateRandomOrientations(double &, double &, double &);
Matrix generateRandomCovarianceMatrix(int);
std::vector<Vector> transform(std::vector<Vector> &, Matrix &);
bool invertMatrix(const Matrix &, Matrix &);
void eigenDecomposition(Matrix, Vector &, Matrix &);
void jacobiRotateMatrix(Matrix &, Matrix &, int, int);
double computeDawsonsIntegral(double);
void track(const std::vector<double> &, const double);
void rhs(const std::vector<double> &, std::vector<double> &, const double);
double Constraint2(const std::vector<double> &, std::vector<double> &, void *);
double Constraint5(const std::vector<double> &, std::vector<double> &, void *);
double Constraint5_2(const std::vector<double> &, std::vector<double> &, void *);
double test_function_integral(double, void *); 
double integrate_wrt_theta(double, void *);
double compute_bivariate_fval(double, void *);
double integrate_wrt_phi(struct IntegrationParams &);
double test_function_integral2(double *, size_t, void *);
double compute_bivariate_fval(double *, size_t, void *);
void display_results(string &, double, double);
double computeTestStatistic(double, double, double, int);
double compute_pvalue(double, chi_squared &);
double compute_aic(int, int, double);
double compute_bic(int, int, double);

double computeConstantTerm(int);
double logLatticeConstant(int);
std::vector<std::vector<int> > updateBins(std::vector<Vector> &, double);
void outputBins(std::vector<std::vector<int> > &, double);
void computeEstimators(struct Parameters &);
void computeResponsibilityGivenMixture(struct Parameters &);
bool gatherData(struct Parameters &, std::vector<Vector> &);
void modelOneComponent(struct Parameters &, std::vector<Vector> &);
void modelMixture(struct Parameters &, std::vector<Vector> &);
void simulateMixtureModel(struct Parameters &);
Vector generateFromSimplex(int);
std::vector<Kent> generateRandomComponents(int);
std::vector<vMF> generateRandomComponents_vMF(int);
Vector generateRandomKappas(int);
Vector generateRandomBetas(Vector &);
void update_tracking_file(int, Mixture &, ostream &);
Mixture inferComponents(std::vector<Vector> &, string &);
Mixture inferComponents(Mixture &, int, ostream &);
void updateInference(Mixture &, Mixture &, int, ostream &, int);
Mixture_vMF inferComponents_vMF(std::vector<Vector> &, string &);
Mixture_vMF jackknife(std::vector<Vector> &, Mixture_vMF &, double, ostream &);
Mixture_vMF inferComponents_vMF(Mixture_vMF &, int, ostream &);
void updateInference_vMF(Mixture_vMF &, Mixture_vMF &, int, ostream &, int);

void TestFunctions(void);
void RunExperiments(int);

Vector sort(Vector &);
Vector sort(Vector &, std::vector<int> &);
void quicksort(Vector &, std::vector<int> &, int, int);
int partition(Vector &, std::vector<int> &, int, int);
std::vector<Vector> flip(std::vector<Vector> &);
double computeMedian(Vector &);
Vector computeMedians(std::vector<Vector> &);
double computeMean(Vector &);
Vector computeMeans(std::vector<Vector> &);
double computeVariance(Vector &);
int minimumIndex(Vector &);
int maximumIndex(Vector &);
void chisquare_hypothesis_testing(
  std::vector<Vector> &, std::vector<struct Estimates> &, Vector &, Vector &
);

#endif

