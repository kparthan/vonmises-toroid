#ifndef SUPPORT_H
#define SUPPORT_H

#include "Header.h"
#include "Mixture_Sine.h"
#include "vMC.h"

struct Parameters
{
  int test;                 // flag to test some modules
  int experiments;          // flag to run some experiments 
  int iterations;           // number of iterations
  double k1,k2,rho;
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

struct TerminationCondition  {
  bool operator() (double min, double max)  {
    return fabs(min - max) <= TOLERANCE;
  }
};

struct EstimatesSine
{
  double mu1,mu2;
  double kappa1,kappa2,lambda,rho;
  double msglen,negloglike,kldiv;
};

struct SufficientStatisticsSine {
  double N;
  double cost1,cost2,sint1,sint2;
  double sint1sint2,sint1cost2,cost1sint2,cost1cost2;
};

// general functions
struct Parameters parseCommandLineInput (int, char **); 
void Usage (const char *, options_description &);
bool checkFile(string &);
void writeToFile(const char *, std::vector<Vector> &, int);
void writeToFile(const char *, std::vector<Vector> &);
void writeToFile(string &file_name, std::vector<Vector> &);
void writeToFile(string &file_name, std::vector<std::vector<int> > &);
string extractName(string &);
void print(ostream &, const Vector &, int);
void print(string &, struct EstimatesSine &);
void print(struct SufficientStatisticsSine &);
void check_and_create_directory(string &);

double scale_to_aom(double &);
std::vector<Vector> scale_to_aom(std::vector<Vector> &);
int sign(double);
double normalize(Vector &, Vector &);
double norm(Vector &);
void cartesian2spherical(Vector &, Vector &);
void spherical2cartesian(Vector &, Vector &);
void toroid2cartesian(Vector &, Vector &);
double computeDotProduct(Vector &, Vector &);
Vector crossProduct(Vector &, Vector &); 
double computeSum(Vector &);
double computeLogSurfaceAreaSphere(int);
double uniform_random();
double computeLogModifiedBesselFirstKind(double, double);
double computeLogRatioBessel(double, double);
double computeRatioBessel(double);

std::vector<Vector> load_data_table(string &, int dim = 2);
Matrix outer_prod(Vector &, Vector &);
Vector prod(Matrix &, Vector &);
Vector prod(Vector &, Matrix &);
double prod_vMv(Vector &, Matrix &);
double prod_xMy(Vector &, Matrix &, Vector &);
double determinant_3d(Matrix &);
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
std::vector<Vector> transform(std::vector<Vector> &, Matrix &);
bool invertMatrix(const Matrix &, Matrix &);
void eigenDecomposition(Matrix, Vector &, Matrix &);
void jacobiRotateMatrix(Matrix &, Matrix &, int, int);
double compute_aic(int, int, double);
double compute_bic(int, int, double);

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
  std::vector<struct EstimatesSine> &, Vector &, Vector &
);

double computeConstantTerm(int);
double logLatticeConstant(int);
std::vector<std::vector<int> > updateBins(std::vector<Vector> &, double);
void outputBins(std::vector<std::vector<int> > &, double);
std::vector<Vector> sample_empirical_distribution(
  int, double, std::vector<std::vector<int> > &
);

void computeEstimators(struct Parameters &);
bool gatherData(struct Parameters &, std::vector<Vector> &);
void modelOneComponent(struct Parameters &, std::vector<Vector> &);
void modelMixture(struct Parameters &, std::vector<Vector> &);
void simulateMixtureModel(struct Parameters &);
Vector generateFromSimplex(int);
std::vector<BVM_Sine> generateRandomComponents(int);
Mixture_Sine inferComponents(std::vector<Vector> &, string &);
Mixture_Sine inferComponents(Mixture_Sine &, int, ostream &);
void updateInference(Mixture_Sine &, Mixture_Sine &, int, ostream &, int);

void TestFunctions(void);
void RunExperiments(struct Parameters &);

double accept_reject_fval_unimodal_marginal_sine(
  double &, double &, double &, double &, double &, double &
);
double accept_reject_fval_bimodal_marginal_sine(
  double &, double &, double &, double &, double &, double &, double &, double &
);
double accept_reject_fval_unimodal_marginal_cosine(
  double &, double &, double &, double &, double &, double &
);
double accept_reject_fval_bimodal_marginal_cosine(
  double &, double &, double &, double &, double &, double &, double &, double &
);

double banerjee_approx(double &);
void computeSufficientStatisticsSineNotParallel(
  std::vector<Vector> &, struct SufficientStatisticsSine &
);
void computeSufficientStatisticsSine(
  std::vector<Vector> &, struct SufficientStatisticsSine &
);
void computeSufficientStatisticsSineNotParallel(
  std::vector<Vector> &, struct SufficientStatisticsSine &, Vector &
);
void computeSufficientStatisticsSine(
  std::vector<Vector> &, struct SufficientStatisticsSine &, Vector &
);
double ConstraintSine(const Vector &, std::vector<double> &, void *);
vMC getConditionalDensitySine(
  double &, double &, double &, double &, double &
);

#endif

