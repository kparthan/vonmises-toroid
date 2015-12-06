#ifndef MIXTURE_SINE_H
#define MIXTURE_SINE_H

#include "BVM_Sine.h"

class Mixture_Sine
{
  private:
    //! ID
    int id;

    //! Sample size
    int N;

    //! Number of components
    int K;

    //! List of components
    std::vector<BVM_Sine> components;

    std::vector<Vector> mix_init_means;
    
    //! Sample (x_i) -- Cartesian coordinates
    //! (on a 3D torus)
    std::vector<Vector> data;

    //! Data weights
    Vector data_weights;

    //! Responsibility matrix (K X N)
    std::vector<Vector> responsibility;

    //! Effective sample size for each component (n_k)
    Vector sample_size;

    //! Weights of the components (a_k)
    Vector weights;

    //! List of message lengths over several iterations
    Vector msglens;

    //! Null model message length
    double null_msglen;

    //! Optimal encoding length
    double Ik,Iw,sum_It,Il,kd_term,part1,part2,minimum_msglen;
    Vector It;

    double negloglike;

    double aic,bic,icl;

  public:
    //! Null constructor
    Mixture_Sine();

    //! Constructor
    Mixture_Sine(int, std::vector<BVM_Sine> &, Vector &);

    //! Constructor
    Mixture_Sine(int, std::vector<Vector> &, Vector &);

    //! Constructor
    Mixture_Sine(int, std::vector<BVM_Sine> &, Vector &, Vector &, 
            std::vector<Vector> &, std::vector<Vector> &, Vector &);

    //! Overloading = operator
    Mixture_Sine operator=(const Mixture_Sine &);

    //! Overloading == operator
    bool operator==(const Mixture_Sine &);

    //! Prepare log file
    string getLogFile();

    //! Gets the list of weights
    Vector getWeights();

    //! Gets the list of components
    std::vector<BVM_Sine> getComponents();

    //! Returns number of components
    int getNumberOfComponents();

    //! Gets the responsibility matrix
    std::vector<Vector> getResponsibilityMatrix();

    //! Gets the sample size
    Vector getSampleSize();

    void setInitialMeans(std::vector<Vector> &);

    //! Initialize parameters
    void initialize();
    void initialize_children_1();
    void initialize_children_2();

    //! Updates the effective sample size
    void updateEffectiveSampleSize();

    //! Update the component weights
    void updateWeights();

    void updateWeights_ML();

    //! Update components
    void updateComponents();

    //! Update the responsibility matrix
    void updateResponsibilityMatrix();

    //! Computes the responsibility matrix
    void computeResponsibilityMatrix(std::vector<Vector> &, string &);
                                          
    //! Probability of a datum
    double log_probability(Vector &);

    //! Computes the negative log likelihood
    double computeNegativeLogLikelihood(std::vector<Vector> &);

    double computeNegativeLogLikelihood(int verbose = 0);

    double compress(std::vector<Vector> &);

    //! Computes the minimum message length
    double computeMinimumMessageLength(int verbose = 0);

    void printIndividualMsgLengths(ostream &);

    //! Gets the minimum message length
    double getMinimumMessageLength();

    //! Gets the first part
    double first_part();

    //! Gets the second part
    double second_part();

    double getNegativeLogLikelihood();

    double getAIC();

    double getBIC();

    double getICL();

    //! Estimate mixture parameters
    double estimateParameters();

    //! EM loop
    void EM();

    void EM(
      ostream &,
      void (Mixture_Sine::*update_weights)(),
      double (Mixture_Sine::*objective_function)(int)
    );

    //! Computes the null model message length
    double computeNullModelMessageLength();

    //! Prints the model parameters
    void printParameters(ostream &, int, double);

    //! Prints the model parameters
    void printParameters(ostream &, int);

    void printParameters(string &);

    //! Prints the model parameters
    void printParameters(ostream &);

    //! Loads the mixture file
    void load(string &);

    //! Loads the mixture file with the corresponding data
    void load(string &, std::vector<Vector> &, Vector &);

    //! Randomly choose a component
    int randomComponent();

    //! Saves the data generated from a component
    void saveComponentData(int, std::vector<Vector> &);

    //! Generate random data from the distribution using mixture proportions
    std::vector<Vector> generate(int, bool);

    //! Splits a component
    Mixture_Sine split(int, ostream &);

    //! Deltes a component
    Mixture_Sine kill(int, ostream &);

    //! Joins two  components
    Mixture_Sine join(int, int, ostream &);

    //! Generate heat maps (for d=3)
    void generateHeatmapData(double);

    //! Get the nearest component
    int getNearestComponent(int);

    //! Computes the approx KL divergence between two mixtures
    double computeKLDivergence(Mixture_Sine &);

    double computeKLDivergence(Mixture_Sine &, std::vector<Vector> &);

    double computeAIC();

    double computeBIC();

    std::vector<std::vector<int> > compute_cluster_indicators();

    double computeICL();
};

#endif

