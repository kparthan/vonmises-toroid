#include "Support.h"
#include "Experiments.h"
#include "BVM_Sine.h"

extern int ESTIMATION,CRITERION;
extern struct stat st;
extern int NUM_THREADS;
extern int ENABLE_DATA_PARALLELISM;

Experiments::Experiments()
{
  iterations = 1;
}

void Experiments::fisher_uncertainty()
{
  int N = 10;

  double kappa1 = 1;
  double kappa2 = 1;
  for (double rho=0.1; rho<0.95; rho+=0.1) {
    double lambda = rho * sqrt(kappa1 * kappa2);
    BVM_Sine bvm_sine(kappa1,kappa2,lambda);
    double log_det_fisher = bvm_sine.computeLogFisherInformation(N);
    double log_vscale = -2.5*logLatticeConstant(5) - 0.5 * log_det_fisher;
    cout << "rho: " << rho << "; log_vscale: " << log_vscale << endl;
  }
}

//void Experiments::simulate_sine(int iterations)
void Experiments::simulate_sine(struct Parameters &parameters)
{
  iterations = parameters.iterations;
  double kappa1 = parameters.k1;
  double kappa2 = parameters.k2;
  double rho = parameters.rho;
  int N = parameters.exp_sample_size;

  cout << "iterations: " << iterations << endl;
  cout << "kappa1: " << kappa1 << endl;
  cout << "kappa2: " << kappa2 << endl;
  cout << "rho: " << rho << endl;

  //int N = 10;

  string n_str = "N_" + boost::lexical_cast<string>(N);
  string parent_dir = "experiments/single_component/bvm_sine/estimates/" + n_str + "/";
  check_and_create_directory(parent_dir);

  std::vector<struct EstimatesSine> all_estimates;
  Vector statistics,pvalues;

  double mu1 = 90; mu1 *= PI/180;
  double mu2 = 90; mu2 *= PI/180;
  //double kappa1 = 100; 
  //double kappa2 = 100;

  ostringstream ssk1;
  ssk1 << fixed << setprecision(0);
  ssk1 << kappa1;
  string kappa1_str = ssk1.str();
  ostringstream ssk2;
  ssk2 << fixed << setprecision(0);
  ssk2 << kappa2;
  string kappa2_str = ssk2.str();

  //double rho = 0.1;
  //while (rho <= 0.95) {
    double lambda = rho * sqrt(kappa1 * kappa2);
    ostringstream ssr;
    ssr << fixed << setprecision(1);
    ssr << rho;
    string rho_str = ssr.str();

    string current_dir = parent_dir + "k1_" + kappa1_str + "_k2_" + kappa2_str 
                         + "_r_" + rho_str + "/";
    check_and_create_directory(current_dir);

    string mu1_est = current_dir + "mu1_est";
    string mu2_est = current_dir + "mu2_est";
    string kappa1_est = current_dir + "kappa1_est";
    string kappa2_est = current_dir + "kappa2_est";
    string lambda_est = current_dir + "lambda_est";
    string rho_file = current_dir + "rho_est";
    ofstream fmu1(mu1_est.c_str(),ios::app);
    ofstream fmu2(mu2_est.c_str(),ios::app);
    ofstream fkappa1(kappa1_est.c_str(),ios::app);
    ofstream fkappa2(kappa2_est.c_str(),ios::app);
    ofstream flambda(lambda_est.c_str(),ios::app);
    ofstream frho(rho_file.c_str(),ios::app);

    string negloglike = current_dir + "negloglike";
    string kldivs = current_dir + "kldivs";
    string msglens = current_dir + "msglens";
    string chisq_stat = current_dir + "chisq_stat";
    string pvalues_file = current_dir + "pvalues";
    ofstream fnlh(negloglike.c_str(),ios::app);
    ofstream fkl(kldivs.c_str(),ios::app);
    ofstream fmsg(msglens.c_str(),ios::app);
    ofstream fchi(chisq_stat.c_str(),ios::app);
    ofstream fpval(pvalues_file.c_str(),ios::app);
    
    cout << "kappa1: " << kappa1 << "; kappa2: " << kappa2 << "; rho: " << rho << endl;
    BVM_Sine bvm_sine(mu1,mu2,kappa1,kappa2,lambda);

    for (int i=0; i<iterations; i++) {
      repeat:
      cout << "Iteration: " << i+1 << endl;
      std::vector<Vector> random_sample = bvm_sine.generate(N);
      bvm_sine.computeAllEstimators(random_sample,all_estimates,0,1);
      Vector msglens(all_estimates.size(),0);

      // ignore PMLE
      //all_estimates[PMLE] = all_estimates[MLE];
      //msglens[PMLE] = HUGE_VAL;
      for (int j=0; j<all_estimates.size(); j++) {
        if (all_estimates[j].kldiv < 0 || all_estimates[j].msglen < 0) goto repeat;
        msglens[j] = all_estimates[j].msglen;
      }
      int min_index = minimumIndex(msglens);
      if (min_index != MML) {  // ignore iteration
        cout << "min_index != MML; repeating iteration ...\n";
        goto repeat;
      }
      chisquare_hypothesis_testing(all_estimates,statistics,pvalues);
      for (int j=0; j<all_estimates.size(); j++) {
        fmu1 << scientific << all_estimates[j].mu1 << "\t";
        fmu2 << scientific << all_estimates[j].mu2 << "\t";
        fkappa1 << scientific << all_estimates[j].kappa1 << "\t";
        fkappa2 << scientific << all_estimates[j].kappa2 << "\t";
        flambda << scientific << all_estimates[j].lambda << "\t";
        frho << scientific << all_estimates[j].rho << "\t";

        fnlh << scientific << all_estimates[j].negloglike << "\t";
        fkl << scientific << all_estimates[j].kldiv << "\t";
        fmsg << scientific << all_estimates[j].msglen << "\t";
        fchi << scientific << statistics[j] << "\t";
        fpval << scientific << pvalues[j] << "\t";
      } // for j ()
      fmu1 << endl; fmu2 << endl; 
      fkappa1 << endl; fkappa2 << endl; flambda << endl; frho << endl; 
      
      fnlh << endl; fkl << endl; fmsg << endl;
      fchi << endl; fpval << endl;
    } // for i ()

    rho += 0.1;

    fmu1.close(); fmu2.close();
    fkappa1.close(); fkappa2.close(); flambda.close(); frho.close();
    fnlh.close(); fkl.close(); fmsg.close(); fchi.close(); fpval.close();
  //} // rho
}

/*void Experiments::traditional_search()
{
}
*/

void Experiments::traditional_search()
{
  string exp_folder = "./experiments/mixture/";

  string data_file = "random_sample_mix.dat";
  std::vector<Vector> data = load_data_table(data_file);

  string all_inferred_mix = exp_folder + "mml_log";
  ofstream mml_log(all_inferred_mix.c_str());

  int K_MAX = 15; 

  Vector data_weights(data.size(),1);
  Mixture_Sine mixture(3,data,data_weights);
  mixture.estimateParameters();

  double msglen_best = mixture.getMinimumMessageLength();
  mml_log << fixed << setw(10) << 3 << "\t\t";
  mml_log << fixed << scientific << mixture.first_part() << "\t\t";
  mml_log << fixed << scientific << mixture.second_part() << "\t\t";
  mml_log << fixed << scientific << msglen_best << endl;
  Mixture_Sine mml_best_mix = mixture;

  for (int k=4; k<=K_MAX; k++) {
    cout << "k: " << k << endl;
    Mixture_Sine mixture(k,data,data_weights);
    mixture.estimateParameters();
    double msglen = mixture.getMinimumMessageLength();
    mml_log << fixed << setw(10) << k << "\t\t";
    mml_log << fixed << scientific << mixture.first_part() << "\t\t";
    mml_log << fixed << scientific << mixture.second_part() << "\t\t";
    mml_log << fixed << scientific << msglen << endl;
    if (msglen < msglen_best) {
      msglen_best = msglen;
      mml_best_mix = mixture;
    }
  } // for()
  mml_log.close();

  /* save inferred mixtures */
  string inferred_mix_file = exp_folder + "best_mml_mixture";
  mml_best_mix.printParameters(inferred_mix_file);
}

