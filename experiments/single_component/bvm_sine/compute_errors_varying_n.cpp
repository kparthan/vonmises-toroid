#include <iostream>
#include <fstream>
#include <sstream>
#include <vector> 
#include <cstdlib>
#include <iomanip>
#include <sys/stat.h>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::program_options;
using namespace boost::filesystem;

typedef std::vector<double> Vector;

#define NUM_METHODS 4 
#define MLE 0 
#define MAP 1
#define MAP_TRANSFORM 2 
#define MML 3

#define MU1 M_PI/2.0
#define MU2 M_PI/2.0

struct stat st = {0};
string kappa1_str,kappa2_str,rho_str,n_str,errors_folder;

std::vector<Vector> load_table(string &file_name, int D)
{
  std::vector<Vector> sample;
  ifstream file(file_name.c_str());
  string line;
  Vector numbers(D,0);
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
    sample.push_back(numbers);
  }
  file.close();
  return sample;
}

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

void quicksort(Vector &list, std::vector<int> &index, int left, int right)
{
	if(left < right)
	{
		int pivotNewIndex = partition(list,index,left,right);
		quicksort(list,index,left,pivotNewIndex-1);
		quicksort(list,index,pivotNewIndex+1,right);
	}
}

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

int minimumIndex(int map_index, Vector &values)
{
  Vector new_list(3,0);
  new_list[MLE] = values[MLE];
  if (map_index == MAP) new_list[MAP] = values[MAP];
  else new_list[MAP] = values[MAP_TRANSFORM];
  new_list[MAP+1] = values[MML];

  int min_index = 0;
  double min_val = new_list[0];
  for (int i=1; i<new_list.size(); i++) { 
    if (new_list[i] <= min_val) {
      min_index = i;
      min_val = new_list[i];
    } // if()
  } // for()
  return min_index;
}

double computeMean(Vector &list)
{
  double sum = 0;
  for (int i=0; i<list.size(); i++) {
    sum += list[i];
  }
  return sum / (double)list.size();
}

// column means
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

/* S **********************************************************************************/ 

Vector computeBiasSquared(std::vector<Vector> &p_est, double p)
{
  Vector p_est_means = computeMeans(p_est);
  Vector biassq(NUM_METHODS,0);
  for (int i=0; i<NUM_METHODS; i++) {
    double diff = p_est_means[i] - p;
    biassq[i] = diff * diff;
  }
  return biassq;
}

Vector computeVariance(std::vector<Vector> &p_est, double p)
{
  Vector p_est_means = computeMeans(p_est);
  Vector variance(NUM_METHODS,0);
  for (int i=0; i<p_est.size(); i++) {
    for (int j=0; j<NUM_METHODS; j++) {
      double diff = p_est_means[j] - p_est[i][j];
      variance[j] += (diff * diff);
    }
  }
  for (int j=0; j<NUM_METHODS; j++) {
    variance[j] /= p_est.size();
  }
  return variance;
}

Vector computeMeanSquaredError(std::vector<Vector> &p_est, double p)
{
  Vector p_est_means = computeMeans(p_est);
  Vector mse(NUM_METHODS,0);
  for (int i=0; i<p_est.size(); i++) {
    for (int j=0; j<NUM_METHODS; j++) {
      double diff = p - p_est[i][j];
      mse[j] += (diff * diff);
    }
  }
  for (int j=0; j<NUM_METHODS; j++) {
    mse[j] /= p_est.size();
  }
  return mse;
}

void compute_mu1_errors(int N) 
{
  //string n_str = "./estimates/N_" + boost::lexical_cast<string>(N) + "/";
  std::vector<Vector> mu1_est;
  Vector biassq_est,variance_est,mse_est;

  string biassq_file = errors_folder + "biassq_mu1";
  string variance_file = errors_folder + "variance_mu1";
  string mse_file = errors_folder + "mse_mu1";
  ofstream biassq(biassq_file.c_str(),ios::app);
  ofstream variance(variance_file.c_str(),ios::app);
  ofstream mse(mse_file.c_str(),ios::app);

  string current_dir = n_str + 
                       + "k1_" + kappa1_str 
                       + "_k2_" + kappa2_str 
                       + "_r_" + rho_str + "/";
  //cout << "current_dir: " << current_dir << endl;

  string mu1_file = current_dir + "mu1_est";
  mu1_est = load_table(mu1_file,NUM_METHODS);
  biassq_est = computeBiasSquared(mu1_est,MU1);
  variance_est = computeVariance(mu1_est,MU1);
  mse_est = computeMeanSquaredError(mu1_est,MU1);

  biassq << fixed << setw(10) << setprecision(0) << N << "\t";
  variance << fixed << setw(10) << setprecision(0) << N << "\t";
  mse << fixed << setw(10) << setprecision(0) << N << "\t";
  for(int i=0; i<NUM_METHODS; i++) {
    biassq << scientific << setprecision(6) << biassq_est[i] << "\t";
    variance << scientific << setprecision(6) << variance_est[i] << "\t";
    mse << scientific << setprecision(6) << mse_est[i] << "\t";
  }
  biassq << endl;
  variance << endl;
  mse << endl;

  biassq.close();
  variance.close();
  mse.close();
}

void compute_mu2_errors(int N) 
{
  std::vector<Vector> mu2_est;
  Vector biassq_est,variance_est,mse_est;

  string biassq_file = errors_folder + "biassq_mu2";
  string variance_file = errors_folder + "variance_mu2";
  string mse_file = errors_folder + "mse_mu2";
  ofstream biassq(biassq_file.c_str(),ios::app);
  ofstream variance(variance_file.c_str(),ios::app);
  ofstream mse(mse_file.c_str(),ios::app);

  string current_dir = n_str + 
                       + "k1_" + kappa1_str 
                       + "_k2_" + kappa2_str 
                       + "_r_" + rho_str + "/";

  string mu2_file = current_dir + "mu2_est";
  mu2_est = load_table(mu2_file,NUM_METHODS);
  biassq_est = computeBiasSquared(mu2_est,MU2);
  variance_est = computeVariance(mu2_est,MU2);
  mse_est = computeMeanSquaredError(mu2_est,MU2);

  biassq << fixed << setw(10) << setprecision(0) << N << "\t";
  variance << fixed << setw(10) << setprecision(0) << N << "\t";
  mse << fixed << setw(10) << setprecision(0) << N << "\t";
  for(int i=0; i<NUM_METHODS; i++) {
    biassq << scientific << setprecision(6) << biassq_est[i] << "\t";
    variance << scientific << setprecision(6) << variance_est[i] << "\t";
    mse << scientific << setprecision(6) << mse_est[i] << "\t";
  }
  biassq << endl;
  variance << endl;
  mse << endl;

  biassq.close();
  variance.close();
  mse.close();
}

void compute_kappa1_errors(
  int N, double kappa1
) {
  std::vector<Vector> kappas_est;
  Vector biassq_est,variance_est,mse_est;

  string biassq_file = errors_folder + "biassq_kappa1";
  string variance_file = errors_folder + "variance_kappa1";
  string mse_file = errors_folder + "mse_kappa1";
  ofstream biassq(biassq_file.c_str(),ios::app);
  ofstream variance(variance_file.c_str(),ios::app);
  ofstream mse(mse_file.c_str(),ios::app);

  string current_dir = n_str + 
                       + "k1_" + kappa1_str 
                       + "_k2_" + kappa2_str 
                       + "_r_" + rho_str + "/";

  string kappas_file = current_dir + "kappa1_est";
  kappas_est = load_table(kappas_file,NUM_METHODS);
  biassq_est = computeBiasSquared(kappas_est,kappa1);
  variance_est = computeVariance(kappas_est,kappa1);
  mse_est = computeMeanSquaredError(kappas_est,kappa1);

  biassq << fixed << setw(10) << setprecision(0) << N << "\t";
  variance << fixed << setw(10) << setprecision(0) << N << "\t";
  mse << fixed << setw(10) << setprecision(0) << N << "\t";
  for(int i=0; i<NUM_METHODS; i++) {
    biassq << scientific << setprecision(6) << biassq_est[i] << "\t";
    variance << scientific << setprecision(6) << variance_est[i] << "\t";
    mse << scientific << setprecision(6) << mse_est[i] << "\t";
  }
  biassq << endl;
  variance << endl;
  mse << endl;

  biassq.close();
  variance.close();
  mse.close();
}

void compute_kappa2_errors(
  int N, double kappa2
) {
  std::vector<Vector> kappas_est;
  Vector biassq_est,variance_est,mse_est;

  string biassq_file = errors_folder + "biassq_kappa2";
  string variance_file = errors_folder + "variance_kappa2";
  string mse_file = errors_folder + "mse_kappa2";
  ofstream biassq(biassq_file.c_str(),ios::app);
  ofstream variance(variance_file.c_str(),ios::app);
  ofstream mse(mse_file.c_str(),ios::app);

  string current_dir = n_str + 
                       + "k1_" + kappa1_str 
                       + "_k2_" + kappa2_str 
                       + "_r_" + rho_str + "/";

  string kappas_file = current_dir + "kappa2_est";
  kappas_est = load_table(kappas_file,NUM_METHODS);
  biassq_est = computeBiasSquared(kappas_est,kappa2);
  variance_est = computeVariance(kappas_est,kappa2);
  mse_est = computeMeanSquaredError(kappas_est,kappa2);

  biassq << fixed << setw(10) << setprecision(0) << N << "\t";
  variance << fixed << setw(10) << setprecision(0) << N << "\t";
  mse << fixed << setw(10) << setprecision(0) << N << "\t";
  for(int i=0; i<NUM_METHODS; i++) {
    biassq << scientific << setprecision(6) << biassq_est[i] << "\t";
    variance << scientific << setprecision(6) << variance_est[i] << "\t";
    mse << scientific << setprecision(6) << mse_est[i] << "\t";
  }
  biassq << endl;
  variance << endl;
  mse << endl;

  biassq.close();
  variance.close();
  mse.close();
}

void compute_lambda_errors(
  int N, double kappa1, double kappa2, double rho
) {
  double lambda = rho * sqrt(kappa1 * kappa2);

  std::vector<Vector> kappas_est;
  Vector biassq_est,variance_est,mse_est;

  string biassq_file = errors_folder + "biassq_lambda";
  string variance_file = errors_folder + "variance_lambda";
  string mse_file = errors_folder + "mse_lambda";
  ofstream biassq(biassq_file.c_str(),ios::app);
  ofstream variance(variance_file.c_str(),ios::app);
  ofstream mse(mse_file.c_str(),ios::app);

  string current_dir = n_str + 
                       + "k1_" + kappa1_str 
                       + "_k2_" + kappa2_str 
                       + "_r_" + rho_str + "/";

  string kappas_file = current_dir + "lambda_est";
  kappas_est = load_table(kappas_file,NUM_METHODS);
  biassq_est = computeBiasSquared(kappas_est,lambda);
  variance_est = computeVariance(kappas_est,lambda);
  mse_est = computeMeanSquaredError(kappas_est,lambda);

  biassq << fixed << setw(10) << setprecision(0) << N << "\t";
  variance << fixed << setw(10) << setprecision(0) << N << "\t";
  mse << fixed << setw(10) << setprecision(0) << N << "\t";
  for(int i=0; i<NUM_METHODS; i++) {
    biassq << scientific << setprecision(6) << biassq_est[i] << "\t";
    variance << scientific << setprecision(6) << variance_est[i] << "\t";
    mse << scientific << setprecision(6) << mse_est[i] << "\t";
  }
  biassq << endl;
  variance << endl;
  mse << endl;

  biassq.close();
  variance.close();
  mse.close();
}

void compute_rho_errors(
  int N, double rho
) {
  std::vector<Vector> rho_est;
  Vector biassq_est,variance_est,mse_est;

  string biassq_file = errors_folder + "biassq_rho";
  string variance_file = errors_folder + "variance_rho";
  string mse_file = errors_folder + "mse_rho";
  ofstream biassq(biassq_file.c_str(),ios::app);
  ofstream variance(variance_file.c_str(),ios::app);
  ofstream mse(mse_file.c_str(),ios::app);

  string current_dir = n_str + 
                       + "k1_" + kappa1_str 
                       + "_k2_" + kappa2_str 
                       + "_r_" + rho_str + "/";

  string rho_file = current_dir + "rho_est";
  rho_est = load_table(rho_file,NUM_METHODS);
  biassq_est = computeBiasSquared(rho_est,rho);
  variance_est = computeVariance(rho_est,rho);
  mse_est = computeMeanSquaredError(rho_est,rho);

  biassq << fixed << setw(10) << setprecision(0) << N << "\t";
  variance << fixed << setw(10) << setprecision(0) << N << "\t";
  mse << fixed << setw(10) << setprecision(0) << N << "\t";
  for(int i=0; i<NUM_METHODS; i++) {
    biassq << scientific << setprecision(6) << biassq_est[i] << "\t";
    variance << scientific << setprecision(6) << variance_est[i] << "\t";
    mse << scientific << setprecision(6) << mse_est[i] << "\t";
  }
  biassq << endl;
  variance << endl;
  mse << endl;

  biassq.close();
  variance.close();
  mse.close();
}

void combine(std::vector<string> &files, string &output_file)
{
  std::vector<std::vector<Vector> > all_tables;

  for (int i=0; i<files.size(); i++) {
    std::vector<Vector> table = load_table(files[i],NUM_METHODS+1);
    all_tables.push_back(table);
  }

  int num_sample_sizes = all_tables[0].size();
  ofstream out(output_file.c_str());

  for (int i=0; i<num_sample_sizes; i++) {
    Vector sum(NUM_METHODS,0);
    for (int j=0; j<all_tables.size(); j++) {
      for (int k=0; k<NUM_METHODS; k++) {
        sum[k] += all_tables[j][i][k+1];
      } // k
    } // j
    out << fixed << setw(10) << setprecision(0) << all_tables[0][i][0] << "\t\t";
    for (int k=0; k<NUM_METHODS; k++) {
      out << scientific << setprecision(6) << sum[k] << "\t\t";
    } // k
    out << endl;
  } // i

  out.close();
}

void compute_all_errors() 
{
  std::vector<string> biassq_files,variance_files,mse_files;
  string biassq_file,variance_file,mse_file,output_file;

  biassq_file = errors_folder + "biassq_mu1";
  biassq_files.push_back(biassq_file);
  biassq_file = errors_folder + "biassq_mu2";
  biassq_files.push_back(biassq_file);
  biassq_file = errors_folder + "biassq_kappa1";
  biassq_files.push_back(biassq_file);
  biassq_file = errors_folder + "biassq_kappa2";
  biassq_files.push_back(biassq_file);
  biassq_file = errors_folder + "biassq_lambda";
  biassq_files.push_back(biassq_file);
  output_file = errors_folder + "biassq_all";
  combine(biassq_files,output_file);

  variance_file = errors_folder + "variance_mu1";
  variance_files.push_back(variance_file);
  variance_file = errors_folder + "variance_mu2";
  variance_files.push_back(variance_file);
  variance_file = errors_folder + "variance_kappa1";
  variance_files.push_back(variance_file);
  variance_file = errors_folder + "variance_kappa2";
  variance_files.push_back(variance_file);
  variance_file = errors_folder + "variance_lambda";
  variance_files.push_back(variance_file);
  output_file = errors_folder + "variance_all";
  combine(variance_files,output_file);

  mse_file = errors_folder + "mse_mu1";
  mse_files.push_back(mse_file);
  mse_file = errors_folder + "mse_mu2";
  mse_files.push_back(mse_file);
  mse_file = errors_folder + "mse_kappa1";
  mse_files.push_back(mse_file);
  mse_file = errors_folder + "mse_kappa2";
  mse_files.push_back(mse_file);
  mse_file = errors_folder + "mse_lambda";
  mse_files.push_back(mse_file);
  output_file = errors_folder + "mse_all";
  combine(mse_files,output_file);
}

//void computeWins(
//  int map_index, ostream &out, std::vector<Vector> &values
//) {
//  std::vector<int> wins(4,0);
//  int min_index;
//
//  for (int i=0; i<values.size(); i++) {
//    min_index = minimumIndex(map_index,values[i]);
//    wins[min_index]++;
//  }
//  double percent_wins;
//  for (int j=0; j<wins.size(); j++) {
//    percent_wins = wins[j] * 100.0 / values.size();
//    out << fixed << setw(10) << setprecision(2) << percent_wins;
//  }
//  out << endl;
//}
//
//void process_kldivs(int N)
//{
//  int num_maps,map_index;
//  if (NUM_METHODS == 5) {
//    num_maps = 2;
//  } else {
//    num_maps = 3;
//  }
//
//  string n_str = "./estimates/N_" + boost::lexical_cast<string>(N);
//  for (int i=1; i<=num_maps; i++) {
//    if (i == 1) map_index = 2;
//    else if (i == 2) map_index = 4;
//    else if (i == 3) map_index = 5;
//    string wins_kldivs_file = errors_folder + common + "wins_kldivs_map"
//                              + boost::lexical_cast<string>(i) + ".dat";
//    ofstream out(wins_kldivs_file.c_str(),ios::app);
//    out << fixed << setw(10) << setprecision(0) << N << "\t\t";
//    string kldivs_file = n_str + "_prior" + prior_str + "/"
//                         "k_" + kappa_str + "_e_" + rho_str + "/kldivs";
//    std::vector<Vector> kldivs = load_table(kldivs_file,NUM_METHODS);
//    computeWins(map_index,out,kldivs);
//    out.close();
//  }
//}

void process_estimates(double kappa1, double kappa2, double rho)
{
  for (int N=10; N<=50; N+=10) {
    n_str = "./estimates/N_" + boost::lexical_cast<string>(N) + "/";

    // mu1 errors
    compute_mu1_errors(N);
    // mu2 errors
    compute_mu2_errors(N);

    // kappa errors
    compute_kappa1_errors(N,kappa1);

    // kappa2 errors
    compute_kappa2_errors(N,kappa2);

    // lambda errors
    compute_lambda_errors(N,kappa1,kappa2,rho);

    // rho errors
    compute_rho_errors(N,rho);

    // combined errors
    compute_all_errors();

    // process kldivs
    process_kldivs(N);
  } // for(N)
}

//void common_plot(
//  string &data_file, string &script_file, string &plot_file, string &ylabel
//) {
//  ofstream out(script_file.c_str());
//  out << "set terminal postscript eps enhanced color\n\n";
//  out << "set output \"" << plot_file << "\"\n\n";
//  out << "set style data linespoints\n";
//  out << "set style fill solid 1.0 noborder\n";
//  out << "set xlabel \"Sample size\\n\"\n";
//  out << "set ylabel \"" << ylabel << "\"\n";
//  out << "set xr [10:50]\n";
//  out << "set xlabel font \"Times-Roman, 25\"\n";
//  out << "set ylabel font \"Times-Roman, 25\"\n";
//  out << "set xtics font \"Times-Roman, 20\"\n";
//  out << "set ytics font \"Times-Roman, 20\"\n";
//  out << "set xtics nomirror\n";
//  out << "set ytics nomirror\n";
//  out << "set border 2 back\n";
//  if (NUM_METHODS == 5) {
//    out << "plot \"" << data_file << "\" using 1:2 t \"MOMENT\" lc rgb \"red\", \\\n"
//        << "\"\" using 1:3 t \"MLE\" lc rgb \"blue\", \\\n"
//        << "\"\" using 1:4 t \"MAP1\" lc rgb \"dark-green\", \\\n"
//        << "\"\" using 1:5 t \"MML\" lc rgb \"dark-magenta\", \\\n"
//        << "\"\" using 1:6 t \"MAP2\" lc rgb \"black\"\n";
//  } else if (NUM_METHODS == 6) {
//    out << "plot \"" << data_file << "\" using 1:2 t \"MOMENT\" lc rgb \"red\", \\\n"
//        << "\"\" using 1:7 t \"MAP3 = MLE\" lc rgb \"blue\", \\\n"
//        << "\"\" using 1:4 t \"MAP1\" lc rgb \"dark-green\", \\\n"
//        << "\"\" using 1:5 t \"MML\" lc rgb \"dark-magenta\", \\\n"
//        << "\"\" using 1:6 t \"MAP2\" lc rgb \"black\"\n";
//  }
//  out.close();
//  string cmd = "gnuplot -persist " + script_file;
//  if(system(cmd.c_str()));
//}
//
//void plot_mu1_errors(
//  string &kappa_str, string &errors_folder
//) {
//  string data_file,script_file,plot_file,ylabel;
//
//  data_file = errors_folder + common + "biassq_mu1";
//  script_file = errors_folder + common + "biassq_mu1.p";
//  plot_file = errors_folder + common + "biassq_mu1.eps";
//  ylabel = "Bias squared";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "variance_mu1";
//  script_file = errors_folder + common + "variance_mu1.p";
//  plot_file = errors_folder + common + "variance_mu1.eps";
//  ylabel = "Variance";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "mse_mu1";
//  script_file = errors_folder + common + "mse_mu1.p";
//  plot_file = errors_folder + common + "mse_mu1.eps";
//  ylabel = "Mean squared error";
//  common_plot(data_file,script_file,plot_file,ylabel);
//}
//
//void plot_mu2_errors(
//  string &kappa_str, string &errors_folder
//) {
//  string data_file,script_file,plot_file,ylabel;
//
//  data_file = errors_folder + common + "biassq_mu2";
//  script_file = errors_folder + common + "biassq_mu2.p";
//  plot_file = errors_folder + common + "biassq_mu2.eps";
//  ylabel = "Bias squared";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "variance_mu2";
//  script_file = errors_folder + common + "variance_mu2.p";
//  plot_file = errors_folder + common + "variance_mu2.eps";
//  ylabel = "Variance";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "mse_mu2";
//  script_file = errors_folder + common + "mse_mu2.p";
//  plot_file = errors_folder + common + "mse_mu2.eps";
//  ylabel = "Mean squared error";
//  common_plot(data_file,script_file,plot_file,ylabel);
//}
//
//void plot_eta_errors(
//  string &kappa_str, string &errors_folder
//) {
//  string data_file,script_file,plot_file,ylabel;
//
//  data_file = errors_folder + common + "biassq_eta";
//  script_file = errors_folder + common + "biassq_eta.p";
//  plot_file = errors_folder + common + "biassq_eta.eps";
//  ylabel = "Bias squared";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "variance_eta";
//  script_file = errors_folder + common + "variance_eta.p";
//  plot_file = errors_folder + common + "variance_eta.eps";
//  ylabel = "Variance";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "mse_eta";
//  script_file = errors_folder + common + "mse_eta.p";
//  plot_file = errors_folder + common + "mse_eta.eps";
//  ylabel = "Mean squared error";
//  common_plot(data_file,script_file,plot_file,ylabel);
//}
//
//void plot_kappa_errors(
//  string &kappa_str, string &errors_folder
//) {
//  string data_file,script_file,plot_file,ylabel;
//
//  data_file = errors_folder + common + "biassq_kappa";
//  script_file = errors_folder + common + "biassq_kappa.p";
//  plot_file = errors_folder + common + "biassq_kappa.eps";
//  ylabel = "Bias squared";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "variance_kappa";
//  script_file = errors_folder + common + "variance_kappa.p";
//  plot_file = errors_folder + common + "variance_kappa.eps";
//  ylabel = "Variance";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "mse_kappa";
//  script_file = errors_folder + common + "mse_kappa.p";
//  plot_file = errors_folder + common + "mse_kappa.eps";
//  ylabel = "Mean squared error";
//  common_plot(data_file,script_file,plot_file,ylabel);
//}
//
//void plot_kappa2_errors(
//  string &kappa2_str, string &errors_folder
//) {
//  string data_file,script_file,plot_file,ylabel;
//
//  data_file = errors_folder + common + "biassq_kappa2";
//  script_file = errors_folder + common + "biassq_kappa2.p";
//  plot_file = errors_folder + common + "biassq_kappa2.eps";
//  ylabel = "Bias squared";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "variance_kappa2";
//  script_file = errors_folder + common + "variance_kappa2.p";
//  plot_file = errors_folder + common + "variance_kappa2.eps";
//  ylabel = "Variance";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "mse_kappa2";
//  script_file = errors_folder + common + "mse_kappa2.p";
//  plot_file = errors_folder + common + "mse_kappa2.eps";
//  ylabel = "Mean squared error";
//  common_plot(data_file,script_file,plot_file,ylabel);
//}
//
//void plot_rho_errors(
//  string &rho_str, string &errors_folder
//) {
//  string data_file,script_file,plot_file,ylabel;
//
//  data_file = errors_folder + common + "biassq_rho";
//  script_file = errors_folder + common + "biassq_rho.p";
//  plot_file = errors_folder + common + "biassq_rho.eps";
//  ylabel = "Bias squared";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "variance_rho";
//  script_file = errors_folder + common + "variance_rho.p";
//  plot_file = errors_folder + common + "variance_rho.eps";
//  ylabel = "Variance";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "mse_rho";
//  script_file = errors_folder + common + "mse_rho.p";
//  plot_file = errors_folder + common + "mse_rho.eps";
//  ylabel = "Mean squared error";
//  common_plot(data_file,script_file,plot_file,ylabel);
//}
//
//void plot_all_errors(
//  string &all_str, string &errors_folder
//) {
//  string data_file,script_file,plot_file,ylabel;
//
//  data_file = errors_folder + common + "biassq_all";
//  script_file = errors_folder + common + "biassq_all.p";
//  plot_file = errors_folder + common + "biassq_all.eps";
//  ylabel = "Bias squared";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "variance_all";
//  script_file = errors_folder + common + "variance_all.p";
//  plot_file = errors_folder + common + "variance_all.eps";
//  ylabel = "Variance";
//  common_plot(data_file,script_file,plot_file,ylabel);
//
//  data_file = errors_folder + common + "mse_all";
//  script_file = errors_folder + common + "mse_all.p";
//  plot_file = errors_folder + common + "mse_all.eps";
//  ylabel = "Mean squared error";
//  common_plot(data_file,script_file,plot_file,ylabel);
//}
//
//void plot_errors()
//{
//  // mu1 errors
//  plot_mu1_errors(kappa_str,errors_folder);
//  // mu2 errors
//  plot_mu2_errors(kappa_str,errors_folder);
//  // eta errors
//  plot_eta_errors(kappa_str,errors_folder);
//
//  // kappa errors
//  plot_kappa_errors(kappa_str,errors_folder);
//  // kappa2 errors
//  plot_kappa2_errors(kappa_str,errors_folder);
//  // rho errors
//  plot_rho_errors(kappa_str,errors_folder);
//
//  // all errors
//  plot_all_errors(kappa_str,errors_folder);
//}
//
//void plot_script_kldivs_wins(string &kldivs_folder, int num_map)
//{
//  string map = boost::lexical_cast<string>(num_map);
//  string wins_kldivs_file = kldivs_folder + common + "wins_kldivs_map" + map + ".dat";
//  string script_file = kldivs_folder + common + "wins_kldivs_map" + map + ".p";
//  string plot_file = kldivs_folder + common + "wins_kldivs_map" + map + ".eps";
//  ofstream out(script_file.c_str());
//  out << "set terminal postscript eps enhanced color\n\n";
//  out << "set output \"" << plot_file << "\"\n\n";
//  out << "set grid y\n";
//  out << "set style data histograms\n";
//  out << "set style fill solid 1.0 noborder\n";
//  out << "set ytics 10 nomirror\n";
//  out << "set yrange [:100]\n";
//  out << "set xlabel font \"Times-Roman, 25\"\n";
//  out << "set ylabel font \"Times-Roman, 25\"\n";
//  out << "set xtics font \"Times-Roman, 20\"\n";
//  out << "set ytics font \"Times-Roman, 20\"\n";
//  out << "set xlabel \"Sample size\\n\"\n";
//  out << "set ylabel \"\% of wins\"\n";
//  out << "set ytics 10\n\n"; 
//  out << "set xtics nomirror\n";
//  out << "set border 2 back\n";
//  if (num_map != 3) {
//    out << "plot \"" << wins_kldivs_file << "\" using 2 t \"MOMENT\" lc rgb \"red\", \\\n"
//        << "\"\" using 3 t \"MLE\" lc rgb \"blue\", \\\n"
//        << "\"\" using 4 t \"MAP" << num_map << "\" lc rgb \"dark-green\", \\\n"
//        << "\"\" using 5:xtic(1) t \"MML\" lc rgb \"dark-magenta\"";
//  } else {
//    out << "plot \"" << wins_kldivs_file << "\" using 2 t \"MOMENT\" lc rgb \"red\", \\\n"
//        << "\"\" using 4 t \"MAP3 = MLE\" lc rgb \"dark-green\", \\\n"
//        << "\"\" using 5:xtic(1) t \"MML\" lc rgb \"dark-magenta\"";
//  }
//  out.close();
//  string cmd = "gnuplot -persist " + script_file;
//  if(system(cmd.c_str()));
//}
//
//void boxplot_test_stats_fixed_kappa(
//  string &errors_folder, int option
//) {
//  string stats_file,script_file,plot_file,ylabel;
//
//  if (option == 1) {
//    script_file = errors_folder + "boxplot_test_statistics.p";
//    plot_file = errors_folder + common + "boxplot_test_statistics.eps";
//    ylabel = "Test statistic";
//  } else if (option == 2) {
//    script_file = errors_folder + "boxplot_pvalues.p";
//    plot_file = errors_folder + common + "boxplot_pvalues.eps";
//    ylabel = "p-values";
//  }
//
//  ofstream out(script_file.c_str());
//  out << "set terminal post eps enhanced color\n";
//  out << "set output \"" << plot_file << "\"\n\n";
//  out << "box_width=0.12\n";
//  out << "set style fill solid 0.25 #noborder\n";
//  out << "set style boxplot outliers pointtype 7\n";
//  out << "set style data boxplot\n";
//  out << "set boxwidth box_width #relative\n";
//  out << "set pointsize 0.5\n";
//  if (option == 1) {
//    out << "set key top right\n";
//  } else if (option == 2) {
//    out << "set key bottom right\n";
//  }
//  out << "set border 2\n";
//  out << "set xtics nomirror\n";
//  out << "set ytics nomirror\n";
//  out << "set xlabel \"Sample size\\n\"\n";
//  out << "set ylabel \"" << ylabel << "\"\n";
//  out << "set xlabel font \"Times-Roman, 25\"\n";
//  out << "set ylabel font \"Times-Roman, 25\"\n";
//  out << "set xtics font \"Times-Roman, 20\"\n";
//  out << "set ytics font \"Times-Roman, 20\"\n";
//  out << "set xtics (\"10\" 1, \"15\" 2, \"20\" 3, "
//      << "\"25\" 4, \"30\" 5, \"35\" 6, \"40\" 7, \"45\" 8, \"50\" 9) scale 0.0\n";
//  out << "d_width=0.5*box_width\n\n";
//
//  int col = 1;
//  for (int N=10; N<=50; N+=5) {
//    string n_str = boost::lexical_cast<string>(N);
//    if (option == 1) {
//      stats_file = "./estimates/N_" + n_str + "_prior" + prior_str + "/k_" + kappa_str + "_e_" + rho_str + "/chisq_stat";
//    } else if (option == 2) {
//      stats_file = "./estimates/N_" + n_str + "_prior" + prior_str + "/k_" + kappa_str + "_e_" + rho_str + "/pvalues";
//    }
//    if (col == 1) {
//      out << "plot ";
//    }
//    if (col < 9) {
//      out << "\"" << stats_file << "\" using ((" << col << ")-3*d_width):1 notitle lt 1 lc rgb \"red\", \\\n"
//          << "\"" << stats_file << "\" using ((" << col << ")-d_width):3 notitle lt 1 lc rgb \"dark-green\", \\\n"
//          << "\"" << stats_file << "\" using ((" << col << ")+d_width):4 notitle lt 1 lc rgb \"dark-magenta\", \\\n"
//          << "\"" << stats_file << "\" using ((" << col << ")+3*d_width):5 notitle lt 1 lc rgb \"black\", \\\n";
//    } else if (col == 9) {
//      out << "\"" << stats_file << "\" using ((" << col << ")-3*d_width):1 title \"MOMENT\" lt 1 lc rgb \"red\", \\\n"
//          << "\"" << stats_file << "\" using ((" << col << ")-d_width):3 title \"MAP1\" lt 1 lc rgb \"dark-green\", \\\n"
//          << "\"" << stats_file << "\" using ((" << col << ")+d_width):4 title \"MML\" lt 1 lc rgb \"dark-magenta\", \\\n"
//          << "\"" << stats_file << "\" using ((" << col << ")+3*d_width):5 title \"MAP2\" lt 1 lc rgb \"black\"\n";
//    }
//    col++;
//  }
//  out.close();
//
//  string cmd = "gnuplot -persist " + script_file;
//  if(system(cmd.c_str()));
//}
//
//void plot_kldivs()
//{
//  int num_maps,map_index;
//  if (NUM_METHODS == 5) {
//    num_maps = 2;
//  } else {
//    num_maps = 3;
//  }
//
//  for (int i=1; i<=num_maps; i++) {
//    plot_script_kldivs_wins(errors_folder,i);
//  }
//
//  boxplot_test_stats_fixed_kappa(errors_folder,1);
//  boxplot_test_stats_fixed_kappa(errors_folder,2);
//}
//
/* E **********************************************************************************/ 

/* S **********************************************************************************/ 
void check_and_create_directory(string &directory)
{
  if (stat(directory.c_str(), &st) == -1) {
    mkdir(directory.c_str(), 0700);
  }
}

void create_required_folders(
  double kappa1, double kappa2, double rho
) {
  ostringstream ssk1;
  ssk1 << fixed << setprecision(0);
  ssk1 << kappa1;
  kappa1_str = ssk1.str();

  ostringstream ssk2;
  ssk2 << fixed << setprecision(0);
  ssk2 << kappa2;
  kappa2_str = ssk2.str();

  ostringstream ssr;
  ssr << fixed << setprecision(1);
  ssr << rho;
  rho_str = ssr.str();

  errors_folder = "./analysis/k1_" + kappa1_str + "_k2_" + kappa2_str 
                         + "_r_" + rho_str + "/";
  check_and_create_directory(errors_folder);

  /*for (int N=10; N<=50; N+=10) {
    string n_str = "N_" + boost::lexical_cast<string>(N) + "/";
    string n_dir = errors_folder + n_str;
    check_and_create_directory(n_dir);
  }*/
}

/* E **********************************************************************************/

struct Parameters
{
  double kappa1,kappa2,rho;
};

struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;

  options_description desc("Allowed options");
  desc.add_options()
       ("k1",value<double>(&parameters.kappa1),"kappa1")
       ("k2",value<double>(&parameters.kappa2),"kappa2")
       ("rho",value<double>(&parameters.rho),"rho")
  ;
  variables_map vm;
  store(command_line_parser(argc,argv).options(desc).run(),vm);
  notify(vm);

  return parameters;
}

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);
  double kappa1 = parameters.kappa1;
  double kappa2 = parameters.kappa2;
  double rho = parameters.rho;

  //for (int kappa1=1; kappa1<=100; kappa1*=10) {
    //for (int kappa2=1; kappa2<=100; kappa2*=10) {
      //for (double rho=0; rho<1; rho+=0.1) {
        create_required_folders(kappa1,kappa2,rho);

        process_estimates(kappa1,kappa2,rho);

        //plot_errors();

        //plot_kldivs();
      //} // rho()
    //} // kappa2()
  //} // kappa1()
}

