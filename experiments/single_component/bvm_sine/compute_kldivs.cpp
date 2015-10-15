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

#define MOMENT 0 
#define MLE 1
#define MAP 2
#define MML 3 
#define MAP_ECCENTRICITY_TRANSFORM 4
#define MAP_UNIFORM_TRANSFORM 5

#define INIT_KAPPA 1
#define MAX_KAPPA 100
#define KAPPA_INCREMENT 10

struct stat st = {0};
int NUM_METHODS;
string common,common_kappa;

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
  Vector new_list(4,0);
  new_list[MOMENT] = values[MOMENT];
  if (map_index == 5) new_list[MLE] = values[map_index];
  else new_list[MLE] = values[MLE];
  new_list[MAP] = values[map_index];
  new_list[MML] = values[MML];

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
void plot_script_kldivs_wins(string &kldivs_folder, int num_map)
{
  string map = boost::lexical_cast<string>(num_map);
  string wins_kldivs_file = kldivs_folder + "wins_kldivs_map" + map + ".dat";
  string script_file = kldivs_folder + "wins_kldivs_map" + map + ".p";
  string plot_file = kldivs_folder + common_kappa + "wins_kldivs_map" + map + ".eps";
  ofstream out(script_file.c_str());
  out << "set terminal postscript eps enhanced color\n\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "set grid y\n";
  out << "set style data histograms\n";
  out << "set style fill solid 1.0 noborder\n";
  out << "set ytics 10 nomirror\n";
  out << "set yrange [:100]\n";
  out << "set xlabel \"eccentricity\\n\"\n";
  out << "set ylabel \"\% of wins\"\n";
  out << "set ytics 10\n\n"; 
  out << "set xtics nomirror\n";
  out << "set ytics nomirror\n";
  out << "set border 2 back\n";
  out << "set xlabel font \"Times-Roman, 25\"\n";
  out << "set ylabel font \"Times-Roman, 25\"\n";
  out << "set xtics font \"Times-Roman, 20\"\n";
  out << "set ytics font \"Times-Roman, 20\"\n";
  if (num_map != 3) {
    out << "plot \"" << wins_kldivs_file << "\" using 2 t \"MOMENT\" lc rgb \"red\", \\\n"
        << "\"\" using 3 t \"MLE\" lc rgb \"blue\", \\\n"
        << "\"\" using 4 t \"MAP" << num_map << "\" lc rgb \"dark-green\", \\\n"
        << "\"\" using 5:xtic(1) t \"MML\" lc rgb \"dark-magenta\"";
  } else {
    out << "plot \"" << wins_kldivs_file << "\" using 2 t \"MOMENT\" lc rgb \"red\", \\\n"
        << "\"\" using 4 t \"MAP3 = MLE\" lc rgb \"dark-green\", \\\n"
        << "\"\" using 5:xtic(1) t \"MML\" lc rgb \"dark-magenta\"";
  }
  out.close();
  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

void boxplot_kldivs_fixed_kappa(string &n_str, string &kappa_str, string &kldivs_folder)
{
  string kldivs_file;
  string script_file = kldivs_folder + "boxplot_kldivs.p";
  string plot_file = kldivs_folder + common_kappa + "boxplot_kldivs.eps";

  ofstream out(script_file.c_str());
  out << "set terminal post eps color enhanced\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "box_width=0.12\n";
  out << "set style fill solid 0.25 noborder\n";
  out << "set style boxplot outliers pointtype 7\n";
  out << "set style data boxplot\n";
  out << "set boxwidth box_width #relative\n";
  out << "set pointsize 0.5\n";
  out << "unset key\n";
  out << "set border 2\n";
  out << "set xtics nomirror\n";
  out << "set ytics nomirror\n";
  out << "set xlabel \"eccentricity\\n\"\n";
  out << "set ylabel \"KL-divergence\"\n";
  out << "set xlabel font \"Times-Roman, 25\"\n";
  out << "set ylabel font \"Times-Roman, 25\"\n";
  out << "set xtics font \"Times-Roman, 20\"\n";
  out << "set ytics font \"Times-Roman, 20\"\n";
  out << "set xtics (\"0.1\" 1, \"0.5\" 2, \"0.9\" 3) scale 0.0\n";
  out << "d_width=0.5*box_width\n\n";

  // e = 0.1
  kldivs_file = n_str + "k_" + kappa_str + "_e_0.1/kldivs";
  if (NUM_METHODS == 6) {
    out << "plot \"" << kldivs_file << "\" using ((1)-5*d_width):1 lt 1 lc rgb \"red\", \\\n"
        << "\"" << kldivs_file << "\" using ((1)-3*d_width):2 lt 1 lc rgb \"blue\", \\\n"
        << "\"" << kldivs_file << "\" using ((1)-d_width):3 lt 1 lc rgb \"dark-green\", \\\n"
        << "\"" << kldivs_file << "\" using ((1)+d_width):4 lt 1 lc rgb \"dark-magenta\", \\\n"
        << "\"" << kldivs_file << "\" using ((1)+3*d_width):5 lt 1 lc rgb \"pink\", \\\n"
        << "\"" << kldivs_file << "\" using ((1)+5*d_width):6 lt 1 lc rgb \"purple\", \\\n";
  } else if (NUM_METHODS == 5) {
    out << "plot \"" << kldivs_file << "\" using ((1)-4*d_width):1 t \"MOMENT\" lt 1 lc rgb \"red\", \\\n"
        << "\"" << kldivs_file << "\" using ((1)-2*d_width):2 lt 1 lc rgb \"blue\", \\\n"
        << "\"" << kldivs_file << "\" using (1):3 lt 1 lc rgb \"dark-green\", \\\n"
        << "\"" << kldivs_file << "\" using ((1)+2*d_width):4 lt 1 lc rgb \"dark-magenta\", \\\n"
        << "\"" << kldivs_file << "\" using ((1)+4*d_width):5 lt 1 lc rgb \"pink\", \\\n";
  }

  // e = 0.5
  kldivs_file = n_str + "k_" + kappa_str + "_e_0.5/kldivs";
  if (NUM_METHODS == 6) {
    out << "\"" << kldivs_file << "\" using ((2)-5*d_width):1 lt 1 lc rgb \"red\", \\\n"
        << "\"" << kldivs_file << "\" using ((2)-3*d_width):2 lt 1 lc rgb \"blue\", \\\n"
        << "\"" << kldivs_file << "\" using ((2)-d_width):3 lt 1 lc rgb \"dark-green\", \\\n"
        << "\"" << kldivs_file << "\" using ((2)+d_width):4 lt 1 lc rgb \"dark-magenta\", \\\n"
        << "\"" << kldivs_file << "\" using ((2)+3*d_width):5 lt 1 lc rgb \"pink\", \\\n"
        << "\"" << kldivs_file << "\" using ((2)+5*d_width):6 lt 1 lc rgb \"purple\", \\\n";
  } else if (NUM_METHODS == 5) {
    out << "\"" << kldivs_file << "\" using ((2)-4*d_width):1 lt 1 lc rgb \"red\", \\\n"
        << "\"" << kldivs_file << "\" using ((2)-2*d_width):2 lt 1 lc rgb \"blue\", \\\n"
        << "\"" << kldivs_file << "\" using (2):3 lt 1 lc rgb \"dark-green\", \\\n"
        << "\"" << kldivs_file << "\" using ((2)+2*d_width):4 lt 1 lc rgb \"dark-magenta\", \\\n"
        << "\"" << kldivs_file << "\" using ((2)+4*d_width):5 lt 1 lc rgb \"pink\", \\\n";
  }

  // e = 0.9
  kldivs_file = n_str + "k_" + kappa_str + "_e_0.9/kldivs";
  if (NUM_METHODS == 6) {
    out << "\"" << kldivs_file << "\" using ((3)-5*d_width):1 lt 1 lc rgb \"red\", \\\n"
        << "\"" << kldivs_file << "\" using ((3)-3*d_width):2 lt 1 lc rgb \"blue\", \\\n"
        << "\"" << kldivs_file << "\" using ((3)-d_width):3 lt 1 lc rgb \"dark-green\", \\\n"
        << "\"" << kldivs_file << "\" using ((3)+d_width):4 lt 1 lc rgb \"dark-magenta\", \\\n"
        << "\"" << kldivs_file << "\" using ((3)+3*d_width):5 lt 1 lc rgb \"pink\", \\\n"
        << "\"" << kldivs_file << "\" using ((3)+5*d_width):6 lt 1 lc rgb \"purple\"\n";
  } else if (NUM_METHODS == 5) {
    out << "\"" << kldivs_file << "\" using ((3)-4*d_width):1 lt 1 lc rgb \"red\", \\\n"
        << "\"" << kldivs_file << "\" using ((3)-2*d_width):2 lt 1 lc rgb \"blue\", \\\n"
        << "\"" << kldivs_file << "\" using (3):3 lt 1 lc rgb \"dark-green\", \\\n"
        << "\"" << kldivs_file << "\" using ((3)+2*d_width):4 lt 1 lc rgb \"dark-magenta\", \\\n"
        << "\"" << kldivs_file << "\" using ((3)+4*d_width):5 lt 1 lc rgb \"pink\"\n";
  }
  out.close();

  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

void computeWins(
  int map_index, ostream &out, std::vector<Vector> &values
) {
  std::vector<int> wins(4,0);
  int min_index;

  for (int i=0; i<values.size(); i++) {
    min_index = minimumIndex(map_index,values[i]);
    wins[min_index]++;
  }
  double percent_wins;
  for (int j=0; j<wins.size(); j++) {
    percent_wins = wins[j] * 100.0 / values.size();
    out << fixed << setw(10) << setprecision(2) << percent_wins;
  }
  out << endl;
}

void tabulate_kappa_kldivs(
  int map_index, string &output_file, string &n_str, string &kappa_str
) {
  std::vector<Vector> kldivs;
  string kldivs_file;
  ofstream out(output_file.c_str());
  double ecc = 0.1;
  while (ecc < 0.95) {
    out << fixed << setw(10) << setprecision(1) << ecc;
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    string ecc_str = sse.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";
    kldivs_file = current_dir + "kldivs";
    kldivs = load_table(kldivs_file,NUM_METHODS);
    computeWins(map_index,out,kldivs);
    ecc += 0.1;
  } // while()
  out.close();
}

void plot_avg_kldivs(string &kldivs_folder)
{
  string avg_kldivs_file = kldivs_folder + "kldivs_avg.dat";
  string script_file = kldivs_folder + "kldivs_avg.p";
  string plot_file = kldivs_folder + common_kappa + "kldivs_avg.eps";
  ofstream out(script_file.c_str());
  out << "set terminal postscript eps enhanced color\n\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "set style data linespoints\n";
  out << "set style fill solid 1.0 noborder\n";
  out << "set ylabel \"Average KL-divergence\"\n";
  out << "set xlabel font \"Times-Roman, 25\"\n";
  out << "set ylabel font \"Times-Roman, 25\"\n";
  out << "set xtics font \"Times-Roman, 20\"\n";
  out << "set ytics font \"Times-Roman, 20\"\n";
  out << "set xtics nomirror\n";
  out << "set ytics nomirror\n";
  out << "set border 2 back\n";
  if (NUM_METHODS == 6) {
    out << "plot \"" << avg_kldivs_file << "\" using 1:2 t \"MOMENT\" lc rgb \"red\", \\\n"
        << "\"\" using 1:7 t \"MAP3 = MLE\" lc rgb \"blue\", \\\n"
        << "\"\" using 1:4 t \"MAP1\" lc rgb \"dark-green\", \\\n"
        << "\"\" using 1:5 t \"MML\" lc rgb \"dark-magenta\", \\\n"
        << "\"\" using 1:6 t \"MAP2\" lc rgb \"black\"\n";
  } else if (NUM_METHODS == 5) {
    out << "plot \"" << avg_kldivs_file << "\" using 1:2 t \"MOMENT\" lc rgb \"red\", \\\n"
        << "\"\" using 1:3 t \"MLE\" lc rgb \"blue\", \\\n"
        << "\"\" using 1:4 t \"MAP1\" lc rgb \"dark-green\", \\\n"
        << "\"\" using 1:5 t \"MML\" lc rgb \"dark-magenta\", \\\n"
        << "\"\" using 1:6 t \"MAP2\" lc rgb \"black\"\n";
  }
  out.close();
  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

void compute_average_kldivs_fixed_kappa(
  string &n_str, string &kappa_str, string &kldivs_folder
) {
  std::vector<Vector> kldivs;
  string kldivs_file;
  string output_file = kldivs_folder + "kldivs_avg.dat";
  ofstream out(output_file.c_str());
  double ecc = 0.1;
  while (ecc < 0.95) {
    out << fixed << setw(10) << setprecision(1) << ecc << "\t\t";
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    string ecc_str = sse.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";
    kldivs_file = current_dir + "kldivs";
    kldivs = load_table(kldivs_file,NUM_METHODS);
    Vector avg_kldivs = computeMeans(kldivs);
    for (int i=0; i<NUM_METHODS; i++) {
      out << fixed << scientific << setprecision(6) << avg_kldivs[i] << "\t\t";
    }
    out << endl;
    ecc += 0.1;
  } // while()
  out.close();
}

void compute_kldivs_diff(
  std::vector<Vector> &kldivs, 
  Vector &kldivs_wins, 
  Vector &kldivs_losses
) {
  double diff;
  kldivs_wins = Vector(NUM_METHODS,0);
  kldivs_losses = Vector(NUM_METHODS,0);
  for (int i=0; i<kldivs.size(); i++) {
    for (int j=0; j<NUM_METHODS; j++) {
      if (j != MML) {
        diff = kldivs[i][j] - kldivs[i][MML];
        if (diff < 0) {  // MML loses ...
          kldivs_wins[j] -= diff;
        } else {  // MML wins ...
          kldivs_losses[j] -= diff;
        } // if ()
      } // j ()
    }
  } // i ()
}

double compute_kldivs_diff_fixed_kappa(string &n_str, string &kappa_str, string &kldivs_folder)
{
  double max = 0;
  Vector kldivs_wins,kldivs_losses;
  string output = kldivs_folder + "kldivs_diff.dat";
  ofstream out(output.c_str());
  double ecc = 0.1;
  while (ecc < 0.95) {
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    string ecc_str = sse.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";
    string kldivs_file = current_dir + "kldivs";
    std::vector<Vector> kldivs_table = load_table(kldivs_file,NUM_METHODS);
    compute_kldivs_diff(kldivs_table,kldivs_wins,kldivs_losses);
    out << fixed << setw(10) << setprecision(1) << ecc << "\t\t";
    for (int i=0; i<NUM_METHODS; i++) {
      if (i != MML) {
        out << fixed << scientific << setprecision(6) << kldivs_wins[i] << "\t\t";
        if (kldivs_wins[i] > max) max = kldivs_wins[i];
      }
    }
    for (int i=0; i<NUM_METHODS; i++) {
      if (i != MML) {
        out << fixed << scientific << setprecision(6) << kldivs_losses[i] << "\t\t";
        if (fabs(kldivs_losses[i]) > max) max = fabs(kldivs_losses[i]);
      }
    }
    out << endl;
    ecc += 0.1;
  } // while()
  out.close();
  return max;
}

void plot_kldivs_diff(string &kldivs_folder, double max)
{
  double ytics;
  if (max >= 1) ytics = 1;
  else if (max > 0.5)  ytics = 0.5;
  else if (max > 0.1) ytics = 0.1;
  else ytics = 0.05;
  string kldivs_diff_file = kldivs_folder + "kldivs_diff.dat";
  string script_file = kldivs_folder + "kldivs_diff.p";
  string plot_file = kldivs_folder + common_kappa + "kldivs_diff.eps";
  ofstream out(script_file.c_str());
  out << "unset xlabel\n";
  out << "set terminal postscript eps enhanced color\n\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "set ylabel \"Difference in KL-divergence\" offset 0,-7,0\n";
  out << "set multiplot layout 2,1\n";
  out << "set lmargin 7\n";
  out << "set rmargin 5\n";
  out << "set ytics " << ytics << "\n";
  out << "unset xtics\n";
  out << "set bmargin 0\n";
  out << "set grid y\n";
  out << "set style data histograms\n";
  out << "set style fill solid 1.0 noborder\n";
  out << "set xlabel font \"Times-Roman, 25\"\n";
  out << "set ylabel font \"Times-Roman, 25\"\n";
  out << "set xtics font \"Times-Roman, 20\"\n";
  out << "set ytics font \"Times-Roman, 20\"\n";
  out << "set yr [0:" << max+0.1 << "]\n";
  out << "plot \"" << kldivs_diff_file << "\" using 2 t \"MOMENT\" lc rgb \"red\", \\\n"
      << "\"\" using 3 t \"MLE\" lc rgb \"blue\", \\\n"
      << "\"\" using 4 t \"MAP1\" lc rgb \"dark-green\", \\\n";
  if (NUM_METHODS == 6) {
    out << "\"\" using 5 t \"MAP2\" lc rgb \"pink\", \\\n"
        << "\"\" using 6 t \"MAP3\" lc rgb \"purple\"\n";
  } else if (NUM_METHODS == 5) {
    out << "\"\" using 5 t \"MAP2\" lc rgb \"pink\"\n";
  }
  out << "unset ylabel\n";
  out << "set xlabel \"eccentricity\\n\"\n";
  out << "set tmargin 0\n";
  out << "set bmargin at screen 0.1\n";
  out << "set xtics nomirror\n";
  out << "set yr [" << -(max+0.1) << ":0]\n";
  int resume_index = NUM_METHODS;
  out << "plot \"" << kldivs_diff_file << "\" using " << resume_index+1 << "notitle lc rgb \"red\", \\\n"
      << "\"\" using " << resume_index+2 << " notitle lc rgb \"blue\", \\\n"
      << "\"\" using " << resume_index+3 << " notitle lc rgb \"dark-green\", \\\n";
  if (NUM_METHODS == 6) {
    out << "\"\" using " << resume_index+4 << " notitle lc rgb \"pink\", \\\n"
        << "\"\" using " << resume_index+5 << ":xtic(1) notitle lc rgb \"purple\", \\\n";
  } else if (NUM_METHODS == 5) {
    out << "\"\" using " << resume_index+4 << ":xtic(1) notitle lc rgb \"pink\", \\\n"; 
  }
  out << "0 lt 1 lw 10 lc rgb \"black\" notitle\n";
  out << "unset multiplot";
  out.close();
  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

void boxplot_test_stats_fixed_kappa(
  string &n_str, string &kappa_str, string &kldivs_folder, int option
) {
  string stats_file,script_file,plot_file,ylabel;

  if (option == 1) {
    script_file = kldivs_folder + "boxplot_test_statistics.p";
    plot_file = kldivs_folder + common_kappa + "boxplot_test_statistics.eps";
    ylabel = "Test statistic";
  } else if (option == 2) {
    script_file = kldivs_folder + "boxplot_pvalues.p";
    plot_file = kldivs_folder + common_kappa + "boxplot_pvalues.eps";
    ylabel = "p-values";
  }

  ofstream out(script_file.c_str());
  out << "set terminal post eps enhanced color\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "box_width=0.12\n";
  out << "set style fill solid 0.25 #noborder\n";
  out << "set style boxplot outliers pointtype 7\n";
  out << "set style data boxplot\n";
  out << "set boxwidth box_width #relative\n";
  out << "set pointsize 0.5\n";
  if (option == 1) {
    out << "set key top right opaque\n";
  } else if (option == 2) {
    out << "set key bottom right opaque\n";
  }
  out << "set border 2\n";
  out << "set xtics nomirror\n";
  out << "set ytics nomirror\n";
  out << "set xlabel \"eccentricity\\n\"\n";
  out << "set ylabel \"" << ylabel << "\"\n";
  out << "set xlabel font \"Times-Roman, 25\"\n";
  out << "set ylabel font \"Times-Roman, 25\"\n";
  out << "set xtics font \"Times-Roman, 20\"\n";
  out << "set ytics font \"Times-Roman, 20\"\n";
  out << "set xtics (\"0.1\" 1, \"0.2\" 2, \"0.3\" 3, \"0.4\" 4, "
      << "\"0.5\" 5, \"0.6\" 6, \"0.7\" 7, \"0.8\" 8, \"0.9\" 9) scale 0.0\n";
  out << "d_width=0.5*box_width\n\n";

  int col = 1;
  for (double ecc=0.1; ecc<0.95; ecc+=0.1) {
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    string ecc_str = sse.str();
    if (option == 1) {
      stats_file = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/chisq_stat";
    } else if (option == 2) {
      stats_file = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/pvalues";
    }
    if (col == 1) {
      out << "plot ";
    }
    if (col < 9) {
      out << "\"" << stats_file << "\" using ((" << col << ")-3*d_width):1 notitle lt 1 lc rgb \"red\", \\\n"
          << "\"" << stats_file << "\" using ((" << col << ")-d_width):3 notitle lt 1 lc rgb \"dark-green\", \\\n"
          << "\"" << stats_file << "\" using ((" << col << ")+d_width):4 notitle lt 1 lc rgb \"dark-magenta\", \\\n"
          << "\"" << stats_file << "\" using ((" << col << ")+3*d_width):5 notitle lt 1 lc rgb \"black\", \\\n";
    } else if (col == 9) {
      out << "\"" << stats_file << "\" using ((" << col << ")-3*d_width):1 title \"MOMENT\" lt 1 lc rgb \"red\", \\\n"
          << "\"" << stats_file << "\" using ((" << col << ")-d_width):3 title \"MAP1\" lt 1 lc rgb \"dark-green\", \\\n"
          << "\"" << stats_file << "\" using ((" << col << ")+d_width):4 title \"MML\" lt 1 lc rgb \"dark-magenta\", \\\n"
          << "\"" << stats_file << "\" using ((" << col << ")+3*d_width):5 title \"MAP2\" lt 1 lc rgb \"black\"\n";
    }
    col++;
  }
  out.close();

  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

void process_kldivs(string &n_str)
{
  string kldivs_folder;
  string wins_kldivs_file;

  int num_maps,map_index;
  if (NUM_METHODS == 5) {
    num_maps = 2;
  } else {
    num_maps = 3;
  }

  string all_ecc = n_str + "fixed_kappa/";
  string kappa_str,kappa_folder;
  double kappa = INIT_KAPPA;
  while (kappa < MAX_KAPPA + 1) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    kappa_str = ssk.str();
    kappa_folder = all_ecc + "kappa_" + kappa_str + "/";
    kldivs_folder = kappa_folder + "kldivs/";

    common_kappa = common + "_k" + kappa_str + "_";
 
    // kldivs wins
    for (int i=1; i<=num_maps; i++) {
      if (i == 1) map_index = 2;
      else if (i == 2) map_index = 4;
      else if (i == 3) map_index = 5;
      wins_kldivs_file = kldivs_folder + "wins_kldivs_map"
                         + boost::lexical_cast<string>(i) + ".dat";
      tabulate_kappa_kldivs(map_index,wins_kldivs_file,n_str,kappa_str);
      plot_script_kldivs_wins(kldivs_folder,i);
    }

    // boxplot kldivs
    boxplot_kldivs_fixed_kappa(n_str,kappa_str,kldivs_folder);

    // average kldivs
    compute_average_kldivs_fixed_kappa(n_str,kappa_str,kldivs_folder);
    plot_avg_kldivs(kldivs_folder);

    // differences
    double max = compute_kldivs_diff_fixed_kappa(n_str,kappa_str,kldivs_folder);
    plot_kldivs_diff(kldivs_folder,max);

    // boxplot hypothesis testing
    boxplot_test_stats_fixed_kappa(n_str,kappa_str,kldivs_folder,1);
    boxplot_test_stats_fixed_kappa(n_str,kappa_str,kldivs_folder,2);

    kappa *= KAPPA_INCREMENT;
  } // while(kappa)
}
/* E **********************************************************************************/ 

/* S **********************************************************************************/ 
void check_and_create_directory(string &directory)
{
  if (stat(directory.c_str(), &st) == -1) {
    mkdir(directory.c_str(), 0700);
  }
}

void create_required_folders(string &n_str)
{
  string kldivs_folder;

  string all_ecc = n_str + "fixed_kappa/";
  check_and_create_directory(all_ecc);
  string kappa_str,kappa_folder;
  double kappa = INIT_KAPPA;
  while (kappa < MAX_KAPPA + 1) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    kappa_str = ssk.str();
    kappa_folder = all_ecc + "kappa_" + kappa_str + "/";
    check_and_create_directory(kappa_folder);

    kldivs_folder = kappa_folder + "kldivs/";
    check_and_create_directory(kldivs_folder);

    //kappa += 10;
    kappa *= KAPPA_INCREMENT;
  } // while(kappa)
}
/* E **********************************************************************************/

struct Parameters
{
  int N;
  int prior;
};

struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;

  options_description desc("Allowed options");
  desc.add_options()
       ("n",value<int>(&parameters.N),"sample size")
       ("prior",value<int>(&parameters.prior),"vMF Kappa (2D/3D) prior")
  ;
  variables_map vm;
  store(command_line_parser(argc,argv).options(desc).run(),vm);
  notify(vm);

  return parameters;
}

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  string n_str = "./estimates/N_" + boost::lexical_cast<string>(parameters.N) 
                + "_prior" + boost::lexical_cast<string>(parameters.prior) + "/";

  common = "n" + boost::lexical_cast<string>(parameters.N) 
           + "_p" + boost::lexical_cast<string>(parameters.prior);

  if (parameters.prior == 2) {
    NUM_METHODS = 6;
  } else if (parameters.prior == 3) {
    NUM_METHODS = 5;
  }

  create_required_folders(n_str);

  process_kldivs(n_str);
}

