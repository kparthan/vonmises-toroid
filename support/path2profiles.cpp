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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/random.hpp>

using namespace std;
using namespace boost::program_options;
using namespace boost::filesystem;
using namespace boost::math;

typedef std::vector<double> Vector;

int main(int argc, char **argv)
{
  string file_name = "../data/class_b_profiles.txt";
  ifstream infile(file_name.c_str());
  file_name = "path2profiles.txt";
  ofstream outfile(file_name.c_str());

  string line;
  std::vector<string> domains;
  string scop_path = "/home/parthan/Research/SCOP/pdbstyle-1.75B/";
  while(getline(infile,line)) {
    string scop_folder = line.substr(2,2);
    string file_path = scop_path + scop_folder + "/" + line;
    //outfile << line << "\t" << line.size() << "\t" << scop_folder << endl;
    outfile << file_path << endl;
  }
  outfile.close();
  infile.close();
}

