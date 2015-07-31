#include "Structure.h"
#include "Support.h"

/*!
 *  \brief Null constructor module.
 */
Structure::Structure()
{}

/*!
 *  \brief This is a constructor function used to instantiate the object 
 *  \param angle_pairs a reference to a std::vector<Vector> 
 *  \param name a reference to a string
 */
Structure::Structure(std::vector<Vector> &angle_pairs, string &name) :
                     angle_pairs(angle_pairs), name(name)
{}

/*!
 *  \brief Loads the profile from an existing file.
 *  \param identifier a reference to a string
 */
void Structure::load(string &file_name)
{
  name = extractName(file_name);
  read_profile(file_name);
}

/*!
 *  \brief Loads the profile from an existing file.
 *  \param path_to_file a reference to path object 
 */
void Structure::load(path &path_to_file)
{
  string file_name = path_to_file.string();
  name = extractName(file_name);
  read_profile(file_name);
}

/*!
 *  \brief This functions reads the profile from a regular file.
 *  (each line contains the unit coordinate list)
 *  \param file_name a reference to a string
 */
void Structure::read_profile(string &file_name)
{
  //cout << "Reading " << file_name << " ..." << endl;
  ifstream profile(file_name.c_str());
  string line;
  angle_pairs.clear();

  while(getline(profile,line)) {
    boost::char_separator<char> sep(",() \t");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    Vector values; 
    BOOST_FOREACH (const string& t, tokens) {
      istringstream iss(t);
      double x;
      iss >> x;
      values.push_back(x);
    }
    angle_pairs.push_back(values);
  }
  profile.close();
}

/*!
 *  \brief This function returns the list of unit coordinates
 */
std::vector<Vector> Structure::getAnglePairs()
{
  return angle_pairs;
}

