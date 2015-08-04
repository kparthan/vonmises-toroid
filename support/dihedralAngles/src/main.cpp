#include<iostream>
#include<fstream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cassert>
#include<string>
#include<stdint.h>
#include<map>
#include<vector>
#include<sstream>
#include<time.h>
#include<algorithm>

#include "pdbClass.h"

using namespace std;

vector<string> parseList(char *fname) {
   ifstream infile(fname, ios::in);
   assert (infile);
   
   vector<string> list;
   char buff[10001];
   while (!infile.eof()) {
      infile.getline(buff,10000);
      if (infile.eof() == true) break;
      list.push_back(buff);
      //cout << list[list.size()-1] << endl;
   }
   infile.close();
   return list;
}

int main( int argc , char *argv[] ) {
   if (argc !=2 ) {
      cerr << "Usage: phipsiomega <file containing a list of PDB files>\n" ;
      exit(1);
   }
   vector<string> pdblist = parseList(argv[1]);

   for (size_t i = 0; i < pdblist.size(); i++) {
      PDB_t pdbobj(pdblist[i].c_str());
      size_t modelidx = 0;
      size_t nChains = pdbobj.getnChains(modelidx);
      for (size_t chainidx = 0; chainidx < nChains; chainidx++) {
         size_t nResidues = pdbobj.getnResidues(modelidx,chainidx);
         for (size_t residx = 1; residx < nResidues-1; residx++) {
            double *phipsiomega = pdbobj.getPhiPsiOmega(modelidx,chainidx,residx);
            /*cout << fixed 
               << setw(10) << setprecision(3)
               << phipsiomega[0] 
               << setw(10) << setprecision(3)
               << phipsiomega[1] 
               << setw(10) << setprecision(3)
               << phipsiomega[2] << endl;*/
            double phi,psi;
            if (phipsiomega[0] < 0) phi = phipsiomega[0] + 360;
            else phi = phipsiomega[0];
            if (phipsiomega[1] < 0) psi = phipsiomega[1] + 360;
            else psi = phipsiomega[1];
            phi *= M_PI/180; psi *= M_PI/180;
            cout << fixed 
               << setw(10) << setprecision(3) << phi
               << setw(10) << setprecision(3) << psi << endl;
         }
      }
   }
}
