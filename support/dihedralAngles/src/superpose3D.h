/*******************************************************************************
 * Copyright (c) 2014,  Arun S Konagurthu, Monash University
 *                                       <http://www.csse.monash.edu.au/~karun>
 * All rights reserved. 
 *
 * This code is released under the terms and conditions of the GNU GENERAL
 * PUBLIC LICENSE V.3. Refer file COPYING in the main directory.
 ******************************************************************************/
#ifndef SUPERPOSE3D_H
#define SUPERPOSE3D_H
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <cassert>
#include <limits>
using namespace std;


// Sufficient statistics class SPECIFIC to least squares 3D superposition
class SuffStatClass {
  public: 
  double sum_xmxm, sum_xmym, sum_xmzm;
  double sum_ymym, sum_ymzm;
  double sum_zmzm;
  double sum_xpxp, sum_xpyp, sum_xpzp;
  double sum_ypyp, sum_ypzp;
  double sum_zpzp;
  double sum_xmyp, sum_xmzp;
  double sum_ymxp, sum_ymzp;
  double sum_zmxp, sum_zmyp;

  // these two store the suff statistics to compute the centers of mass
  // of the two sets
  double sum_Ax, sum_Ay, sum_Az;
  double sum_Bx, sum_By, sum_Bz;

  size_t nVecs;
  bool setflag;
  SuffStatClass();
  void printSufficientStatistics();
};

/* A superposition class to find the optimal least-square 
 * transformation of two corresponding vectors sets.
 * Also allows updates using precomputed sufficient statistics */
class Superpose3DClass {
 private:
   SuffStatClass stats;   // specialized sufficient Statistics class
   double rotacenterA[3]; // movingset rotational center
   double rotacenterB[3]; // fixedset  rotational center

   size_t spA; // start point index in vecSetA (supplied moving coords set)
   size_t spB; // start point index in vecSetB (supplied fixed  coords set)

   /* These are used when updating superposition using sufficient statistics */
   SuffStatClass stats1, stats2;
   double cm1A[3], cm1B[3]; // First part centers of mass (two sets)
   double cm2A[3], cm2B[3]; // Second part centers of mass (two sets)
   double mu1_xm, mu1_ym, mu1_zm;
   double mu1_xp, mu1_yp, mu1_zp;
   double mu2_xm, mu2_ym, mu2_zm;
   double mu2_xp, mu2_yp, mu2_zp;
   SuffStatClass updated_stats;

   double rmsd;
   double eulerRotMat[3][3];

   double quat[4][4]; // A 4x4 quaternion matrices
   double eigenValues[4];
   double eigenVectors[4][4];

   double INFINITYVAL;
   size_t nVecs;
   bool success_flag;

   void computeRotationalCenter(vector<vector<double> >&, vector<vector<double> >&);
   void assignRotationalCenter(vector<double>, vector<double>);
   
   void updateCenterOfMasses_addition();
   void updateCenterOfMasses_subtraction();

   void computeQuaternionMatrix(vector<vector<double> >&, vector<vector<double> >&);
   void computeSufficientStatistics(vector<vector<double> >&, vector<vector<double> >&);
   void computeQuaternionMatrixUsingSufficientStats();
   void updateSufficientStatistics_addition();
   void updateSufficientStatistics_subtraction();

   void eigsrt();
   bool diagonalize();
   void computeRotationMatrix();
   bool computeLSqFit(vector<vector<double> >&, vector<vector<double> >&);
   bool computeLSqFit(vector<vector<double> >&, vector<vector<double> >&,vector<double>&rcA,vector<double>&rcB);
   bool updateLSqFit_addition();
   bool updateLSqFit_subtraction();

   //debug
   void printRotationalCenters();
   void printEigens();
   void printQuaternionMatrix();
   void printRotationMatrix();

 public:
   Superpose3DClass();
   Superpose3DClass(vector<vector<double> >&, vector<vector<double> >&);
   Superpose3DClass(vector<vector<double> >&, vector<vector<double> >&, size_t, size_t, size_t);
   Superpose3DClass(SuffStatClass&, SuffStatClass&,char );
   Superpose3DClass(vector<double>&, vector<double>&);
   Superpose3DClass(SuffStatClass&, vector<double>&, vector<double>&, char);
   double getRMSD();	
   vector<double> getDeviations(size_t);
   void transformVector( vector<double> & );
   void transformVectors( vector<vector<double> > & );
   void copyRotationMatrixInto(double [][3]);
   void copyRotationalCentersInto(double [3], double[3]);
   SuffStatClass getSufficientStatistics();
};
#endif
