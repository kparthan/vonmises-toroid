/*******************************************************************************
 * Copyright (c) 2014,  Arun S Konagurthu, Monash University
 *                                       <http://www.csse.monash.edu.au/~karun>
 * All rights reserved. 
 *
 * This code is released under the terms and conditions of the GNU GENERAL
 * PUBLIC LICENSE V.3. Refer file COPYING in the main directory.
 ******************************************************************************/
#include "superpose3D.h"

#define UPDATE(k,l,n) ((n)*(k)*(l))
#define DELTA_PLUS(x,y)  ((vecSetA[x+spA][y]-rotacenterA[y])+(vecSetB[x+spB][y]-rotacenterB[y]))
#define DELTA_MINUS(x,y) ((vecSetA[x+spA][y]-rotacenterA[y])-(vecSetB[x+spB][y]-rotacenterB[y]))
#define ROTATE(quat,i,j,k,l)  g=quat[i][j];h=quat[k][l];\
            quat[i][j]=g-s*(h+g*tau);\
            quat[k][l]=h+s*(g-h*tau);

#define SQR(x) ((x)*(x))

/* constructor of sufficient statistics class */
SuffStatClass::SuffStatClass(){
  setflag = false;
  sum_xmxm= sum_xmym= sum_xmzm=
  sum_ymym= sum_ymzm= sum_zmzm= sum_xpxp= sum_xpyp= sum_xpzp= sum_ypyp= 
  sum_ypzp= sum_zpzp= sum_xmyp= sum_xmzp= sum_ymxp= sum_ymzp=
  sum_zmxp= sum_zmyp= sum_Ax= sum_Ay= sum_Az= sum_Bx= sum_By= sum_Bz=
  numeric_limits<double>::infinity();
  nVecs = ~(0LU);
}

/* Print sufficient statistics -- used for debugging */
void SuffStatClass::printSufficientStatistics() {
   cout << "********************SUFFICIENT STATS***********************\n";
   cout << "sum_xmxm = " << fixed << setw(8) << setprecision(1) << sum_xmxm << endl;
   cout << "sum_ymym = " << fixed << setw(8) << setprecision(1) << sum_ymym << endl;
   cout << "sum_zmzm = " << fixed << setw(8) << setprecision(1) << sum_zmzm << endl;

   cout << "sum_xpxp = " << fixed << setw(8) << setprecision(1) << sum_xpxp << endl;
   cout << "sum_ypyp = " << fixed << setw(8) << setprecision(1) << sum_ypyp << endl;
   cout << "sum_zpzp = " << fixed << setw(8) << setprecision(1) << sum_zpzp << endl;

   cout << "sum_xmyp = " << fixed << setw(8) << setprecision(1) << sum_xmyp << endl;
   cout << "sum_xmzp = " << fixed << setw(8) << setprecision(1) << sum_xmzp << endl;
   cout << "sum_ymxp = " << fixed << setw(8) << setprecision(1) << sum_ymxp << endl;
   cout << "sum_ymzp = " << fixed << setw(8) << setprecision(1) << sum_ymzp << endl;
   cout << "sum_zmxp = " << fixed << setw(8) << setprecision(1) << sum_zmxp << endl;
   cout << "sum_zmyp = " << fixed << setw(8) << setprecision(1) << sum_zmyp << endl;
   cout << "****************************END***********************\n";
}


/* Compute the rotational center of the two fragments. 
   HERE IT IS THE BARRY CENTERS OF THE TWO FRAGMENTS.  */
void Superpose3DClass::computeRotationalCenter(
  vector<vector<double> > &vecSetA,
  vector<vector<double> > &vecSetB
) {
   double temp = 0;
   rotacenterA[0] = rotacenterA[1] = rotacenterA[2] = INFINITYVAL;
   rotacenterB[0] = rotacenterB[1] = rotacenterB[2] = INFINITYVAL;

   for( int i = 0; i < 3; i++ ) {
      temp = 0;
      for( size_t j = 0; j < nVecs; j++ ) {
         temp += vecSetA[j+spA][i];
      }
      rotacenterA[i] = temp/nVecs;
   }
   for( size_t i = 0; i < 3; i++ ) {
      temp = 0;
      for( size_t j = 0; j < nVecs; j++ ) {
         temp += vecSetB[j+spB][i];
      }
      rotacenterB[i] = temp/nVecs;
   }
}

void Superpose3DClass::assignRotationalCenter(
  vector<double> rcA,
  vector<double> rcB
) {
   double temp = 0;
   rotacenterA[0] = rcA[0]; rotacenterA[1] = rcA[1]; rotacenterA[2] = rcA[2];
   rotacenterB[0] = rcB[0]; rotacenterB[1] = rcB[1]; rotacenterB[2] = rcB[2];
}


/* Compute the square-symmetric 4x4 Quaternion matrix. */
void Superpose3DClass::computeQuaternionMatrix(
  vector<vector<double> > &vecSetA,
  vector<vector<double> > &vecSetB
) {
   //initialize quaternion matrix
   quat[0][0] = quat[0][1] = quat[0][2] = quat [0][3] =
   quat[1][0] = quat[1][1] = quat[1][2] = quat [1][3] =
   quat[2][0] = quat[2][1] = quat[2][2] = quat [2][3] =
   quat[3][0] = quat[3][1] = quat[3][2] = quat [3][3] = 0.0;
   
   // Filling upper triangle of the quaternion matrix
   for( size_t i = 0; i < nVecs  ; i++ )
   {
     double xm = (DELTA_MINUS(i,0));
     double ym = (DELTA_MINUS(i,1));
     double zm = (DELTA_MINUS(i,2));

     double xp = (DELTA_PLUS(i,0));
     double yp = (DELTA_PLUS(i,1));
     double zp = (DELTA_PLUS(i,2));

      //Diags =Sum of squared cyclic coordinate differences
      quat[0][0] += ((xm*xm)+(ym*ym)+(zm*zm));
      quat[1][1] += ((yp*yp)+(zp*zp)+(xm*xm));
      quat[2][2] += ((xp*xp)+(zp*zp)+(ym*ym));
      quat[3][3] += ((xp*xp)+(yp*yp)+(zm*zm));
      // Cross differences
      quat[0][1] += (yp*zm-ym*zp);
      quat[0][2] += (xm*zp-xp*zm);
      quat[0][3] += (xp*ym-xm*yp);
      quat[1][2] += (xm*ym-xp*yp);
      quat[1][3] += (xm*zm-xp*zp);
      quat[2][3] += (ym*zm-yp*zp);
   }
   // Fill the rest by transposing it onto itself
   quat[1][0] = quat[0][1]; 
   quat[2][0] = quat[0][2]; 
   quat[2][1] = quat[1][2];  
   quat[3][0] = quat[0][3]; 
   quat[3][1] = quat[1][3]; 
   quat[3][2] = quat[2][3];
}


void Superpose3DClass::computeSufficientStatistics(
vector<vector<double> > &vecSetA,
vector<vector<double> > &vecSetB
) {
  stats.sum_xmxm  = stats.sum_ymym  = stats.sum_zmzm  =
  stats.sum_xpxp  = stats.sum_ypyp  = stats.sum_zpzp  =
  stats.sum_xmyp  = stats.sum_xmzp  = stats.sum_ymxp  =
  stats.sum_ymzp  = stats.sum_zmxp  = stats.sum_zmyp  = 
  stats.sum_xmym  = stats.sum_xmzm  = stats.sum_ymzm  =
  stats.sum_xpyp  = stats.sum_xpzp  = stats.sum_ypzp  = 0;

  for( size_t i = 0; i < nVecs; i++ ) {
  // sum_xm, sum_ym, sum_zm;
		double t_sum_xm = DELTA_MINUS(i,0);
		double t_sum_ym = DELTA_MINUS(i,1);
		double t_sum_zm = DELTA_MINUS(i,2);
		double t_sum_xp = DELTA_PLUS(i,0);
		double t_sum_yp = DELTA_PLUS(i,1);
		double t_sum_zp = DELTA_PLUS(i,2);
  // sum_xmxm, sum_ymym, sum_zmzm;
      stats.sum_xmxm += t_sum_xm*t_sum_xm;
      stats.sum_ymym += t_sum_ym*t_sum_ym;
      stats.sum_zmzm += t_sum_zm*t_sum_zm;
  // sum_xpxp, sum_ypyp, sum_zpzp;
      stats.sum_xpxp += t_sum_xp*t_sum_xp;
      stats.sum_ypyp += t_sum_yp*t_sum_yp;
      stats.sum_zpzp += t_sum_zp*t_sum_zp;
  // sum_xmyp, sum_xmzp;
      stats.sum_xmyp += t_sum_xm*t_sum_yp;
      stats.sum_xmzp += t_sum_xm*t_sum_zp;
  // sum_ymxp, sum_ymzp;
      stats.sum_ymxp += t_sum_ym*t_sum_xp;
      stats.sum_ymzp += t_sum_ym*t_sum_zp;
  // sum_zmxp, sum_zmyp;
      stats.sum_zmxp += t_sum_zm*t_sum_xp;
      stats.sum_zmyp += t_sum_zm*t_sum_yp;
  //sum_xmym, sum_xmzm, sum_ymzm;    
      stats.sum_xmym += t_sum_xm*t_sum_ym;
      stats.sum_xmzm += t_sum_xm*t_sum_zm;
      stats.sum_ymzm += t_sum_ym*t_sum_zm;
  //sum_xpyp, sum_xpzp, sum_ypzp;    
      stats.sum_xpyp += t_sum_xp*t_sum_yp;
      stats.sum_xpzp += t_sum_xp*t_sum_zp;
      stats.sum_ypzp += t_sum_yp*t_sum_zp;
	}

  stats.sum_Ax = rotacenterA[0]*nVecs;
  stats.sum_Ay = rotacenterA[1]*nVecs;
  stats.sum_Az = rotacenterA[2]*nVecs;

  stats.sum_Bx = rotacenterB[0]*nVecs;
  stats.sum_By = rotacenterB[1]*nVecs;
  stats.sum_Bz = rotacenterB[2]*nVecs;

  stats.nVecs = nVecs;
  stats.setflag = true;
}

/* Compute the Barry centers of the two fragments. */
void Superpose3DClass::updateCenterOfMasses_addition() {
  /*
  cout << "nPoints1 = " << stats1.nVecs << endl;
  cout << "nPoints2 = " << stats2.nVecs << endl;
  cout << "nPoints = " << nVecs << endl;
  */
  cm1A[0] = stats1.sum_Ax/stats1.nVecs;
  cm1A[1] = stats1.sum_Ay/stats1.nVecs;
  cm1A[2] = stats1.sum_Az/stats1.nVecs;

  cm1B[0] = stats1.sum_Bx/stats1.nVecs;
  cm1B[1] = stats1.sum_By/stats1.nVecs;
  cm1B[2] = stats1.sum_Bz/stats1.nVecs;

  cm2A[0] = stats2.sum_Ax/stats2.nVecs;
  cm2A[1] = stats2.sum_Ay/stats2.nVecs;
  cm2A[2] = stats2.sum_Az/stats2.nVecs;

  cm2B[0] = stats2.sum_Bx/stats2.nVecs;
  cm2B[1] = stats2.sum_By/stats2.nVecs;
  cm2B[2] = stats2.sum_Bz/stats2.nVecs;

  rotacenterA[0] = (stats1.sum_Ax+stats2.sum_Ax)/nVecs;
  rotacenterA[1] = (stats1.sum_Ay+stats2.sum_Ay)/nVecs;
  rotacenterA[2] = (stats1.sum_Az+stats2.sum_Az)/nVecs;

  rotacenterB[0] = (stats1.sum_Bx+stats2.sum_Bx)/nVecs;
  rotacenterB[1] = (stats1.sum_By+stats2.sum_By)/nVecs;
  rotacenterB[2] = (stats1.sum_Bz+stats2.sum_Bz)/nVecs;
  /*
  cout << "cm1A "<< cm1A[0] << " "<< cm1A[1]<< " "<< cm1A[2] << endl;
  cout << "cm1B "<< cm1B[0] << " "<< cm1B[1]<< " "<< cm1B[2] << endl;
  cout << "cm2A "<< cm2A[0] << " "<< cm2A[1]<< " "<< cm2A[2] << endl;
  cout << "cm2B "<< cm2B[0] << " "<< cm2B[1]<< " "<< cm2B[2] << endl;
  cout << "rotacenterA   "<< rotacenterA[0] << " "<< rotacenterA[1]<< " "<< rotacenterA[2] << endl;
  cout << "rotacenterB   "<< rotacenterB[0] << " "<< rotacenterB[1]<< " "<< rotacenterB[2] << endl;
  */
}

/* Compute the Barry centers of the two fragments. */
void Superpose3DClass::updateCenterOfMasses_subtraction() {
  /*
  cout << "nPoints1 = " << stats1.nVecs << endl;
  cout << "nPoints2 = " << stats2.nVecs << endl;
  cout << "nPoints = " << nVecs << endl;
  //*/
  cm1A[0] = stats1.sum_Ax/stats1.nVecs;
  cm1A[1] = stats1.sum_Ay/stats1.nVecs;
  cm1A[2] = stats1.sum_Az/stats1.nVecs;

  cm1B[0] = stats1.sum_Bx/stats1.nVecs;
  cm1B[1] = stats1.sum_By/stats1.nVecs;
  cm1B[2] = stats1.sum_Bz/stats1.nVecs;

  cm2A[0] = stats2.sum_Ax/stats2.nVecs;
  cm2A[1] = stats2.sum_Ay/stats2.nVecs;
  cm2A[2] = stats2.sum_Az/stats2.nVecs;

  cm2B[0] = stats2.sum_Bx/stats2.nVecs;
  cm2B[1] = stats2.sum_By/stats2.nVecs;
  cm2B[2] = stats2.sum_Bz/stats2.nVecs;

  rotacenterA[0] = (stats1.sum_Ax-stats2.sum_Ax)/nVecs;
  rotacenterA[1] = (stats1.sum_Ay-stats2.sum_Ay)/nVecs;
  rotacenterA[2] = (stats1.sum_Az-stats2.sum_Az)/nVecs;

  rotacenterB[0] = (stats1.sum_Bx-stats2.sum_Bx)/nVecs;
  rotacenterB[1] = (stats1.sum_By-stats2.sum_By)/nVecs;
  rotacenterB[2] = (stats1.sum_Bz-stats2.sum_Bz)/nVecs;
  /*
  cout << "cm1A "<< cm1A[0] << " "<< cm1A[1]<< " "<< cm1A[2] << endl;
  cout << "cm1B "<< cm1B[0] << " "<< cm1B[1]<< " "<< cm1B[2] << endl;
  cout << "cm2A "<< cm2A[0] << " "<< cm2A[1]<< " "<< cm2A[2] << endl;
  cout << "cm2B "<< cm2B[0] << " "<< cm2B[1]<< " "<< cm2B[2] << endl;
  cout << "rotacenterA   "<< rotacenterA[0] << " "<< rotacenterA[1]<< " "<< rotacenterA[2] << endl;
  cout << "rotacenterB   "<< rotacenterB[0] << " "<< rotacenterB[1]<< " "<< rotacenterB[2] << endl;
  */
}



/* Compute the square-symmetric 4x4 Quaternion matrix, using sufficient stats.*/
void Superpose3DClass::computeQuaternionMatrixUsingSufficientStats() {
	//initialize quaternion matrix
  for (size_t i = 0; i < 4; i++) {
    for (size_t j = 0; j < 4; j++) {
      quat[i][j] = 0;
    }
  }
	
	// Filling upper triangle of the quaternion matrix
  quat[0][0] = stats.sum_xmxm + stats.sum_ymym + stats.sum_zmzm;
  quat[1][1] = stats.sum_xmxm + stats.sum_ypyp + stats.sum_zpzp;
  quat[2][2] = stats.sum_xpxp + stats.sum_ymym + stats.sum_zpzp;
  quat[3][3] = stats.sum_xpxp + stats.sum_ypyp + stats.sum_zmzm;

  quat[0][1] = stats.sum_zmyp - stats.sum_ymzp;
  quat[0][2] = stats.sum_xmzp - stats.sum_zmxp;
  quat[0][3] = stats.sum_ymxp - stats.sum_xmyp;

  quat[1][2] = stats.sum_xmym - stats.sum_xpyp;
  quat[1][3] = stats.sum_xmzm - stats.sum_xpzp;
  
  quat[2][3] = stats.sum_ymzm - stats.sum_ypzp;
  // Fill the rest by transposing it onto itself
  quat[1][0] = quat[0][1]; quat[2][0] = quat[0][2]; quat[2][1] = quat[1][2];
  quat[3][0] = quat[0][3]; quat[3][1] = quat[1][3]; quat[3][2] = quat[2][3];

  
  //printQuaternionMatrix();
}

void Superpose3DClass::updateSufficientStatistics_addition() {
  stats.sum_xmxm = stats.sum_ymym = stats.sum_zmzm =
  stats.sum_xpxp = stats.sum_ypyp = stats.sum_zpzp =
  stats.sum_xmyp = stats.sum_xmzp = stats.sum_ymxp =
  stats.sum_ymzp = stats.sum_zmxp = stats.sum_zmyp = 
  stats.sum_xmym = stats.sum_xmzm = stats.sum_ymzm =
  stats.sum_xpyp = stats.sum_xpzp = stats.sum_ypzp = 0;

  stats.sum_Ax = stats1.sum_Ax + stats2.sum_Ax;
  stats.sum_Ay = stats1.sum_Ay + stats2.sum_Ay;
  stats.sum_Az = stats1.sum_Az + stats2.sum_Az;

  stats.sum_Bx = stats1.sum_Bx + stats2.sum_Bx;
  stats.sum_By = stats1.sum_By + stats2.sum_By;
  stats.sum_Bz = stats1.sum_Bz + stats2.sum_Bz;

  stats.nVecs = stats1.nVecs + stats2.nVecs;
  double mu1A[3]; //correction to center of mass
  mu1A[0] = cm1A[0] - rotacenterA[0];
  mu1A[1] = cm1A[1] - rotacenterA[1];
  mu1A[2] = cm1A[2] - rotacenterA[2];
  double mu1B[3]; //correction to center of mass
  mu1B[0] = cm1B[0] - rotacenterB[0];
  mu1B[1] = cm1B[1] - rotacenterB[1];
  mu1B[2] = cm1B[2] - rotacenterB[2];

  mu1_xm = mu1A[0] - mu1B[0];
  mu1_ym = mu1A[1] - mu1B[1];
  mu1_zm = mu1A[2] - mu1B[2];

  mu1_xp = mu1A[0] + mu1B[0];
  mu1_yp = mu1A[1] + mu1B[1];
  mu1_zp = mu1A[2] + mu1B[2];

  double mu2A[3]; //correction to center of mass
  mu2A[0] = cm2A[0] - rotacenterA[0];
  mu2A[1] = cm2A[1] - rotacenterA[1];
  mu2A[2] = cm2A[2] - rotacenterA[2];
  double mu2B[3]; //correction to center of mass
  mu2B[0] = cm2B[0] - rotacenterB[0];
  mu2B[1] = cm2B[1] - rotacenterB[1];
  mu2B[2] = cm2B[2] - rotacenterB[2];

  mu2_xm = mu2A[0] - mu2B[0];
  mu2_ym = mu2A[1] - mu2B[1];
  mu2_zm = mu2A[2] - mu2B[2];

  mu2_xp = mu2A[0] + mu2B[0];
  mu2_yp = mu2A[1] + mu2B[1];
  mu2_zp = mu2A[2] + mu2B[2];

  /*
  cout << "######\n";
  cout << cm1A[0] << " "<< cm1A[1]<< " "<< cm1A[2] << endl;
  cout << cm1B[0] << " "<< cm1B[1]<< " "<< cm1B[2] << endl;
  cout << cm2A[0] << " "<< cm2A[1]<< " "<< cm2A[2] << endl;
  cout << cm2B[0] << " "<< cm2B[1]<< " "<< cm2B[2] << endl;
  cout << rotacenterA[0] << " "<< rotacenterA[1]<< " "<< rotacenterA[2] << endl;
  cout << rotacenterB[0] << " "<< rotacenterB[1]<< " "<< rotacenterB[2] << endl;
  cout << mu1_xm << " "<< mu1_ym << " "<< mu1_zm << endl;
  cout << "------\n";
*/

    // sum_xmxm, sum_ymym, sum_zmzm;
    stats.sum_xmxm =  stats1.sum_xmxm + UPDATE(mu1_xm,mu1_xm,stats1.nVecs)
             + stats2.sum_xmxm + UPDATE(mu2_xm,mu2_xm,stats2.nVecs);
    stats.sum_ymym =  stats1.sum_ymym + UPDATE(mu1_ym,mu1_ym,stats1.nVecs)
             + stats2.sum_ymym + UPDATE(mu2_ym,mu2_ym,stats2.nVecs);
    stats.sum_zmzm =  stats1.sum_zmzm + UPDATE(mu1_zm,mu1_zm,stats1.nVecs)
             + stats2.sum_zmzm + UPDATE(mu2_zm,mu2_zm,stats2.nVecs);
    // sum_xpxp, sum_ypyp, sum_zpzp;
    stats.sum_xpxp =  stats1.sum_xpxp + UPDATE(mu1_xp,mu1_xp,stats1.nVecs)
             + stats2.sum_xpxp + UPDATE(mu2_xp,mu2_xp,stats2.nVecs);
    stats.sum_ypyp =  stats1.sum_ypyp + UPDATE(mu1_yp,mu1_yp,stats1.nVecs)
             + stats2.sum_ypyp + UPDATE(mu2_yp,mu2_yp,stats2.nVecs);
    stats.sum_zpzp =  stats1.sum_zpzp + UPDATE(mu1_zp,mu1_zp,stats1.nVecs)
             + stats2.sum_zpzp + UPDATE(mu2_zp,mu2_zp,stats2.nVecs);
    // sum_xmyp, sum_xmzp;
    stats.sum_xmyp =  stats1.sum_xmyp + UPDATE(mu1_xm,mu1_yp,stats1.nVecs)
             + stats2.sum_xmyp + UPDATE(mu2_xm,mu2_yp,stats2.nVecs);
    stats.sum_xmzp =  stats1.sum_xmzp + UPDATE(mu1_xm,mu1_zp,stats1.nVecs)
             + stats2.sum_xmzp + UPDATE(mu2_xm,mu2_zp,stats2.nVecs);
    // sum_ymxp, sum_ymzp;
    stats.sum_ymxp =  stats1.sum_ymxp + UPDATE(mu1_ym,mu1_xp,stats1.nVecs)
             + stats2.sum_ymxp + UPDATE(mu2_ym,mu2_xp,stats2.nVecs);
    stats.sum_ymzp =  stats1.sum_ymzp + UPDATE(mu1_ym,mu1_zp,stats1.nVecs)
             + stats2.sum_ymzp + UPDATE(mu2_ym,mu2_zp,stats2.nVecs);
    // sum_zmxp, sum_zmyp;
    stats.sum_zmxp =  stats1.sum_zmxp + UPDATE(mu1_zm,mu1_xp,stats1.nVecs)
             + stats2.sum_zmxp + UPDATE(mu2_zm,mu2_xp,stats2.nVecs);
    stats.sum_zmyp =  stats1.sum_zmyp + UPDATE(mu1_zm,mu1_yp,stats1.nVecs)
             + stats2.sum_zmyp + UPDATE(mu2_zm,mu2_yp,stats2.nVecs);
    //sum_xmym, sum_xmzm, sum_ymzm;    
    stats.sum_xmym =  stats1.sum_xmym + UPDATE(mu1_xm,mu1_ym,stats1.nVecs)
             + stats2.sum_xmym + UPDATE(mu2_xm,mu2_ym,stats2.nVecs);
    stats.sum_xmzm =  stats1.sum_xmzm + UPDATE(mu1_xm,mu1_zm,stats1.nVecs)
             + stats2.sum_xmzm + UPDATE(mu2_xm,mu2_zm,stats2.nVecs);
    stats.sum_ymzm =  stats1.sum_ymzm + UPDATE(mu1_ym,mu1_zm,stats1.nVecs)
             + stats2.sum_ymzm + UPDATE(mu2_ym,mu2_zm,stats2.nVecs);
    //sum_xpyp, sum_xpzp, sum_ypzp;    
    stats.sum_xpyp =  stats1.sum_xpyp + UPDATE(mu1_xp,mu1_yp,stats1.nVecs)
             + stats2.sum_xpyp + UPDATE(mu2_xp,mu2_yp,stats2.nVecs);
    stats.sum_xpzp =  stats1.sum_xpzp + UPDATE(mu1_xp,mu1_zp,stats1.nVecs)
             + stats2.sum_xpzp + UPDATE(mu2_xp,mu2_zp,stats2.nVecs);
    stats.sum_ypzp =  stats1.sum_ypzp + UPDATE(mu1_yp,mu1_zp,stats1.nVecs)
             + stats2.sum_ypzp + UPDATE(mu2_yp,mu2_zp,stats2.nVecs);

  stats.setflag = true;
}

void Superpose3DClass::updateSufficientStatistics_subtraction() {
  stats.sum_xmxm = stats.sum_ymym = stats.sum_zmzm =
  stats.sum_xpxp = stats.sum_ypyp = stats.sum_zpzp =
  stats.sum_xmyp = stats.sum_xmzp = stats.sum_ymxp =
  stats.sum_ymzp = stats.sum_zmxp = stats.sum_zmyp = 
  stats.sum_xmym = stats.sum_xmzm = stats.sum_ymzm =
  stats.sum_xpyp = stats.sum_xpzp = stats.sum_ypzp = 0;

  stats.sum_Ax = stats1.sum_Ax - stats2.sum_Ax;
  stats.sum_Ay = stats1.sum_Ay - stats2.sum_Ay;
  stats.sum_Az = stats1.sum_Az - stats2.sum_Az;

  stats.sum_Bx = stats1.sum_Bx - stats2.sum_Bx;
  stats.sum_By = stats1.sum_By - stats2.sum_By;
  stats.sum_Bz = stats1.sum_Bz - stats2.sum_Bz;

  stats.nVecs = stats1.nVecs - stats2.nVecs;


  double mu1A[3]; //correction to center of mass
  mu1A[0] = -cm1A[0] + cm2A[0];
  mu1A[1] = -cm1A[1] + cm2A[1];
  mu1A[2] = -cm1A[2] + cm2A[2];
  double mu1B[3]; //correction to center of mass
  mu1B[0] = -cm1B[0] + cm2B[0];
  mu1B[1] = -cm1B[1] + cm2B[1];
  mu1B[2] = -cm1B[2] + cm2B[2];

  mu1_xm = mu1A[0] - mu1B[0];
  mu1_ym = mu1A[1] - mu1B[1];
  mu1_zm = mu1A[2] - mu1B[2];

  mu1_xp = mu1A[0] + mu1B[0];
  mu1_yp = mu1A[1] + mu1B[1];
  mu1_zp = mu1A[2] + mu1B[2];

  double mu2A[3]; //correction to center of mass
  mu2A[0] = rotacenterA[0] - cm1A[0];
  mu2A[1] = rotacenterA[1] - cm1A[1];
  mu2A[2] = rotacenterA[2] - cm1A[2];
  double mu2B[3]; //correction to center of mass
  mu2B[0] = rotacenterB[0] - cm1B[0];
  mu2B[1] = rotacenterB[1] - cm1B[1];
  mu2B[2] = rotacenterB[2] - cm1B[2];

  mu2_xm = mu2A[0] - mu2B[0];
  mu2_ym = mu2A[1] - mu2B[1];
  mu2_zm = mu2A[2] - mu2B[2];

  mu2_xp = mu2A[0] + mu2B[0];
  mu2_yp = mu2A[1] + mu2B[1];
  mu2_zp = mu2A[2] + mu2B[2];

  /*
  cout << "######\n";
  cout << cm1A[0] << " "<< cm1A[1]<< " "<< cm1A[2] << endl;
  cout << cm1B[0] << " "<< cm1B[1]<< " "<< cm1B[2] << endl;
  cout << cm2A[0] << " "<< cm2A[1]<< " "<< cm2A[2] << endl;
  cout << cm2B[0] << " "<< cm2B[1]<< " "<< cm2B[2] << endl;
  cout << rotacenterA[0] << " "<< rotacenterA[1]<< " "<< rotacenterA[2] << endl;
  cout << rotacenterB[0] << " "<< rotacenterB[1]<< " "<< rotacenterB[2] << endl;
  cout << mu1_xm << " "<< mu1_ym << " "<< mu1_zm << endl;
  cout << "------\n";
*/

    // sum_xmxm, sum_ymym, sum_zmzm;
    stats.sum_xmxm =  stats1.sum_xmxm - stats2.sum_xmxm
       - UPDATE(mu1_xm,mu1_xm,stats2.nVecs) - UPDATE(mu2_xm,mu2_xm,stats.nVecs);
    stats.sum_ymym =  stats1.sum_ymym - stats2.sum_ymym
       - UPDATE(mu1_ym,mu1_ym,stats2.nVecs) - UPDATE(mu2_ym,mu2_ym,stats.nVecs);
    stats.sum_zmzm =  stats1.sum_zmzm - stats2.sum_zmzm
       - UPDATE(mu1_zm,mu1_zm,stats2.nVecs) - UPDATE(mu2_zm,mu2_zm,stats.nVecs);
    // sum_xpxp, sum_ypyp, sum_zpzp;
    stats.sum_xpxp =  stats1.sum_xpxp - stats2.sum_xpxp
       - UPDATE(mu1_xp,mu1_xp,stats2.nVecs) - UPDATE(mu2_xp,mu2_xp,stats.nVecs);
    stats.sum_ypyp =  stats1.sum_ypyp - stats2.sum_ypyp
       - UPDATE(mu1_yp,mu1_yp,stats2.nVecs) - UPDATE(mu2_yp,mu2_yp,stats.nVecs);
    stats.sum_zpzp =  stats1.sum_zpzp - stats2.sum_zpzp
       - UPDATE(mu1_zp,mu1_zp,stats2.nVecs) - UPDATE(mu2_zp,mu2_zp,stats.nVecs);

    // sum_xmyp, sum_xmzp;
    stats.sum_xmyp =  stats1.sum_xmyp - stats2.sum_xmyp
       - UPDATE(mu1_xm,mu1_yp,stats2.nVecs) - UPDATE(mu2_xm,mu2_yp,stats.nVecs);
    stats.sum_xmzp =  stats1.sum_xmzp - stats2.sum_xmzp
       - UPDATE(mu1_xm,mu1_zp,stats2.nVecs) - UPDATE(mu2_xm,mu2_zp,stats.nVecs);
    // sum_ymxp, sum_ymzp;
    stats.sum_ymxp =  stats1.sum_ymxp - stats2.sum_ymxp
       - UPDATE(mu1_ym,mu1_xp,stats2.nVecs) - UPDATE(mu2_ym,mu2_xp,stats.nVecs);
    stats.sum_ymzp =  stats1.sum_ymzp - stats2.sum_ymzp
       - UPDATE(mu1_ym,mu1_zp,stats2.nVecs) - UPDATE(mu2_ym,mu2_zp,stats.nVecs);
    // sum_zmxp, sum_zmyp;
    stats.sum_zmxp =  stats1.sum_zmxp - stats2.sum_zmxp
       - UPDATE(mu1_zm,mu1_xp,stats2.nVecs) - UPDATE(mu2_zm,mu2_xp,stats.nVecs);
    stats.sum_zmyp =  stats1.sum_zmyp - stats2.sum_zmyp
       - UPDATE(mu1_zm,mu1_yp,stats2.nVecs) - UPDATE(mu2_zm,mu2_yp,stats.nVecs);
    //sum_xmym, sum_xmzm, sum_ymzm;    
    stats.sum_xmym =  stats1.sum_xmym - stats2.sum_xmym
       - UPDATE(mu1_xm,mu1_ym,stats2.nVecs) - UPDATE(mu2_xm,mu2_ym,stats.nVecs);
    stats.sum_xmzm =  stats1.sum_xmzm - stats2.sum_xmzm
       - UPDATE(mu1_xm,mu1_zm,stats2.nVecs) - UPDATE(mu2_xm,mu2_zm,stats.nVecs);
    stats.sum_ymzm =  stats1.sum_ymzm - stats2.sum_ymzm
       - UPDATE(mu1_ym,mu1_zm,stats2.nVecs) - UPDATE(mu2_ym,mu2_zm,stats.nVecs);
    //sum_xpyp, sum_xpzp, sum_ypzp;    
    stats.sum_xpyp =  stats1.sum_xpyp - stats2.sum_xpyp
       - UPDATE(mu1_xp,mu1_yp,stats2.nVecs) - UPDATE(mu2_xp,mu2_yp,stats.nVecs);
    stats.sum_xpzp =  stats1.sum_xpzp - stats2.sum_xpzp
       - UPDATE(mu1_xp,mu1_zp,stats2.nVecs) - UPDATE(mu2_xp,mu2_zp,stats.nVecs);
    stats.sum_ypzp =  stats1.sum_ypzp - stats2.sum_ypzp
       - UPDATE(mu1_yp,mu1_zp,stats2.nVecs) - UPDATE(mu2_yp,mu2_zp,stats.nVecs);

  stats.setflag = true;
}




/* sort in descending order eigen-values and corresponding -vectors */
void Superpose3DClass::eigsrt()
{
   int i, j, k;
   double p;

   for( i = 0; i < 4; i++ ) {
      p = eigenValues[k=i];
      for( j = i+1; j < 4; j++ ) {
         if( eigenValues[j] >= p ) {
            p = eigenValues[k=j];
         }
      }
      if(k != i) {
         eigenValues[k] = eigenValues[i];
         eigenValues[i] = p;
         for( j = 0; j < 4; j++ ) {
            p = eigenVectors[j][i];
            eigenVectors[j][i] = eigenVectors[j][k];
            eigenVectors[j][k] = p;
         }
      }
   }
}

/* diagonolize the quaternion matrix on success compute the RMSD
 * from the leading eigenvalue. Returns a true on success and false on
 * failure. On failure, the rmsd is set to INFINITYVAL.
   Uses routine from numerical recipies with small variations */
bool Superpose3DClass::diagonalize() {
   int  iq , ip;
   double tresh , theta , tau , t , sm , s , h , g , c;
   double b[4], z[4];
   int NROT;

   /* Initialize eigenVectors to an Identity matrix */
   // set everything to zero.
   for (size_t i = 0; i < 4; i++) eigenValues[i] = 0;
   for (size_t i = 0; i < 4; i++) {
    for (size_t j = 0; j < 4; j++) {
      eigenVectors[i][j] = 0;
    }
   }
   // set diagonal elements to 1.
   eigenVectors[0][0] = eigenVectors[1][1] = 
   eigenVectors[2][2] = eigenVectors[3][3] = 1;
   
   // set b and eigenValues to diagonal of quat  
   for( ip = 0; ip < 4; ip++ ) {
      b[ip] = eigenValues[ip] = quat[ip][ip];
      z[ip] = 0;
   }
   NROT = 0;
   for( int i = 0; i <= 50; i++ ) {
      // calculating sum of off-diagonal elements
      sm = 0;
      sm += (fabs(quat[0][1]) + fabs(quat[0][2]) 
         + fabs(quat[0][3]) + fabs(quat[1][2]) 
         + fabs(quat[1][3]) + fabs(quat[2][3])) ; 

      if( sm == 0 ) {
        //printEigens();
         eigsrt();
         //compute RMSD and return;
         rmsd = sqrt(fabs(eigenValues[3]/nVecs));
         return true;
      }

      if( i < 4 ) tresh = 0.2 * sm / ( 4 * 4 );
      else tresh = 0;
      
      for( ip = 0; ip < 3; ip++ ) {
         for( iq = ip + 1; iq < 4; iq++ ) {
            g = 100*fabs( quat[ip][iq] );
            // After four sweeps skip rotation if the off-diagonal element
            // is small
            if( i > 4 && fabs(eigenValues[ip])+g == fabs(eigenValues[ip])
             && fabs( eigenValues[iq] )+g == fabs(eigenValues[iq]) ) {
               quat[ip][iq] = 0;
            }
            else if( fabs(quat[ip][iq]) > tresh ) {
               h = eigenValues[iq]-eigenValues[ip];

               if(fabs(h)+g == fabs(h)) {
                  // t equals half theta
                  t = (quat[ip][iq])/h;
               }
               else {
                  theta = 0.5*h/(quat[ip][iq]);
                  t =  1/(fabs(theta)+ sqrt(1+theta*theta));
                  if(theta < 0) t = -t;
               }
               c = 1/sqrt(1+t*t );
               s = t*c;
               tau = s/(1+c);
               h = t * quat[ip][iq];
               z[ ip ] -= h;
               z[ iq ] += h;
               eigenValues[ ip ] -= h;
               eigenValues[ iq ] += h;
               quat[ip][iq] = 0;
               for( int j = 0; j <= ip - 1; j++ ){
                  ROTATE(quat,j,ip,j,iq) 
               }
               for( int j=ip+1; j<= iq -1; j++ ) {
                  ROTATE(quat,ip,j,j,iq) 
               }
               for( int j = iq+1; j < 4; j++ ) {
                  ROTATE(quat,ip,j,iq,j) 
               }
               for( int j = 0; j < 4; j++ ) {
                     ROTATE(eigenVectors,j,ip,j,iq) 
               }
               ++NROT;
            }
         }
      }

      for( ip = 0; ip < 4; ip++ ) {
         b[ip] += z[ ip ];
         eigenValues[ ip ] = b[ip];
         z[ip] = 0; 
      }
   }

   //At this stage diagonalizing failed. 
   rmsd = INFINITYVAL;
   eulerRotMat[0][0] = INFINITYVAL;
   return false;
}

void Superpose3DClass::computeRotationMatrix() {
  double q1 = eigenVectors[0][3];
  double q2 = eigenVectors[1][3];
  double q3 = eigenVectors[2][3];
  double q4 = eigenVectors[3][3];


   eulerRotMat[0][0] = SQR(q1)+SQR(q2)-SQR(q3)-SQR(q4);
   eulerRotMat[1][1] = SQR(q1)+SQR(q3)-SQR(q2)-SQR(q4);
   eulerRotMat[2][2] = SQR(q1)+SQR(q4)-SQR(q2)-SQR(q3);

   eulerRotMat[0][1] = 2*(q2*q3 - q1*q4);
   eulerRotMat[1][0] = 2*(q2*q3 + q1*q4);

   eulerRotMat[0][2] = 2*(q2*q4 + q1*q3);
   eulerRotMat[2][0] = 2*(q2*q4 - q1*q3);

   eulerRotMat[1][2] = 2*(q3*q4 - q1*q2);
   eulerRotMat[2][1] = 2*(q3*q4 + q1*q2);
}

bool Superpose3DClass::computeLSqFit(
  vector<vector<double> > &vecSetA,
  vector<vector<double> > &vecSetB
) {
  computeRotationalCenter(vecSetA,vecSetB);
  //printRotationalCenters();
  computeQuaternionMatrix(vecSetA,vecSetB);
  //printQuaternionMatrix();
  return diagonalize(); // rmsd is computed inside this routine.
}

bool Superpose3DClass::computeLSqFit(
  vector<vector<double> > &vecSetA,
  vector<vector<double> > &vecSetB,
  vector<double>&rcA,
  vector<double>&rcB
) {
  assignRotationalCenter(rcA,rcB);
  //printRotationalCenters();
  computeQuaternionMatrix(vecSetA,vecSetB);
  //printQuaternionMatrix();
  return diagonalize(); // rmsd is computed inside this routine.
}



bool Superpose3DClass::updateLSqFit_addition() {
	updateCenterOfMasses_addition();
  //printRotationalCenters();
  updateSufficientStatistics_addition();
  computeQuaternionMatrixUsingSufficientStats();
  //printQuaternionMatrix();
	return diagonalize(); // rmsd is computed inside this routine.
}

bool Superpose3DClass::updateLSqFit_subtraction() {
	updateCenterOfMasses_subtraction();
  //printRotationalCenters();
  updateSufficientStatistics_subtraction();
  computeQuaternionMatrixUsingSufficientStats();
  //printQuaternionMatrix();
	return diagonalize(); // rmsd is computed inside this routine.
}



/* Main interface. 
 * Performs superpositions from scratch using supplied vector sets 
 * vecSetA and vecSetB */
Superpose3DClass::Superpose3DClass(
  vector<vector<double> > &vecSetA, //moving
  vector<vector<double> > &vecSetB //fixed
) {
   INFINITYVAL = numeric_limits<double>::infinity();
   success_flag = 0;
   assert(vecSetA.size() == vecSetB.size());
   assert(vecSetA.size() >=3);

   //set spA and spB indexes for each of the supplied vector sets 
   spA = 0;
   spB = 0;

   //set nVecs variable
   nVecs = vecSetA.size(); 
   //sanity check on spA and spB
   assert(spA+nVecs == vecSetA.size());
   assert(spB+nVecs == vecSetB.size());


   // compute fit
   bool stat = computeLSqFit(vecSetA,vecSetB);
   //compute rotation matrix
   computeRotationMatrix();
   computeSufficientStatistics(vecSetA,vecSetB);
   //printRotationMatrix();
   //printEigens();
   success_flag = 1;
}


/* Overloaded constructor allowing superposition on subsets (size N) 
 * of moving and fixed coordinates vectors. */ 
Superpose3DClass::Superpose3DClass(
  vector<vector<double> > &vecSetA, //moving
  vector<vector<double> > &vecSetB, //fixed
  size_t spA,// start index in the moving set vectors
  size_t spB,// start index in the fixed  set vectors
  size_t N   // size of vectors  
) {
   success_flag = 0;
   INFINITYVAL = numeric_limits<double>::infinity();

   assert(N >= 3);
   assert(vecSetA.size() >= N);
   assert(vecSetB.size() >= N);

   //set nVecs variable
   nVecs = N;
   //set spA and spB indexes for each of the supplied vector sets 
   this->spA = spA;
   this->spB = spB;
   //sanity check on spA and spB
   assert(spA+N <= vecSetA.size());
   assert(spB+N <= vecSetB.size());

   // compute fit
   bool stat = computeLSqFit(vecSetA,vecSetB);
   //compute rotation matrix
   computeRotationMatrix();
   computeSufficientStatistics(vecSetA,vecSetB);
   //printRotationMatrix();
   //printEigens();
   success_flag = 1;
}


/* overloaded constructor allowing superposition using
   sufficient statistics of previous two partial superpositions 
   under ADDITION '+' or SUBTRACTION '-' */
Superpose3DClass::Superpose3DClass(
  SuffStatClass &stats1, 
  SuffStatClass &stats2, 
  char operation // can be only '+' or '-'
) {
  this->stats1 = stats1;
  this->stats2 = stats2;
  assert (operation == '+' || operation == '-');

  if (operation == '+') {
     this->nVecs = stats1.nVecs+stats2.nVecs;
     updateLSqFit_addition();
  }
  else {
     this->nVecs = stats1.nVecs-stats2.nVecs;
     updateLSqFit_subtraction();
  }
  computeRotationMatrix();
  success_flag = 1;
}

/* Overloaded constructor. 
   Accepts a pair of corresponding coords and simply computes 
   sufficient statistics and exits 
*/
Superpose3DClass::Superpose3DClass(
  vector<double>  &Acoord, //moving
  vector<double>  &Bcoord //fixed
) {
  spA = 0;
  spB = 0;
  vector<vector<double> > vecSetA;
  vector<vector<double> > vecSetB;
   success_flag = 0;
   vecSetA.push_back(Acoord);
   vecSetB.push_back(Bcoord);
   assert(vecSetA.size() == vecSetB.size() );

   nVecs = vecSetA.size(); 
   assert(nVecs == 1);
   INFINITYVAL = numeric_limits<double>::infinity();

   computeRotationalCenter(vecSetA,vecSetB);
   computeSufficientStatistics(vecSetA,vecSetB);
   success_flag = 0;
}


/* overloaded constructor allowing extension of a superposition by ONE PAIR
   building on the sufficient statistics of previous superposition */
Superpose3DClass::Superpose3DClass(
  SuffStatClass &stats1, 
  vector<double> &addedPointA, // moving point
  vector<double> &addedPointB,  // fixed point
  char operation
) {
  spA = 0;
  spB = 0;
  this->stats1 = stats1;
  Superpose3DClass supobj(addedPointA, addedPointB);  
  SuffStatClass stats2 = supobj.getSufficientStatistics();

  this->stats2 = stats2;
  assert (operation == '+' || operation == '-');

  if (operation == '+') {
     this->nVecs = stats1.nVecs+stats2.nVecs;
     updateLSqFit_addition();
  }
  else {
     this->nVecs = stats1.nVecs-stats2.nVecs;
     updateLSqFit_subtraction();
  }
  computeRotationMatrix();
  success_flag = 1;
}

Superpose3DClass::Superpose3DClass() {
}


double Superpose3DClass::getRMSD() {
   assert(success_flag == 1);
   return rmsd;
}


/* Transform a set of vectors using previously computed transform.  */
void Superpose3DClass::transformVectors(vector<vector<double> > &v) {
   assert(success_flag==1);
   for ( size_t i = 0; i < v.size(); i++ ) {
       double s[3];
       s[0] = v[i][0]; 
       s[1] = v[i][1]; 
       s[2] = v[i][2];

       //subtract moving set cm 
       s[0] -= rotacenterA[0];
       s[1] -= rotacenterA[1];
       s[2] -= rotacenterA[2];

      //rotate s using computed rotmat
       v[i][0] = v[i][1] = v[i][2] = 0;

       double t[3];
       t[0] = s[0]*eulerRotMat[0][0] + s[1]*eulerRotMat[0][1] + s[2]*eulerRotMat[0][2];
       t[1] = s[0]*eulerRotMat[1][0] + s[1]*eulerRotMat[1][1] + s[2]*eulerRotMat[1][2]; 
       t[2] = s[0]*eulerRotMat[2][0] + s[1]*eulerRotMat[2][1] + s[2]*eulerRotMat[2][2]; 

      //add fixed set cm
      t[0] += rotacenterB[0];
      t[1] += rotacenterB[1];
      t[2] += rotacenterB[2];

      v[i][0] = t[0];
      v[i][1] = t[1];
      v[i][2] = t[2];
   }
}

/* Transform a vector using previously computed transform.  */
void Superpose3DClass::transformVector(vector<double> &v) {
   assert(success_flag==1);
       double s[3];
       s[0] = v[0]; 
       s[1] = v[1]; 
       s[2] = v[2];

       //subtract moving set cm 
       s[0] -= rotacenterA[0];
       s[1] -= rotacenterA[1];
       s[2] -= rotacenterA[2];

      //rotate s using computed rotmat
       v[0] = v[1] = v[2] = 0;

       double t[3];
       t[0] = (s[0]*eulerRotMat[0][0]) + (s[1]*eulerRotMat[0][1]) 
    + (s[2]*eulerRotMat[0][2]);
       t[1] = (s[0]*eulerRotMat[1][0]) + (s[1]*eulerRotMat[1][1]) 
    + (s[2]*eulerRotMat[1][2]); 
       t[2] = (s[0]*eulerRotMat[2][0] + s[1]*eulerRotMat[2][1]) 
    + (s[2]*eulerRotMat[2][2]); 

      //add fixed set cm
      t[0] += rotacenterB[0];
      t[1] += rotacenterB[1];
      t[2] += rotacenterB[2];


      v[0] = t[0];
      v[1] = t[1];
      v[2] = t[2];
}

/* Copy the computed rotation matrix into the passed 
   argument.
*/
void Superpose3DClass::copyRotationMatrixInto(double rotmat[][3]) {
  assert(success_flag==1);
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      rotmat[i][j] = eulerRotMat[i][j];
    }
  }
}

/* Copy the computed rotational centers into the passed 
   arguments.
*/
void Superpose3DClass::copyRotationalCentersInto(
 double center_moving[3],
 double center_fixed[3]
) {
  assert(success_flag==1);
  for (size_t i = 0; i < 3; i++) {
    center_moving[i] = rotacenterA[i];
    center_fixed[i]  = rotacenterB[i];
  }
}



/* DEBUGGING */
void Superpose3DClass::printRotationalCenters() {
   cout << "Print Centers of Mass of moving set:\n";
   for (size_t  i = 0;  i < 3; i++) {
     cout << fixed << setw(8) << setprecision(3) << rotacenterA[i];
   }
   cout << endl;

   cout << "Print Centers of Mass of fixed set:\n";
   for (size_t  i = 0; i < 3; i++) {
     cout << fixed << setw(8) << setprecision(3) << rotacenterB[i];
   }
   cout << endl << "---------------------\n";
}

void Superpose3DClass::printEigens() {
   cout << "Print Eigen values\n";
   for (size_t  i = 0; i < 4; i++) {
     cout << fixed << setw(16) << setprecision(6) << eigenValues[i];
   }
   cout << endl << "-----\n";

   cout << "Print Eigen vectors\n";
   for( int i = 0; i < 4; i ++ ) {
      for( int j =0; j< 4; j++ ) {
        cout << fixed << setw(16) << setprecision(6) << eigenVectors[i][j];
      }
      cout << endl;
   }
   cout << "Quaternion norms: (should be ~1)" << endl;
   for( int i = 0; i < 4; i ++ ) {
        double q1 = eigenVectors[0][i];
        double q2 = eigenVectors[1][i];
        double q3 = eigenVectors[2][i];
        double q4 = eigenVectors[3][i];
        double qnorm =  SQR(q1)+SQR(q2) + SQR(q3) + SQR(q4);
        qnorm = sqrt(qnorm);
        cout << fixed << setw(16) << setprecision(6) << qnorm;
   }
   cout << endl << "---------------------\n";
}

void Superpose3DClass::printRotationMatrix() {
   cout << "Print Rotation matrix\n";
   for( int i = 0; i < 3; i++ ) {
      for( int j = 0; j < 3; j++ ) {
         cout << eulerRotMat[i][j] << " ";
      }
      cout << endl;
   }
}

void Superpose3DClass::printQuaternionMatrix() {
	cout << "Print Quaternion matrix\n";
  cout << fixed;
	for( int i = 0; i < 4; i++ ) {
		for( int j = 0; j < 4; j++ ) {
			cout << setw(10) << setprecision(3) << quat[i][j] << " ";
		}
		cout << endl;
	}
}


SuffStatClass Superpose3DClass::getSufficientStatistics() {
  assert(stats.setflag);
  return stats;
}
