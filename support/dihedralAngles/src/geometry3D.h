/* This file contains some useful functions for 3D Geometry manipulation */ 
#ifndef GEOMETRY3D_H
#define GEOMETRY3D_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <vector>
#include <assert.h>
using namespace std ;

/* This function computes direction cosines of a vector given its 
 * direction ratios. */ 
void computeDirectionCosines( double dratios[3], double dcosines[3] ); 


/* This functions computes the unit vector along the normal to a 
 * plane formed by two lines A and B */
void computeNormal( double dratiosA[3], double dratiosB[3], 
		double dcosnormal[3] ) ;


/* This function computes the cross product of two vectors A and B */
void computeCrossProduct( double dratiosA[3], double dratiosB[3], 
		double dratnormal[3] ) ;


/* This function computes the dot product of two vectors A and B */
void computeDotProduct( double dratiosA[3], double dratiosB[3], 
		double& dotproduct ); 


/* This function computes the box product (scalar triple product) of 
 * three vectors A, B, C */
void computeBoxProduct( double dratiosA[3], double dratiosB[3], 
		double dratiosC[3], double& boxproduct );


/* This function computes the rotation matrix for rotation through an 
 * angle theta about a line whose direction cosines are given by dcosines 
*/
void computeRotationMatrix( double dcosines[3], double theta, 
		double rotation_matrix[3][3] ) ;


/* This function rotates initial_vector to final_vector using the 
 * rotation matrix */
void rotateVector( double rotation_matrix[3][3], double initial_vector[3], 
		double final_vector[3] ); 


/* This function rotates a point about the origin using the rotation matrix */
void rotatePoint( double rotation_matrix[3][3], double xq, double yq, 
		double zq, double& newx, double& newy, double& newz ); 

/* This function performs tetrahedral fixing, i.e given coordinates x1 of 
 * the central atom and coordinates x2 and x3 of two attached atoms, and the 
 * tetrahedral angle theta, computes the coordinates x4 and x5 of the other 
 * two atoms attached to the central atom */
void tetrahedralFix( double x1[3], double x2[3], double x3[3], 
		double theta, double bond_length1, double bond_length2, 
		double x4[3], double x5[3] );

/* This function computes the dihedral angle between the planes 
 * ABC and BCD defined by the points A, B, C and D */
void computeDihedralAngle( const double A[3], const double B[3], 
	const double C[3], const double D[3], double& dihedral_angle);
  

/*compute angle between 2 vecotrs A and B*/
double computeAngle( double A[3], double B[3] ) ;



/* This function transforms an array of points X to a new coordinate system 
 * EX where the point NAT1 is origin, NAT2 lies on the x axis, and NAT3 lies 
 * on the xy plane */
void coordinateTransformer( double X[4][500],double EX[4][500], int N, 
		int NAT1, int NAT2, int NAT3);

/* This function finds the norm of difference between 2 vectors. 
 * Let A and B be 2 vectors, it returns the norm of the vector A-B */
double normAminusB( const double A[3] , const double B[3] ) ;

/* This function finds the norm of a position vector.*/
double vectorNorm( const double A[3] ) ;


/* overloading the above with std::vector arguments*/
double normAminusB( vector<double> A , vector<double> B ) ;

void projectPoint2Line( const double A[3], const double B[3], const double V[3], 
		double P[3] ) ;

/* distance between a point and a line */
double distPoint2Line( const double A[3], const double B[3], const double P[3] );

void projectPoint2Plane( double A[3], double X0[3], double vect1[3], 
		double vect2[3], double P[3] ) ;

/* Calculate the line of shortest distance between two 3D non-planar vectors.
 * Let the first vector be denoted by the points A1 and A2.
 * Let the second vector be denoted by the points B1 and B2.
 * Assume that the shortest line intersects the first vector at A and 
 * second at B.
 * L1 = A2-A1 and L2 = B2 - B1 
 * LP = A - B // perpendicular
 * where A = A1 + \lambda_1( A2 - A1 )
 * and   B = B1 + \lambda_2( B2 - B1 )
 * 
 * solving the equations L1.LP = 0  and L2.LP = 0 
 * gives \lambda_1 and \lambda_2
 *
 * This function takes A1, A2, B1, B2 and calculates points A and 
 * B of shortest distance between 2 vectors. 
 */

void computeMutalPerpendicularPoints( 
  double A1[3], 
  double A2[3], 
  double B1[3], 
  double B2[3], 
  double A[3], 
  double B[3] 
) ;


/* determined whether four given points A, B, C, D are coplanar or not
 */
bool areCoplanar( double A[3], double B[3], double C[3], double D[3] ) ;

/* determine whether the vecotrs AB and CD defined by points A, B, C, D are parallel/antiparallel or not
 */
bool areParallel( double A[3], double B[3], double C[3], double D[3] ) ;
bool areParallel2( double A[3], double B[3], double C[3], double D[3] ) ;



/* Determine orientation angle between two non-coplanar vectors.
 * Specifically, given are:
 * P1  a position vector
 * V1  direction cosines of some vector
 * P2  another position vector
 * V2  direction cosines of another vector,
 * we then want to find the orientation angle between vectors:
 * P1+\lambda_1{V1} and 
 * P2+\lambda_2{V2}
 *
 * The orientation angle will be in the range (-180, 180].
 */
	
double getOrientationAngle( double P1[3], double V1[3], double P2[3], double V2[3] ) ;

/* This functions determines whether a point C lies in between points A and B*/
bool isCbetweenAandB( double A[3], double B[3], double C[3] ) ;

/* This function takes a mx3 matrix A (std::vector type) and decomposes
 * it into UDV, where U, V are left and right orthogonal transformation 
 * matrices, and D is a diagonal matrix of singular values.
*/
int fsvd(vector<vector<double> >&, double [], double [][3]) ;
#endif
           
