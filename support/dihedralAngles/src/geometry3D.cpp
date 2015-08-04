/* This file contains some useful functions for 3D Geometry manipulation */ 
#include "geometry3D.h"
#define SQ(x) (x*x)

/* This function computes direction cosines of a vector given its 
 * direction ratios. */ 
void computeDirectionCosines( 
  double dratios[3], 
  double dcosines[3] 
) {
	double sqsum = 
		1/(sqrt( SQ(dratios[0]) + SQ(dratios[1]) + SQ(dratios[2])) ); 
	for( int i = 0 ; i < 3 ; i++ ) dcosines[i] = dratios[i] * sqsum ; 
}


/* This function computes the cross product of two vectors A and B */
void computeCrossProduct( 
  double dratiosA[3], 
  double dratiosB[3], 
  double dratnormal[3] 
) { 
	dratnormal[0] = (dratiosA[1]*dratiosB[2])-(dratiosA[2]*dratiosB[1]);
	dratnormal[1] = (dratiosA[2]*dratiosB[0])-(dratiosA[0]*dratiosB[2]);
	dratnormal[2] = (dratiosA[0]*dratiosB[1])-(dratiosA[1]*dratiosB[0]);
}


/* This functions computes the unit vector along the normal to a plane 
 * formed by two lines A and B */
void computeNormal( 
  double dratiosA[3], 
  double dratiosB[3], 
  double dcosnormal[3] 
) {
	double dratnormal[3] ; 
	computeCrossProduct(dratiosA,dratiosB,dratnormal) ; 
	computeDirectionCosines( dratnormal, dcosnormal ) ; 
} 


/* This function computes the dot product of two vectors A and B */
void computeDotProduct( 
  double dratiosA[3], 
  double dratiosB[3], 
  double& dotproduct 
) {
	dotproduct = 0.0 ; 
	for( int i = 0 ; i < 3 ; i++ ) dotproduct += dratiosA[i]*dratiosB[i];
}

/* This function computes the box product ( scalar triple product ) of 
 * three vectors A,B,C */
void computeBoxProduct( 
  double dratiosA[3], 
  double dratiosB[3], 
  double dratiosC[3], 
  double& boxproduct ) {
	double cross[3] ; 
	computeCrossProduct(dratiosB,dratiosC,cross); 
	computeDotProduct(dratiosA,cross,boxproduct); 
}

/* This function computes the rotation matrix for rotation through an angle 
 * theta about a line whose direction cosines are given by dcosines */
void computeRotationMatrix( 
  double dcosines[3], 
  double theta, 
  double rotation_matrix[3][3] 
) {
	double cost = cos(theta) ; 
	double sint = sin(theta) ; 
	double a1 = dcosines[0] * sint ; 
	double a2 = dcosines[1] * sint ; 
	double a3 = dcosines[2] * sint ; 
	double b1 = dcosines[1] * dcosines[2] * ( 1 - cost ) ; 
	double b2 = dcosines[2] * dcosines[0] * ( 1 - cost ) ; 
	double b3 = dcosines[0] * dcosines[1] * ( 1 - cost ) ; 
	for( int i = 0 ; i < 3 ; i++ ) {
		rotation_matrix[i][i] = cost 
			+ ( dcosines[i] * dcosines[i] *(1-cost) ) ; 
	}
	rotation_matrix[0][1] = b3 - a3 ; 
	rotation_matrix[1][0] = b3 + a3 ; 
	rotation_matrix[2][0] = b2 - a2 ; 
	rotation_matrix[0][2] = b2 + a2 ; 
	rotation_matrix[1][2] = b1 - a1 ; 
	rotation_matrix[2][1] = b1 + a1 ; 
}

/* This function rotates initial_vector to final_vector using the rotation 
 * matrix */
void rotateVector( 
  double rotation_matrix[3][3], 
  double initial_vector[3], 
  double final_vector[3] ) {
	for( int i = 0 ; i < 3 ; i++ ) { 
		final_vector[i] = 0.0 ; 
		for( int j = 0 ; j < 3 ; j++ ) {
			final_vector[i] += (initial_vector[j] 
					* rotation_matrix[i][j] ) ;
		}
	}	
}

/* This function rotates a point about the origin using the rotation matrix */
void rotatePoint( 
  double rotation_matrix[3][3], 
  double xq, 
  double yq, 
  double zq, 
  double& newx, 
  double& newy, 
  double& newz 
) {
	double initial_vector[3], final_vector[3] ; 
	initial_vector[0] = xq ; 
	initial_vector[1] = yq ; 
	initial_vector[2] = zq ; 
	rotateVector( rotation_matrix, initial_vector, final_vector) ;
	newx = final_vector[0] ; 
	newy = final_vector[1] ; 
	newz = final_vector[2] ; 
} 

/* This function performs tetrahedral fixing, i.e given coordinates x1 of 
 * the central atom and coordinates x2 and x3 of two attached atoms, and 
 * the tetrahedral angle theta, computes the coordinates x4 and x5 of the 
 * other two atoms attached to the central atom */
void tetrahedralFix( 
  double x1[3], 
  double x2[3], 
  double x3[3], 
  double theta, 
  double bond_length1, 
  double bond_length2, 
  double x4[3], 
  double x5[3] 
) {
	theta = theta * (M_PI/180) ; 
	double angle1 = theta/2.0 ;
	double v1[3], v2[3], v3[3], v4[3]; 
	double wrat[3], wcos[3], u1[3], u2[3], wfinal[3] ;
	int i = 0 ; 
	for( i = 0 ; i < 3 ; i++ ) { 
		v1[i] = x2[i] - x1[i] ; v2[i] = x3[i] - x1[i] ; 
	} 
	computeDirectionCosines( v1, u1 ) ; 
	computeDirectionCosines( v2, u2 ) ; 
	for( i = 0 ; i < 3 ; i++ ) {
		wrat[i] = - 0.5 * ( u1[i] + u2[i] ) ; 
	}

	computeNormal( v1, v2, v3 ) ; 
	computeNormal( v3, wrat, v4 ) ; 
	computeDirectionCosines( wrat, wcos ) ; 

	double angle2 = -(angle1) ; 
	double rotation_matrix[3][3] ; 
	computeRotationMatrix( v4, angle1, rotation_matrix ) ;          

	rotateVector( rotation_matrix, wcos, wfinal ) ; 
	for( i = 0 ; i < 3 ; i++ ) {
		x4[i] = x1[i] + ( wfinal[i] * bond_length1 ) ; 
	}
	computeRotationMatrix( v4, angle2, rotation_matrix ) ; 
	rotateVector( rotation_matrix, wcos, wfinal ) ; 

	for( i = 0 ; i < 3 ; i++ ) {
		x5[i] = x1[i] + ( wfinal[i] * bond_length2 ) ; 
	}
}

/* This function computes the dihedral angle between the planes ABC and 
 * BCD defined by the points A, B, C and D */
void computeDihedralAngle( 
  const double A[3], 
  const double B[3], 
  const double C[3], 
  const double D[3], 
  double& dihedral_angle
) {
	int i ; 
	double sign = 0.0 ; 
	double AB[3], BC[3], CD[3], normalABC[3], normalBCD[3] ; 
	for( i = 0 ; i < 3 ; i++ ) 
	 { 
	   AB[i] = B[i] - A[i] ; 
	   BC[i] = C[i] - B[i] ; 
	   CD[i] = D[i] - C[i] ; 
	 }
	computeNormal( AB, BC, normalABC ) ; 
	computeNormal( BC, CD, normalBCD ) ;
	computeDotProduct( normalABC, normalBCD, dihedral_angle ) ; 
	if( dihedral_angle > 1 ) {
		dihedral_angle = 1 ;
	}
	else if( dihedral_angle < -1 ) {
		dihedral_angle = -1 ;
	}
	dihedral_angle = acos( dihedral_angle ) ; 
	computeBoxProduct( BC, normalABC, normalBCD, sign ) ; 
	dihedral_angle = ( sign > 0 ) ? dihedral_angle : -dihedral_angle ; 
	dihedral_angle *= 180/M_PI ; 
}
  
/* This function transforms an array of points X to a new coordinate system 
 * EX where the point NAT1 is origin, NAT2 lies on the x axis, and NAT3 
 * lies on the xy plane */
void coordinateTransformer( 
  double X[4][500],
  double EX[4][500], 
  int N, 
  int NAT1, 
  int NAT2, 
  int NAT3
){ 
	int i,j,k ; 
	double AX[4][500] ; 
	double T[4][4] ; 
	double Tee[4][4] ; 
	double Tphi[4][4] ; 
	double en[4] ; 
	double sq[4] ; 
	for( k = 1 ; k <= 3 ; k++ ) {
		sq[k] = X[k][NAT1] ; 
		for( j = 1 ; j <= N ; j++ ) { 
			AX[k][j] = X[k][j] ; 
			X[k][j] -= sq[k] ;
		}
	}

	for( i = 1 ; i <= 3 ; i++ ) sq[i] = X[i][NAT2] * X[i][NAT2] ; 

	sq[1] = sqrt( sq[1] + sq[2] + sq[3] ) ; 
	T[1][1] = X[1][NAT2] / sq[1] ; 
	sq[2] = sqrt( sq[2] + sq[3] ) ; 
	sq[3] = 1.0 / sq[2] ; 
	sq[1] = 1.0 - ( T[1][1] * T[1][1] ) ; 
	sq[1] = sqrt( sq[1] ) ; 
	double cssq = 0.5 * ( 1.0 + T[1][1] ) ; 
	double snsq = 1 - cssq ;  

	for( k = 2 ; k <= 3 ; k++ ) {
		en[k] = X[k][NAT2] * sq[3] ;
		T[1][k] = en[k] * sq[1] ; 
		T[k][1] = - T[1][k] ; 
		sq[k] = en[k] * en[k] ; 
	}

	sq[1] = ( sq[3] - sq[2] ) * snsq ; 
	T[3][3] = cssq - sq[1] ; 
	T[2][2] = cssq + sq[1] ; 
	T[2][3] = snsq * en[2] * en[3] ; 
	T[2][3] *= -2.0 ; 
	T[3][2] = T[2][3] ; 

	for( i = 1 ; i <= 3 ; i++ ) { 
		en[i] = 0 ; 
		for( k =1 ; k <= 3 ; k++ ) 
			en[i] += T[i][k] * X[k][NAT3] ; 
	}
	sq[2] = ( en[2] * en[2] ) + ( en[3] * en[3] ) ; 
	sq[2] = 1.0 / sqrt( sq[2] ) ; 
	Tphi[2][2] = en[2] * sq[2] ; 
	Tphi[3][3] = Tphi[2][2] ; 
	Tphi[2][3] = en[3] * sq[2] ; 
	Tphi[3][2] = -Tphi[2][3] ; 
	Tphi[1][1] = 1 ; 
	Tphi[1][2] = 0 ; 
	Tphi[1][3] = 0 ; 
	Tphi[2][1] = 0 ; 
	Tphi[3][1] = 0 ; 
	for( i = 1 ; i <= 3 ; i++ ) {
		for( j = 1 ; j <= 3 ; j++ ) {
			Tee[i][j] = 0 ; 
			for( k = 1 ; k <= 3 ; k++ ) 
				Tee[i][j] += Tphi[i][k] * T[k][j] ; 
	   	}
	}

	for( j = 1 ; j <= N ; j++ ) {  
	   for( i = 1 ; i <= 3 ; i++ ) {
		EX[i][j] = 0 ; 
		for( k = 1 ; k <= 3 ; k++ ) 
			EX[i][j] += Tee[i][k] * X[k][j] ;
	   }
	}

}
           
double normAminusB( const double A[3] , const double B[3] )
{
	double dratios[3] ;
	dratios[0] = A[0] - B[0] ;
	dratios[1] = A[1] - B[1] ;
	dratios[2] = A[2] - B[2] ;

	return sqrt( SQ(dratios[0]) + SQ(dratios[1]) + SQ(dratios[2]) )  ;
}
 
double normAminusB( vector<double> A , vector<double> B )
{
	double dratios[3] ;
	dratios[0] = A[0] - B[0] ;
	dratios[1] = A[1] - B[1] ;
	dratios[2] = A[2] - B[2] ;

	return sqrt( SQ(dratios[0]) + SQ(dratios[1]) + SQ(dratios[2]) )  ;
}

double vectorNorm( const double A[3] ) {
	static double origin[3] = {0,0,0} ;
	return normAminusB( A, origin ) ;
}



double computeAngle( double A[3], double B[3] ) {

	double origin[3] = {0.0,0.0,0.0} ;
	double dp  = 0.0 ;
       	computeDotProduct( A, B, dp) ;

	double normA = normAminusB( A, origin ) ;
	double normB = normAminusB( B, origin ) ;
	//cout << dp << " " << normA << " " << normB << endl ;
	double cost = (dp/ (normA*normB)) ;
	  if( cost > 1 ) 
	  {
		  cost = 1 ;
	  }
	  else if( cost < -1 ) 
	  {
		  cost = -1 ;
	  }
	  double theta = acos( cost ) ; 

	  //find sign 
	  double normalAB[3] ;
	  double sign = 0.0 ;
	  computeNormal( A, B, normalAB ) ; 

	  computeBoxProduct( normalAB, A, B, sign ) ; 
	  theta = ( sign > 0 ) ? theta : -theta; 
	  //theta *= 180/M_PI ; 

	return ( theta )	;
}
/* This function computes the dihedral angle between the planes ABC and BCD defined by vecotrs 
   AB, BC and CD 
*/
void compute_dihedral_angle( double AB[3], double BC[3], double CD[3],  double& dihedral_angle)
{
  double sign = 0.0 ; 

  double normalABC[3], normalBCD[3] ; 
  computeNormal( AB, BC, normalABC ) ; 
  computeNormal( BC, CD, normalBCD ) ;
  computeDotProduct( normalABC, normalBCD, dihedral_angle ) ; 
  if( dihedral_angle > 1 ) 
  {
	  //cout << setprecision(13) <<  (double)dihedral_angle << endl ; 
	  dihedral_angle = 1 ;
  	  //exit(0);
  }
  else if( dihedral_angle < -1 ) 
  {
	  //cout << "*" << setprecision(13) <<  (double)dihedral_angle << endl ; 
	  dihedral_angle = -1 ;
	  //exit(0);
  }
  dihedral_angle = acos( dihedral_angle ) ; 
  computeBoxProduct( BC, normalABC, normalBCD, sign ) ; 
  dihedral_angle = ( sign > 0 ) ? dihedral_angle : -dihedral_angle ; 
  dihedral_angle *= 180/M_PI ; 
}
 
/* 
 * This function calculates the projection P of a point A
 * onto a line represented by a vector V, passing
 * through another point B. 
 *        P             B
 *    ----*-------------*-----------> V
 *        |
 *        |
 *        * A
 *
 * The general solution for P is as follows:
 * V.AP = 0 
 * where P = B + lambda*V
 * Then solving for lambda (in respective components), 
 * we get:
 * lambda = sum_i=[0,3) (a[i]-b[i])*v[i] /sum_i=[0,3) v[i]^2
 * */
void projectPoint2Line( 
  const double A[3],
  const double B[3],
  const double V[3],
  double P[3]
) {
	static double origin[3] = {0,0,0} ;
	double normV =  normAminusB(V,origin) ;

	double lambda = 0 ;
	for( int i = 0 ; i < 3 ; i++ ) {
		lambda += ((A[i]-B[i])*V[i]) ;
	}

	lambda /= (normV*normV) ;


	for( int i = 0 ; i < 3 ; i++ ) {
		P[i] = B[i] +lambda*(V[i]) ;
	}
}

/*
 * This function calculates the shortest distance of a point P
 * to a line passing through two points A and B. 
 *      A   x               B
 *     -*---+---------------*------>
 *          |
 *          | d
 *          |
 *          * P
 * */

double distPoint2Line(
  const double A[3],
  const double B[3],
  const double P[3]
) {
	//unit vec passing through A in the direction of B
	double V[3] ;
	V[0] = B[0]-A[0] ;
	V[1] = B[1]-A[1] ;
	V[2] = B[2]-A[2] ;


	double X[3] ;
	projectPoint2Line( P, A, V, X ) ;
	return normAminusB(P,X) ;
}



/* This function calculates the projection P of a point A onto a plane
 * The function needs a point X0 through which the plane passes 
 * and two vectors in the plane */
void projectPoint2Plane( 
  double A[3], 
  double X0[3], 
  double vect1[3], 
  double vect2[3], 
  double P[3] 
) {

	double normal[3] ;
	computeNormal( vect1, vect2, normal) ;

	/* Rx+Sy+Tz+U = 0 defines a plane with a non-zero normal vector N (R,S,T) passing through a point X0 = (x0,y0,z0).
	 * N.(X-X0) = 0 ;
	 * =>(R,S,T).([x-x0],[y-y0],[z-z0]) = 0
	 * =>Rx+Sy+Tz - (R*X0 + S*Y0 + T*Z0) = 0 ;
	 * Therefore U = - (R*X0 + S*Y0 + T*Z0) 
	 */

	double plane[4] ;  // (R,S,T,U)
	plane[0] =  normal[0] ;
	plane[1] =  normal[1] ;
	plane[2] =  normal[2] ;
	plane[3] =  -( (normal[0]*X0[0]) + (normal[1]*X0[1]) + (normal[2]*X0[2]) );

	/*
	 * Projection distance (point-plane distance)
	 * D = (R*X0+S*Y0+T*Z0+U) / sqrt( R**2 + S**2 + T**2 ) 
	 */

	double sqrsum = 1/ ( sqrt( (plane[0]*plane[0]) + (plane[1]*plane[1]) + (plane[2]*plane[2]) ) );  // denomenator
	double distance = (plane[0]*A[0]) + (plane[1]*A[1]) + (plane[2]*A[2]) + plane[3] ; // numerator
	distance *= sqrsum ;

	/* Note the distance is signed. Positive distance indicates that
	 * the point A -- that is to be projected -- is on the same side of the plane as the normal vector
	 * and negative if on the other side.
	 */
	if( distance >= 0 ) {
		P[0] = A[0]-normal[0]*distance ;
		P[1] = A[1]-normal[1]*distance ;
		P[2] = A[2]-normal[2]*distance ;
	}
	else{
		P[0] = A[0]+normal[0]*distance ;
		P[1] = A[1]+normal[1]*distance ;
		P[2] = A[2]+normal[2]*distance ;
	}
}

/* The routine below calculates the points associated with the
 * mutual perpendicular between two 3D non-planar vectors.
 * Note that the mutual perpendicular is also the line associated
 * with the shortest distance between 2 points.
 *
 * The general method to 
 * calculate the line of shortest distance between two 3D non-planar vectors
 * is as follows:
 *
 * *  Let the first vector be denoted by the points A1 and A2.
 * *  Let the second vector be denoted by the points B1 and B2.
 * *  Assume that the shortest line intersects the first vector at A and 
 *    second at B.
 *
 * * Then L1 = A2-A1 and L2 = B2 - B1 
 * * and  LP = A - B  is the mutual perpendicular
 * where A = A1 + \lambda_1( A2 - A1 )
 * and   B = B1 + \lambda_2( B2 - B1 )
 * 
 * solving the equations L1.LP = 0  and L2.LP = 0 gives 
 * \lambda_1 and \lambda_2
 *
 * This function takes as arguments points A1, A2, B1, B2, where 
 * A2-A1 forms one non-planar vector and
 * B2-B1 forms another non-planar vector.
 * 
 *
 * The function computes points A and B of shortest distance between 
 * 2 vectors where B-A or A-B is the mutual perpendicular. 
 */

/*First, the funtion below is just a subroutine used by 
 * computeMutalPerpendicularPoints */
double C_value( double m[3], double n[3], double o[3], double p[3] ) {
	double val = 0 ;
	val += (m[0]-n[0]) * (o[0]-p[0]) ;
	val += (m[1]-n[1]) * (o[1]-p[1]) ;
	val += (m[2]-n[2]) * (o[2]-p[2]) ;
	return val ;
}
void computeMutalPerpendicularPoints( 
  double A1[3], 
  double A2[3], 
  double B1[3], 
  double B2[3], 
  double A[3], 
  double B[3] ) {
	
	double lambda1, lambda2 ;
	/* coefficients of the equation C1 + lambda_1 C2 - lambda_2 C3 */
	double C1 = C_value( A2, A1, A1, B1 ) ; 
	double C2 = C_value( A2, A1, A2, A1 ) ; 
	double C3 = C_value( A2, A1, B2, B1 ) ; 

	/* coefficients of the equation C4 + lambda_1 C5 - lambda_2 C6 */
	double C4 = C_value( B2, B1, A1, B1 ) ; 
	double C5 = C_value( B2, B1, A2, A1 ) ; 
	double C6 = C_value( B2, B1, B2, B1 ) ; 

	lambda1 = ( (C1*C6) - (C3*C4) )/ ( (C5*C3) - (C2*C6) ) ;
	lambda2 = ( (C1*C5) - (C2*C4) )/ ( (C5*C3) - (C2*C6) ) ;

	/* Calculate A as A1 + lambda1 (A2-A1) */
	A[0] = A1[0] + lambda1* ( A2[0] - A1[0] ) ;
	A[1] = A1[1] + lambda1* ( A2[1] - A1[1] ) ;
	A[2] = A1[2] + lambda1* ( A2[2] - A1[2] ) ;

	/* Calculate B as B1 + lambda2 (B2-B1) */
	B[0] = B1[0] + lambda2* ( B2[0] - B1[0] ) ;
	B[1] = B1[1] + lambda2* ( B2[1] - B1[1] ) ;
	B[2] = B1[2] + lambda2* ( B2[2] - B1[2] ) ;
}

/* Determine whether or not four given points A, B, C, D are coplanar */
bool areCoplanar( double A[3], double B[3], double C[3], double D[3] ) {
	double AB[3], AC[3], AD[3] ;

	AB[0] = B[0] - A[0] ; AB[1] = B[1] - A[1] ; AB[2] = B[2] - A[2] ;
	AC[0] = C[0] - A[0] ; AC[1] = C[1] - A[1] ; AC[2] = C[2] - A[2] ;
	AD[0] = D[0] - A[0] ; AD[1] = D[1] - A[1] ; AD[2] = D[2] - A[2] ;

	double volume ; 
	computeBoxProduct(AB, AC, AD, volume ) ;
	if( volume == 0 ) return 1 ;
	else return 0 ;
}

/* Determine whether or not vectors B-A and D-C defined by points A, B, C, D 
 * are parallel/antiparallel */
bool areParallel2( double A[3], double B[3], double C[3], double D[3] ) {
	double AB[3], CD[3] ; 
	AB[0] = B[0] - A[0] ; AB[1] = B[1] - A[1] ; AB[2] = B[2] - A[2] ;
	CD[0] = D[0] - C[0] ; CD[1] = D[1] - C[1] ; CD[2] = D[2] - C[2] ;

	double normal[3] ; 
	computeCrossProduct( AB, CD, normal ) ;
	if ( normal[0] == 0 && normal[1] == 0 && normal[2] == 0 ) return 1 ;
	
	return 0;
}
bool areParallel( double A[3], double B[3], double C[3], double D[3] ) {
	double AB[3], CD[3] ; 
	AB[0] = B[0] - A[0] ; AB[1] = B[1] - A[1] ; AB[2] = B[2] - A[2] ;
	CD[0] = D[0] - C[0] ; CD[1] = D[1] - C[1] ; CD[2] = D[2] - C[2] ;

	double dp ; 
	computeDotProduct( AB, CD, dp) ;
	double normAB = normAminusB(A,B) ;
	double normCD = normAminusB(C,D) ;
	assert( normAB > 0 ) ;
	assert( normCD > 0 ) ;
	double cost = dp/(normAB*normCD) ;
	if ( cost > 1- 0.00001 ) return 1 ;
	else {
		cout << cost << endl ;
		return 0;
	}
}

/* Determine orientation angle between two non-coplanar vectors.
 * Specifically, given are:
 * P1  a position vector
 * A1  direction cosines of some vector
 * P2  another position vector
 * A2  direction cosines of another vector,
 * we then want to find the orientation angle between vectors:
 * P1+\lambda_1{A1} and 
 * P2+\lambda_2{A2}
 *
 * The orientation angle will be in the range (-180, 180].
 */
double getOrientationAngle( double P1[3], double V1[3], double P2[3], double V2[3] ) {
	static double A1[3], A2[3], B1[3], B2[3] ;
	static double X[3], Y[3] ;
	static double Xp[3], Yp[3] ;
	static double origin[3] = {0,0,0} ;
	const double lambda = 10 ;


	A1[0] = P1[0] ; 
	A1[1] = P1[1] ; 
	A1[2] = P1[2] ;
	
	A2[0] = P1[0]+lambda*V1[0] ; 
	A2[1] = P1[1]+lambda*V1[1] ; 
	A2[2] = P1[2]+lambda*V1[2] ;

	B1[0] = P2[0] ; 
	B1[1] = P2[1] ; 
	B1[2] = P2[2] ;
	
	B2[0] = P2[0]+lambda*V2[0] ; 
	B2[1] = P2[1]+lambda*V2[1] ; 
	B2[2] = P2[2]+lambda*V2[2] ;
	assert( areParallel(A1,A2, origin, V1) ) ;
	assert( areParallel(B1,B2, origin, V2) ) ;

	computeMutalPerpendicularPoints(A1,A2,B1,B2,X,Y) ;
	// find points Xp Yp such that vector X-Xp is parallel to axisA
	// and Yp-Y is parallel to axisB
	Xp[0] = X[0]+lambda*V1[0] ;
	Xp[1] = X[1]+lambda*V1[1] ;
	Xp[2] = X[2]+lambda*V1[2] ;

	Yp[0] = Y[0]+lambda*V2[0] ;
	Yp[1] = Y[1]+lambda*V2[1] ;
	Yp[2] = Y[2]+lambda*V2[2] ;


	assert( areParallel(X,Xp, origin, V1) ) ;
	assert( areParallel(Y,Yp, origin, V2) ) ;

	double omega = -999.99 ;
	computeDihedralAngle( Xp, X, Y, Yp, omega) ;
	return omega ;

}

/* Given points A, B, and C, this routine, checks whether C lies 
 * between A and B. Returns True or False depending on the result. */
bool isCbetweenAandB( double A[3], double B[3], double C[3] ) {
	double BminusA[3], CminusA[3] ;
	BminusA[0] = B[0] - A[0] ;
	BminusA[1] = B[1] - A[1] ;
	BminusA[2] = B[2] - A[2] ;

	CminusA[0] = C[0] - A[0] ;
	CminusA[1] = C[1] - A[1] ;
	CminusA[2] = C[2] - A[2] ;

	double normal[3] ;
	double origin[3] = {0,0,0} ;
	computeCrossProduct( CminusA, BminusA, normal) ;
	double norm = normAminusB( normal, origin ) ;
	
	if ( norm > 0.001 ) return false;
	
	double dp ;
	computeDotProduct(CminusA,BminusA,dp) ;

	if(dp < 0) return false ;

	norm = normAminusB(A,B) ;		          
	if(dp > norm*norm) return false ;

	return true ;
}



/* The svd code has been taken from the following link:
 * http://www.public.iastate.edu/~dicook/JSS/paper/code/svd.c 
 * Changes have been made to comply with my code.
 * The following header was embedding in the code and is produced /in toto/:
 * --------------------------------------------------------------------
 * svdcomp - SVD decomposition routine. 
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is 
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to fsvd is as follows:
 *   a = mx3 matrix to be decomposed, gets overwritten with u
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
 *   -----------------------------------------------------------------
*/

#define SIGN_SVD(a,b) ( (b >= 0.0) ? (fabs(a)) : (-fabs(a)) ) 
#define MAX_SVD(x,y)  ( (x) > (y) ? (x) : (y) ) 
 
static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}


int fsvd(vector<vector<double> > &a, double w[], double v[][3])
{
    const int m = a.size() ;
    const int n = 3 ;    
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        //fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }
  
    rv1 = (double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs((double)a[k][i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k][i] = (double)((double)a[k][i]/scale);
                    s += ((double)a[k][i] * (double)a[k][i]);
                }
                f = (double)a[i][i];
                g = -SIGN_SVD(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (double)(f - g);
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += (double)(f * (double)a[k][i]);
                    }
                }
                for (k = i; k < m; k++) 
                    a[k][i] = (double)((double)a[k][i]*scale);
            }
        }
        w[i] = (double)(scale * g);
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs((double)a[i][k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i][k] = (double)((double)a[i][k]/scale);
                    s += ((double)a[i][k] * (double)a[i][k]);
                }
                f = (double)a[i][l];
                g = -SIGN_SVD(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (double)(f - g);
                for (k = l; k < n; k++) 
                    rv1[k] = (double)a[i][k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += ((double)a[j][k] * (double)a[i][k]);
                        for (k = l; k < n; k++) 
                            a[j][k] += (double)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] = (double)((double)a[i][k]*scale);
            }
        }
        anorm = MAX_SVD(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v[j][i] = (double)(((double)a[i][j] / (double)a[i][l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += ((double)a[i][k] * (double)v[k][j]);
                    for (k = l; k < n; k++) 
                        v[k][j] += (double)(s * (double)v[k][i]);
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
  
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i][j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += ((double)a[k][i] * (double)a[k][j]);
                    f = (s / (double)a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += (double)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] = (double)((double)a[j][i]*g);
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (double)h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (double)(y * c + z * s);
                            a[j][i] = (double)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k] = (double)(-z);
                    for (j = 0; j < n; j++) 
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                //fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN_SVD(g, f)))-h))/x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (double)(x * c + z * s);
                    v[jj][i] = (double)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (double)z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (double)(y * c + z * s);
                    a[jj][i] = (double)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (double)x;
        }
    }
    free((void*) rv1);
    return(1);
}
