#include "StochasticParameters.h"

StochasticParameters::StochasticParameters(void)
{
	seed = (long)(time (NULL)) * _getpid();
	rng = gsl_rng_alloc (gsl_rng_mt19937);     // pick random number generator
	 
	gsl_rng_env_setup();
	gsl_rng_set(rng, seed);
	MAXIMUM_NUMBER = gsl_rng_max(rng);
}

double StochasticParameters::giveRandomNumber(){
	return gsl_rng_get(rng)/MAXIMUM_NUMBER;
}


void StochasticParameters::GenerateParameters(double *&ParametersSet,bool Gender, double tmpCovariance[COVARIANCE_MATRIX_DIMENSION][COVARIANCE_MATRIX_DIMENSION], double tmpUtilities[UTILITY_ESTIMATES_DIMENSION1][DISTRIBUTION_ESTIMATE_DETAILS_DIMENSION], double tmpRisk[RISK_ESTIMATES_DIMENSION1][DISTRIBUTION_ESTIMATE_DETAILS_DIMENSION], double Hazard[COVARIANCE_MATRIX_DIMENSION][HAZARD_ESTIMATES_DIMENSION], double Costs[COST_ESTIMATES_DIMENSION1][DISTRIBUTION_ESTIMATE_DETAILS_DIMENSION], bool runningMode)
{
	double myRND; 

	//Define local variables
	double Coeff[COVARIANCE_MATRIX_DIMENSION]     = {0};
	double Lamdba = 0;
	double Gamma  = 0;



	if(runningMode == DETERMINISTIC){		
		//omrPTHR
		ParametersSet[omrPTHR] = tmpRisk[omrPTHR][DIST_MEAN];
		//omrRTHR
		ParametersSet[omrRTHR] = tmpRisk[omrRTHR][DIST_MEAN];
		//rrr
		ParametersSet[rrr]     = tmpRisk[rrr][DIST_MEAN];
		//rrNP1
		ParametersSet[rrNP1] = gsl_sf_exp(Hazard[NP1][DIST_MEAN]);
		//rrNP2
		ParametersSet[rrNP2] = gsl_sf_exp(Hazard[NP2][DIST_MEAN]);
		//lambda
		ParametersSet[lambda]= gsl_sf_exp(Hazard[consta][DIST_MEAN] + DEFAULT_AGE * Hazard[Age][DIST_MEAN] - (int)(Gender) * Hazard[Male][DIST_MEAN]);
		//gamma
		ParametersSet[gamma] = gsl_sf_exp(Hazard[lngamma][DIST_MEAN]);
		//cPrimary
		ParametersSet[cPrimary] = Costs[uSuccessP][DIST_MEAN];
		//cRevision
		ParametersSet[cRevision]= Costs[uRev][DIST_MEAN];
		//cSuccess
		ParametersSet[cSuccess] = Costs[uSuccessR][DIST_MEAN];
		//uSuccessPrim
		ParametersSet[uSuccessPrim] = tmpUtilities[uSuccessP][DIST_MEAN];
		//uSuccessHavingRev
		ParametersSet[uRevision]    = tmpUtilities[uRev][DIST_MEAN];
		//uRevisionDuring
		ParametersSet[uSuccess2nd]  = tmpUtilities[uSuccessR][DIST_MEAN];	
	}
	else
	{
		//CALCULATE COEFFICENTS BASED ON COVARIANCE MATRIX
		EstimateCoefficients(tmpCovariance,Coeff,Hazard);

		//PROBABILISTIC OR EVPI	OR EVPPI	
		//omrPTHR  
		myRND = this->giveRandomNumber();
		ParametersSet[omrPTHR] = gsl_cdf_beta_Qinv(myRND, tmpRisk[STD][DIST_ALPHA], tmpRisk[STD][DIST_BETA]);
		//omrRTHR
		myRND = this->giveRandomNumber();
		ParametersSet[omrRTHR] = gsl_cdf_beta_Qinv(myRND, tmpRisk[TxNP1][DIST_ALPHA], tmpRisk[TxNP1][DIST_BETA]);
		//rrr
		myRND = this->giveRandomNumber();
		ParametersSet[rrr] = gsl_cdf_beta_Qinv(myRND, tmpRisk[TxNP2][DIST_ALPHA], tmpRisk[TxNP2][DIST_BETA]);
		//rrNP1
		ParametersSet[rrNP1]  = gsl_sf_exp(Coeff[NP1]);
		//rrNP1
		ParametersSet[rrNP2]  = gsl_sf_exp(Coeff[NP2]);
		//lambda
		ParametersSet[lambda] = gsl_sf_exp(Coeff[consta] + DEFAULT_AGE * Coeff[Age] - (int)(Gender) * Coeff[Male]);      
		//gamma
		ParametersSet[gamma]  = gsl_sf_exp(Coeff[lngamma]);
		//cPrimary      
		ParametersSet[cPrimary] = 0;
		//cRevision
		myRND = this->giveRandomNumber();
		ParametersSet[cRevision] = gsl_cdf_gamma_Qinv(myRND, Costs[uRev][DIST_ALPHA], Costs[uRev][DIST_BETA]);
		//cSuccess
		ParametersSet[cSuccess] = 0;
		//uSuccessPrim
		myRND = this->giveRandomNumber();
		ParametersSet[uSuccessPrim]      = gsl_cdf_beta_Qinv(myRND, tmpUtilities[uSuccessP][DIST_ALPHA], tmpUtilities[uSuccessP][DIST_BETA]);
		//uSuccessHavingRev
		myRND = this->giveRandomNumber();
		ParametersSet[uRevision]         = gsl_cdf_beta_Qinv(myRND, tmpUtilities[uRev][DIST_ALPHA], tmpUtilities[uRev][DIST_BETA]);
		//uRevisionDuring
		myRND = this->giveRandomNumber();
		ParametersSet[uSuccess2nd]   =     gsl_cdf_beta_Qinv(myRND, tmpUtilities[uSuccessR][DIST_ALPHA], tmpUtilities[uSuccessR][DIST_BETA]);
	}
	return;
}

void StochasticParameters::EstimateCoefficients(double Covariance[COVARIANCE_MATRIX_DIMENSION][COVARIANCE_MATRIX_DIMENSION], double (&Coeff)[COVARIANCE_MATRIX_DIMENSION], double Hazard[COVARIANCE_MATRIX_DIMENSION][HAZARD_ESTIMATES_DIMENSION])
{
	double myRND; 
	double choleskyMatrix[COVARIANCE_MATRIX_DIMENSION][COVARIANCE_MATRIX_DIMENSION]={{0}};
	double z_matrix[COVARIANCE_MATRIX_DIMENSION] = {{0}};
	bool successful_decomposition_process = choleskyDecomposition(Covariance,choleskyMatrix);
	if(successful_decomposition_process == true){
		for(int i=0;i<COVARIANCE_MATRIX_DIMENSION;i++){
			myRND = this->giveRandomNumber();
			z_matrix[i] = gsl_cdf_gaussian_Qinv(myRND,1);
		}
	}
	else{
		exit(-1);
	}

	for(int col = 0; col<COVARIANCE_MATRIX_DIMENSION; col++){
		for(int row = 0; row <COVARIANCE_MATRIX_DIMENSION;row++){
			//Cholesky matrix here is upper triangle not the lower triangle then we transpose it
			Coeff[col] += choleskyMatrix[row][col]*z_matrix[row];
		}
	}

	for(int coefficient = 0; coefficient <COVARIANCE_MATRIX_DIMENSION;coefficient++){
		Coeff[coefficient] += Hazard[coefficient][DIST_MEAN];
	}
}
						 
												 
bool StochasticParameters::choleskyDecomposition(double orig[COVARIANCE_MATRIX_DIMENSION][COVARIANCE_MATRIX_DIMENSION], double (&chol)[COVARIANCE_MATRIX_DIMENSION][COVARIANCE_MATRIX_DIMENSION])
{
	/* file: choesky.c 
    Take the cholesky decomposition in the manner described in FA Graybill
   (1976).

	Do the augmented cholesky decomposition as described in FA Graybill
	(1976) Theory and Application of the Linear Model. The original matrix
	must be symmetric positive definite. The augmentation matrix, or
	series of column vectors, are multiplied by C^-t, where C is the
	upper triangular cholesky matrix, ie C^t * C = M and M is the original
	matrix. Returns with a value of 0 if M is a non-positive definite 
	matrix. Returns with a value of 1 with succesful completion.

	Arguments:

	orig (input) double n x n array. The matrix to take the Cholesky
	      decomposition of.
	n    (input) integer. Number of rows and columns in orig.
	aug  (input) double n x mcol array. The matrix for the augmented
	      part of the decomposition.
	mcol (input) integer. Number of columns in aug.
	chol (output) double n x n array. Holds the upper triangular matrix
	      C on output. The lower triangular portion remains unchanged.
	      This maybe the same as orig, in which case the upper triangular
	      portion of orig is overwritten.
	cholaug (output) double n x mcol array. Holds the product C^-t * aug.
	         May be the same as aug, in which case aug is over written.
	ofs (input) integer. The index of the first element in the matrices.
	     Normally this is 0, but commonly is 1 (but may be any integer).
     */

   int i, j, k;//, l;
   bool retval = true;
   int ofs    = 0;
   int n = COVARIANCE_MATRIX_DIMENSION;


   

   for (i=ofs; i<n+ofs; i++) {
	   chol[i][i] = orig[i][i];
	   for (k=ofs; k<i; k++)
		   chol[i][i] -= chol[k][i]*chol[k][i];
	   if (chol[i][i] <= 0) {
		   retval = false;
		   return retval;
	   }
	   chol[i][i] = sqrt(chol[i][i]);

       for (j=i+1; j<n+ofs; j++) {
		   chol[i][j] = orig[i][j];
		   for (k=ofs; k<i; k++)
			   chol[i][j] -= chol[k][i]*chol[k][j];
		   chol[i][j] /= chol[i][i];
	   }
   }

   return retval;
}

StochasticParameters::~StochasticParameters(void)
{
	gsl_rng_free(rng);
}
