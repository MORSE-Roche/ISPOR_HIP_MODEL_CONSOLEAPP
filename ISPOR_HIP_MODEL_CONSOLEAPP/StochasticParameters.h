#pragma once
#include "constantVariables.h"


class StochasticParameters
{
private:
	long seed;
	gsl_rng * rng;  // random number generator
	double MAXIMUM_NUMBER;

	double giveRandomNumber();
	void EstimateCoefficients(double Covariance[COVARIANCE_MATRIX_DIMENSION][COVARIANCE_MATRIX_DIMENSION], double (&Coeff)[COVARIANCE_MATRIX_DIMENSION], double Hazard[COVARIANCE_MATRIX_DIMENSION][HAZARD_ESTIMATES_DIMENSION]);

public:
	void GenerateParameters(double *&Parameters,
						bool Gender,
							   double tmpCovariance[COVARIANCE_MATRIX_DIMENSION][COVARIANCE_MATRIX_DIMENSION],
							   double tmpUtilities[UTILITY_ESTIMATES_DIMENSION1][DISTRIBUTION_ESTIMATE_DETAILS_DIMENSION],
							   double tmpRisk[RISK_ESTIMATES_DIMENSION1][DISTRIBUTION_ESTIMATE_DETAILS_DIMENSION],
							   double Hazard[COVARIANCE_MATRIX_DIMENSION][HAZARD_ESTIMATES_DIMENSION],
							   double Costs[COST_ESTIMATES_DIMENSION1][DISTRIBUTION_ESTIMATE_DETAILS_DIMENSION],
							   bool runningMode);

	
	bool choleskyDecomposition(double orig[COVARIANCE_MATRIX_DIMENSION][COVARIANCE_MATRIX_DIMENSION], double (&chol)[COVARIANCE_MATRIX_DIMENSION][COVARIANCE_MATRIX_DIMENSION]);
	
	StochasticParameters(void);
	~StochasticParameters(void);
};

