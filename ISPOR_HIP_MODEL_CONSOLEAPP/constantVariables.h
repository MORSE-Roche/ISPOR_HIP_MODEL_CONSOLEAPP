#pragma once
#include "omp.h"
#include <string>
#include <memory.h>
#include <process.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_exp.h>
using namespace std;

//Define maximum number of physical computational cores avialable in the machine
static const int OMP_MAX_NO_THREADS = omp_get_max_threads();

//
const long unsigned int TIME_HORIZON = 60;
const long unsigned int CohortSize = 1000;
const long unsigned int NumHealthStates = 5;
const long unsigned int NUMBER_OF_MODEL_PARAMETERS = 13;
const long unsigned int WTP = 20000;
const long unsigned int CEAC_WTP_VECTOR_SIZE = 57; //How many times WTP to calculate the NMB
const long unsigned int EFFECTIVE_TECHNOLOGY_LIFE_YEARS   = 10; //We assume that the technology is adopted for 10 years
const long unsigned int EFFECTIVE_POPULATION_SIZE_ADOPTING_TECHNOLOGY = 40000; 

//Define Treatment arms
const unsigned int NUMBER_OF_TREATMENT_ARMS = 3;
const unsigned int     MAX_NMB = NUMBER_OF_TREATMENT_ARMS +1;

const unsigned int STD   = 0;
const unsigned int TxNP1 = 1;
const unsigned int TxNP2 = 2;


//Define Genders
const bool FEMALE = false;
const bool MALE   = true;
const int NUMBER_OF_GENDERS = 2;


//Define AGE_INDECIES
const unsigned int DEFAULT_AGE     =60;
const unsigned int BaseCaseAge     =44;
const unsigned int AgeIntervals    =10;
const unsigned int AGE_INDEX_35_44 = 0;
const unsigned int AGE_INDEX_45_54 = 1;
const unsigned int AGE_INDEX_55_64 = 2;
const unsigned int AGE_INDEX_65_74 = 3;
const unsigned int AGE_INDEX_75_84 = 4;
const unsigned int AGE_INDEX_85    = 5;
const unsigned int AGE_DATA_DIMENSION1 = 6;
const unsigned int Mortality_RATES_FEMALE = 0;
const unsigned int Mortality_RATES_MALE = 1;


                 

//Define Running modes
const bool DETERMINISTIC = true; // Deterministic
const bool PSA           = false;// Probabilitic Sensitivity Analysis
const bool WITH_EVPPI            = true; // with    EVPPI
const bool WITHOUT_EVPPI         = false;// eithout EVPPI

//Health States
const unsigned int PrimaryTHR = 0;
const unsigned int SuccessP   = 1;
const unsigned int RevisionTHR= 2;
const unsigned int SuccessR   = 3;
const unsigned int Death      = 4;


//Define Cohort Model parameters
const unsigned int omrPTHR           = 0;
const unsigned int omrRTHR           = 1;
const unsigned int rrr               = 2;
const unsigned int rrNP1             = 3;
const unsigned int rrNP2             = 4;
const unsigned int lambda            = 5;
const unsigned int gamma             = 6;
const unsigned int cPrimary          = 7;
const unsigned int cRevision         = 8;
const unsigned int cSuccess          = 9;
const unsigned int uSuccessPrim      = 10;
const unsigned int uRevision         = 11;
const unsigned int uSuccess2nd       = 12;

//Define Regression Parameters
const unsigned int lngamma = 0;
const unsigned int consta  = 1;
const unsigned int Age     = 2;
const unsigned int Male    = 3;
const unsigned int NP1     = 4;
const unsigned int NP2     = 5;
const unsigned int COVARIANCE_MATRIX_DIMENSION = 6;
const unsigned int HAZARD_ESTIMATES_DIMENSION = 3;
const unsigned int HAZARD_RATIO =2;


//Define  distribution estimate
const unsigned int DISTRIBUTION_ESTIMATE_DETAILS_DIMENSION = 4;
const unsigned int DIST_MEAN    = 0;
const unsigned int DIST_SE      = 1;
const unsigned int DIST_ALPHA   = 2;
const unsigned int DIST_BETA    = 3;

//Define Utility Matrix dimensions
const unsigned int uSuccessP = 0;
const unsigned int uRev      = 1;
const unsigned int uSuccessR = 2;
const unsigned int UTILITY_ESTIMATES_DIMENSION1 = 3;

//Define Risk Matrix dimensions
const unsigned int RISK_ESTIMATES_DIMENSION1 = 3;

//Define Cost Matrix dimensions
const unsigned int COST_ESTIMATES_DIMENSION1 = 3;


//Define Raws of Heaslthstate input matrix
const unsigned int CHECK_RAW_SUM                 =0;
const unsigned int INCLUDE_FOR_LYGs              =1;
const unsigned int UTILIY_VALUES_PER_HEALTH_STATE =2;
const unsigned int SURGERY_COSTS_PER_HEALTH_STATE =3;
const unsigned int INPUT_PER_HEALTH_STATE_MATRIX_DIMENSION1 =4;


//Define Outcome RESULTS matrix indicies
const unsigned int Outcome_Result_QALY        = 0;
const unsigned int Outcome_Result_COST        = 1;
const unsigned int OUTCOME_RESULTS_DIMENSION1 = 2;

//Define Discount rates
const unsigned int Discountrate_Effect       = 0;
const unsigned int Discountrate_Cost         = 1;
const unsigned int Discountrate_DIMENSION1   = 2;


