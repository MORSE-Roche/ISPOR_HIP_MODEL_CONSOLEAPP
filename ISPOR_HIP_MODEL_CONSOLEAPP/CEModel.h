#pragma once
#include "ModelOutcome.h"
#include "ModelInput.h"
#include "StochasticParameters.h"

class CEModel
{
private:
	bool runningmode;
	bool EVPPImode;
	unsigned int Number_of_Trials;
	double DiscountRates[Discountrate_DIMENSION1][TIME_HORIZON];
	
	double*** Cohort;
	long unsigned int NMB_SIZE1, NMB_SIZE2;

	double** NMB_RESULTS; // Net Monetry Benefits
	double* VoI_RESULTS;
	double EVPI;
	
	

	void reset(unsigned int myTx);

	double RevisionRisk(unsigned int cycle, 
						double RelativeRisk, 
						double Lambda, 
						double Gamma);

	double CalculateMortalityRisk(bool Gender, 
								unsigned int cycle);
public:
	double Cumulative_Effective_Population;
	double** VoI_per_WTP;
	//Model Parameters
	StochasticParameters myStochasticParametersGenerator;
	//Model Outcome
	ModelOutcome         myModelOutcome;
	ModelInput           myModelInput;
	
	//Functions
	unsigned int getNumberOfTrials();
	void ModelTreatmentArms(double Parameters[NUMBER_OF_MODEL_PARAMETERS], 
							bool Gender,
							unsigned int runID);
	CEModel(bool Mode,
			bool EVPPI,
			unsigned int number_of_PSA_trials,
			double DiscEff, 
			double DiscCost);


	double CalculateNMB(long unsigned int myWTP);
	void   CalculatePopulationEVPI(void);
	void   setEVPI(double evpi);

	void CleanUp(void);
	void completeReset(void);
	~CEModel(void);
};
