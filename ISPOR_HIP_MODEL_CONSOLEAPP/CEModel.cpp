#pragma once
#include "CEModel.h"
#include <iostream>


using namespace std;

CEModel::CEModel(bool Mode,
				 bool EVPPI,
				 unsigned int number_of_PSA_trials, 
				 double DiscEff, 
				 double DiscCost){

	
	omp_set_num_threads(OMP_MAX_NO_THREADS);

	this->runningmode = Mode;
	this->EVPPImode   = EVPPI;
	this->EVPI = 0;
	//Initialize the Model Outcome sub-class
	if(runningmode == DETERMINISTIC ){
		this->Number_of_Trials = 1;
	}
	else{
		this->Number_of_Trials = number_of_PSA_trials;
	}
	this->myModelOutcome.Initialize(this->Number_of_Trials);
	
	for(int cycle = 0; cycle < TIME_HORIZON; cycle++){
			DiscountRates[Discountrate_Effect][cycle] = 1 /pow((1+ DiscEff ),cycle);
			DiscountRates[Discountrate_Cost][cycle]   = 1 /pow((1+ DiscCost),cycle);
	}

	//Initialize the NMB_RESULTS array
	NMB_RESULTS = new double*[this->Number_of_Trials];
	this -> NMB_SIZE1 = NUMBER_OF_TREATMENT_ARMS + 1; // NMB for each treatment arm as well as MaxNMB
	
	for(unsigned int RunID = 0; RunID < this->Number_of_Trials; RunID++){
			NMB_RESULTS[RunID] = new double[this->NMB_SIZE1];
			memset(NMB_RESULTS[RunID],0,this->NMB_SIZE1*sizeof(double));
	}	
	
	
	//Initialize Patient cohort
	Cohort = new double**[NUMBER_OF_TREATMENT_ARMS];
	for( int treatment = 0; treatment < NUMBER_OF_TREATMENT_ARMS; treatment++){
		Cohort[treatment] = new double*[TIME_HORIZON];
		for( int cycle = 0; cycle < TIME_HORIZON; cycle++){
			Cohort[treatment][cycle] = new double[NumHealthStates];
			memset(Cohort[treatment][cycle], 0, NumHealthStates*sizeof(double));
		}
	}

	//Calculating the cumulative effective population size
	double Discountrate_effPopulation= 0.06; //We consider the annual discount rate of effective population of 6% adopting the technology
	this->Cumulative_Effective_Population=0;
	for(unsigned int year=1; year <= EFFECTIVE_TECHNOLOGY_LIFE_YEARS; year++){
		this->Cumulative_Effective_Population += EFFECTIVE_POPULATION_SIZE_ADOPTING_TECHNOLOGY/(pow(1+Discountrate_effPopulation,year-1));
	}


	//Initialize VoI_per_WTP 
	VoI_per_WTP = new double*[CEAC_WTP_VECTOR_SIZE];
	unsigned int myWTP=0;
	unsigned int increment = 0;
	for(int i = 0; i < CEAC_WTP_VECTOR_SIZE ; i++){
		VoI_per_WTP[i] = new double[4];// WTP, EVPI, popEVPI, EVPPI
		myWTP += increment;
		if(i <= 9){
			increment = 100;
		}
		else if(i > 9 && i <= 17){
			increment = 500;
		}
		else if(i > 17 && i <= 42){
			increment = 1000;
		}
		else if(i > 42 && i <= 57){
			increment = 5000;
		}
		VoI_per_WTP[i][0] = myWTP;
		VoI_per_WTP[i][1] = 0;
		VoI_per_WTP[i][2] = 0;
		VoI_per_WTP[i][3] = 0;
	}

	this->VoI_RESULTS = new double[this->NMB_SIZE1];
	memset(this->VoI_RESULTS,0,sizeof(double)*this->NMB_SIZE1);
}

unsigned int CEModel::getNumberOfTrials(){
	return this->Number_of_Trials;
}

void CEModel::setEVPI(double evpi){
	this->EVPI = evpi;
}

double CEModel::CalculateMortalityRisk(bool Gender, unsigned int cycle)
{

	long unsigned int i = (long unsigned)(ceil((double)((DEFAULT_AGE -1 + cycle - BaseCaseAge)/AgeIntervals)))+1;
	if(i>=5){i=5;}

	if(Gender == FEMALE){
		return this->myModelInput.MortalityRates[i][Mortality_RATES_FEMALE];
	}
	else{
		return this->myModelInput.MortalityRates[i][Mortality_RATES_MALE];
	}
}

void CEModel::ModelTreatmentArms(double Parameters[NUMBER_OF_MODEL_PARAMETERS], bool Gender, unsigned int trialID)
{	
	// Treatment effect
	double txEffect_Array[NUMBER_OF_TREATMENT_ARMS];
	txEffect_Array[STD] = 1;
	txEffect_Array[TxNP1] = Parameters[rrNP1];
	txEffect_Array[TxNP2] = Parameters[rrNP2];
	
	for( int myTx = 0 ; myTx < NUMBER_OF_TREATMENT_ARMS ; myTx ++){
		//reset the Cohort matrix
		double txEffect = txEffect_Array[myTx];
		this->reset(myTx);
		
		//PrimaryTHR
		Cohort[myTx][0][PrimaryTHR]        =  CohortSize;	
		for(unsigned int cycle = 1 ; cycle < TIME_HORIZON; cycle++){
			//SuccessP
			Cohort[myTx][cycle][SuccessP]   =	Cohort[myTx][cycle - 1][PrimaryTHR]  * (1 - Parameters[omrPTHR]) + 
													Cohort[myTx][cycle - 1][SuccessP]    * (1 - RevisionRisk(cycle - 1,  txEffect, Parameters[lambda], Parameters[gamma]) - CalculateMortalityRisk(Gender, cycle));
			//RevisionTHR
			Cohort[myTx][cycle][RevisionTHR]=	Cohort[myTx][cycle - 1][SuccessP]    * RevisionRisk(cycle - 1, txEffect, Parameters[lambda], Parameters[gamma]) + 
													Cohort[myTx][cycle - 1][SuccessR]    * Parameters[rrr];
			//SuccessR
			Cohort[myTx][cycle][SuccessR]   =	Cohort[myTx][cycle - 1][RevisionTHR] * (1 - CalculateMortalityRisk(Gender, cycle) - Parameters[omrRTHR]) + 
													Cohort[myTx][cycle - 1][SuccessR]    * (1 - CalculateMortalityRisk(Gender, cycle) - Parameters[rrr]);
			//Death
			Cohort[myTx][cycle][Death]      =	Cohort[myTx][cycle - 1][PrimaryTHR]  * Parameters[omrPTHR] + 
													Cohort[myTx][cycle - 1][SuccessP]    * CalculateMortalityRisk(Gender,cycle) + 
													Cohort[myTx][cycle - 1][RevisionTHR] * Parameters[omrRTHR] + 
													Cohort[myTx][cycle - 1][RevisionTHR] * CalculateMortalityRisk(Gender,cycle) + 
													Cohort[myTx][cycle - 1][SuccessR]    * CalculateMortalityRisk(Gender, cycle) + Cohort[myTx][cycle - 1][Death];
		}
			
		//ADD COST FOR THE FIRST PROCEDURE
		myModelOutcome.Costs[myTx][0] = this->myModelInput.TxCosts[myTx]*CohortSize;
		for(unsigned int cycle = 1; cycle < TIME_HORIZON; cycle++){
			//Don't look at the Death health state since LYG,QALY and cost is zero for this health state
			for(unsigned int HealthState =0; HealthState<NumHealthStates-1;HealthState++){
				myModelOutcome.LYGs[myTx][cycle]  += Cohort[myTx][cycle][HealthState] ;//* DiscountRates[Discountrate_Effect][cycle];
				myModelOutcome.QALYs[myTx][cycle] += Cohort[myTx][cycle][HealthState] * this->myModelInput.Input_per_HS[UTILIY_VALUES_PER_HEALTH_STATE][HealthState] * DiscountRates[Discountrate_Effect][cycle];
				myModelOutcome.Costs[myTx][cycle] += Cohort[myTx][cycle][HealthState] * this->myModelInput.Input_per_HS[SURGERY_COSTS_PER_HEALTH_STATE][HealthState] * DiscountRates[Discountrate_Cost][cycle];
			}
			myModelOutcome.Res[trialID][myTx][Outcome_Result_COST] += myModelOutcome.Costs[myTx][cycle];
			myModelOutcome.Res[trialID][myTx][Outcome_Result_QALY] += myModelOutcome.QALYs[myTx][cycle];
		}

		myModelOutcome.Res[trialID][myTx][Outcome_Result_COST] = (myModelOutcome.Res[trialID][myTx][Outcome_Result_COST]+myModelOutcome.Costs[myTx][0])/CohortSize;
		myModelOutcome.Res[trialID][myTx][Outcome_Result_QALY] /= CohortSize;
	}
}


double CEModel::RevisionRisk(unsigned int cycle, double RelativeRisk, double Lambda, double Gamma)
{
	return (1 - gsl_sf_exp(Lambda * RelativeRisk * (pow(cycle, Gamma) - pow((cycle + 1),  Gamma))));
}

void CEModel::reset(unsigned int myTx){

	//Resetting the memory
		for( int cycle = 0; cycle < TIME_HORIZON; cycle++){
			memset(this->Cohort[myTx][cycle], 0, sizeof(double)*NumHealthStates);
		}
		memset(this->myModelOutcome.LYGs[myTx] , 0 , sizeof(double)*TIME_HORIZON);
		memset(this->myModelOutcome.QALYs[myTx], 0 , sizeof(double)*TIME_HORIZON);
		memset(this->myModelOutcome.Costs[myTx], 0 , sizeof(double)*TIME_HORIZON);
}

void CEModel::completeReset(void){
	//Resetting the memory
	for( int treatment = 0; treatment < NUMBER_OF_TREATMENT_ARMS; treatment++){
		for( int cycle = 0; cycle < TIME_HORIZON; cycle++){
			memset(this->Cohort[treatment][cycle], 0, sizeof(double)*NumHealthStates);
		}
	}
	
	for(int treatment = 0; treatment < NUMBER_OF_TREATMENT_ARMS; treatment++){
		memset(myModelOutcome.LYGs[treatment],0,sizeof(double)*TIME_HORIZON);
		memset(myModelOutcome.QALYs[treatment],0,sizeof(double)*TIME_HORIZON);
		memset(myModelOutcome.Costs[treatment],0,sizeof(double)*TIME_HORIZON);			
	}
}

double CEModel::CalculateNMB(long unsigned int myWTP){

	memset(this->VoI_RESULTS,0,sizeof(double)*this->NMB_SIZE1);
	double max =0;

	//CALCULATE NMB PER TREATMENT per trial
	for(unsigned int runID=0; runID < this->Number_of_Trials;runID++){
		max =0;
		for(unsigned int treatment=0; treatment < NUMBER_OF_TREATMENT_ARMS;treatment++){
			NMB_RESULTS[runID][treatment] = this->myModelOutcome.Res[runID][treatment][Outcome_Result_QALY]*myWTP - this->myModelOutcome.Res[runID][treatment][Outcome_Result_COST];
			VoI_RESULTS[treatment] += NMB_RESULTS[runID][treatment];
			//find the maximum NMB per trial
			if(max<= NMB_RESULTS[runID][treatment]){max=NMB_RESULTS[runID][treatment];}
		}
		NMB_RESULTS[runID][NUMBER_OF_TREATMENT_ARMS]=max;
		VoI_RESULTS[NUMBER_OF_TREATMENT_ARMS] += max;
	}


	//CALCULATE EVPI
	max =0; //Max of the Mean NMBs
	for(unsigned int treatment=0; treatment < NUMBER_OF_TREATMENT_ARMS;treatment++){
		VoI_RESULTS[treatment] /= this->Number_of_Trials;
		if(max<=VoI_RESULTS[treatment]){max = VoI_RESULTS[treatment];}
	}

	VoI_RESULTS[NUMBER_OF_TREATMENT_ARMS] /= this->Number_of_Trials; //Mean of the Max NMB 
	return VoI_RESULTS[NUMBER_OF_TREATMENT_ARMS]-max;
}

void CEModel::CalculatePopulationEVPI(void){
	for(int i = 0; i < CEAC_WTP_VECTOR_SIZE ; i++){
		long unsigned int myWTP = (long unsigned int)(VoI_per_WTP[i][0]);
		double myEVPI = CalculateNMB(myWTP);
		VoI_per_WTP[i][1] = myEVPI;
		VoI_per_WTP[i][2] = myEVPI*this->Cumulative_Effective_Population;
	}
}

void CEModel::CleanUp(void){
	
	myModelOutcome.CleanUp();

	//free Cohort multi-dimentional array
	for(unsigned int i=0;i< NUMBER_OF_TREATMENT_ARMS;i++){
		for(unsigned int j=0;j<TIME_HORIZON;j++){
			delete [] Cohort[i][j];
		}
		delete [] Cohort[i];
	}
	delete [] Cohort;

	//free VoI_per_WTP array
	for(unsigned int i=0 ; i < CEAC_WTP_VECTOR_SIZE; i++){
		delete []  VoI_per_WTP[i];
	}
	delete [] VoI_per_WTP;

	//free  VoI_RESULTS array
	delete [] VoI_RESULTS;


	//free NMB_RESULTS array
	for(unsigned int i=0 ; i < this->Number_of_Trials; i++){
		delete []  NMB_RESULTS[i];
	}
	delete [] NMB_RESULTS;
}

CEModel::~CEModel(void)
{
}		


