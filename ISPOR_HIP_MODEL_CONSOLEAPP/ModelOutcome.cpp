#pragma once
#include "ModelOutcome.h"

using namespace std; 

ModelOutcome::ModelOutcome(void)
{
	LYGs  = new double*[NUMBER_OF_TREATMENT_ARMS];
	QALYs = new double*[NUMBER_OF_TREATMENT_ARMS];
	Costs = new double*[NUMBER_OF_TREATMENT_ARMS];
	for(unsigned int treatment = 0; treatment < NUMBER_OF_TREATMENT_ARMS; treatment++){
		 LYGs[treatment] = new double[TIME_HORIZON];
		QALYs[treatment] = new double[TIME_HORIZON];
		Costs[treatment] = new double[TIME_HORIZON];	
		//Initializing the memory
		memset(LYGs[treatment] ,0,sizeof(double)*TIME_HORIZON);
		memset(QALYs[treatment],0,sizeof(double)*TIME_HORIZON);
		memset(Costs[treatment],0,sizeof(double)*TIME_HORIZON);			
	}
}
void ModelOutcome::Initialize(unsigned int NumberOfTrials)
{
	this->myNumberOfTrials =NumberOfTrials;
	
	Res = new double**[NumberOfTrials];
	for(unsigned int trial = 0; trial < NumberOfTrials; trial++){
		Res[trial] = new double*[NUMBER_OF_TREATMENT_ARMS];
		for(unsigned int treatment = 0; treatment < NUMBER_OF_TREATMENT_ARMS; treatment++){
			Res[trial][treatment] = new double[OUTCOME_RESULTS_DIMENSION1]; //Average Cost and QALYs for each treatment arm
			memset(Res[trial][treatment],0,sizeof(double)*OUTCOME_RESULTS_DIMENSION1);
		}
	}
}

void ModelOutcome::Reset(unsigned int NumberOfTrials){
	for(unsigned int trial = 0; trial < NumberOfTrials; trial++){
		for(unsigned int treatment = 0; treatment < NUMBER_OF_TREATMENT_ARMS; treatment++){
			memset(Res[trial][treatment],0,sizeof(double)*OUTCOME_RESULTS_DIMENSION1);
		}
	}
}

void ModelOutcome::CleanUp(void){
	for(unsigned int i=0;i<NUMBER_OF_TREATMENT_ARMS;i++){
		delete []  LYGs[i];
		delete [] QALYs[i];
		delete [] Costs[i];
	}
	delete []  LYGs;
	delete [] QALYs;
	delete [] Costs;

	for(unsigned int i=0;i< this->myNumberOfTrials;i++){
		for(unsigned int j=0;j<NUMBER_OF_TREATMENT_ARMS;j++){
			delete [] Res[i][j];
		}
		delete [] Res[i];
	}
	delete [] Res;
}

ModelOutcome::~ModelOutcome(void)
{

}
