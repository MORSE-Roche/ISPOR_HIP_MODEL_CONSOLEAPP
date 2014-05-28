#pragma once
#include "CEModel.h"
#include <iostream>
using namespace std;

int main(void){
	
	
	bool runningMode = PSA;
	bool EVPPImode = WITH_EVPPI;
	int number_of_PSA_trials =500;
	int patient_cohort_size =1000;
	unsigned int cRatio= 2200;
	double Discount_Eff = 0.035;
	double Discount_Cost= 0.035;
	unsigned int innerloopSize = 500;


	string ParameterNames[NUMBER_OF_MODEL_PARAMETERS];
	ParameterNames[omrPTHR] = "omrPTHR";
	ParameterNames[omrRTHR] = "omrRTHR";       
	ParameterNames[rrr]     = "rrr";
	ParameterNames[rrNP1]   = "rrNP1";
	ParameterNames[rrNP2]   = "rrNP2";
	ParameterNames[lambda]  = "lambda";
	ParameterNames[gamma]   = "gamma";
	ParameterNames[cPrimary]= "cPrimary";
	ParameterNames[cRevision]="cRevision";          
	ParameterNames[cSuccess] = "cSuccess";
	ParameterNames[uSuccessPrim] ="uSuccessPrim";
	ParameterNames[uRevision]    ="uRevision";
	ParameterNames[uSuccess2nd]  ="uSuccess2nd";
	
	bool Gender = FEMALE;									   
										   
	CEModel myCEModel =	CEModel(runningMode, EVPPImode, number_of_PSA_trials, Discount_Eff, Discount_Cost);
	//Initialize the EVPPI Parameter
	double* Parameters_Original = new double[NUMBER_OF_MODEL_PARAMETERS];
	memset(Parameters_Original, 0, sizeof(double)*NUMBER_OF_MODEL_PARAMETERS);
	double* Parameters = new double[NUMBER_OF_MODEL_PARAMETERS];
	memset(Parameters, 0, sizeof(double)*NUMBER_OF_MODEL_PARAMETERS);
	myCEModel.myStochasticParametersGenerator.GenerateParameters(Parameters_Original, Gender,  myCEModel.myModelInput.CovMat, myCEModel.myModelInput.Utilities, myCEModel.myModelInput.Risk, myCEModel.myModelInput.Hazard, myCEModel.myModelInput.Costs, DETERMINISTIC);
	
	//CALL MARKOV MODEL BASED - CALL MODEL PER TREATMENT DONE IN THE MODEL FUCTION
	myCEModel.ModelTreatmentArms(Parameters_Original, Gender, 0);
	
	//Calculate EVPI by Calculating NMB FOR EACH TREATMENT AND THE MAXIMUM NMB BASED UPON THE RESULTS AND THE RANGE OF WTP
	myCEModel.setEVPI(myCEModel.CalculateNMB(WTP));
	cout << "DETERMINISTIC results\n";
	cout << "Treatment Arm\tTotal Cost\tTotal QALY\n";
	for(unsigned int myTx =0 ;myTx < NUMBER_OF_TREATMENT_ARMS;myTx++){
		cout << myTx << "\t"<< myCEModel.myModelOutcome.Res[0][myTx][Outcome_Result_COST] << "\t" << myCEModel.myModelOutcome.Res[0][myTx][Outcome_Result_QALY] << endl;
	}


	//Reset the Model
	myCEModel.completeReset();
	
	cout << "\n\n\n";	
	if( runningMode == PSA){
		cout << "Probability Sensitivity Analysis\n";		
		for(unsigned int trialID = 0 ; trialID < myCEModel.getNumberOfTrials() ; trialID++){
			myCEModel.myStochasticParametersGenerator.GenerateParameters(Parameters, Gender,  myCEModel.myModelInput.CovMat, myCEModel.myModelInput.Utilities, myCEModel.myModelInput.Risk, myCEModel.myModelInput.Hazard, myCEModel.myModelInput.Costs, PSA);			
			
			//CALL MARKOV MODEL BASED - CALL MODEL PER TREATMENT DONE IN THE MODEL FUCTION
			myCEModel.ModelTreatmentArms(Parameters, Gender, trialID);
	
			//Reset the Model
			myCEModel.completeReset();	
		}
		//Calculate EVPI by Calculating NMB FOR EACH TREATMENT AND THE MAXIMUM NMB BASED UPON THE RESULTS AND THE RANGE OF WTP
		myCEModel.setEVPI(myCEModel.CalculateNMB(cRatio));
		cout << "cRatio = " << cRatio << endl; 
		cout << "EVPI = " << myCEModel.CalculateNMB(cRatio) << endl;
				
		//Calculate Population EVPI 
		myCEModel.CalculatePopulationEVPI();
		cout << "WTP\tEVPI\tPopulationEVPI\n";
		for(int i =0 ; i< CEAC_WTP_VECTOR_SIZE ;i++){
			for( int j =0 ; j<3 ; j++){
				cout << myCEModel.VoI_per_WTP[i][j] << "\t";
			}
			cout << "\n";
		}
		cout << "\n\n\n";
		if(	EVPPImode == WITH_EVPPI){
			//WITH EVPPI
			cout << "EVPPI Analysis\n";
			cout << "cRatio = " << cRatio << endl;
			//LOOP OVER ALL THE PARAMETER
			double Parameter_of_Interest =0;
			for(int myParameter = 0; myParameter < NUMBER_OF_MODEL_PARAMETERS; myParameter++){
				for(unsigned int innerloop = 0; innerloop <  innerloopSize; innerloop++){
					myCEModel.myStochasticParametersGenerator.GenerateParameters(Parameters , Gender, myCEModel.myModelInput.CovMat, myCEModel.myModelInput.Utilities, myCEModel.myModelInput.Risk, myCEModel.myModelInput.Hazard, myCEModel.myModelInput.Costs, PSA);
					Parameter_of_Interest = Parameters[myParameter];	
					#pragma omp parallel for num_threads(OMP_MAX_NO_THREADS)
					for( int trialID = 0 ; (unsigned int)(trialID) < myCEModel.getNumberOfTrials() ; trialID++){
						//CALL GENERATE INPUTS FUNCTION
						myCEModel.myStochasticParametersGenerator.GenerateParameters(Parameters , Gender, myCEModel.myModelInput.CovMat, myCEModel.myModelInput.Utilities, myCEModel.myModelInput.Risk, myCEModel.myModelInput.Hazard, myCEModel.myModelInput.Costs, PSA);
						
						//REPLACE THE PARAMETER OF INTEREST WITH THE CONSTANT VALUE FROM THE EVPPI parameter vector
						Parameters[myParameter] = Parameter_of_Interest;
						
						//CALL MARKOV MODEL BASED - CALL MODEL PER TREATMENT DONE IN THE MODEL FUCTION
						myCEModel.ModelTreatmentArms(Parameters, Gender, trialID);
					}
				}

				//Calculate Population EVPI 
				myCEModel.CalculatePopulationEVPI();		 
				double tmp_evppi = myCEModel.CalculateNMB(cRatio);
				tmp_evppi *= myCEModel.Cumulative_Effective_Population;
				cout << "EVPPI for the parameter " << ParameterNames[myParameter] << " = " << tmp_evppi << endl;
			
				//Reset the Model
				myCEModel.completeReset();
			}
		}
	}
	//Free the allocated memory
	delete [] Parameters_Original;
	delete [] Parameters;
	myCEModel.CleanUp();
	return 0;
}