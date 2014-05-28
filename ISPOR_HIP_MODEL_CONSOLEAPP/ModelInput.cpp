#include "ModelInput.h"


ModelInput::ModelInput(void)
{
	MortalityRates[AGE_INDEX_35_44][Mortality_RATES_FEMALE] = 0.00099;	
	MortalityRates[AGE_INDEX_45_54][Mortality_RATES_FEMALE] = 0.00260;
	MortalityRates[AGE_INDEX_55_64][Mortality_RATES_FEMALE] = 0.00670;
	MortalityRates[AGE_INDEX_65_74][Mortality_RATES_FEMALE] = 0.01930;
	MortalityRates[AGE_INDEX_75_84][Mortality_RATES_FEMALE] = 0.05350;
	MortalityRates[AGE_INDEX_85][Mortality_RATES_FEMALE]    = 0.15480;
	
	MortalityRates[AGE_INDEX_35_44][Mortality_RATES_MALE]   = 0.00151;
	MortalityRates[AGE_INDEX_45_54][Mortality_RATES_MALE]   = 0.00393;
	MortalityRates[AGE_INDEX_55_64][Mortality_RATES_MALE]   = 0.01090;
	MortalityRates[AGE_INDEX_65_74][Mortality_RATES_MALE]   = 0.03160;
	MortalityRates[AGE_INDEX_75_84][Mortality_RATES_MALE]   = 0.08010;  
	MortalityRates[AGE_INDEX_85][Mortality_RATES_MALE]      = 0.18790;
//////////////////////////////////////////////////////////////////////////////////
	Risk[omrPTHR][DIST_ALPHA]= 2;
	Risk[omrPTHR][DIST_BETA] = 98;
	Risk[omrPTHR][DIST_MEAN] = Risk[omrPTHR][DIST_ALPHA]/(Risk[omrPTHR][DIST_ALPHA]+Risk[omrPTHR][DIST_BETA]); 
	Risk[omrPTHR][DIST_SE]   = sqrt(Risk[omrPTHR][DIST_ALPHA]*Risk[omrPTHR][DIST_BETA]/(pow((double)(Risk[omrPTHR][DIST_ALPHA]+Risk[omrPTHR][DIST_BETA]),(int)(2))*(Risk[omrPTHR][DIST_ALPHA]+Risk[omrPTHR][DIST_BETA]+1)));
	

	Risk[omrRTHR][DIST_ALPHA]= 2;
	Risk[omrRTHR][DIST_BETA] = 98;
	Risk[omrRTHR][DIST_MEAN] = Risk[omrRTHR][DIST_ALPHA]/(Risk[omrRTHR][DIST_ALPHA]+Risk[omrRTHR][DIST_BETA]); 
	Risk[omrRTHR][DIST_SE]   = sqrt(Risk[omrRTHR][DIST_ALPHA]*Risk[omrRTHR][DIST_BETA]/(pow((double)(Risk[omrRTHR][DIST_ALPHA]+Risk[omrRTHR][DIST_BETA]),(int)(2))*(Risk[omrRTHR][DIST_ALPHA]+Risk[omrRTHR][DIST_BETA]+1)));
	

	Risk[rrr][DIST_ALPHA]    = 4;
	Risk[rrr][DIST_BETA]     = 96;
	Risk[rrr][DIST_MEAN]     = Risk[rrr][DIST_ALPHA]/(Risk[rrr][DIST_ALPHA]+Risk[rrr][DIST_BETA]);
	Risk[rrr][DIST_SE]       = sqrt(Risk[rrr][DIST_ALPHA]*Risk[rrr][DIST_BETA]/(pow((double)(Risk[rrr][DIST_ALPHA]+Risk[rrr][DIST_BETA]),(int)(2))*(Risk[rrr][DIST_ALPHA]+Risk[rrr][DIST_BETA]+1)));
//////////////////////////////////////////////////////////////////////////////////
	Hazard[lngamma][DIST_MEAN]= 0.3740968;
	Hazard[lngamma][DIST_SE]  = 0.0474501;
	Hazard[lngamma][HAZARD_RATIO]        = exp(Hazard[lngamma][DIST_MEAN]);

	Hazard[consta][DIST_MEAN] = -5.490935;
	Hazard[consta][DIST_SE]   =  0.207892;
	Hazard[consta][HAZARD_RATIO]         = exp(Hazard[consta][DIST_MEAN]);

	Hazard[Age][DIST_MEAN]    = -0.0367022;
	Hazard[Age][DIST_SE]      =  0.0052112;
	Hazard[Age][HAZARD_RATIO]            = exp(Hazard[Age][DIST_MEAN]);

	Hazard[Male][DIST_MEAN]   = 0.768536;
	Hazard[Male][DIST_SE]     = 0.109066;
	Hazard[Male][HAZARD_RATIO]           = exp(Hazard[Male][DIST_MEAN]);

	Hazard[NP1][DIST_MEAN]    = -1.344474;
	Hazard[NP1][DIST_SE]      = 0.3825815;
	Hazard[NP1][HAZARD_RATIO]            = exp(Hazard[NP1][DIST_MEAN]);

	Hazard[NP2][DIST_MEAN]    = -1.6687;
	Hazard[NP2][DIST_SE]      = 0.517653;
	Hazard[NP2][HAZARD_RATIO]            = exp(Hazard[NP2][DIST_MEAN]);
//////////////////////////////////////////////////////////////////////////////////
	Utilities[uSuccessP][DIST_MEAN] = 0.85;
	Utilities[uSuccessP][DIST_SE]   = 0.030;
	Utilities[uSuccessP][DIST_ALPHA]= Utilities[uSuccessP][DIST_MEAN]*(Utilities[uSuccessP][DIST_MEAN]*(1-Utilities[uSuccessP][DIST_MEAN])/(pow(Utilities[uSuccessP][DIST_SE],2))-1);
	Utilities[uSuccessP][DIST_BETA] = Utilities[uSuccessP][DIST_MEAN]*(1-Utilities[uSuccessP][DIST_MEAN])/(pow(Utilities[uSuccessP][DIST_SE],2))-1-Utilities[uSuccessP][DIST_ALPHA];

	Utilities[uRev][DIST_MEAN] = 0.30;
	Utilities[uRev][DIST_SE]   = 0.030;
	Utilities[uRev][DIST_ALPHA]= Utilities[uRev][DIST_MEAN]*(Utilities[uRev][DIST_MEAN]*(1-Utilities[uRev][DIST_MEAN])/(pow(Utilities[uRev][DIST_SE],2))-1);
	Utilities[uRev][DIST_BETA] = Utilities[uRev][DIST_MEAN]*(1-Utilities[uRev][DIST_MEAN])/(pow(Utilities[uRev][DIST_SE],2))-1-Utilities[uRev][DIST_ALPHA];

	Utilities[uSuccessR][DIST_MEAN] = 0.75;
	Utilities[uSuccessR][DIST_SE]   = 0.040;
	Utilities[uSuccessR][DIST_ALPHA]= Utilities[uSuccessR][DIST_MEAN]*(Utilities[uSuccessR][DIST_MEAN]*(1-Utilities[uSuccessR][DIST_MEAN])/(pow(Utilities[uSuccessR][DIST_SE],2))-1);
	Utilities[uSuccessR][DIST_BETA] = Utilities[uSuccessR][DIST_MEAN]*(1-Utilities[uSuccessR][DIST_MEAN])/(pow(Utilities[uSuccessR][DIST_SE],2))-1-Utilities[uSuccessR][DIST_ALPHA];
//////////////////////////////////////////////////////////////////////////////////
	CovMat[lngamma][lngamma] = 0.002251511990010;
	CovMat[lngamma][consta]  =-0.005691000000000; 
    CovMat[lngamma][Age]     = 0.000000028000000;
	CovMat[lngamma][Male]    = 0.000005100000000;
    CovMat[lngamma][NP1]     = 0.000259000000000;
	CovMat[lngamma][NP2]     = 0.000351000000000;

    CovMat[consta][lngamma]  =-0.005691000000000;
	CovMat[consta][consta]   = 0.043219083664000;
    CovMat[consta][Age]      =-0.000783000000000;
	CovMat[consta][Male]     =-0.007247000000000;
	CovMat[consta][NP1]      =-0.000642000000000;
	CovMat[consta][NP2]      =-0.000537000000000;

	CovMat[Age][lngamma]     = 0.000000028000000;
	CovMat[Age][consta]	     =-0.000783000000000;
	CovMat[Age][Age]         = 0.000027156605440;
	CovMat[Age][Male]        = 0.000033000000000;
	CovMat[Age][NP1]         =-0.000111000000000;
	CovMat[Age][NP2]         =-0.000299000000000;

	CovMat[Male][lngamma]    = 0.000005100000000;
	CovMat[Male][consta]     =-0.007247000000000;
	CovMat[Male][Age]        = 0.000033000000000;
	CovMat[Male][Male]       = 0.011895392356000;
	CovMat[Male][NP1]        = 0.000184000000000;
	CovMat[Male][NP2]        = 0.000098000000000;

	CovMat[NP1][lngamma]     = 0.000259000000000;
	CovMat[NP1][consta]      =-0.000642000000000;
	CovMat[NP1][Age]         =-0.000111000000000;
	CovMat[NP1][Male]        = 0.000184000000000;
	CovMat[NP1][NP1]         = 0.146368604142250;
	CovMat[NP1][NP2]         = 0.000354680000000;

	CovMat[NP2][lngamma]     = 0.000351000000000;
	CovMat[NP2][consta]      =-0.000537000000000;
	CovMat[NP2][Age]         =-0.000299000000000;
	CovMat[NP2][Male]        = 0.000098000000000;
	CovMat[NP2][NP1]         = 0.000354680000000;
	CovMat[NP2][NP2]         = 0.267964628409000;
//////////////////////////////////////////////////////////////////////////////////
	TxCosts[STD]   = 394;
	TxCosts[TxNP1] = 579;
	TxCosts[TxNP2] = 788;
//////////////////////////////////////////////////////////////////////////////////
	Costs[STD][DIST_MEAN] = 0;
	Costs[STD][DIST_SE]   = 0;
	Costs[STD][DIST_ALPHA]= 0;
	Costs[STD][DIST_BETA] = 0;
	
	Costs[TxNP1][DIST_MEAN] = 5294;
	Costs[TxNP1][DIST_SE]   = 0;
	Costs[TxNP1][DIST_ALPHA]= 12.67;
	Costs[TxNP1][DIST_BETA] = 417.67;

	Costs[TxNP2][DIST_MEAN] = 0;
	Costs[TxNP2][DIST_SE]   = 0;
	Costs[TxNP2][DIST_ALPHA]= 0;
	Costs[TxNP2][DIST_BETA] = 0;
//////////////////////////////////////////////////////////////////////////////////
	Input_per_HS[0][0]=1;
	Input_per_HS[0][1]=1;
	Input_per_HS[0][2]=1;
	Input_per_HS[0][3]=1;
	Input_per_HS[0][4]=1;
	
	Input_per_HS[1][0]=1;
	Input_per_HS[1][1]=1;
	Input_per_HS[1][2]=1;
	Input_per_HS[1][3]=1;
	Input_per_HS[1][4]=0;

	Input_per_HS[UTILIY_VALUES_PER_HEALTH_STATE][PrimaryTHR]=0;
	Input_per_HS[UTILIY_VALUES_PER_HEALTH_STATE][SuccessP]=0.85;
	Input_per_HS[UTILIY_VALUES_PER_HEALTH_STATE][RevisionTHR]=0.3;
	Input_per_HS[UTILIY_VALUES_PER_HEALTH_STATE][SuccessR]=0.75;
	Input_per_HS[UTILIY_VALUES_PER_HEALTH_STATE][Death]=0;

	Input_per_HS[SURGERY_COSTS_PER_HEALTH_STATE][PrimaryTHR]=0;
	Input_per_HS[SURGERY_COSTS_PER_HEALTH_STATE][SuccessP]=0;
	Input_per_HS[SURGERY_COSTS_PER_HEALTH_STATE][RevisionTHR]=5294;
	Input_per_HS[SURGERY_COSTS_PER_HEALTH_STATE][SuccessR]=0;
	Input_per_HS[SURGERY_COSTS_PER_HEALTH_STATE][Death]=0;	
}

ModelInput::~ModelInput(void)
{
}
