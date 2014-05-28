#pragma once
#include "StochasticParameters.h"

class ModelOutcome
{
private:
	unsigned int myNumberOfTrials;

public:

	double**    LYGs;
	double**   QALYs;
	double**   Costs;
	double***    Res;

	void Initialize(unsigned int NumberOfTrials);
	void Reset(unsigned int NumberOfTrials);
	void CleanUp(void);
    ModelOutcome(void);
	~ModelOutcome(void);
};


