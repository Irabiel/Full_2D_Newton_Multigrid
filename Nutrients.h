#ifndef NUTRIENTS_H_
#define NUTRIENTS_H_

#include "tools.h"
#include "Constants.h"
#include "Array.h"

template<typename T>
class Array2D;

struct Testprob
{
    
#if defined(COUPLED)
	Testprob()
	{
		nB = B0;
        nA = A0;
	}

	Testprob(double _B, double _A): nB(_B), nA(_A) {}

	double nB;
    double nA;
#else
	Testprob()
	{
		nB = B0;
	}

	Testprob(double _B): nB(_B) {}

	double nB;
#endif

};

#endif
