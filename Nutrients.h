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
		U = 0.0;
        V = 0.0;
	}

	Testprob(double _u, double _v): U(_u), V(_v) {}

	double U;
    double V;
#else
	Testprob()
	{
		Approx = 0.0;
	}

	Testprob(double _v): Approx(_v) {}

	double Approx;
#endif

};

#endif
