#ifndef _SIMULATION_H_clTabCtrl
#define _SIMULATION_H_

#include <stdio.h>

#include "Array.h"
#include "2D_MG.h"

struct OutputFiles;

class Test_prob
{
public:
    Test_prob(char* fname);
    ~Test_prob();
    void simulation(Multigrid& mg);
    void writeToFile(int OutputID);
	void passvalue(Multigrid& mg);
private:
    param P;
	ApproxArray2D* D;
};

#endif // _SIMULATION_H_
