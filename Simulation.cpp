#include "Simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


Test_prob::Test_prob(char* fname)
{
    P = ReadParameters(fname);
	D = new ApproxArray2D(P.BoxX+1,P.BoxY+1);
};

Test_prob::~Test_prob()
{
    delete D;
};

void Test_prob::passvalue(Multigrid& mg)
{
	for (int i =0; i < P.BoxX; i++)
    {
        for (int j =0; j < P.BoxX; j++)
        {
#if defined(COUPLED)
            D[0](i,j).U = mg.Sol[0][0](i,j).U;
            D[0](i,j).V = mg.Sol[0][0](i,j).V;
#else
            D[0](i,j).Approx = mg.Sol[0][0](i,j).Approx;
#endif
        }
    }
};

void Test_prob::simulation(Multigrid& mg)
{
    double cycle = 0;
    int outputID = 0;
	double er;
//    printf("Consol output records the following: Output step; time; num events; total weight; CellType1 Num; CellType2 Num; Total Cell Num\n");

    do {
		mg.solve(3, 20);
		passvalue(mg);//*******************
		writeToFile(outputID);
		//cout<< "Cycle = " <<cycle<<", Error = " << mg.evalUpdateErr() << ", Residual = " << mg.evalUpdateRes() <<"\n";
//		printf("%d %.4e\n", outputID, mg.ResVcycle());
		fflush(stdout);
		outputID++;
		cycle++;
    } while (cycle<P.Cycles);
};

void Test_prob::writeToFile(int OutputID)
{
    char fileName[500];
    FILE* out;
    strcpy(fileName,P.DirName);
    mkDir(fileName);
    sprintf(fileName,"%s/%d.txt",P.DirName,OutputID);
    out = fopen(fileName, "w");
    D->Output(out);
    fclose(out);
};