#include <iostream>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "2D_MG.h"
#include "Simulation.h"


//
// #####################################
// ######## instruction ################
// #####################################
//
// Compile options:
// 1. Compile the code 
//    g++ *.cpp -O3
//    
// 2. Compile the code for coupled Newton Multigrid
//    g++ *.cpp -O3 -D COUPLED
//
// #####################################
//
// To run the code:
//    ./a.out in1.txt 	 	for mac or linux
//	  ./a.exe in1.txt 		for window 
// in1.txt is the name of the input file, 100 the number of files you would like the program to output.
// All the parameters are included in the input file, and can be changed.
// The name of the input file can also be changed.
//

using namespace std;


int main(int argc, char* argv[])
{   
#if defined(COUPLED)
    printf("Applying 2D Newton Multgrid to solve the non linear coupled Dirichlet system, \n-L(nB) + (nB)^2 + (nA)^2 = f \n-L(nA) + (nA)^2 + (nB)^2 = f \nWhere L is the laplacian operator \n");
#else
    printf("Applying 2D Newton Multgrid to solve the non linear Dirichlet system, \n-L(nB) + (nB)^2 = f \nWhere L is the laplacian operator \n");
#endif
    Test_prob Mod(argv[1]);
	Multigrid mg(argv[1]);
	Mod.simulation(mg);

    return 0;
}