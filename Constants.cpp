#include "Constants.h"

// time
double t_max = 20.0; // total simulation time
double initial_dt = 0.00001;	// initial time step
double OutputTime = 100*initial_dt;		// how often to output
double t0 = 0.0;	// time at start of simulation
double UpdateTime = 100*initial_dt;

// non-physical constants
int BoxX = 150;	// number of boxes per side of the grid
int BoxY = 150;
int maxLevelsMG = 1;
double BoxLength = 4;
int FilterLen = 5;

// directory name
char DirName[500] = "";

// Restart information
int restart = 0;
int restartIndex = 0;
char restartDir[500] = "";


// nutrient constants

int minIter = 500;
int maxIter = 20000;

// Boundary condition

double B0 = 0.0;
double A0 = 0.0;
