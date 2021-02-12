#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// default values for all variables

// time
extern double t_max;
extern double initial_dt;
extern double OutputTime;
extern double t0;
extern double UpdateTime;

// non-physical constants
extern int BoxX;
extern int BoxY;
extern int levels;
extern double B0;
extern double A0;

// directory name
extern char DirName[500];

// Restart information
extern int restart;
extern int restartIndex;
extern char restartDir[500];


#endif /* CONSTANTS_H_ */
