#ifndef __2D_MG_H__
#define __2D_MG_H__

#include <cmath>
#include <iostream>
#include "Array.h"
#include "tools.h"
using namespace std;

struct Tensor2D;
struct Tensor2D2;
struct Tensor4D;
struct Tensor4D2;
class Multigrid;

struct Tensor2D
{
    Tensor2D(double v);
    Tensor2D(){resetZero();};
    double& operator()(int i, int j){return val[i][j];};
    const double& operator()(int i, int j) const{return val[i][j];};
    void Interpolate(const Tensor2D2& u, Array2D<Tensor2D>* v, int d1, int d2, double a);
    void resetZero();
    void checknan();
    double val[3][3];
};

struct Tensor2D2
{
    Tensor2D2(double v);
    Tensor2D2(){resetZero();};
    double& operator()(int i, int j){return val[i][j];};
    const double& operator()(int i, int j) const{return val[i][j];};
    void resetZero();
    void addMultiplyCopyVale(const Tensor2D & u, double m, int i, int j);
    void checknan();
    double val[7][7];
};

struct Tensor4D
{
    Tensor4D(double v);
    Tensor4D(){resetZero();};
    void resetZero();
    double& operator()(int i, int j, int l, int m){return val[i][j][l][m];};
    void Interpolate(const Tensor4D2& u, Array2D<Tensor2D>* v, int d1, int d2, double a);
    const double& operator()(int i, int j, int l, int m) const{return val[i][j][l][m];};
    void print(){
        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
                {
                    for (int l=0; l<2; l++)
                    {
                        for (int m=0; m<2; m++)
                        {
                            cout<<val[i][j][l][m]<<" ";
                        }
                        cout<<endl;
                    }
                    cout<<endl;
                }
        }
    };
    double val[3][3][2][2];
};

struct Tensor4D2
{
    Tensor4D2(double v);
    Tensor4D2(){resetZero();};
    void resetZero();
    void addMultiplyCopyVale(const Tensor4D & u, double m, int i, int j);
    double& operator()(int i, int j, int l, int m){return val[i][j][l][m];};
    const double& operator()(int i, int j, int l, int m) const{return val[i][j][l][m];};
    void checknan();
    double val[7][7][2][2];
};

class Multigrid
{
public:
    Multigrid(char* fname);
    ~Multigrid();
    void setRestrInterp();
    void setF();
    void setOperators();
    void setOperatorsNextLevel();
    void resetZero();
    void relaxationAtLevel(int l);
    void relaxationElementwise(int l, int i, int j, int nx, int ny);
    void residueRestriction(int l);
    void restrictionElementwise1(int l, int i, int j, int nx, int ny);
    void restrictionElementwise2(int l, int i, int j, int nx, int ny);
    void interpolation(int l);
    void interpolationElementwise(int l, int i, int j, int nx, int ny);
    void updateSol();
    void updateSolElementwise(int i, int j);
    void setup();
    void Vcycle(int mu1, int mu2);
    void solve(int mu1, int mu2);
    void reset();
    double evalUpdateErr();
    double evalUpdateRes();
    ApproxArray2D** Sol; // solution vector 
    ApproxArray2D** SolRes; // risdual vector 
    double ResVcycle();
private:
	param P;
    // private variables.
    ApproxArray2D** SolResTmp; // for red and black GS
    ApproxArray2D** SolErr; // error vector
#if defined(COUPLED)
    Array2D<Tensor4D>** SolStencil; // A matrix stencil
#else
    Array2D<Tensor2D>** SolStencil; // A matrix stencil
#endif
    
    Array2D<Tensor2D>** SolRestrict; // Restriction operator 
    Array2D<Tensor2D>** SolInterpolate; // interpolating operator 
    double hx;
    double hy;
    
    double fB(double u) {return u*u;};
    double dfB(double u) {return 2*u;};
    
    double fA(double v) {return v*v;};
    double dfA(double v) {return 2*v;};
    
    double fij(double hx, double hy, double i, double j) {return (-2*( (i*hx)*((i*hx) - 1) + (j*hy)*((j*hy) - 1) ) + 2*(i*hx)*((i*hx) - 1)*(j*hy)*((j*hy) - 1)*(i*hx)*((i*hx) - 1)*(j*hy)*((j*hy) - 1));};
    
    // private member functions.
    int ipow(int base, int exp);
};


#endif
