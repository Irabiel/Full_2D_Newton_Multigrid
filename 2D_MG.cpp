#include "2D_MG.h"
#include <cmath>

#include <omp.h>

Tensor2D::Tensor2D(double v)
{
    for (int i1 = 0;i1<3;i1++)
        for (int j1 = 0;j1<3;j1++)
            val[i1][j1] = v;
};

void Tensor2D::checknan()
{
    for (int i1 = 0;i1<3;i1++)
        for (int j1 = 0;j1<3;j1++)
            if (isnan(val[i1][j1]))
                cout<<i1<< " " <<j1<<endl;
};

void Tensor2D::Interpolate(const Tensor2D2& u, Array2D<Tensor2D>* v, int d1, int d2, double a)
{
    for (int i=0; i<3; i++)
        for (int j = 0;j<3;j++)
            {
                (*this)(i,j) = 0;
                for (int i3=0; i3<3; i3++)
                    for (int j3 = 0;j3<3;j3++)
                        (*this)(i,j) += u(2*i+1+i3-1,2*j+1+j3-1)*v[0](d1+i-1,d2+j-1)(i3,j3)*a;
            }
}

void Tensor2D::resetZero()
{
    for (int i1 = 0;i1<3;i1++)
        for (int j1 = 0;j1<3;j1++)
            val[i1][j1] = 0;
};

Tensor2D2::Tensor2D2(double v)
{
    for (int i1 = 0;i1<7;i1++)
        for (int j1 = 0;j1<7;j1++)
            val[i1][j1] = v;
};

void Tensor2D2::checknan()
{
    for (int i1 = 0;i1<7;i1++)
        for (int j1 = 0;j1<7;j1++)
        {        
            if (isnan(val[i1][j1]))
                cout<<i1<<" " << j1<<endl;
        }
};

void Tensor2D2::resetZero()
{
    for (int i1 = 0;i1<7;i1++)
        for (int j1 = 0;j1<7;j1++)
            val[i1][j1] = 0;
};

void Tensor2D2::addMultiplyCopyVale(const Tensor2D & u, double m, int i, int j)
{
    for (int i1 = -1;i1<2;i1++)
        for (int j1 = -1;j1<2;j1++)
            val[i+i1][j+j1] += u(i1+1,j1+1)*m;
};


Tensor4D::Tensor4D(double v)
{
    for (int i1 = 0;i1<3;i1++)
        for (int i2 = 0;i2<3;i2++)
            for (int i3 = 0;i3<2;i3++)
                for (int i4 = 0;i4<2;i4++)
                    val[i1][i2][i3][i4] = v;
};

void Tensor4D::resetZero()
{
    for (int i1 = 0;i1<3;i1++)
        for (int i2 = 0;i2<3;i2++)
            for (int i3 = 0;i3<2;i3++)
                for (int i4 = 0;i4<2;i4++)
                    val[i1][i2][i3][i4] = 0;
};

void Tensor4D::Interpolate(const Tensor4D2& u, Array2D<Tensor2D>* v, int d1, int d2, double a)
{
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            for (int i1=0; i1<2; i1++)
                for (int i2=0; i2<2; i2++)
                {
                    (*this)(i,j,i1,i2) = 0;
                    for (int i3=0; i3<3; i3++)
                        for (int i4=0; i4<3; i4++)
                            (*this)(i,j,i1,i2) += u(2*i+1+i3-1, 2*j+1+i4-1,i1,i2)*v[0](d1+i-1,d2+j-1)(i3,i4)*a;
                }
};

Tensor4D2::Tensor4D2(double v)
{
    for (int i1 = 0;i1<7;i1++)
        for (int i2 = 0;i2<7;i2++)
            for (int i4 = 0;i4<2;i4++)
                for (int i5 = 0;i5<2;i5++)
                    val[i1][i2][i4][i5] = v;
};

void Tensor4D2::resetZero()
{
    for (int i1 = 0;i1<7;i1++)
    {
        for (int i2 = 0;i2<7;i2++)
        {
            for (int i4 = 0;i4<2;i4++)
            {
                for (int i5 = 0;i5<2;i5++)
                {
                    val[i1][i2][i4][i5] = 0;
                }
            }
        }
    }
};
    
void Tensor4D2::addMultiplyCopyVale(const Tensor4D & u, double m, int i, int j)
{
    for (int i1 = -1;i1<2;i1++)
    {
        for (int i2 = -1;i2<2;i2++)
        {
            for (int i4 = 0;i4<2;i4++)
            {
                for (int i5 = 0;i5<2;i5++)
                {
                    val[i+i1][j+i2][i4][i5] += u(i1+1,i2+1,i4,i5)*m;
                }
            }
        }
    }
};

void Tensor4D2::checknan()
{
    for (int i1 = 0;i1<7;i1++)
    {
        for (int i2 = 0;i2<7;i2++)
        {
            for (int i4 = 0;i4<2;i4++)
            {
                for (int i5 = 0;i5<2;i5++)
                {
                    if (isnan(val[i1][i2][i4][i5]))
                        cout<<i1<<" "<<i2<<" "<<i4<<" "<<i5<<endl;
                }
            }
        }
    
    }
}

Multigrid::Multigrid(char* fname)
{
	P = ReadParameters(fname);
    int nx = P.BoxX;
    int ny = P.BoxY;
    hx = 1.0/ipow(nx,2);
    hy = 1.0/ipow(ny,2);
    
    Sol = new ApproxArray2D*[P.levels];
    SolRes = new ApproxArray2D*[P.levels];
    SolResTmp = new ApproxArray2D*[P.levels];
    SolErr = new ApproxArray2D*[P.levels];
#if defined(COUPLED)
    SolStencil = new Array2D<Tensor4D>*[P.levels];
#else
    SolStencil = new Array2D<Tensor2D>*[P.levels];
#endif
    SolRestrict = new Array2D<Tensor2D>*[P.levels];
    SolInterpolate = new Array2D<Tensor2D>*[P.levels];
    for (int i=0; i< P.levels; i++)
    {
        Sol[i] = new ApproxArray2D(nx+1,ny+1);
        SolRes[i] = new ApproxArray2D(nx+1,ny+1);
        SolResTmp[i] = new ApproxArray2D(nx+1,ny+1);
        SolErr[i] = new ApproxArray2D(nx+1,ny+1);
#if defined(COUPLED)
        SolStencil[i] = new Array2D<Tensor4D>(nx+1,ny+1);
#else
        SolStencil[i] = new Array2D<Tensor2D>(nx+1,ny+1);
#endif
        SolRestrict[i] = new Array2D<Tensor2D>(nx+1,ny+1);
        SolInterpolate[i] = new Array2D<Tensor2D>(nx+1,ny+1);
        nx = nx/2;
        ny = ny/2;
    }

    
};

Multigrid::~Multigrid()
{
    for (int i=0; i< P.levels; i++)
    {
        delete Sol[i];
        delete SolRes[i];
        delete SolResTmp[i];
        delete SolErr[i];
        delete SolStencil[i];
        delete SolRestrict[i];
        delete SolInterpolate[i];
    }
    delete []Sol;
}

void Multigrid::setRestrInterp()
{
    for (int l=1; l< P.levels; l++)
    {
        int nx = P.BoxX/ipow(2,l);
        int ny = P.BoxY/ipow(2,l);
        int i, j, i1, j1;
//#pragma omp parallel for default(shared) private(i, i1)
        for ( i=1; i<nx;i++)
        {
            for ( j=1; j<ny;j++)
            {
                for ( i1=0; i1<3;i1++)
                {
                    for ( j1=0; j1<3;j1++)
                    {
                        (*SolRestrict[l])(i,j)(i1,j1) = 0.25/(1.0+(i1+1)%2)/(1.0+(j1+1)%2);
                        (*SolInterpolate[l-1])(2*i+i1-1,2*j+j1-1)(2-i1,2-j1) = 1.0/(1.0+(i1+1)%2)/(1.0+(j1+1)%2);
                    }
                }
            }
        }
    }
}

void Multigrid::setF()
{
    int nx = P.BoxX;
    int ny = P.BoxY;
    int oddEvenIndex = 0;
    double hx = 1.0/nx, hy = 1.0/ny;
    
    int i, j, i1, j1, i2, j2;
//#pragma omp parallel for default(shared) private(i, i1, i2, oddEvenIndex)
    for ( i=0;i<=nx;i++)
    {
        for ( j=0;j<=ny;j++)
        {
#if defined(COUPLED)
            if (i == 0 || i == nx || j == 0 || j == ny)
            {
                SolRes[0][0](i,j).nB = Sol[0][0](i,j).nB - SolStencil[0][0](i,j)(1,1,0,0)*Sol[0][0](i,j).nB;
                SolRes[0][0](i,j).nA = Sol[0][0](i,j).nA - SolStencil[0][0](i,j)(1,1,1,1)*Sol[0][0](i,j).nA;
            }
            else  
            {              
                SolRes[0][0](i,j).nB = fij(hx,hy,i,j) - (SolStencil[0][0](i,j)(0,1,0,0)*Sol[0][0](i-1,j).nB + (2.0/(hx*hx) + 2.0/(hy*hy))*Sol[0][0](i,j).nB + fA(Sol[0][0](i,j).nA) + fB(Sol[0][0](i,j).nB) + SolStencil[0][0](i,j)(2,1,0,0)*Sol[0][0](i+1,j).nB + SolStencil[0][0](i,j)(1,0,0,0)*Sol[0][0](i,j-1).nB + SolStencil[0][0](i,j)(1,2,0,0)*Sol[0][0](i,j+1).nB);
            
                SolRes[0][0](i,j).nA = fij(hx,hy,i,j) - (SolStencil[0][0](i,j)(0,1,1,1)*Sol[0][0](i-1,j).nA + fB(Sol[0][0](i,j).nB) + fA(Sol[0][0](i,j).nA)+(2.0/(hx*hx) + 2.0/(hy*hy))*Sol[0][0](i,j).nA + SolStencil[0][0](i,j)(2,1,1,1)*Sol[0][0](i+1,j).nA + SolStencil[0][0](i,j)(1,0,1,1)*Sol[0][0](i,j-1).nA + SolStencil[0][0](i,j)(1,2,1,1)*Sol[0][0](i,j+1).nA);
            }
#else
            if (i == 0 || i == nx || j == 0 || j == ny)
                SolRes[0][0](i,j).nB = Sol[0][0](i,j).nB - SolStencil[0][0](i,j)(1,1)*Sol[0][0](i,j).nB;
            else
            {
                SolRes[0][0](i,j).nB = fij(hx,hy,i,j) - (SolStencil[0][0](i,j)(0,1)*Sol[0][0](i-1,j).nB + (2.0/(hx*hx) + 2.0/(hy*hy))*Sol[0][0](i,j).nB + fB(Sol[0][0](i,j).nB ) + SolStencil[0][0](i,j)(2,1)*Sol[0][0](i+1,j).nB + SolStencil[0][0](i,j)(1,0)*Sol[0][0](i,j-1).nB + SolStencil[0][0](i,j)(1,2)*Sol[0][0](i,j+1).nB);
            }
#endif

        }
    }
}


void Multigrid::resetZero()
{    
    for (int l=0;l< P.levels;l++)
    {
        int nx = P.BoxX/ipow(2,l);
        int ny = P.BoxY/ipow(2,l);
        int i, j;
//#pragma omp parallel for default(shared) private(i)
        for ( i=0; i<=nx; i++)
        {
            for ( j=0; j<=ny; j++)
            {
#if defined(COUPLED)
                SolErr[l][0](i,j).nB = 0;
                SolErr[l][0](i,j).nA = 0;
#else
                SolErr[l][0](i,j).nB = 0;
#endif
            }
        }
    }
}

//To do Newton Multigrid change this code 
void Multigrid::setOperators()
{
    int nx = P.BoxX;
    int ny = P.BoxY;
    int oddEvenIndex = 0;
	double hx = 1.0 / nx;
    double hy = 1.0 / ny;
    
    int i, j;
    
    for (i = 0; i <= nx; i++)
    {
        for (j = 0; j <= ny; j++)
        {
#if defined(COUPLED)
            if (i == 0 || i == nx || j ==0 || j == ny)
            {
                for (int k=0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        if (k*l == 1)
                        {
                            if (i == 0 || i == nx)
                            {
                                SolStencil[0][0](i,j)(k,l,1,1) = 1.0/(hx*hx);
                                SolStencil[0][0](i,j)(k,l,0,0) = 1.0/(hx*hx);
                            }
                            else
                            {
                                SolStencil[0][0](i,j)(k,l,1,1) = 1.0/(hy*hy);
                                SolStencil[0][0](i,j)(k,l,0,0) = 1.0/(hy*hy);
                            }
                            SolStencil[0][0](i,j)(k,l,0,1) = 0;
                            SolStencil[0][0](i,j)(k,l,1,0) = 0;
                        }
                        else
                        {
                            SolStencil[0][0](i,j)(k,l,1,1) = 0;
                            SolStencil[0][0](i,j)(k,l,0,0) = 0;
                            SolStencil[0][0](i,j)(k,l,0,1) = 0;
                            SolStencil[0][0](i,j)(k,l,1,0) = 0;
                        }
                    }
                }
            }
            else
            {
                for (int k=0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        if (k*l == 1)
                        {
                            SolStencil[0][0](i,j)(k,l,0,0) = 2.0/(hx*hx) + 2.0/(hy*hy) + dfB(Sol[0][0](i,j).nB);
                            SolStencil[0][0](i,j)(k,l,0,1) = dfA(Sol[0][0](i,j).nA);
                            SolStencil[0][0](i,j)(k,l,1,0) = dfB(Sol[0][0](i,j).nB);
                            SolStencil[0][0](i,j)(k,l,1,1) = 2.0/(hx*hx) + 2.0/(hy*hy) + dfA(Sol[0][0](i,j).nA);
                        }
                        else
                        {
                            if (k == 1 && l!= 1)
                            {
                                SolStencil[0][0](i,j)(k,l,1,1) = -1.0/(hy*hy);
                                SolStencil[0][0](i,j)(k,l,0,0) = -1.0/(hy*hy);
                                SolStencil[0][0](i,j)(k,l,0,1) = 0;
                                SolStencil[0][0](i,j)(k,l,1,0) = 0;
                            }
                            else if (l == 1 && k!= 1)
                            {
                                SolStencil[0][0](i,j)(k,l,1,1) = -1.0/(hx*hx);
                                SolStencil[0][0](i,j)(k,l,0,0) = -1.0/(hx*hx);
                                SolStencil[0][0](i,j)(k,l,0,1) = 0;
                                SolStencil[0][0](i,j)(k,l,1,0) = 0;
                            }
                            else
                            {
                                SolStencil[0][0](i,j)(k,l,1,1) = 0;
                                SolStencil[0][0](i,j)(k,l,0,0) = 0;
                                SolStencil[0][0](i,j)(k,l,0,1) = 0;
                                SolStencil[0][0](i,j)(k,l,1,0) = 0;
                            }
                        }
                    }
                }
            }
#else
            for (int l = 0; l <3; l++)
            {
                for (int k=0; k<3;k++)
                {
                    SolStencil[0][0](i,j)(l,k) = 0;
                }
            }
            
            if (i == 0 || i == nx)
            {
                SolStencil[0][0](i,j)(1,1) = 1.0/(hx*hx);
            }            
            else if (j == 0 || j == ny)
            {
                SolStencil[0][0](i,j)(1,1) = 1.0/(hy*hy);
            }
            else
            {   
                SolStencil[0][0](i,j)(0,1) = -1.0/(hx*hx);
                SolStencil[0][0](i,j)(2,1) = -1.0/(hx*hx);
                SolStencil[0][0](i,j)(1,1) =  2.0/(hx*hx) + 2.0/(hy*hy)+dfB(Sol[0][0](i,j).nB);
                SolStencil[0][0](i,j)(1,0) = -1.0/(hy*hy);
                SolStencil[0][0](i,j)(1,2) = -1.0/(hy*hy);
            }
#endif
        }
    }
    

}

void Multigrid::setOperatorsNextLevel()
{
    for (int l=1; l< P.levels; l++)
    {
        int nx = P.BoxX/ipow(2,l);
        int ny = P.BoxY/ipow(2,l);
        
        int i, j, i1, j1;
#if defined(COUPLED)
        Tensor4D2 RstctA; // RA^h = R * A^h
#else
        Tensor2D2 RstctA; // RA^h = R * A^h
#endif        
//#pragma omp parallel for default(shared) private(i, i1)
        for ( i=1; i<nx; i++)
        {
            for ( j=1; j<ny; j++)
            {
                RstctA.resetZero();
                for ( i1=0; i1<3; i1++)
                {
                    for ( j1=0; j1<3; j1++)
                    {
                        RstctA.addMultiplyCopyVale((*SolStencil[l-1])(i*2+i1-1,j*2+j1-1),(*SolRestrict[l])(i,j)(i1,j1), 2+i1,2+j1);
                    }
                    (*SolStencil[l])(i,j).Interpolate(RstctA, SolRestrict[l], i, j, 4); // A^2h = RA^h * I
                }
            }
        }
    }

}

int Multigrid::ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp % 2)
            result *= base;
        exp /= 2;
        base *= base;
    }
    return result;
}

//red and black gauss seidel
void Multigrid::relaxationAtLevel(int l)
{
    int d = ipow(2,l);
    int nx = P.BoxX/d;
    int ny = P.BoxY/d;
    
    int i, j;
//#pragma omp parallel for default(shared) private(i,j)
    for (i=1; i<nx; i+=2)
    {
        for (j=1; j<ny; j++)
        {
            relaxationElementwise(l, i, j, nx, ny);
        }
    }
//#pragma omp parallel for default(shared) private(i,j)
    for (i=2; i<nx; i+=2)
    {
        for (j=1; j<ny; j++)
        {
            relaxationElementwise(l, i, j, nx, ny);
        }
    }
}

void Multigrid::relaxationElementwise(int l,int i, int j, int nx, int ny)
{
    
#if defined(COUPLED)
    double rB = SolRes[l][0](i,j).nB, trB = SolRes[l][0](i,j).nB;
    double rv = SolRes[l][0](i,j).nA, trv = SolRes[l][0](i,j).nA;
    for (int i1=0; i1<3; i1++)
    {
        for (int j1=0; j1<3; j1++)
        {
            if (i1*j1!=1)
            {
                rB = rB-(*SolStencil[l])(i,j)(i1,j1,0,0)*SolErr[l][0](i+i1-1,j+j1-1).nB-(*SolStencil[l])(i,j)(i1,j1,0,1)*SolErr[l][0](i+i1-1,j+j1-1).nA;
                rv = rv-(*SolStencil[l])(i,j)(i1,j1,1,0)*SolErr[l][0](i+i1-1,j+j1-1).nB-(*SolStencil[l])(i,j)(i1,j1,1,1)*SolErr[l][0](i+i1-1,j+j1-1).nA;
            }

        }
    }
    
    double d1 = det2x2((*SolStencil[l])(i,j)(1,1,0,0),(*SolStencil[l])(i,j)(1,1,0,1),(*SolStencil[l])(i,j)(1,1,1,0),(*SolStencil[l])(i,j)(1,1,1,1));
    double nu = det2x2(rB,rv,(*SolStencil[l])(i,j)(1,1,0,1),(*SolStencil[l])(i,j)(1,1,0,0));
    double nv = det2x2(rv,rB,(*SolStencil[l])(i,j)(1,1,1,0),(*SolStencil[l])(i,j)(1,1,1,1));
    
    if (abs(nu/d1)<10)
        SolErr[l][0](i,j).nB = nu/d1;
    if (abs(nv/d1)<10)
        SolErr[l][0](i,j).nA = nv/d1;
#else
    double rB = SolRes[l][0](i,j).nB, trB = SolRes[l][0](i,j).nB;
    for (int i1=0; i1<3; i1++)
    {
        for (int j1=0; j1<3; j1++)
        {
            if (i1*j1!=1)
            {
                rB = rB-(*SolStencil[l])(i,j)(i1,j1)*SolErr[l][0](i+i1-1,j+j1-1).nB;
            }
        }
    }
    SolErr[l][0](i,j).nB = rB/SolStencil[l][0](i,j)(1,1);
#endif

}

void Multigrid::residueRestriction(int l)
{
    int nx = P.BoxX/ipow(2,l);
    int ny = P.BoxY/ipow(2,l);
    int i, j;
//#pragma omp parallel for default(shared) private(i)
    for ( i=1; i<nx; i++)
    {
        for ( j=1; j<ny; j++)
        {
            restrictionElementwise1(l, i, j, nx, ny);
        }
    }   
    nx = P.BoxX/ipow(2,l+1);
    ny = P.BoxY/ipow(2,l+1);
//#pragma omp parallel for default(shared) private(i)
    for ( i=1; i<nx; i++)
    {
        for ( j=1; j<ny; j++)
        {
            restrictionElementwise2(l+1, i, j, nx, ny);
        }
    }
}

//Creating residual on level l, r = r - Ae
void Multigrid::restrictionElementwise1(int l,int i, int j, int nx, int ny)
{
#if defined(COUPLED)
    SolResTmp[l][0](i,j).nB = SolRes[l][0](i,j).nB;
    SolResTmp[l][0](i,j).nA = SolRes[l][0](i,j).nA;
    for (int i1=0;i1<3;i1++)
    {
        for (int j1=0;j1<3;j1++)
        {            
            SolResTmp[l][0](i,j).nB -= (*SolStencil[l])(i,j)(i1,j1,0,0)*SolErr[l][0](i+i1-1,j+j1-1).nB+(*SolStencil[l])(i,j)(i1,j1,0,1)*SolErr[l][0](i+i1-1,j+j1-1).nA;
            SolResTmp[l][0](i,j).nA -= (*SolStencil[l])(i,j)(i1,j1,1,0)*SolErr[l][0](i+i1-1,j+j1-1).nB+(*SolStencil[l])(i,j)(i1,j1,1,1)*SolErr[l][0](i+i1-1,j+j1-1).nA;
        }
    }
#else
    SolResTmp[l][0](i,j).nB = SolRes[l][0](i,j).nB;
    for (int i1=0;i1<3;i1++)
    {
        for (int j1=0;j1<3;j1++)
        {
            SolResTmp[l][0](i,j).nB -= SolStencil[l][0](i,j)(i1,j1)*SolErr[l][0](i+i1-1,j+j1-1).nB;
        }
    }
#endif
}

//interpolating residual to the next level r^{l+1} = R r^{l}
void Multigrid::restrictionElementwise2(int l,int i, int j, int nx, int ny)
{
#if defined(COUPLED)
    SolRes[l][0](i,j).nB = 0;
    SolRes[l][0](i,j).nA = 0;
    for (int i1=0;i1<3;i1++)
    {
        for (int j1=0; j1<3; j1++)
        {
            SolRes[l][0](i,j).nB += SolRestrict[l][0](i,j)(i1,j1)*SolResTmp[l-1][0](2*i+i1-1,2*j+j1-1).nB;
            SolRes[l][0](i,j).nA += SolRestrict[l][0](i,j)(i1,j1)*SolResTmp[l-1][0](2*i+i1-1,2*j+j1-1).nA;
        }
    }
#else
    SolRes[l][0](i,j).nB = 0;
    for (int i1=0;i1<3;i1++)
    {
        for (int j1=0;j1<3;j1++)
        {
            SolRes[l][0](i,j).nB += SolRestrict[l][0](i,j)(i1,j1)*SolResTmp[l-1][0](2*i+i1-1,2*j+j1-1).nB;
        }
    }
#endif

}

void Multigrid::interpolation(int l)
{
    int nx = P.BoxX/ipow(2,l);
    int ny = P.BoxY/ipow(2,l);
    int i, j;
//#pragma omp parallel for default(shared) private(i)
    for ( i=1; i<nx; i++)
    {
        for ( j=1; j<ny; j++)
        {
            interpolationElementwise(l, i, j, nx, ny);
        }
    }
}

// interpolating error. e^{l-1} = I e^l
void Multigrid::interpolationElementwise(int l, int i, int j, int nx, int ny)
{
    
#if defined(COUPLED)
    int i0 = i%2;
    int j0 = j%2;
    if (i0==0)
    {
        if (j0 == 0)
        {
            SolErr[l][0](i,j).nB += SolInterpolate[l][0](i,j)(1,1)*SolErr[l+1][0](i/2,j/2).nB;
            SolErr[l][0](i,j).nA += SolInterpolate[l][0](i,j)(1,1)*SolErr[l+1][0](i/2,j/2).nA;
        }
        else
        {
            SolErr[l][0](i,j).nB += SolInterpolate[l][0](i,j)(1,0)*SolErr[l+1][0](i/2,(j-1)/2).nB+SolInterpolate[l][0](i,j)(1,2)*SolErr[l+1][0](i/2,(j+1)/2).nB; 
            SolErr[l][0](i,j).nA += SolInterpolate[l][0](i,j)(1,0)*SolErr[l+1][0](i/2,(j-1)/2).nA+SolInterpolate[l][0](i,j)(1,2)*SolErr[l+1][0](i/2,(j+1)/2).nA; 
        }
    }
    else
    {
        if (j0 == 0)
        {
            SolErr[l][0](i,j).nB += SolInterpolate[l][0](i,j)(0,1)*SolErr[l+1][0]((i-1)/2,j/2).nB+SolInterpolate[l][0](i,j)(2,1)*SolErr[l+1][0]((i+1)/2,j/2).nB; 
            SolErr[l][0](i,j).nA += SolInterpolate[l][0](i,j)(0,1)*SolErr[l+1][0]((i-1)/2,j/2).nA+SolInterpolate[l][0](i,j)(2,1)*SolErr[l+1][0]((i+1)/2,j/2).nA; 
        }
        else
        {
            SolErr[l][0](i,j).nB += SolInterpolate[l][0](i,j)(0,0)*SolErr[l+1][0]((i-1)/2,(j-1)/2).nB+SolInterpolate[l][0](i,j)(0,2)*SolErr[l+1][0]((i+1)/2,(j-1)/2).nB+SolInterpolate[l][0](i,j)(2,0)*SolErr[l+1][0]((i-1)/2,(j+1)/2).nB+SolInterpolate[l][0](i,j)(2,2)*SolErr[l+1][0]((i+1)/2,(j+1)/2).nB; 
            SolErr[l][0](i,j).nA += SolInterpolate[l][0](i,j)(0,0)*SolErr[l+1][0]((i-1)/2,(j-1)/2).nA+SolInterpolate[l][0](i,j)(0,2)*SolErr[l+1][0]((i+1)/2,(j-1)/2).nB+SolInterpolate[l][0](i,j)(2,0)*SolErr[l+1][0]((i-1)/2,(j+1)/2).nA+SolInterpolate[l][0](i,j)(2,2)*SolErr[l+1][0]((i+1)/2,(j+1)/2).nA; 
        }
    }
#else
    int i0 = i%2;
    int j0 = j%2;
    if (i0==0)
    {
        if (j0 == 0)
        {
            SolErr[l][0](i,j).nB += SolInterpolate[l][0](i,j)(1,1)*SolErr[l+1][0](i/2,j/2).nB;
        }
        else
        {
            SolErr[l][0](i,j).nB += SolInterpolate[l][0](i,j)(1,0)*SolErr[l+1][0](i/2,(j-1)/2).nB+SolInterpolate[l][0](i,j)(1,2)*SolErr[l+1][0](i/2,(j+1)/2).nB; 
        }
    }
    else
    {
        if (j0 == 0)
        {
            SolErr[l][0](i,j).nB += SolInterpolate[l][0](i,j)(0,1)*SolErr[l+1][0]((i-1)/2,j/2).nB+SolInterpolate[l][0](i,j)(2,1)*SolErr[l+1][0]((i+1)/2,j/2).nB; 
        }
        else
        {
            SolErr[l][0](i,j).nB += SolInterpolate[l][0](i,j)(0,0)*SolErr[l+1][0]((i-1)/2,(j-1)/2).nB+SolInterpolate[l][0](i,j)(0,2)*SolErr[l+1][0]((i+1)/2,(j-1)/2).nB+SolInterpolate[l][0](i,j)(2,0)*SolErr[l+1][0]((i-1)/2,(j+1)/2).nB+SolInterpolate[l][0](i,j)(2,2)*SolErr[l+1][0]((i+1)/2,(j+1)/2).nB; 
        }
    }
#endif
}

void Multigrid::updateSol()
{
    int i, j;
//#pragma omp parallel for default(shared) private(i)
    for ( i=0;i<=P.BoxX;i++)
    {
        for ( j=0;j<=P.BoxY;j++)
        {
            updateSolElementwise(i, j);
        }
    }
}

// correction step v = v + e
void Multigrid::updateSolElementwise(int i, int j)
{
#if defined(COUPLED)
    Sol[0][0](i,j).nB += SolErr[0][0](i,j).nB;
    Sol[0][0](i,j).nA += SolErr[0][0](i,j).nA;
#else
    Sol[0][0](i, j).nB += SolErr[0][0](i, j).nB;
#endif

}

double Multigrid::evalUpdateErr()
{
#if defined(COUPLED)
    double er2U = 0;
    double er2V = 0;
    double erinfU = 0;
    double erinfV = 0;
    int i,j;
//#pragma omp parallel for default(shared) private(i) reduction(+:er2) reduction(max : erinf)
    for ( i=0; i<=P.BoxX; i++)
        for (j=0; j <=P.BoxY; j++)
        {
            //er2 += max(SolErr[0][0](i).nB,0.0)*BoxLength*BoxLength;
            er2U += max(SolErr[0][0](i,j).nB,0.0);
            er2V += max(SolErr[0][0](i,j).nA,0.0);
            erinfU = max(erinfU, abs(SolErr[0][0](i,j).nA));
            erinfV = max(erinfV, abs(SolErr[0][0](i,j).nA));
        }
    //cout<<" "<<erinf<<" " << er2 << "\n";
    return max(erinfU,erinfV);
#else
    double er2 = 0;
    double erinf = 0;
    int i, j;
//#pragma omp parallel for default(shared) private(i) reduction(+:er2) reduction(max : erinf)
    for ( i=0; i<=P.BoxX; i++)
    {
        for ( j=0; j<=P.BoxY; j++)
        {
            //er2 += max(SolErr[0][0](i).nB,0.0)*BoxLength*BoxLength;
            er2 += max(SolErr[0][0](i,j).nB,0.0);
            erinf = max(erinf, abs(SolErr[0][0](i,j).nB));
        }
    }
    //cout<<" "<<erinf<<" " << er2 << "\n";
    return erinf;
#endif

}

double Multigrid::evalUpdateRes()
{
#if defined(COUPLED)
    double res2U = 0;
    double res2V = 0;
    double resinfU = 0;
    double resinfV = 0;
    int i,j;
//#pragma omp parallel for default(shared) private(i) reduction(+:res2) reduction(max : resinf)
    for ( i=0; i<=P.BoxX; i++)
        for (j=0; j <=P.BoxY; j++)
        {
            //res2 += max(SolRes[0][0](i).nB,0.0)*BoxLength*BoxLength;
            res2U += max(SolRes[0][0](i,j).nB,0.0);
            res2V += max(SolRes[0][0](i,j).nA,0.0);
            resinfU = max(resinfU, abs(SolRes[0][0](i,j).nB));
            resinfV = max(resinfV, abs(SolRes[0][0](i,j).nA));
        }
    //cout<<" "<<resinf<<" "<<res2<<" "<<endl;
    return max(resinfU,resinfV);
#else
    double res2 = 0;
    double resinf = 0;
    int i, j;
//#pragma omp parallel for default(shared) private(i) reduction(+:res2) reduction(max : resinf)
    for ( i=0; i<=P.BoxX; i++)
    {
        for ( j=0; j<=P.BoxY; j++)
        {
            //res2 += max(SolRes[0][0](i).nB,0.0)*BoxLength*BoxLength;
            res2 += max(SolRes[0][0](i,j).nB,0.0);
            resinf = max(resinf, abs(SolRes[0][0](i,j).nB));
        }                                                                   
    }
    //cout<<" "<<resinf<<" "<<res2<<" "<<endl;
    return resinf;
#endif

}

void Multigrid::setup()
{
    resetZero();
}

void Multigrid::reset()
{
    updateSol();
    evalUpdateRes();
    
    resetZero();
}

void Multigrid::Vcycle(int mu1, int mu2)
{
    setRestrInterp();
    setF();
    setOperators();
    setOperatorsNextLevel();
    for (int l=0; l< P.levels; l++)
    {
        for (int j=0;j<mu1; j++)
        {
            relaxationAtLevel(l); // solving r = Ae
        }
        if (l< P.levels-1)
            residueRestriction(l);
    }
    for (int j=0;j<mu2; j++)
    {
        relaxationAtLevel(P.levels-1);
    }
    for (int l= P.levels-1; l>=0; l--)
    {
        if (l< P.levels-1)
            interpolation(l);
        for (int j=0;j<mu1; j++)
        {
            relaxationAtLevel(l);
        }
    }
}

//Need to work on this
double Multigrid::ResVcycle()
{
#if defined(COUPLED)
    int i,j, l=0;
    double rB=0, rV=0;
    for (i=1; i<P.BoxX; i++)
        for (j=0; j <=P.BoxY; j++)
        {
            double ru = SolRes[l][0](i,j).nB;
            double rv = SolRes[l][0](i,j).nA;
            for (int i1=0; i1<3; i1++)
                for (int j1=0; j1<3; j1++)
                {
                    ru = ru-(*SolStencil[l])(i,j)(i1,j1,0,0)*SolErr[l][0](i+i1-1,j+j1-1).nB-(*SolStencil[l])(i,j)(i1,j1,0,1)*SolErr[l][0](i+i1-1,j+j1-1).nA;
                    rv = rv-(*SolStencil[l])(i,j)(i1,j1,1,0)*SolErr[l][0](i+i1-1,j+j1-1).nB-(*SolStencil[l])(i,j)(i1,j1,1,1)*SolErr[l][0](i+i1-1,j+j1-1).nA;
                }  
            rB = max(rB,abs(rB));
            rV = max(rV,abs(rv));
        }
    
	//cout << "V-cycle Colony residue: "<< rB << " " << rV <<"\n";
    return max(rB,rV);
#else
    int i, j, l;
    double rB=0;
	for (i=1; i<P.BoxX; i++)
	{
        for (j=1; j<P.BoxY; j++)
        {
            double ru = SolRes[0][0](i,j).nB;
            for (int i1=0; i1<3; i1++)
            {
                for (int j1=0; j1<3; j1++)
                {
                    ru = ru-(*SolStencil[l])(i,j)(i1,j1)*SolErr[l][0](i+i1-1,j+j1-1).nB;
                }
            }
            rB = max(rB,abs(ru));
            //rB = rB;
        }
    }
	//cout << "V-cycle Colony residue: "<<rB<< " " << "Level: " << l <<"\n";
    return rB;
#endif

}

void Multigrid::solve(int mu1, int mu2)
{
    setup();
    //OutputFiles Files;
    double er = 1; 
	double rB = 1;
    double tol = 1e-14;
    int max_It = 10;
    int it = 1;
	
    while (it<max_It && er>tol)
    {
        Vcycle(mu1, mu2);
        rB = ResVcycle();
        er = evalUpdateErr();
        cout<< "Cycle = " <<it<<" Error = " << er << " Residual = " << rB <<  "\n";
        reset();
        it++;
    }
    //cout << it << " " << ResVcycle() << " " << rB << " " << evalUpdateRes() << endl;
    
}