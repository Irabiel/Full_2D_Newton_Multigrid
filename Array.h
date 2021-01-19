#ifndef ARRAY_H_
#define ARRAY_H_
#include <assert.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <memory.h>

#include "tools.h"
#include "Nutrients.h"


/*************************************
 * A bunch of functions for array manipulations
 ************************************/
template<typename T>
class Array2D
{
public:
	IntCoord2D m_Size;

	Array2D()
	{
		Data = NULL;
	}

	Array2D(int xSize, int ySize)
			:m_Size(xSize, ySize)
	{
		Data = new T [m_Size.x*m_Size.y];	// 1D array of doubles
	}

	const IntCoord2D& Size() { return m_Size; }

	~Array2D(void)
	{
		if (Data!=NULL)
			delete[] Data;
	}

	T Get(const IntCoord2D& coord)
	{
		return Get(coord.x, coord.y);
	}

	T Get(int xi, int yi)
	{
		return Data[Index(xi, yi)];
	}

	T& At(const IntCoord2D& coord)
	{
		return At(coord.x, coord.y);
	}

	T& At(int xi, int yi)
	{
		return Data[Index(xi, yi)];
	}
    
    T& operator()(int a, int b)
    {
        return At(a, b);
    }
    
    const T& operator()(int a, int b) const
    {
        return At(a, b);
    }

	const T& At(const IntCoord2D& coord) const
	{
		return At(coord.x, coord.y);
	}

	const T& At(int xi, int yi) const
	{
		return Data[Index(xi, yi)];
	}

	void Set(int xi, int yi, T value)
	{
		const int index = Index(xi, yi);
		MyAssert(m_Size.x*m_Size.y>index,"Invalid index into array");
		Data[index] = value;
	}

	void SetData(T* DataArray, int Length)
	{
		MyAssert(m_Size.x*m_Size.y==Length,"Incorrect array size");
		memcpy(Data,DataArray,Length*sizeof(T));
	}

	void Initialize(T val)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			for (int yi = 0; yi < m_Size.y; yi++)
			{
				Data[xi*m_Size.y+yi] = val;
			}
		}
	}

	void Output(char* Mname)//FILE* FID
	{
		FILE* FID = fopen(Mname,"w");
		Output(FID);
		fclose(FID);
	}

	void Output(FILE* FID);
    
    void Append(FILE* FID);
    
    void Input(char* Fname)
    {
        std::ifstream infile("data.txt");
        Input(infile);
    }
    
    void Input(std::ifstream& infile);

	double linear_interp(const int x0, const int y0, const double dx, const double dy);

private:
	T* Data;

	int Index(const IntCoord2D& coord) { return Index(coord.x, coord.y); }
	int Index(const int x, const int y) { return x*m_Size.y+y; }
	int Index(const int x, const int y) const { return x*m_Size.y+y; }
};

template <>
inline void Array2D<double>::Output(FILE* FID)
{
	rewind(FID);
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi]);
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline void Array2D<int>::Output(FILE* FID)
{
	rewind(FID);
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%d\t", Data[xi*m_Size.y + yi]);
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline double Array2D<double>::linear_interp(const int x0, const int y0, const double dx, const double dy)
{
	float x = dx<0? 1+dx: dx;
	float y = dy<0? 1+dy: dy;

	return Data[Index(x0, y0)]*(1-x)*(1-y) + Data[Index(x0+1, y0)]*x*(1-y) + Data[Index(x0, y0+1)]*(1-x)*y + Data[Index(x0+1, y0+1)]*x*y;
};

template <>
inline void Array2D<DoubleCoord>::Output(FILE* FID)
{
	rewind(FID);
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].x);
		}
		fprintf(FID, "\n");
	}
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].y);
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline void Array2D<Testprob>::Output(FILE* FID)
{
#if defined(COUPLED)
    rewind(FID);
    for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].nB);
		}
		fprintf(FID, "\n");
	}
    for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].nA);
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
#else
	rewind(FID);
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].nB);
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
#endif

};

typedef Array2D<double> DoubleArray2D;
typedef Array2D<DoubleCoord> CoordArray2D;
typedef Array2D<Testprob> ApproxArray2D;


#endif /* ARRAY_H_ */
