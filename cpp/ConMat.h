#ifndef CONMAT_INC
#define CONMAT_INC
#include <vector>
#include "Neuron.h"

typedef std::vector<dvector*> ddvector;

class ConMat
{
public:
	ddvector* W;
	ConMat(int num_);
	ConMat(int n, int m);
	~ConMat();

	void add(int i, int j, double alpha);
	
	// TODO: 
	// M(i,j) = alpha
	// 	a = M(i,j)
	
};


typedef std::pair<int, double> NeuroConn;
typedef std::vector<NeuroConn> cvector;
typedef std::vector<cvector> ccvector;

class ConMat2
{
public:
	ddvector W;
	ConMat2(size_t n);
	ConMat2(size_t n, size_t m);
	~ConMat2();

	void add(size_t i, size_t j, double alpha);
};




#endif
