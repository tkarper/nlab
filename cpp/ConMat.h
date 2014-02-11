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




#endif