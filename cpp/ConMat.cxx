#include "ConMat.h"
#include <iostream>


ConMat::ConMat(int num)
{
	ConMat(num,num);
}


ConMat::~ConMat()
{
}


ConMat::ConMat(int n, int m)
{		
	W = new ddvector(n);
	for(int i=0; i<n; i++)
	{
		W->at(i) = new dvector();
	}	
}


void ConMat::add(int i, int j, double alpha)
{
	if(static_cast<size_t>(i) >= W->size())
	{
		std::cout << "index exceeds size of ConMat"<< std::endl;
		return;
	}
	
	W->at(i)->push_back(j);
	W->at(i)->push_back(alpha);		
}
