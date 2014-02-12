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






ConMat2::ConMat2(size_t n)
{
	ConMat(n, n);
}


ConMat2::ConMat2(size_t n, size_t m)
{
	W.resize(n);
	for(int i=0; i<n; i++)
	{
		W[i].reserve(m);
	}	
}


ConMat2::~ConMat2()
{
}


void ConMat2::add(size_t i, size_t j, double alpha)
{
	if(i >= W.size() || j >= W[i].size())
	{
		std::cout << "index exceeds size of ConMat"<< std::endl;
		return;
	}
	
	W[i].push_back(NeuroConn(j, alpha));
}
