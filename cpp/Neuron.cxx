#include "Neuron.h"
#include <cmath>
#include <iostream>

Neuron::Neuron()
{
	con 	= new nvector();
	weight	= new dvector();
	sp = 0.0;
	s = 0.0;
	I = 0.0;
}

Neuron::~Neuron()
{
	std::cout << "CALLED DESTRUCTOR" << std::endl;
}



void Neuron::connect(Neuron* n1, double val)
{
	con->push_back(n1);
	weight->push_back(val);
}


Neuron_IF::Neuron_IF()
{
}

Neuron_IF::~Neuron_IF()
{
	std::cout << "CALLED IF's DESTRUCTOR" << std::endl;
}



void Neuron_IF::step(double dt, double inp)
{
	double ip=0.0;
	Neuron * nr;
	for(int n =0; n< con->size();n++)
	{
		nr = con->at(n);
		ip += weight->at(n)*(nr->sp);
	}
	
	s = (1.0 - dt)*sp + dt*fmax(ip + I + inp,0.0);
}


void Neuron_IF::step(double dt)
{
	step(dt,0.0);
}





