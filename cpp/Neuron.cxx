#include "Neuron.h"
#include <cmath>
#include <iostream>

Neuron::Neuron() :
Vsyn(0),
gsyn(1)
{
	con 	= new nvector();
	weight	= new dvector();
	sp = 0.0;
	s = 0.0;
	I = 0.0;
}

Neuron::~Neuron()
{

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

}



void Neuron_IF::step(double t, double dt, double inp)
{
	double ip=0.0;
	Neuron * nr;
	for(size_t n=0; n<con->size(); n++)
	{
		nr = con->at(n);
		ip += weight->at(n)*(nr->sp);
	}
	
	s = (1.0 - dt)*sp + dt*fmax(ip + I + inp,0.0);
}




