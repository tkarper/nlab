#include <cmath>
#include <iostream>
#include "Neuron_TIF.h"


double Neuron_TIF::C1 = 1.0;
double Neuron_TIF::C2 = 1.0;


void Neuron_TIF::step(double t, double dt, double inp)
{
	double ip=0.0;
	Neuron * nr;
	for(size_t n=0; n<con->size(); n++)
	{
		nr = con->at(n);
		ip += weight->at(n)*(nr->sp);
	}
	
	s = sp + dt*(-sp + fmax(ip + I + inp - d,0.0));
	d = dp + dt*(-C1*dp + C2*sp);
}
