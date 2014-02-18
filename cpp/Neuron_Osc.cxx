#include <iostream>
#include <cmath>
#include "Neuron_Osc.h"


Neuron_Osc::Neuron_Osc(double p, double d, double s) :
period(p),
duration(d),
strength(s)
{
}


void Neuron_Osc::step(double t, double dt, double input) {
	if (fmod(t, period) < duration)
		s = strength;
	else
		s = 0;
}
