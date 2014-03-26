#ifndef NEURON_TIF_INC
#define NEURON_TIF_INC

#include <vector>
#include "Neuron.h"

/*
		"TIRED" INTEGRATE-AND-FIRE NEURON MODEL
*/
class Neuron_TIF: public Neuron
{
public:
	double d, dp;	// Tiredness variable
	double C1, C2;	// Dampening and input constants

	void step(double t, double dt, double input);
	void update() { Neuron::update(); dp = d; }
};


#endif
