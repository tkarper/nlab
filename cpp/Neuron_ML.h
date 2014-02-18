#ifndef NEURON_ML_INC
#define NEURON_ML_INC


#include "Neuron.h"


// Morris-Lecar
class Neuron_ML: public Neuron
{
public:
	// Synaptic potential
	double V, Vp;
	// "Recovery variable" -- The probability that the K+ channel is conducting
	double n, np;
	// Model parameters, common to all neurons of this class
	static double phi, gCa, V3, V4, ECa, EK, EL, gK, gL, V1, V2, CM, alpha, beta, VT;
	
	Neuron_ML();
	void step(double t, double dt, double input_);
	// Switch both s, V and n with the previous values
	void update() { Neuron::update(); np = n; Vp = V; }
};




#endif
