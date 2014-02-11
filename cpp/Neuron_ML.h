#ifndef NEURON_ML_INC
#define NEURON_ML_INC


#include "Neuron.h"
#include <cmath>


// Morris-Lecar
class Neuron_ML: public Neuron
{
public:
	// "Recovery variable" -- The probability that the K+ channel is conducting
	double n, np;
	// Model parameters, common to all neurons of this class
	static double phi, gCa, V3, V4, ECa, EK, EL, gK, gL, V1, V2, CM;
	
	Neuron_ML();
	void step(double dt, double input_);
	void step(double dt);
	// Switch both s and n with the previous values
	void update() { Neuron::update(); np = n; }

private:
	// Returns the synaptic potential, given a presynaptic potential of 'input'
	double synPot(double input);

	// Functions appearing in the model
	double mInf(double V) { return (1+tanh((V-V1)/V2))/2; }
	double tauN(double V) { return 1 / cosh((V-V3)/(2*V4)); }
	double nInf(double V) { return (1 + tanh((V-V3)/V4))/2; }
};




#endif
