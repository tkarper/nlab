#ifndef NEURON_OSC_INC
#define NEURON_OSC_INC


#include "Neuron.h"


// On/off oscillator
class Neuron_Osc: public Neuron
{
public:
	double period;		// Period of signal activation
	double duration;	// Duration of each signal
	double strength;	// Strength of signal once activated
	
	Neuron_Osc(double p, double d, double s);
	void step(double t, double dt, double input);
};




#endif
