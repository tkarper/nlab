#ifndef NEURON_STEL_INC
#define NEURON_STEL_INC


#include "Neuron.h"


// Hodgkin-Huxley by Acker et al. (2003) for stellate cells
class Neuron_Stel: public Neuron
{
public:
	double V, VP;	// Synaptic potential
	// m=activation, h=deactivation
	double mNa, mNaP;
	double hNa, hNaP;
	double n, nP;
	double mNap, mNapP;
	
	// Model parameters, common to all neurons of this class
	static double VL, VNa, VK,
		gNap, gL, gNa, gk,
		CM;		// Time scaling factor
	double a_r, a_d;	// Rise/decay rate of synaptic conductance
	
	Neuron_Stel();
	void step(double t, double dt, double input_);
	// Switch both s, V and n with the previous values
	void update() { Neuron::update(); VP=V; mNaP=mNa; hNaP=hNa; nP=n; mNapP=mNap; }
};



class Neuron_IntN : public Neuron_Stel
{
public:
	Neuron_IntN();
};




#endif
