#ifndef NEURON_ACK_INC
#define NEURON_ACK_INC


#include "Neuron.h"


// Hodgkin-Huxley by Acker et al. (2003)
class Neuron_Ack: public Neuron
{
public:
	double V, VP;	// Synaptic potential
	// m=activation, h=deactivation
	double mNa, mNaP;
	double hNa, hNaP;
	double n, nP;
	double mNap, mNapP;
//	double mKs, mKsP;
	double mhf, mhfP;
	double mhs, mhsP;
	double msyn, msynP;
	
	// Model parameters, common to all neurons of this class
	static double VL, VNa, VK, Vh, Vsyn,
		gNap, gh, gsyn, gL, gNa, gk,
		CM;		// Time scaling factor
	
	Neuron_Ack();
	void step(double t, double dt, double input_);
	// Switch both s, V and n with the previous values
	void update() { Neuron::update(); VP=V; mNaP=mNa; hNaP=hNa; nP=n; mNapP=mNap; mhfP=mhf; mhsP=mhs; msynP=msyn; }
};




#endif
