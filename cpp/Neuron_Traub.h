#ifndef NEURON_TRAUB_INC
#define NEURON_TRAUB_INC


#include "Neuron.h"


// Conductance model by Traub (from Ermentrout ea 2001)
class Neuron_Traub: public Neuron
{
public:
	double V, Vp;	// Synaptic potential
	// m=activation, h=deactivation
	double m, mp,
		n, np,
		h, hp,
		w, wp;
	
	// Model parameters, common to all neurons of this class
	double EK, ENa, EL, gL, gK, gM, gNa, CM, alpha_s, beta_s;
	
	Neuron_Traub();
	void step(double t, double dt, double input_);
	// Switch both s, V and n with the previous values
	void update() { Neuron::update(); std::swap(V,Vp); mp=m; np=n; hp=h; wp=w; }
};




#endif
