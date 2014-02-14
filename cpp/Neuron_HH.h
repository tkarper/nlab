#ifndef NEURON_HH_INC
#define NEURON_HH_INC


#include "Neuron.h"


// Hodgkin-Huxley
class Neuron_HH: public Neuron
{
public:
	double V, Vp;	// Synaptic potential
	double n, np;	// K+ channel activation
	double m, mp;	// Na+ channel activation
	double h, hp;	// Na+ channel deactivation
	// Model parameters, common to all neurons of this class
	static double phi, gNa, gK, gL, ENa, EK, EL, CM, alpha, beta, VT, VM;
	
	Neuron_HH();
	void step(double dt, double input_);
	void step(double dt);
	// Switch both s, V and n with the previous values
	void update() { Neuron::update(); np = n; Vp = V; mp=m; hp=h;}
};




#endif
