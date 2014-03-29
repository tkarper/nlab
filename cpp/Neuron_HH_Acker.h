#ifndef NEURON_HH_INC
#define NEURON_HH_INC


#include "Neuron.h"


// Hodgkin-Huxley by Acker et al. (2003)
class Neuron_HH_Acker: public Neuron
{
public:
	double V, VP;	// Synaptic potential
	// m=activation, h=deactivation
	double mNa, mNaP;
	double hNa, hNaP;
	double n, nP;
	double mNap, mNapP;
	double mKs, mKsP;
	
	double n, np;	// K+ channel activation
	double m, mp;	// Na+ channel activation
	double h, hp;	// Na+ channel deactivation
	// Model parameters, common to all neurons of this class
	static double phi, gNa, gK, gL, ENa, EK, EL, CM, alpha, beta, VT, VM;
	
	Neuron_HH();
	void step(double t, double dt, double input_);
	// Switch both s, V and n with the previous values
	void update() { Neuron::update(); np = n; Vp = V; mp=m; hp=h;}
};




#endif
