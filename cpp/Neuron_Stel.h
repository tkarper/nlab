#ifndef NEURON_STEL_INC
#define NEURON_STEL_INC


#include "Neuron.h"


// Hodgkin-Huxley by Acker et al. (2003) for stellate cells
class Neuron_Stel: public Neuron
{
public:
	double V, VP;	// Membrane potential
	// m=activation, h=deactivation
	double mNa, mNaP;
	double hNa, hNaP;
	double n, nP;
	double mNap, mNapP;
//	double mDep, mDepP;		// Degree of depression
//	double mFac, mFacP;		// Degree of facilitation
	double mSyn, mSynP;		// Synaptic potential (non-depressed, non-facilitated)
	
	// Model parameters, common to all neurons of this class
	static double VL, VNa, VK,
		gNap, gL, gNa, gK,
		CM;		// Time scaling factor
//	double a_r, a_d;	// Rise/decay rate of synaptic conductance, measured in mM^(-1)ms^(-1)
//	double mDep0, mFac0;	// Resting values for depression and facilitation probabilities
//	double mDepRise, mDepDecay, mFacRise, mFacDecay;
	
	Neuron_Stel();
	void step(double t, double dt, double input_);
	// Switch both s, V and n with the previous values
	void update() { Neuron::update(); VP=V; mNaP=mNa; hNaP=hNa; nP=n; mNapP=mNap; mDepP=mDep; mFacP=mFac; mSynP = mSyn; }
};



class Neuron_Pyr : public Neuron_Stel
{
public:
	Neuron_Pyr();
};




#endif
