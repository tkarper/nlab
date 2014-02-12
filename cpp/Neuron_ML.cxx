#include "Neuron_ML.h"


// Set model parameters (could be done in Python script)
double Neuron_ML::phi = 0.04;
double Neuron_ML::gCa = 4.4;
double Neuron_ML::V3 = 2;
double Neuron_ML::V4 = 30;
double Neuron_ML::ECa = 120;
double Neuron_ML::EK = -84;
double Neuron_ML::EL = -60;
double Neuron_ML::gK = 8;
double Neuron_ML::gL = 2;
double Neuron_ML::V1 = -1.2;
double Neuron_ML::V2 = 18;
double Neuron_ML::CM = 20;



Neuron_ML::Neuron_ML() :
n(0), np(0)
{
}


void Neuron_ML::step(double dt, double inp)
{
	// Compute input (synaptic potentials) coming from other neurons
	double ip = 0.0;
	for(int n=0; n<con->size(); n++)
	{
		Neuron* nr = con->at(n);
		ip += weight->at(n) * synPot(nr->sp);
	}
	
	s = sp - dt*((gL*(sp-EL) + gK*np*(sp-EK) + gCa*mInf(sp)*(sp-ECa)))/CM;
	n = np + dt*phi*(nInf(sp)-np)/tauN(sp);
}


void Neuron_ML::step(double dt)
{
	step(dt, 0.0);
}
