#include <cmath>
#include <iostream>
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
double Neuron_ML::VT = -50;
double Neuron_ML::alpha = 1;
double Neuron_ML::beta = 1;



// Functions appearing in the model
double mInf(double V) { return (1+tanh((V-Neuron_ML::V1)/Neuron_ML::V2))/2; }
double tauN(double V) { return 1 / cosh((V-Neuron_ML::V3)/(2*Neuron_ML::V4)); }
double nInf(double V) { return (1 + tanh((V-Neuron_ML::V3)/Neuron_ML::V4))/2; }


Neuron_ML::Neuron_ML() :
V(EL), Vp(EL),
n(0), np(0)
{
}


void Neuron_ML::step(double dt, double inp)
{
	// Compute input (synaptic potentials) coming from other neurons
	double ip = 0.0;
	for(size_t i=0; i<con->size(); i++)
	{
		Neuron* nr = con->at(i);
		ip += weight->at(i) * nr->sp;
	}
	
	V = Vp - dt*((gL*(Vp-EL) + gK*np*(Vp-EK) + gCa*mInf(Vp)*(Vp-ECa)))/CM;
	n = np + dt*phi*(nInf(Vp)-np)/tauN(Vp);
	s = sp + dt*(alpha*(1-sp)*fmax(Vp - VT, 0) - beta*sp);
}


void Neuron_ML::step(double dt)
{
	step(dt, 0.0);
}
