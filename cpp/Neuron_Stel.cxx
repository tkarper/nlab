#include <cmath>
#include <iostream>
#include "Neuron_Stel.h"


// Set model parameters (could be done in Python script)
double Neuron_Stel::VL	= -65;	
double Neuron_Stel::VNa	= 55;
double Neuron_Stel::VK	= -90;
double Neuron_Stel::gNa	= 52;
double Neuron_Stel::gNap	= 0.5;
double Neuron_Stel::gL	= 0.5;
double Neuron_Stel::gk	= 11;
double Neuron_Stel::CM	= 1.5;



Neuron_Stel::Neuron_Stel() :
V(VL), VP(VL), mNa(0), mNaP(0), hNa(0), hNaP(0), n(0), nP(0), mNap(0), mNapP(0)
{
	Vsyn = 0;
	gsyn = 0.006;
	a_r = 1.1;
	a_d = 0.19;
}


void Neuron_Stel::step(double t, double dt, double input)
{
	// Compute input (synaptic potentials) coming from other neurons
	double synCurrent = 0;
	for(size_t i=0; i<con->size(); i++)
	{
		Neuron* nr = con->at(i);
		synCurrent += weight->at(i) * nr->sp * nr->gsyn * (VP - nr->Vsyn);
	}
	
	// Functions appearing in the model
	//const double tol = 1e-10;
	double alpha_mNa = 0.1*(VP+23)/(1 - exp(-(VP+23)/10));
	double alpha_n = 0.01*(VP+27)/(1 - exp(-(VP+27)/10));
	double alpha_mNap = 1/(0.15*(1 + exp(-(VP+38)/6.5)));
	double alpha_hNa = 0.07*exp(-(VP+37)/20);
	double beta_mNa = 4*exp(-(VP+48)/18);
	double beta_hNa = 1/(1+exp(-(VP+7)/10));
	double beta_n = 0.125*exp(-(VP+37)/80);
	double beta_mNap = exp(-(VP+38)/6.5)/(0.15*(1+exp(-(VP+38)/6.5)));

	V = VP + dt/CM*(I+input - (
		(gNa*pow(mNaP,3)*hNaP + gNap*mNap)*(VP-VNa) + 
		gk*pow(nP,4)*(VP-VK) +
		gL*(VP-VL) + 
		synCurrent));
	mNa = mNaP + dt*(alpha_mNa*(1-mNaP) - beta_mNa*mNaP);
	hNa = hNaP + dt*(alpha_hNa*(1-hNaP) - beta_hNa*hNaP);
	mNap = mNapP + dt*(alpha_mNap*(1-mNapP) - beta_mNap*mNapP);
	n = nP + dt*(alpha_n*(1-nP) - beta_n*nP);
	s = sp + dt*((VP>0 ? a_r : 0.0)*(1-sp) - a_d*sp);		// What should the transmitter pulse threshold be???
}



Neuron_IntN::Neuron_IntN()
{
	Vsyn = -75;
	a_r = 5;
	a_d = 0.18;
}
