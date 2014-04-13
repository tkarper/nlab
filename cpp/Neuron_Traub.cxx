#include <cstdio>
#include <cmath>
#include <iostream>
#include "Neuron_Traub.h"



Neuron_Traub::Neuron_Traub() :
V(-62), Vp(-62),
m(0), mp(0),
n(0), np(0),
h(0), hp(0),
w(0), wp(0),
EK(-100),
ENa(50),
EL(-67),
gL(0.2),
gK(80),
gM(0),		// No M-current
gNa(100),
CM(1),
alpha_s(10),	// Excitatory
beta_s(5)		// Excitatory
{
	Esyn = 0;	// Excitatory
}


void Neuron_Traub::step(double t, double dt, double input)
{
	// Compute input (synaptic potentials) coming from other neurons
	double synPot = 0;
	for(size_t i=0; i<con->size(); i++)
	{
		Neuron* nr = con->at(i);
		synPot += weight->at(i) * nr->sp * nr->gsyn * (Vp - nr->Esyn);
	}
	
	// Functions appearing in the model
	//const double tol = 1e-10;
	double alpha_m = 0.32*(54+Vp) / (1-exp(-(Vp+54)/4)),
		alpha_h = 0.128*exp(-(50+Vp)/18),
		alpha_n = 0.032*(Vp+52) / (1-exp(-(Vp+52)/5)),
		beta_m = 0.28*(Vp+27) / (exp((Vp+27)/5) - 1),
		beta_h = 4 / (1+exp(-(Vp+27)/5)),
		beta_n = 0.5*exp(-(57+Vp)/40),
		wInf = 1 / (1+exp(-(Vp+35)/10)),
		tauW = 100 / (3.3*exp((Vp+35)/20) + exp(-(Vp+35)/20));

	double VNew = Vp + dt/CM*(I+input 
		- (Vp-ENa) * gNa*hp*pow(mp,3)
		- (Vp-EK) * (gK*pow(np,4) + gM*wp)
		- (Vp-EL) * gL
		- (Vp-Esyn) * gsyn*synPot);
	if (VNew > 0 && V <= Vp && Vp > VNew)
		isFiring = true;
	else
		isFiring = false;
	V = VNew;
	
	m = mp + dt*(alpha_m*(1-mp) - beta_m*mp);
	h = hp + dt*(alpha_h*(1-hp) - beta_h*hp);
	n = np + dt*(alpha_n*(1-np) - beta_n*np);
	w = wp + dt*(wInf - wp)/tauW;
	s = sp + dt*(alpha_s*(1-sp)/(1+exp(-(Vp+10)/10)) - beta_s*sp);
}
