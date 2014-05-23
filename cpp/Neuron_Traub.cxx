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
EK(-85),	// From Destexe ea 1999; EK=-100 in Ermentrout ea 2001
ENa(50),
EL(-67),
gL(0.2),
gK(80),
gM(0),		// No M-current
gNa(100),
CM(1),
alpha_s(10),	// Excitatory
beta_s(5),		// Excitatory
VT(-67),	// (Ermentrout ea 2001: VT=-67)
VS(0)
{
	Esyn = 0;	// Excitatory
}


Neuron_Traub_IN::Neuron_Traub_IN()
{
	Esyn = -80;
	alpha_s = 2;
	beta_s = 0.1;
}


double linSmooth(double x, double theta) {
	const double tol = 1e-8;
	return std::abs(x)<tol ? theta : x / (1 - exp(-x/theta));
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
//	double alpha_m = 0.32*(Vp-VT-13) / (1-exp(-(Vp-VT-13)/4)),
//		alpha_h = 0.128*exp(-(Vp-VT-VS-17)/18),
//		alpha_n = 0.032*(Vp-VT-15) / (1-exp(-(Vp-VT-15)/5)),
//		beta_m = -0.28*(Vp-VT-40) / (1-exp((Vp-VT-40)/5)),
//		beta_h = 4 / (1+exp(-(Vp-VT-VS-40)/5)),
//		beta_n = 0.5*exp(-(Vp-VT-10)/40);
	double alpha_m = 0.32*linSmooth(Vp-VT-13, 4),
		alpha_h = 0.128*exp(-(Vp-VT-VS-17)/18),
		alpha_n = 0.032*linSmooth(Vp-VT-15, 5),
		beta_m = -0.28*linSmooth(Vp-VT-40, -5),
		beta_h = 4 / (1+exp(-(Vp-VT-VS-40)/5)),
		beta_n = 0.5*exp(-(Vp-VT-10)/40);
	
	// Ermentrout ea 2001
//	double 	wInf = 1 / (1+exp(-(Vp+35)/10)),
//		tauW = 100 / (3.3*exp((Vp+35)/20) + exp(-(Vp+35)/20));
//	w = wp + dt*(wInf - wp)/tauW;
	
	// Destexhe and Paré 1999
//	double alpha_w = 0.0001*(Vp+30) / (1-exp(-(Vp+30)/9)),
//		beta_w = -0.0001*(Vp+30) / (1-exp((Vp+30)/9));
//	double alpha_w = 0.0001*linSmooth(Vp+30, 9),
//		beta_w = -0.0001*linSmooth(Vp+30, -9);
//	w = wp + dt*(alpha_w*(1-wp) - beta_w*wp);
	double wTau = 1/(alpha_w+beta_w),
		wInf = 1 / (1 + exp(-(-30-Vp)/-9));
	w = wp + dt*(wInf - w)/wTau;

	double VNew = Vp + dt/CM*(I+input 
		- (Vp-ENa) * gNa*hp*pow(mp,3)
		- (Vp-EK) * (gK*pow(np,4) + gM*wp)
		- (Vp-EL) * gL
		- synPot);
	if (VNew > 0 && V <= Vp && Vp > VNew)
		isFiring = true;
	else
		isFiring = false;
	V = VNew;
	
	m = mp + dt*(alpha_m*(1-mp) - beta_m*mp);
	h = hp + dt*(alpha_h*(1-hp) - beta_h*hp);
	n = np + dt*(alpha_n*(1-np) - beta_n*np);
	s = sp + dt*(alpha_s*(1-sp)/(1+exp(-(Vp+10))) - beta_s*sp);		// In Ermentrout ea 2001: (...)exp(-(Vp+10)/10))
}
