#include <cmath>
#include <iostream>
#include "Neuron_HH.h"


// Set model parameters (could be done in Python script)
double Neuron_HH::phi = 1;		// Temperature factor
double Neuron_HH::gNa = 120;
double Neuron_HH::gK = 36;
double Neuron_HH::gL = 0.3;
double Neuron_HH::ENa = 50;
double Neuron_HH::EK = -77;
double Neuron_HH::EL = -54.4;
double Neuron_HH::CM = 1;
double Neuron_HH::alpha = 1;
double Neuron_HH::beta = 1;
double Neuron_HH::VT = -50;
double Neuron_HH::VM = 40;


// Functions appearing in the model
const double tol = 1e-10;
double alpha_n(double V) { return (fabs(V+55)<tol) ? 0.1 : 0.01*(V+55)/(1 - exp(-(V+55)/10)); }
double alpha_m(double V) { return (fabs(V+40)<tol) ? 1 : 0.1*(V+40)/(1 - exp(-(V+40)/10)); }
double alpha_h(double V) { return 0.07*exp(-(V+65)/20); }
double beta_n(double V) { return 0.125*exp(-(V+65)/80); }
double beta_m(double V) { return 4*exp(-(V+65)/18); }
double beta_h(double V) { return 1/(1+exp(-(V+35)/10)); }
int heaviside(double V) { return V > 0 ? 1 : 0; }


Neuron_HH::Neuron_HH() :
V(-65), Vp(-65),
n(0.317), np(0.317), 
m(0.053), mp(0.053),
h(0.6), hp(0.6)
{
}


void Neuron_HH::step(double dt, double inp)
{
	// Compute input (synaptic potentials) coming from other neurons
	double ip = 0.0;
	for(size_t i=0; i<con->size(); i++)
	{
		Neuron* nr = con->at(i);
		ip += weight->at(i) * nr->sp;
	}

	V = Vp + dt*(/*I+inp+ip*/ - (gL*(Vp-EL) + gK*pow(np,4)*(Vp-EK) + gNa*pow(mp,3)*hp*(Vp-ENa)))/CM;
	n = np + dt*phi*(alpha_n(Vp)*(1-np) - beta_n(Vp)*np);
	// Signals are summed up in the Na+ gate variable
	m = mp + dt*(I+inp+ip + phi*(alpha_m(Vp)*(1-mp) - beta_m(Vp)*mp));
	h = hp + dt*phi*(alpha_h(Vp)*(1-hp) - beta_h(Vp)*hp);
//	s = sp + dt*(alpha*(1-sp)*heaviside(Vp - VT) - beta*sp);
	s = sp + dt*phi*(alpha*fmax(Vp-VT,0)/(VM-VT) - beta*sp);
}


void Neuron_HH::step(double dt)
{
	step(dt, 0.0);
}
