#include <cmath>
#include "Neuron_ML.h"


// Set model parameters (could be done in Python script)
double Neuron_ML::phi = 1;	% Temperature factor
double Neuron_ML::gNa = 120;
double Neuron_ML::gK = 36;
double Neuron_ML::gL = 0.3;
double Neuron_ML::ENa = 50;
double Neuron_ML::EK = -77;
double Neuron_ML::EL = -54.4;
double Neuron_ML::CM = 1;
double Neuron_ML::alpha = 1;
double Neuron_ML::beta = 1;
double Neuron_ML::VT = -40;


// Functions appearing in the model
double alpha_n(double V) { return 0.01*(V+55)/(1 - exp(-(V+55)/10)); }
double alpha_m(double V) { return 0.1*(V+40)/(1 - exp(-(V+40)/10)); }
double alpha_h(double V) { return 0.07*exp(-(V+65)/20); }
double beta_n(double V) { return 0.125*exp(-(V+65)/80); }
double beta_m(double V) { return 4*exp(-(V+65)/18); }
double beta_h(double V) { return 1/(1+exp(-(V+35)/10)); }


Neuron_ML::Neuron_ML() :
V(-65), Vp(-65),
n(0.317), np(0.317), 
m(0.053), mp(0.053),
h(0.6), hp(0.6)
{
}


void Neuron_ML::step(double dt, double inp)
{
	// Compute input (synaptic potentials) coming from other neurons
	double ip = 0.0;
	for(int n=0; n<con->size(); n++)
	{
		Neuron* nr = con->at(n);
		ip += weight->at(n) * nr->sp;
	}
	
	V = Vp + dt*(I+inp+ip - ((gL*(Vp-EL) + gK*pow(np,4)*(Vp-EK) + gNa*pow(mp,3)*hp*(Vp-ENa))))/CM;
	n = np + dt*phi*(alpha_n(V)*(1-n) - beta_n(V)*n);
	m = mp + dt*phi*(alpha_m(V)*(1-m) - beta_m(V)*m);
	h = hp + dt*phi*(alpha_h(V)*(1-h) - beta_h(V)*h);
	s = sp + dt*(alpha*(1-sp)*fmax(Vp - VT, 0) - beta*sp);
}


void Neuron_ML::step(double dt)
{
	step(dt, 0.0);
}
