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
//double alpha_n(double V) { return (fabs(V+55)<tol) ? 0.1 : 0.01*(V+55)/(1 - exp(-(V+55)/10)); }
double alpha_mNa(double V) { 0.1*(V+23)/(1 - exp(-(V+23)/10)); }
double alpha_n(double V) { 0.01*(V+27)/(1 - exp(-(V+27)/10)); }
double alpha_mNap(double V) { 1/(0.15*(1 + exp(-(V+38)/6.5))); }
double beta_mNa(double V) { return 4*exp(-(V+48)/18); }
double beta_hNa(double V) { return 1/(1+exp(-(V+7)/10)); }
double beta_n(double V) { return 0.125*exp(-(V+37)/80); }
double beta_mNap(double V) { return exp(-(V+38)/6.5)/(0.15*(1+exp(-(V+38)/6.5))); }

int heaviside(double V) { return V > 0 ? 1 : 0; }


Neuron_HH::Neuron_HH() :
V(-65), VP(-65),
n(0.317), np(0.317), 
m(0.053), mp(0.053),
h(0.6), hp(0.6)
{
}


void Neuron_HH::step(double t, double dt, double input)
{
	// Compute input (synaptic potentials) coming from other neurons
	double ipPos = 0, ipNeg = 0;
	for(size_t i=0; i<con->size(); i++)
	{
		Neuron* nr = con->at(i);
		double signal = weight->at(i) * nr->sp;
		if (signal > 0)
			ipPos += signal;
		else
			ipNeg -= signal;
	}
	if (I > 0)
		ipPos += I;
	else
		ipNeg -= I;
	if (input > 0)
		ipPos += input;
	else
		ipNeg -= input;

	double ip = ipPos - ipNeg;
	V = VP + dt*(ip - (gL*(VP-EL) + (gK*pow(np,4) + gKs*mKs)*(VP-EK) + (gNa*pow(mNaP,3)*hNaP + gh*(0.65*mhfP+0.35*mhsP)*(VP-Eh) + gNap*mNap)*(VP-ENa) + gsyn*msynP*(VP-Esyn)))/CM;
	mNa = mNaP + dt*phi*(alpha_mNa(VP)*(1-mNaP) - beta_mNa(VP)*mNaP);
	hNa = hNaP + dt*phi*(alpha_hNa(VP)*(1-hNaP) - beta_hNa(VP)*hNaP);
	mNap = mNapP + dt*phi*(alpha_mNap(VP)*(1-mNapP) - beta_mNap(VP)*mNapP);
	n = nP + dt*phi*(alpha_n(VP)*(1-nP) - beta_n(VP)*nP);
	
//	s = sp + dt*(alpha*(1-sp)*heaviside(VP - VT) - beta*sp);
	s = sp + dt*phi*(alpha*fmax(VP-VT,0)/(VM-VT) - beta*sp);
}
