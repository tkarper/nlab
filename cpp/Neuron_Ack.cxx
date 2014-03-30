#include <cmath>
#include <iostream>
#include "Neuron_Ack.h"


// Set model parameters (could be done in Python script)
double Neuron_Ack::VL	= -65;	
double Neuron_Ack::VNa	= 55;
double Neuron_Ack::VK	= -90;
double Neuron_Ack::Vh	= -20;
double Neuron_Ack::Vsyn	= -20;	// ???
double Neuron_Ack::gNa	= 52;
double Neuron_Ack::gNap	= 0.5;
double Neuron_Ack::gh	= 1.5;
double Neuron_Ack::gsyn	= 0.006;
double Neuron_Ack::gL	= 0.5;
double Neuron_Ack::gk	= 11;
double Neuron_Ack::CM	= 1.5;



Neuron_Ack::Neuron_Ack() :
V(VL), VP(VL),
mNa(0), mNaP(0), hNa(0), hNaP(0), n(0), nP(0), mNap(0), mNapP(0), mhf(0), mhfP(0), mhs(0), mhsP(0), msyn(0), msynP(0)
{
}


void Neuron_Ack::step(double t, double dt, double input)
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
	double mhfInf = 1/(1+exp((VP+79.2)/9.78));
	double mhsInf = 1/(1+exp((VP+71.3)/7.9));
	double taumhf = 0.51/(exp((VP-1.7)/10) + exp(-(VP+340)/52)) + 1;
	double taumhs = 5.6/(exp((VP-1.7)/14) + exp(-(VP+260)/43)) + 1;

	double ip = ipPos - ipNeg;
	V = VP + dt/CM*(ip - (
		(gNa*pow(mNaP,3)*hNaP + gNap*mNap)*(VP-VNa) + 
		(gk*pow(nP,4) /*+ gKs*mKsP*/)*(VP-VK) + 
		gh*(0.65*mhfP+0.35*mhsP)*(VP-Vh) + 
		gL*(VP-VL) + 
		gsyn*msynP*(VP-Vsyn)));
	mNa = mNaP + dt*(alpha_mNa*(1-mNaP) - beta_mNa*mNaP);
	hNa = hNaP + dt*(alpha_hNa*(1-hNaP) - beta_hNa*hNaP);
	mNap = mNapP + dt*(alpha_mNap*(1-mNapP) - beta_mNap*mNapP);
	n = nP + dt*(alpha_n*(1-nP) - beta_n*nP);
	mhf = mhfP + dt*(mhfInf - mhfP)/taumhf;
	mhs = mhsP + dt*(mhsInf - mhsP)/taumhs;
	msyn = msynP + dt*((VP>-20 ? 1.1 : 0.0)*(1-msynP) - 0.19*msynP);

	s = msyn;
//	s = sp + dt*(alpha*(1-sp)*heaviside(VP - VT) - beta*sp);
//	s = sp + dt*(alpha*fmax(VP-VT,0)/(VM-VT) - beta*sp);
}
