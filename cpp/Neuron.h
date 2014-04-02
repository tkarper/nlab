#ifndef NEURON_INC
#define NEURON_INC
#include <vector>

class Neuron
{
protected:
	double X, Y, Z; 	// position of the Neuron (in the brain)
	
	
public:
	std::vector<Neuron*>* con;		// list of connections
	std::vector<double>* weight;	    // weight of each connection
	
	double s, sp;	// current and previous spike value
	double I;		// external input
	double Vsyn;	// Synaptic reversal potential. Only meaningful for HH-type models
	double gsyn;	// Conductance for synapse
	
	Neuron();
	virtual ~Neuron();

	virtual void step(double t, double dt, double input) {}		// step using time-step dt with extra input
	void step(double t, double dt) { step(t, dt, 0.0); }	// step using time-step dt without extra input
	virtual void update() { sp = s; }					  // goes without saying
	void connect(Neuron* n1, double val);		  // make a new connection to neuron n1 with strength val
	std::vector<Neuron*>*	 getConnections(){return con;}

};


/*
	
		SIMPLE INTEGRATE-AND-FIRE NEURON MODEL
	
*/
class Neuron_IF: public Neuron
{
	
public:
	Neuron_IF();
	~Neuron_IF();
	void step(double t, double dt, double input);
};


typedef std::vector<Neuron*> nvector;
typedef std::vector<double> dvector;

#endif
