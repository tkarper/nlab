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
	
	double s, sp;		// current and previous spike value
	double I;		// external input
	
	Neuron();
	virtual ~Neuron();

	virtual void step(double dt, double input_){}	  // step using time-step dt with extra input
	virtual void step(double dt){}				  	// step using time-step dt with extra input
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
	void step(double dt, double input_);
	void step(double dt);
	
	
};

typedef std::vector<Neuron*> nvector;
typedef std::vector<double> dvector;

#endif
