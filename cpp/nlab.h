#ifndef NLAB_INC
#define NLAB_INC

#include <vector>
#include "Neuron.h"
#include "ConMat.h"


// Create connections
void connect(Neuron* from, nvector* to, double w);						// Create all-to-all connections with weight w.
// void connect(std::vector<Neuron*> from, std::vector<Neuron*> to, double w);						// Create all-to-all connections with weight w.
// void connect(std::vector<Neuron*> from, std::vector<Neuron*> to, std::vector<double> w);		// Create all-to-all connections with weights w.
void connect(nvector* from, nvector* to,  ConMat* M);											// Create connections according to ConMat.

// Update a list of Neurons
void stepNetwork(nvector* nlist, double dt);   // Solve for s
void updateNetwork(nvector* nlist);			   // sp = s





// Create ConMat methods
/*
	Create a ConMat W for a stripcell network with onebump. Headcells must be put in separately. 
	s1 - Neurons connected to h1
	s2 - Neurons connected to h2
	weight - goes without saying
*/
ConMat* stripCell_onebump(nvector* s1, nvector* s2, int l, double w);		
//ConMat gcMat_Strips(Neuron** s1, Neuron** s2, double weight);

#endif
