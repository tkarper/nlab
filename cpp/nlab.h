#ifndef NLAB_INC
#define NLAB_INC

#include <vector>
#include "Neuron.h"
#include "ConMat.h"


// Create connections
void connect_one_to_many(Neuron* from, nvector* to, double w);			// Create all-to-all connections with weight w.
void connect_with_matrix(nvector* from, nvector* to,  ConMat* M);		// Create connections according to ConMat.

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
#ifdef SWIG
%feature("autodoc","1");
%feature("docstring") "Create one-bump connectivity matrix for a strip cell network. 

nvector* left   - Leftward pointing neurons (decreasing index)
nvector* right  - Rightward pointing neurons (increasing index)
int      l      - Size of region NOT inhibited by a Neuron 
double   weight - The weight on each connection (should be NEGATIVE) 	 

";
#endif
ConMat* strip_matrix_OneBump(nvector* left, nvector* right, int l, double w);		


#endif
