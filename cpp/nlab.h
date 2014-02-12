#ifndef NLAB_INC
#define NLAB_INC
#include <Python.h>
#include <vector>
#include "Neuron.h"
#include "ConMat.h"


// Create connections
void connect_one_to_many(Neuron* from, nvector* to, double w);			// Create all-to-all connections with weight w.
void connect_with_matrix(nvector* from, nvector* to, ConMat* M);		// Create connections according to ConMat.
void connect_with_matrix2(nvector* from, nvector* to, ConMat2* M);		// Create connections according to ConMat.

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
num_neuro:	int 		 - Number of neurons in each strip (e.g num of left neurons) 
l:		int              - Size of region NOT inhibited by a Neuron 
weight:	double           - The weight on each connection (should be NEGATIVE) 	 
	
Returns connection Matrix of the form:

M_ij = connection from Strip(j) to Strip(i)
where Strip = [Left, Right]

";
#endif
ConMat* strip_matrix_OneBump(int num_neuro, int l, double w);
ConMat2* strip_matrix_OneBump2(int num_neuro, int l, double w);


#ifdef SWIG
%feature("autodoc","1");
%feature("docstring") "Create a gridcell connectivity matrix from strip cells (Up, Down, Right, Left)

Returning Connection Matrix of form:

GridCell | Up  | Down  |  Right  | Left  |
-----------------------------------------
	0   |
	1   |	
	2   |
		
	";
#endif
ConMat* gridcell_matrix_from_updrl(int num_neuro_in_strip);



#ifdef SWIG
%feature("autodoc","1");
%feature("docstring") "Create a gridcell connectivity matrix from phi(x,y,v)

Returning Connection Matrix of form:

GridCell    | Layer-Up  | Layer-Down  |  Layer-Right  | Layer-Left |
--------------------------------------------------------------------
Layer-Up    |           
Layer-Down  |           
Layer-Right |           
Layer-Left  |
	
	";
#endif


ConMat* gridcell_matrix_from_phi(int num_neuro_in_layer, PyObject* phi_p);











#endif
