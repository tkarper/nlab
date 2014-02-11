#include "nlab.h"
#include <iostream>

/*
	*********	UPDATE ROUTINES	*************
*/




void stepNetwork(nvector* nlist, double dt)
{
	Neuron* nr;
	for(int n= 0;n<nlist->size();n++)
	{
		nr = nlist->at(n);
		nr->step(0.1);
	}
		
}

void updateNetwork(nvector* nlist)
{
	for(int n= 0;n<nlist->size();n++)
		(nlist->at(n))->update();
	
}



/*


	*********** CONNECTION ROUTINES *************

*/
void connect_one_to_many(Neuron* from, nvector* to, double w)	
{
	Neuron* nr;
	for(int j=0;j<to->size();j++)
	{
		nr = to->at(j);
		nr->connect(from, w);	
	}
		
}


void connect_with_matrix(nvector* from, nvector* to,  ConMat* M)
{
	int ix;
	Neuron* frm; 
	Neuron* tm;
	ddvector* A;
	for(int n=0;n<M->W->size(); n++)
		for(int m = 0;m<M->W->at(n)->size()-1; m+=2)
		{
			A = M->W;
			ix = A->at(n)->at(m);
			frm = from->at(ix);
			tm = to->at(n);
			tm->connect(frm, 1.0);
		}
			
	
}





/*


	************ CREATE CONMAT ROUTINES ************

*/

ConMat* stripCell_onebump(nvector* s1, nvector* s2, int l, double w)
{
	int num_neuro = s1->size();
	ConMat* M = new ConMat(2*s1->size(), 2*s2->size());
	
	for(int n=0; n< num_neuro; n++)
	{
		for(int m=n+l;m<n+num_neuro;m++)
		{
			// Remove self-interaction
			if(m == n)
				continue;

			double ix = m;
			if(m>=num_neuro)
				ix = ix - num_neuro;	
			
		 	M->add(n,ix,w);
			M->add(n+num_neuro, ix, w);
			
			double ix2 = 2*n - m;
			if(ix2<0)
				ix2 = ix2 + num_neuro;
			
			M->add(n,ix2+num_neuro, w);
			M->add(n+num_neuro, ix2 + num_neuro, w);
		}
		
		
		M->add(n,n + num_neuro,w);
		M->add(n+num_neuro, n, w);
	}
	
	return M;
} 










