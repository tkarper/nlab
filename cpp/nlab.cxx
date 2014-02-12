#include "nlab.h"
#include <iostream>

/*
	*********	UPDATE ROUTINES	*************
*/




void stepNetwork(nvector* nlist, double dt)
{
	Neuron* nr;
	for(size_t n=0; n<nlist->size(); n++)
	{
		nr = nlist->at(n);
		nr->step(0.1);
	}
		
}

void updateNetwork(nvector* nlist)
{
	for(size_t n=0; n<nlist->size(); n++)
		(nlist->at(n))->update();
}



/*


	*********** CONNECTION ROUTINES *************

*/
void connect_one_to_many(Neuron* from, nvector* to, double w)	
{
	for(size_t j=0; j<to->size(); j++)
	{
		Neuron* nr = to->at(j);
		nr->connect(from, w);	
	}
		
}


void connect_with_matrix(nvector* from, nvector* to,  ConMat* M)
{
	for(size_t n=0; n<M->W->size(); n++)
	{
		for(size_t m=0; m<M->W->at(n)->size()-1; m+=2)
		{
			ddvector* A = M->W;
			int ix = A->at(n)->at(m);
			Neuron* frm = from->at(ix);
			Neuron* tm = to->at(n);
			tm->connect(frm, A->at(n)->at(m+1));
		}
	}
}





/*


	************ CREATE CONMAT ROUTINES ************

*/

ConMat* strip_matrix_OneBump(int num_neuro, int l, double w)
{

	ConMat* M = new ConMat(2*num_neuro, 2*num_neuro);
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


ConMat* gridcell_matrix_from_updrl(int num_in_strip)
{
	ConMat* M = new ConMat(num_in_strip*num_in_strip,num_in_strip);	
	int ix;
	for(int i = 0;i<num_in_strip; i++)
		for(int j=0;j< num_in_strip;j++)
	{
		ix = i*num_in_strip + j;
		M->add(ix, j,0.25);
		M->add(ix, j,0.25);
		M->add(ix, i,0.25);
		M->add(ix, i,0.25);
	}
	return M;
}







