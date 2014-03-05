#include <iostream>
#include <stdexcept>
#include "nlab.h"

/*
	*********	UPDATE ROUTINES	*************
*/




void stepNetwork(nvector* nlist, double t, double dt)
{
	Neuron* nr;
	for(size_t n=0; n<nlist->size(); n++)
	{
		nr = nlist->at(n);
		nr->step(t, dt);
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


void connect_with_matrix2(nvector* from, nvector* to, ConMat2* M)
{
	for(size_t n=0; n<M->W.size(); n++)
	{
		for(size_t m=0; m<M->W[n].size(); m++)
		{
			ccvector& A = M->W;
			int ix = A[n][m].first;
			Neuron* frm = from->at(ix);
			Neuron* tm = to->at(n);
			tm->connect(frm, A[n][m].second);
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

			int ix = m;
			if(m>=num_neuro)
				ix = ix - num_neuro;	
			
		 	M->add(n,ix,w);
			M->add(n+num_neuro, ix, w);
			
			int ix2 = 2*n - m;
			if(ix2<0)
				ix2 = ix2 + num_neuro;
			
			M->add(n, ix2+num_neuro, w);
			M->add(n+num_neuro, ix2+num_neuro, w);
		}
		
		
		M->add(n, n+num_neuro, w);
		M->add(n+num_neuro, n, w);
	}
	
	return M;
}


ConMat2* strip_matrix_OneBump2(int num_neuro, int l, double w)
{
	if (l < 0 || l >= num_neuro)
		throw std::invalid_argument("in strip_matrix_OneBump2(): parameter l must lie between 0 and num_neuro-1");
		
	ConMat2* M = new ConMat2(2*num_neuro, 2*num_neuro);
	for(int n=0; n<num_neuro; n++)
	{
		for(int m=l; m<num_neuro; m++)
		{
			int iRight = (n+m) % num_neuro;
			int iLeft = (n-m+num_neuro) % num_neuro;
			// Avoid self-interaction
			if(iRight != n) {
				M->add(iRight, n, w);
				M->add(iRight+num_neuro, n, w);
			}
			if(iLeft != n) {
				M->add(iLeft, n+num_neuro, w);
				M->add(iLeft+num_neuro, n+num_neuro, w);
			}
		}
		// Add connection between opposing neurons
		M->add(n, n+num_neuro, w);
		M->add(n+num_neuro, n, w);
	}
	return M;
}


ConMat2* strip_matrix_local_OneBump(int num_neuro, int l, double w)
{
	if (l < 0 || l >= num_neuro)
		throw std::invalid_argument("in strip_matrix_OneBump2(): parameter l must lie between 0 and num_neuro-1");
	if (w < 0)
		std::cout << "Warning: strip cell weight 'w' is negative." << std::endl;
	double wIn = -w, wEx = w;
		
	ConMat2* M = new ConMat2(2*num_neuro, 2*num_neuro);
	for(int n=0; n<num_neuro; n++)
	{
		const int nR = n, nL = n+num_neuro;
		for(int m=1; m<=l; m++)
		{
			int iRLeft = (n-m+num_neuro) % num_neuro;
			int iRRight = (n+m) % num_neuro;
			int iLLeft = iRLeft + num_neuro;
			int iLRight = iRRight + num_neuro;
			// Add leftgoing connections
			M->add(iRLeft, nR, wIn);
//			M->add(iLLeft, nR, wIn);
//			M->add(iRLeft, nL, wEx);
			M->add(iLLeft, nL, wEx);
			// Add rightgoing connections
			M->add(iRRight, nR, wEx);
//			M->add(iLRight, nR, wEx);
//			M->add(iRRight, nL, wIn);
			M->add(iLRight, nL, wIn);
		}
		// Add exitatory connection between opposing neurons
		M->add(nL, nR, wEx);
		M->add(nR, nL, wEx);
	}
	return M;
}
 

// Up, Down, Right, Left
ConMat* gridcell_matrix_from_updrl(int num_in_strip)
{
	ConMat* M = new ConMat(num_in_strip*num_in_strip,4*num_in_strip);	
	int ix;
	for(int i = 0;i<num_in_strip; i++)
		for(int j=0;j< num_in_strip;j++)
		{
			ix = i*num_in_strip + j;
			M->add(ix, j                 , 1.0);	// up(j)
			M->add(ix, j +   num_in_strip, 1.0);	// down(j) 
			M->add(ix, i + 2*num_in_strip, 1.0); 	// right(i)
			M->add(ix, i + 3*num_in_strip, 1.0);	// left(i)
		}
	return M;
}


// A convenient wrapper for the process of calling a python function
#ifdef SWIG
%ignore phi(double, double, double, double, double, double, PyObject*);
#endif
double phi(double x1, double y1, double x2, double y2, double u, double v, PyObject* phi_p)
{
	if(!PyCallable_Check(phi_p))
	{
		std::cout << "Python object phi_p is not callable" << std::endl;
	}
	PyObject* pArgs = PyTuple_New(6);
	PyTuple_SetItem(pArgs, 0, PyFloat_FromDouble(x1));
	PyTuple_SetItem(pArgs, 1, PyFloat_FromDouble(y1));
	PyTuple_SetItem(pArgs, 2, PyFloat_FromDouble(x2));
	PyTuple_SetItem(pArgs, 3, PyFloat_FromDouble(y2));
	PyTuple_SetItem(pArgs, 4, PyFloat_FromDouble(u));
	PyTuple_SetItem(pArgs, 5, PyFloat_FromDouble(v));
	
	double val = PyFloat_AsDouble(PyObject_CallObject(phi_p, pArgs));
	Py_DECREF(pArgs);
	
	return val;
}
ConMat* gridcell_matrix_from_phi(int num_neuro, PyObject* phi_p)
{
	
	int N = (int) sqrt(num_neuro);
	
	// Find relevant R and shift l 
	double val = 2.0;
	int mx =1;
	while(fabs(val)>0.001 && mx < 100)
	{
		val = phi(0,0,mx,0,1.0,0.0,phi_p);
		mx++;
	}

	
	val = 2.0;
	int Mx=1;
	while(fabs(val)>0.001 && Mx < 100)
	{
		val = phi(0,0,Mx,0,0.0,0.0,phi_p);
		Mx++;
	}
	
	int Lx = fabs(Mx - mx);  // Shift in x direction
	int Rx  = fmin(Mx, mx);  // Radius in x direction

	val = 2.0;
	int my=1;
	while(fabs(val)>0.001 && my < 100)
	{
		val = phi(0,0,0,my,0.0,1.0,phi_p);
		my++;
	}
	
	val = 2.0;
	int My=1;
	while(fabs(val)>0.001 && My < 100)
	{
		val = phi(0,0,0,My,0.0,0.0,phi_p);
		My++;
	}
	
	int Ly = fabs(My - my);  // Shift in x direction
	int Ry  = fmin(My, my);  // Radius in x direction
	
	std::cout << "R=" << Rx << " with shift l=" << Lx << std::endl;
	std::cout << "R=" << Ry << " with shift l=" << Ly << std::endl;
	

	// Create a matrix holding the values
	double** up = new double*[2*(Ry+Ly)];
	double** down = new double*[2*(Ry+Ly)];
	double** right = new double*[2*(Ry+Ly)];
	double** left = new double*[2*(Ry+Ly)];
	
	for(int i =0;i<2*(Ry+Ly);i++)
	{
		up[i]    = new double[2*(Rx+Lx)];
		down[i]  = new double[2*(Rx+Lx)];
		right[i] = new double[2*(Rx+Lx)];
		left[i]  = new double[2*(Rx+Lx)];
		
		for(int j=0;j<2*(Rx+Lx);j++)
		{
			double x2 = (i-Ry-Ly);
			double y2 = (j-Rx-Lx);
			
			up[i][j]    = phi(0.0, 0.0, x2, y2,  0.0, 1.0,phi_p);
			down[i][j]  = phi(0.0, 0.0, x2, y2,  0.0, -1.0,phi_p);
			right[i][j] = phi(0.0, 0.0, x2, y2,  1.0, 0.0,phi_p);
			left[i][j]  = phi(0.0, 0.0, x2, y2, -1.0, 0.0,phi_p);
			
		}
		
	}
	
	int sz = 2*(Rx+Lx)*(Rx+Lx) + 2*(Ry+Ly)*(Ry+Ly);
	ConMat* M = new ConMat(4*num_neuro, sz);
	
	
	
	for(int i=0;i<N;i++)
	{
		std::cout << (100.0*i)/(N+0.0) << "%" << std::endl;
		for(int j=0;j<N;j++)
		{
			int to_index = i*N + j;
		
			for(int k=1; k< 2*(Ly + Ry);k++)
				for(int l=1; l< 2*(Lx + Rx);l++)
			{
				int j2 = j+l-(Ly+Ry);
				int i2 = i+k-(Lx+Rx);
				if(j2 < 0)
					j2 += N;
				if(j2 >= N)
					j2 -= N;
				if(i2 < 0)
					i2 += N;
				if(i2 >= N)
					i2 -= N;
			 
				int from_index = i2*N + j2;
			
				// Add connections to up
				if(to_index != from_index) 
				M->add(to_index, from_index, up[k][l]);
				M->add(to_index, from_index + num_neuro, down[k][l]);
				M->add(to_index, from_index + 2*num_neuro, right[k][l]);
				M->add(to_index, from_index + 3*num_neuro, left[k][l]);
			
				// Add connection to down
				M->add(to_index+num_neuro, from_index, up[k][l]);
				if(to_index != from_index) 
				M->add(to_index+num_neuro, from_index + num_neuro, down[k][l]);
				M->add(to_index+num_neuro, from_index + 2*num_neuro, right[k][l]);
				M->add(to_index+num_neuro, from_index + 3*num_neuro, left[k][l]);

				// Add connection to right
				M->add(to_index+2*num_neuro, from_index, up[k][l]);
				M->add(to_index+2*num_neuro, from_index + num_neuro, down[k][l]);
				if(to_index != from_index)
				M->add(to_index+2*num_neuro, from_index + 2*num_neuro, right[k][l]);
				M->add(to_index+2*num_neuro, from_index + 3*num_neuro, left[k][l]);

				// Add connections to left
				M->add(to_index+3*num_neuro, from_index, up[k][l]);
				M->add(to_index+3*num_neuro, from_index + num_neuro, down[k][l]);
				M->add(to_index+3*num_neuro, from_index + 2*num_neuro, right[k][l]);
				if(to_index != from_index)
				M->add(to_index+3*num_neuro, from_index + 3*num_neuro, left[k][l]);

			
			}
		
		}
	}
	
	for(int i=0;i< Ry+Ly;i++)
	{
		delete [] up[i];
		delete [] down[i];
		delete [] right[i];
		delete [] left[i];
	}
	delete [] up;
	delete [] down;
	delete [] right;
	delete [] left;

	
	return M;
}
// ConMat* gridcell_matrix_from_phi(int num_neuro, PyObject* phi_p)
// {
// 	
// 	double Xdir [4] = {0.0, 0.0, 1.0,-1.0};
// 	double Ydir [4] = {1.0, -1.0, 0.0,0.0};
// 	int N = (int) sqrt(num_neuro);
// 	
// 	
// 	
// 	ConMat* M = new ConMat(4*num_neuro, N);
// 	
// 	
// 	if(!PyCallable_Check(phi_p))
// 	{
// 		std::cout << "Python object phi_p is not callable" << std::endl;
// 	}
// 	
// 	PyObject* pArgs = PyTuple_New(6);
// 	for(int i = 0;i<N; i++)
// 	{
// 		std::cout << i*100.0/(N+0.0) << "%" << std::endl;
// 		for(int j =0;j<N; j++)
// 		{
// 			int to = i*N + j;
// 			double x1 = j;
// 			double y1 = i;
// 			PyTuple_SetItem(pArgs, 0, PyFloat_FromDouble(x1));
// 			PyTuple_SetItem(pArgs, 1, PyFloat_FromDouble(y1));
// 		
// 			for(int k=0;k<N;k++)
// 				for(int l=0;l<N;l++)
// 			{
// 				int from = k*N + l;
// 				double x2 = l;
// 				double y2 = k;
// 			
// 				PyTuple_SetItem(pArgs, 2, PyFloat_FromDouble(x2));
// 				PyTuple_SetItem(pArgs, 3, PyFloat_FromDouble(y2));
// 			
// 				for(int r = 0;r<4;r++)
// 					for(int s=0;s<4;s++)
// 					{
// 						PyTuple_SetItem(pArgs, 4, PyFloat_FromDouble(Xdir[s]));
// 						PyTuple_SetItem(pArgs, 5, PyFloat_FromDouble(Ydir[s]));
// 	
// 						double val = PyFloat_AsDouble(PyObject_CallObject(phi_p, pArgs));
// 						if(val!=0.0 && (to+r*num_neuro != from+s*num_neuro) )
// 							M->add(to+r*num_neuro, from+s*num_neuro, val);		
// 					} 		
// 			}	
// 		}
// 	}
// 	Py_DECREF(pArgs);
// 	
// 	return M;
// }






