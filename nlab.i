%{
    #define SWIG_FILE_WITH_INIT
%}

%module nlab

%include "typemaps.i"

%{
	#include <numpy/arrayobject.h>
	#include <iostream>
    #include "Neuron.h"
    #include "nlab.h"
    #include "ConMat.h"
%}

%init %{
    import_array();
    %}

%typemap(in) nvector*
{
	PyArrayObject* arr = (PyArrayObject*) PyArray_FROM_O($input);
	PyObject** list = (PyObject**) PyArray_DATA(arr);
	const int n = PyArray_DIM(arr, 0);
	void* vptr = 0;
	Neuron* nr = 0;

	nvector* v = new nvector();
	v->reserve(n);

	for(int i =0;i< n; i++)
	{
		SWIG_ConvertPtr(list[i], &vptr, SWIGTYPE_p_Neuron, 0);
		nr = (Neuron*) vptr;
		v->push_back(nr);
	}
	
	$1 = v;

}

%typemap(freearg) nvector*
{
	delete $1;
}

%include "Neuron.h"
%include "ConMat.h"
%include "nlab.h"

%ignore connect(Neuron* from, nvector* to, double w);
%ignore connect(nvector* from, nvector* to, ConMat* con);

%inline{

PyObject* get_spike_rates(nvector* nr)
{
	  double* S = new double[nr->size()];
	  for(int i =0;i<nr->size();i++)
	  	S[i] = (nr->at(i))->sp;
	
	  npy_intp* npy_size = new npy_intp[1];
      npy_size[0] = (npy_intp) nr->size();
	  PyArrayObject* array =  (PyArrayObject*) PyArray_SimpleNewFromData(1, npy_size, NPY_DOUBLE, (void*) S); 
	  return PyArray_Return(array);
}



nvector* nconvert_p2c(PyObject* parray)
{
	PyArrayObject* arr = (PyArrayObject*) PyArray_FROM_O(parray);
	PyObject** list = (PyObject**) PyArray_DATA(arr);
	int n = PyArray_DIM(arr, 0);
	void* vptr = 0;
	Neuron* nr = 0;

	nvector* carray = new nvector();
	carray->reserve(n);

	for(int i =0;i< n; i++)
	{
		SWIG_ConvertPtr(list[i], &vptr, SWIGTYPE_p_Neuron, 0);
		nr = (Neuron*) vptr;
		carray->push_back(nr);
	}
	
	return carray;
}


void connect_with_matrix(PyObject* from_, PyObject* to_, PyObject* con_)
{
	nvector* to = nconvert_p2c(to_);
	nvector* from = nconvert_p2c(from_);
	
	void* vptr = 0;
	ConMat* con = 0;
	int res = SWIG_ConvertPtr(con_, &vptr, SWIGTYPE_p_ConMat, 0);
	con = (ConMat*) vptr;
	
	connect(from,to,con);
	
}

void connect_one_to_many(PyObject* from_, PyObject* to_, double val)
{
	nvector* to = nconvert_p2c(to_);
	
	void* vptr = 0;
	Neuron* from = 0;
	SWIG_ConvertPtr(from_, &vptr, SWIGTYPE_p_Neuron, 0);
	from = (Neuron*) vptr;
	
	connect(from,to,val);
}



/*
	
		THIS WORKS!!!! :) Please leave as this was the breaking-point
	
*/
	
void mongo(PyObject* inn)
{
	PyArrayObject* arr = (PyArrayObject*) PyArray_FROM_O(inn);
	void* vptr = 0;
	Neuron* nr = 0;
	PyObject** list = (PyObject**) PyArray_DATA(arr);

	for(int i =0;i< PyArray_DIM(arr,0); i++)
	{
		SWIG_ConvertPtr(list[i], &vptr, SWIGTYPE_p_Neuron, 0);
		nr = (Neuron*) vptr;
	}
		
	
}


}