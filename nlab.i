%{
    #define SWIG_FILE_WITH_INIT
%}

%module nlab

%include "typemaps.i"

%{
	#include <numpy/arrayobject.h>
	#include <iostream>
	#include <string>
    #include "Neuron.h"
    #include "Neuron_HH.h"
    #include "Neuron_ML.h"
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
%include "Neuron_ML.h"
%include "Neuron_HH.h"
%include "ConMat.h"
%include "nlab.h"


%inline{


PyObject* get_neuron_entry(nvector* nr, const char* vn)
{
	std::string varname = vn;
	npy_intp* npy_size = new npy_intp[1];
	npy_size[0] = (npy_intp) nr->size();
	PyArrayObject* array =  (PyArrayObject*) PyArray_SimpleNew(1, npy_size, NPY_DOUBLE); 
	double * rr = (double*) PyArray_DATA(array);
	if (varname == "sp")
		for(size_t i=0; i<nr->size(); i++)
			rr[i] = (nr->at(i))->sp;
	else if (varname == "Vp")
		for(size_t i=0; i<nr->size(); i++)
			rr[i] = dynamic_cast<Neuron_HH*>(nr->at(i))->Vp;
	else if (varname == "np")
		for(size_t i=0; i<nr->size(); i++)
			rr[i] = dynamic_cast<Neuron_HH*>(nr->at(i))->np;
	else if (varname == "mp")
		for(size_t i=0; i<nr->size(); i++)
			rr[i] = dynamic_cast<Neuron_HH*>(nr->at(i))->mp;
	else
		throw std::invalid_argument("In get_spike_rates: unrecognized varname " + varname);
	return PyArray_Return(array);
}


PyObject* get_spike_rates(nvector* nr)
{
	return get_neuron_entry(nr, "sp");
}


nvector* nconvert_p2c(PyObject* parray)
{
	PyArrayObject* arr = (PyArrayObject*) PyArray_FROM_O(parray);
	PyObject** list = (PyObject**) PyArray_DATA(arr);
	int n = PyArray_DIM(arr, 0);
	void* vptr = 0;

	nvector* carray = new nvector();
	carray->reserve(n);

	for(int i =0;i< n; i++)
	{
		SWIG_ConvertPtr(list[i], &vptr, SWIGTYPE_p_Neuron, 0);
		Neuron* nr = (Neuron*) vptr;
		carray->push_back(nr);
	}
	
	return carray;
}

/*
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
	SWIG_ConvertPtr(from_, &vptr, SWIGTYPE_p_Neuron, 0);
	Neuron* from = (Neuron*) vptr;
	
	connect(from,to,val);
}


*/
/*
	
		THIS WORKS!!!! :) Please leave as this was the breaking-point
	
*/
	
void mongo(PyObject* inn)
{
	PyArrayObject* arr = (PyArrayObject*) PyArray_FROM_O(inn);
	PyObject** list = (PyObject**) PyArray_DATA(arr);

	for(int i =0;i< PyArray_DIM(arr,0); i++)
	{
		void* vptr = 0;
		SWIG_ConvertPtr(list[i], &vptr, SWIGTYPE_p_Neuron, 0);
		Neuron* nr = (Neuron*) vptr;
	}
}


}
