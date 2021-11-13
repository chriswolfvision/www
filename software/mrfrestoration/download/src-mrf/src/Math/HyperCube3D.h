/***********************************************************************
 * A 5D hyper cube for data storage
 * 
 * Author: Christian Wolf
 * Begin: 11.4.2007
 ***********************************************************************/
 
#ifndef _WOLF_HYPERCUBE3D_H_
#define _WOLF_HYPERCUBE3D_H_

#include <ostream>

using namespace std;

#define NO_DIM				3

template <class T>
class HyperCube3D
{
	public:
	
		HyperCube3D (unsigned int dim0, unsigned int dim1, unsigned int dim2);
		~HyperCube3D();
		
		// Accessors
		T get (unsigned int dim0, unsigned int dim1, unsigned int dim2);		
		void set (unsigned int dim0, unsigned int dim1, unsigned int dim2, T value);
		void add (unsigned int dim0, unsigned int dim1, unsigned int dim2, T value)
			{ set(dim0,dim1,dim2,get(dim0,dim1,dim2)+value); }
		
		void setToValue(T val);
		unsigned int getDim(unsigned dim_ind) { return dimensions[dim_ind]; }
		
		void printFlat(ostream &os);
	
	private:
	
        T ***data;
        T *flat_shortcut;

	unsigned int dimensions[NO_DIM];
};

/***********************************************************************
 * The constructor
 ***********************************************************************/

template <class T>
HyperCube3D<T>::HyperCube3D (unsigned int dim0, unsigned int dim1, unsigned int dim2)
{
	typedef T*   p1T;
	typedef T**  p2T;
	T *C;
	p2T B;
	
	// Allocate the data and the management arrays
	data = new p2T [dim0];
	B = new p1T [dim1*dim0];
	C = new T [dim2*dim1*dim0];
	
	// Set up the pointer structure
	for (unsigned int i=0; i<dim0; ++i)
		data[i] = B+i*dim1;
	for (unsigned int i=0; i<dim0*dim1; ++i)
		B[i]    = C+i*dim2;
		
	flat_shortcut = C;
	dimensions[0] = dim0;
	dimensions[1] = dim1;
	dimensions[2] = dim2;
}

/***********************************************************************
 * The destructor
 ***********************************************************************/

template <class T>
HyperCube3D<T>::~HyperCube3D ()
{
	delete [] data[0][0];
	delete [] data[0];
	delete [] data;
}

/***********************************************************************
 * The accessors
 ***********************************************************************/

template <class T>
T HyperCube3D<T>::get (unsigned int dim0, unsigned int dim1, unsigned int dim2)
{
	return data[dim0][dim1][dim2];
}

template <class T>
void HyperCube3D<T>::set (unsigned int dim0, unsigned int dim1, unsigned int dim2, T value)
{
	data[dim0][dim1][dim2] = value;
}

/***********************************************************************
 * The destructor
 ***********************************************************************/

template <class T>
void HyperCube3D<T>::setToValue (T val)
{
	unsigned int cnt=dimensions[0];
	for (unsigned int i=1; i<NO_DIM; ++i)
		cnt*=dimensions[i];
			
	for (unsigned int i=0; i<cnt; ++i)
		flat_shortcut[i]=val;
}

/***********************************************************************
 * DEBUG:
 * Print all values in a flat manner
 ***********************************************************************/

template <class T>
void HyperCube3D<T>::printFlat(ostream & os)
{
    unsigned int cnt=dimensions[0];
    for (unsigned int i=1; i<NO_DIM; ++i)
        cnt*=dimensions[i];
    
    for (int i=0; i<cnt; ++i)
        os << flat_shortcut[i] << " ";
    os << endl;
}    

#undef NO_DIM

#endif
