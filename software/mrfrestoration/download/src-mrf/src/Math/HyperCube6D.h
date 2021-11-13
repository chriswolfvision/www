/***********************************************************************
 * A 6D hyper cube for data storage
 * 
 * Author: Christian Wolf
 * Begin: 16.5.2007
 ***********************************************************************/
 
#ifndef _WOLF_HYPERCUBE6D_H_
#define _WOLF_HYPERCUBE6D_H_

#define NO_DIM				6

template <class T>
class HyperCube6D
{
	public:
	
		HyperCube6D (unsigned int dim0, unsigned int dim1, unsigned int dim2, 
					 unsigned int dim3, unsigned int dim4, unsigned int dim5);
		~HyperCube6D();
		
		// Accessors
		T get (unsigned int dim0, unsigned int dim1, unsigned int dim2, 
			   unsigned int dim3, unsigned int dim4, unsigned int dim5);
		void set (unsigned int dim0, unsigned int dim1, unsigned int dim2, 
			      unsigned int dim3, unsigned int dim4, unsigned int dim5,
			      T value);
		
		void setToValue(T val);
		unsigned int getDim(unsigned dim_ind) { return dimensions[dim_ind]; }
	
	private:
	
	T ******data;
	T *flat_shortcut;

	unsigned int dimensions[NO_DIM];
};

/***********************************************************************
 * The constructor
 ***********************************************************************/

template <class T>
HyperCube6D<T>::HyperCube6D (
	unsigned int dim0, unsigned int dim1, unsigned int dim2, 
	unsigned int dim3, unsigned int dim4, unsigned int dim5)
{
	typedef T*    p1T;
	typedef T**   p2T;
	typedef T***  p3T;
	typedef T**** p4T;
	typedef T*****p5T;
	p1T F;
	p2T E;
	p3T D;
	p4T C;
	p5T B;
	
	// Allocate the data and the management arrays
	data = new p5T [dim0];
	B = new p4T [dim1*dim0];
	C = new p3T [dim2*dim1*dim0];
	D = new p2T [dim3*dim2*dim1*dim0];
	E = new p1T [dim4*dim3*dim2*dim1*dim0];
	F = new   T [dim5*dim4*dim3*dim2*dim1*dim0];
	
	// Set up the pointer structure
	for (unsigned int i=0; i<dim0; ++i)
		data[i] = B+i*dim1;
	for (unsigned int i=0; i<dim0*dim1; ++i)
		B[i]    = C+i*dim2;
	for (unsigned int i=0; i<dim0*dim1*dim2; ++i)
		C[i] =    D+i*dim3;
	for (unsigned int i=0; i<dim0*dim1*dim2*dim3; ++i)
		D[i] =    E+i*dim4;
	for (unsigned int i=0; i<dim0*dim1*dim2*dim3*dim4; ++i)
		E[i] =    F+i*dim5;
		
	flat_shortcut = F;
	dimensions[0] = dim0;
	dimensions[1] = dim1;
	dimensions[2] = dim2;
	dimensions[3] = dim3;
	dimensions[4] = dim4;
	dimensions[5] = dim5;
}

/***********************************************************************
 * The destructor
 ***********************************************************************/

template <class T>
HyperCube6D<T>::~HyperCube6D ()
{
	delete [] data[0][0][0][0][0];
	delete [] data[0][0][0][0];
	delete [] data[0][0][0];
	delete [] data[0][0];
	delete [] data[0];
	delete [] data;
}

/***********************************************************************
 * The accessors
 ***********************************************************************/

template <class T>
T HyperCube6D<T>::get (
	unsigned int dim0, unsigned int dim1, unsigned int dim2, 
	unsigned int dim3, unsigned int dim4, unsigned int dim5)
{
	return data[dim0][dim1][dim2][dim3][dim4][dim5];
}

template <class T>
void HyperCube6D<T>::set (
	unsigned int dim0, unsigned int dim1, unsigned int dim2, 
	unsigned int dim3, unsigned int dim4, unsigned int dim5,
	T value)
{
	data[dim0][dim1][dim2][dim3][dim4][dim5] = value;
}

/***********************************************************************
 * The destructor
 ***********************************************************************/

template <class T>
void HyperCube6D<T>::setToValue (T val)
{
	unsigned int cnt=dimensions[0];
	for (unsigned int i=1; i<NO_DIM; ++i)
		cnt*=dimensions[i];
			
	for (unsigned int i=0; i<cnt; ++i)
		flat_shortcut[i]=val;
}

#undef NO_DIM

#endif
