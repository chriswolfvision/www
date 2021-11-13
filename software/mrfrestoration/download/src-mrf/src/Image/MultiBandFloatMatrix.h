/**********************************************************
 * MultiBandFloatMatrix.h
 *
 * Christian Wolf, chriswolf@gmx.at
 * Beginn: 22.3.2005
 **********************************************************/
 
#ifndef _WOLF_MULTIBANDFLOATMATRIX_H_
#define _WOLF_MULTIBANDFLOATMATRIX_H_
 
// C++
#include <vector>

// From the main module
#include <CIL.h>

// From this module
#include "FloatMatrix.h"

class MultiBandFloatMatrix 
{
	public:
			
		// Constructors and destructor
		MultiBandFloatMatrix ();
		MultiBandFloatMatrix (int noBands, int xsize, int ysize);
#ifdef HAVE_LIBCFITSIO 		
		MultiBandFloatMatrix (const char *filename);
#endif 		
		~MultiBandFloatMatrix();		
		
		// Accessors
		inline int getNoBands () { return noBands; }
		inline FloatMatrix *&operator[](unsigned index);
		inline float get (int plane, int x, int y);
		inline void  set (int plane, int x, int y, float v);
		inline void add (FloatMatrix *);
		
		// The iterators
		typedef	vector<FloatMatrix *>::iterator iterator;
		typedef	vector<FloatMatrix *>::const_iterator const_iterator;
		iterator begin() { return matrices.begin(); }
		iterator end()	{ return matrices.end(); }
		const_iterator begin() const { return matrices.begin(); }
		const_iterator end()   const { return matrices.end(); }
		
		// returns the pointer to the plane and 
		// removes the plane from the structure
		inline FloatMatrix *unlink (int plane);
		
#ifdef HAVE_LIBCFITSIO 				
		void saveFITS (const char *filename);
#endif		 
	
	private:
	
		int noBands;
		int xsize, ysize;
		vector <FloatMatrix *> matrices;
};

/************************************************************************
 * Accessor
 ************************************************************************/
 
FloatMatrix *&MultiBandFloatMatrix::operator[](unsigned index)
{
	return matrices[index];
}

/************************************************************************
 * Accessor
 ************************************************************************/
 
float MultiBandFloatMatrix::get (int plane, int x, int y)
{
#ifdef PEDANTIC_CHECK_CODE
	if (x<0 || x>=(int)xsize || y<0 || y>=(int)ysize || plane>=noBands)
	{	
		CORE_DUMP_ON_EXCEPTION;
		throw EBounds3DError (plane,x,y,noBands,xsize,ysize);
	}
#endif
	return matrices[plane]->get(x,y);
}

/************************************************************************
 * Accessor
 ************************************************************************/
 
void MultiBandFloatMatrix::set (int plane, int x, int y, float val)
{
#ifdef PEDANTIC_CHECK_CODE
	if (x<0 || x>=(int)xsize || y<0 || y>=(int)ysize || plane>=noBands)
	{	
		CORE_DUMP_ON_EXCEPTION;
		throw EBounds3DError (plane,x,y,noBands,xsize,ysize);
	}
#endif
	return matrices[plane]->set(x,y,val);
}

/************************************************************************
 * Adds a plane
 ************************************************************************/

inline void MultiBandFloatMatrix::add (FloatMatrix *i) 
{
	matrices.push_back(i);
	if (xsize==-1)
		xsize = i->xsize;
	if (ysize==-1)
		ysize = i->ysize;
}

/************************************************************************
 * returns the pointer to the plane and
 * removes the plane from the structure
 ************************************************************************/
 
FloatMatrix *MultiBandFloatMatrix::unlink (int plane)
{
	FloatMatrix *rv = matrices[plane];
	matrices[plane] = NULL;
	return rv;
}

#endif
