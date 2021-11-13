/***************************************************************************
                          FloatPyramid.h  -  description
                             -------------------    
 ***************************************************************************/

#ifndef _WOLF_FLOATPYRAMID_H_ 
#define _WOLF_FLOATPYRAMID_H_
 
// C++
#include <vector>

// From the IMAGE library
#include <FloatMatrix.h>

// From the IMAGE PROCESSING library
#include <FilterMask.h>

enum PyramidType {
	PYRT_FIL33_RED22=0,
	PYRT_FIL53_RED22,
	PYRT_FIL35_RED22,
	PYRT_FIL73_RED22,
	PYRT_FIL37_RED22,
	PYRT_FIL77_RED22,
	PYRT_FIL77RD_RED22,
	PYRT_FIL77LD_RED22,
	PYRT_MAX_RED22
};

class FloatPyramid {

	public:

		// Constructor and destructor	
		FloatPyramid (FloatMatrix *I, int nLevels, PyramidType t);
		~FloatPyramid ();
		
		// Accessors
		inline FloatMatrix *&operator[](unsigned level);			
		inline int levels() { return noLevels; }
		
		// Change the base image and rebuild the pyramid
		void changeImage (FloatMatrix *I);		
		
		void build();
		void collapse();		
		
		// returns the pointer to the level and 
		// removes the level from the structure
		inline FloatMatrix *unlink (int level);
				
	private:	// METHODS				
	
		void filterForRed2x2 (FloatMatrix &parent, FloatMatrix &buf, FilterMask &fm);
		
		// Isotropic and ani-isotropic gaussian filtering,
		// but always the same isotropic reduction
		void gaussian33red22 (vector <FloatMatrix *> &levelVector, int begLevel, int endLevel);
		void gaussian53red22 (vector <FloatMatrix *> &levelVector, int begLevel, int endLevel);
		void gaussian35red22 (vector <FloatMatrix *> &levelVector, int begLevel, int endLevel);
		void gaussian73red22 (vector <FloatMatrix *> &levelVector, int begLevel, int endLevel);
		void gaussian37red22 (vector <FloatMatrix *> &levelVector, int begLevel, int endLevel);
		void gaussianFMred22 (vector <FloatMatrix *> &levelVector, int begLevel, int endLevel,
			FilterMask &fm);
			
		// The maximum pyramid
		void maxRed22 (vector <FloatMatrix *> &vec, int begLevel, int endLevel);
					
	private:	// DATA
		
		// The type of pyramid
		// (determines filter and reduction function)
		PyramidType type;
	
		// The levels
		vector<FloatMatrix *> images;		
		int noLevels;
				
		FloatMatrix *reductionBuffer;	

		// The masks for the gaussian diagonal image reduction		
		// (if chosen)
		FilterMask mask;
};
/************************************************************************
 * Accessor
 ************************************************************************/
 
FloatMatrix *&FloatPyramid::operator[](unsigned index)
{
	return images[index];
}


/************************************************************************
 * returns the pointer to the level and
 * removes the level from the structure
 ************************************************************************/
 
FloatMatrix *FloatPyramid::unlink (int level)
{
#ifdef CHECK_CODE
	if (level<0 || level>=noLevels) 
		ERR_THROW ("FloatPyramid::unlink(): invalid level=" << level 
			<< " in pyramid with " << noLevels << " levels.\n" << endl);
#endif
	FloatMatrix *rv = images[level];
	images[level] = NULL;
	return rv;
}

#endif
