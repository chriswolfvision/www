/***********************************************************************
 * Confusion matrix
 * or: how to evaluate
 *
 * Author: Christian Wolf, christian.wolf@insa-lyon.fr
 * Begin: 14.6.2006
 *
 * Changelog:
 * 
 ***********************************************************************/

#ifndef _WOLF_CONFUSION_MATRIX_H_
#define _WOLF_CONFUSION_MATRIX_H_

// C++
#include <vector>

// The main libarary module
#include <CIL.h>

// From the MATHEMATICS module
#include "Matrix.h"

class Image;
class ConfusionMatrix;
ConfusionMatrix * CreateConfusionMatrix (Image *seg, Image *gt);

class ConfusionMatrix
{

	friend 		
		ostream & operator << (ostream &os, ConfusionMatrix &m);

	public:
	
		// Constructor
		ConfusionMatrix (unsigned int noSegLabels, unsigned int noGtLabels);
		ConfusionMatrix () 				{ mat.setZero(); }
		
		// The following "Constructor" can be found in the 
		// module "IMAGE PROCESSING":
		friend ConfusionMatrix * CreateConfusionMatrix (Image *seg, Image *gt);
		
		void inc(unsigned int reslabel, unsigned int gtlabel);
		
		// Accessors
		float & operator() (unsigned int r, unsigned int c);	
		unsigned int rows() 			{ return mat.rows(); }
		unsigned int columns() 			{ return mat.columns(); }
		unsigned char get_gv_seg(int i)	{ return gv_seg[i]; }
		unsigned char get_gv_gt(int j)	{ return gv_gt[j]; }
				
		float getAccuracy();
		float getRecall(int objectClassNumber);
		float getPrecision(int objectClassNumber);
		float sum()						{ return mat.sum(); }
		
		// Operators for combining multiple matrices
		void operator += (ConfusionMatrix &o);
		void operator /= (float div);
		
		vector<int> * searchMostProbableLabelPermutation();
		void normalize();
		
		void print (ostream &os, bool latex, int width=-1, int precision=-1);

	private: // METHODS
	
		ConfusionMatrix *getPermutatedColumns (vector<int> permutation);
	
	private: // DATA
	
		Matrix<float> mat;
		vector<unsigned char> gv_seg, gv_gt;
};

/**************************************************************
 * Constructor
 **************************************************************/

inline ConfusionMatrix::ConfusionMatrix (unsigned int noSegLabels, unsigned int noGtLabels) 
{
	mat.resize (noSegLabels, noGtLabels);
	gv_seg.resize (noSegLabels);
	gv_gt.resize  (noGtLabels);
	mat.setZero();
}

/**************************************************************
 * Add a new classification result: increase the corresponding
 * element in the matrix
 **************************************************************/

inline void ConfusionMatrix::inc(unsigned int reslabel, unsigned int gtlabel)
{
	mat(reslabel,gtlabel) = mat(reslabel,gtlabel)+1;
}

/**************************************************************
 * Access methods
 **************************************************************/

inline float & ConfusionMatrix::operator() (unsigned int r, unsigned int c)  
{
	return mat(r,c);
}

/**************************************************************
 * Divide the matrix by a scalar
 **************************************************************/

inline void ConfusionMatrix::operator /= (float div) 
{	
	mat /= div;
}
	
ostream & operator << (ostream &os, ConfusionMatrix &m);

#endif

