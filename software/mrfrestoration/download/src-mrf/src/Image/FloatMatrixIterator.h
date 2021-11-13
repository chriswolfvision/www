/*******************************************************************
 * FloatMatrixIterator
 *
 * An iterator for the FloatMatrix class, which can be used to access the
 * pixel values of a _COLUMN_ (due to to data layout it is _not_ a
 * line). The ++ and -- operators move downward and upward in a
 * column, they do not move into the next column and for speed
 * reasons _NO_ boundary checking is performed.
 *
 * If you need to change into another column, use the reset() method
 * to change the position of the iterator.
 *******************************************************************/

#ifndef	_WOLF_FLOATMATRIXITERATOR_H_
#define	_WOLF_FLOATMATRIXITERATOR_H_

#include "Image.h"

class FloatMatrixIterator	
{

	public:

		// Constructor
		FloatMatrixIterator () {};
		FloatMatrixIterator (FloatMatrix &f, int xpos, int ypos)
			:positionR (f.data1[xpos]+ypos)	{}

		// Reset the iterator to a new position
		void reset (FloatMatrix &f, int xpos,	int	ypos) 
		{
			positionR =	f.data1[xpos]+ypos;
		}

		// Access Operators
		float & operator* ()	{ return *positionR; }
		float get()				{ return *positionR; }

		// Move	the	cursor
		void operator++	() { ++positionR; }
		void operator--	() { --positionR; }

	private:

		float *positionR;
};

// The class method of the IMAGE class returning the iterator
inline FloatMatrix::iterator FloatMatrix::iterPos (int x, int y)
{
	return FloatMatrix::iterator (*this, x, y);
}


#endif
