/*******************************************************************
 * ImageIterator
 *
 * An iterator for the image class, which can be used to access the
 * pixel values of a _COLUMN_ (due to to data layout it is _not_ a
 * line). The ++ and -- operators move downward and upward in a
 * column, they do not move into the next column and for speed
 * reasons _NO_ boundary checking is performed.
 *
 * If you need to change into another column, use the reset() method
 * to change the position of the iterator.
 *******************************************************************/

#ifndef	_WOLF_IMAGEITERATOR_H_
#define	_WOLF_IMAGEITERATOR_H_

#include "Image.h"

class ImageIterator	{

	public:

		// Constructor
		ImageIterator () {};
		ImageIterator (Image &image, int xpos, int ypos)
			:positionR (image.R[xpos]+ypos)	{}

		// Reset the iterator to a new position
		void reset (Image &image, int xpos,	int	ypos) {
			positionR =	image.R[xpos]+ypos;
		}

		// Access Operator
		byte & operator* ()	{ return *positionR; }
		byte get()			{ return *positionR; }

		// Move	the	cursor
		void operator++	() { ++positionR; }
		void operator--	() { --positionR; }

	private:

		byte *positionR;
};

// The class method of the IMAGE class returning the iterator
inline Image::iterator Image::iterPos (int x, int y)
{
	return Image::iterator (*this, x, y);
}


#endif

