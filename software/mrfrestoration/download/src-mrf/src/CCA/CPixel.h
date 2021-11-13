#ifndef	_WOLF_CPIXEL_H_
#define	_WOLF_CPIXEL_H_

typedef	int	coord;

// **************************************************************
// A Pixel "old	style" from	the	C-time.	It incorporates	a pointer
// to the next pixel for pixel list, and is	still from the time
// when	I did not use the STL lists	:(
// **************************************************************

class CPixel {

	public:

		CPixel (int	nx,	int	ny)	{ x=nx;	y=ny; }

		coord x,y;
		CPixel *next;

};

#endif

