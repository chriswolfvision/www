#ifndef	_WOLF_PIXEL_H_
#define	_WOLF_PIXEL_H_

class Pixel	{

	public:

		// Constructors
		Pixel () {}
		Pixel (int nx, int ny) { x=nx; y=ny; }
		
		int	x,y;
};

#endif

