/***********************************************************************
 * Some algorithms with calculated component related stuff
 * without actually performing a CCA
 *
 * Author: Christian Wolf
 ***********************************************************************/

// From the IMAGE library
#include <Image.h>

// From this library
#include "CCA.h" 
 
// *************************************************************
// Store per component: The height of the column in every pixel
// of the column
// *************************************************************

void componentColumnHeight (Image &im, Matrix<int> &dst) 
{
	Image::iterator srciter;
	int max, lastval, dstval;

	// ----------------------------------------------
	// Pass 1: downwards.
	// Count the components's height in this column
	// and pull it downwards across the column.
	// ----------------------------------------------
		
	for (int x=0; x<im.xsize; ++x) 
	{
			
		lastval = 0;
		srciter.reset (im, x, 0);				
		for (int y=0; y<im.ysize; ++y) 
		{
		
			// The pixel is set -> increase the current components height
			if (*srciter) 
			{
				++lastval;
				dst.set (x,y,lastval);
			}
			
			// The pixel is not set
			else 
			{
				dst.set (x,y,0);				
				lastval = 0;
			}
		
			++srciter;
		}		
	}	
	
	// ----------------------------------------------
	// Pass 2: upwards.
	// Carry the column's height, which is currently
	// only in the lowest pixel of the column, into
	// all pixels of the column of this component
	// ----------------------------------------------
		
	for (int x=0; x<im.xsize; ++x) 
	{
	
		max = 0;
		srciter.reset (im, x, im.ysize-1);		
		for (int y=im.ysize-1; y>=0; --y) 
		{
		
			// The pixel is set
			if (*srciter) 
			{
				dstval = dst.get (x,y);
				max = dstval > max ? dstval : max;				
			}
			
			// The pixel is not set
			else
				max = 0;
				
			dst.set (x,y,max);
		
			--srciter;
		}		
	}
}

// *************************************************************
// Store per component: The width of the line in every pixel
// of the line
// *************************************************************

void componentLineWidth (Image &im, Matrix<int> &dst) 
{
	int max, lastval, dstval;

	// ----------------------------------------------
	// Pass 1: to the right
	// Count the components's width in this line
	// and pull it to the right across the line
	// ----------------------------------------------
		
	for (int y=0; y<im.ysize; ++y) 
	{	
			
		lastval = 0;
		for (int x=0; x<im.xsize; ++x) 
		{
		
			// The pixel is set
			if (im.R[x][y]) 
			{
				++lastval;
				dst.set (x,y,lastval);
			}
			
			// The pixel is not set
			else 
			{
				dst.set (x,y,0);				
				lastval = 0;
			}		
		}		
	}	
	
	// ----------------------------------------------
	// Pass 2: to the right
	// Carry the line's width, which is currently
	// only in the right'est pixel of the line, into
	// all pixels of the line of this component
	// ----------------------------------------------
		
	for (int y=0; y<im.ysize; ++y) {	
	
		max = 0;
		for (int x=im.xsize-1; x>=0; --x) {
		
			// The pixel is set
			if (im.R[x][y]) {
				dstval = dst.get (x,y);
				max = dstval > max ? dstval : max;				
			}
			
			// The pixel is not set
			else
				max = 0;
				
			dst.set (x,y,max);	
		}		
	}
}

