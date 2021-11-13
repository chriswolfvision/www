/***************************************************************************
                          FloatPyramid.cc  -  description
                             -------------------
 ***************************************************************************/

// C
#include <math.h>

// From the main module
#include <CIL.h>

// From the IMAGE module
#include <FloatMatrix.h>
#include <Image.h>
#include <ImageFuncs.h>

// From the IMAGE PROCESSING module
#include <ImageProc.h>

// From this library
#include "FloatPyramid.h"
#include "FloatPyramid_Parents.h"

// ***************************************************************
// Constructor
// ***************************************************************

FloatPyramid::FloatPyramid (FloatMatrix *I,	int nl, PyramidType t) 
{
	
	noLevels = nl;
	type=t;
	reductionBuffer = NULL;	
	
	for (int i=0;i<noLevels; ++i)
		images.push_back(NULL);		
	
	changeImage (I);
	
	// construct the convolution kernel for the
	// L-diagonal gaussian filter
	if (type==PYRT_FIL77LD_RED22)
	{
		mask << 0.0 << 1.0 <<  1.0 <<  0.0 <<  0.0 << 0.0 << 0.0;
		mask << 1.0 << 3.0 <<  7.0 <<  4.0 <<  1.0 << 0.0 << 0.0;
		mask << 1.0 << 7.0 << 22.0 << 23.0 <<  8.0 << 1.0 << 0.0;
		mask << 0.0 << 4.0 << 23.0 << 42.0 << 23.0 << 4.0 << 0.0;
		mask << 0.0 << 1.0 <<  8.0 << 23.0 << 22.0 << 7.0 << 1.0;
		mask << 0.0 << 0.0 <<  1.0 <<  4.0 <<  7.0 << 3.0 << 1.0;
		mask << 0.0 << 0.0 <<  0.0 <<  0.0 <<  1.0 << 1.0 << 0.0;
		mask.xsize = 7;
		mask.ysize = 7;
		mask.normalize();
	}
	
	// construct the convolution kernel for the
	// R-diagonal gaussian filter
	if (type==PYRT_FIL77RD_RED22)
	{
		mask << 0.0 << 0.0 <<  0.0 <<  0.0 <<  1.0 << 1.0 << 0.0;
		mask << 0.0 << 0.0 <<  1.0 <<  4.0 <<  7.0 << 3.0 << 1.0;
		mask << 0.0 << 1.0 <<  8.0 << 23.0 << 22.0 << 7.0 << 1.0;
		mask << 0.0 << 4.0 << 23.0 << 42.0 << 23.0 << 4.0 << 0.0;
		mask << 1.0 << 7.0 << 22.0 << 23.0 <<  8.0 << 1.0 << 0.0;
		mask << 1.0 << 3.0 <<  7.0 <<  4.0 <<  1.0 << 0.0 << 0.0;
		mask << 0.0 << 1.0 <<  1.0 <<  0.0 <<  0.0 << 0.0 << 0.0;			
		mask.xsize = 7;
		mask.ysize = 7;
		mask.normalize();
	}
}

// ***************************************************************
// Destructor
// ***************************************************************

FloatPyramid::~FloatPyramid () 
{

	// The first level is a _POINTER_ to the input image, the other ones
	// are created.
	for (int i=1; i<noLevels; ++i)
		if (images[i]!=NULL) 
			delete images[i];		
	
	if (reductionBuffer!=NULL)
		delete reductionBuffer;
}

// ***************************************************************
// Change the baseimage of a pyramid and rebuild it.
// ***************************************************************

void FloatPyramid::changeImage (FloatMatrix *X) 
{
	if (X==NULL)
		return;

	if ((images[0]!=NULL) && 
		((images[0]->xsize != X->xsize) ||
		 (images[0]->ysize != X->ysize))) 
	{
		ERR_THROW ("Internal ERROR in FloatPyramid::changeImage()!!!!\n"
				"Changing size of a pyramid not yet implemented!!!\n");	
	}

	// Due to the inheritance properties of the pyramid,
	// this method can be called multiple times.	
	// (E.g. when a TextPyramid has been created)
	if (images[0]!=X) 
	{
    	
      	// The first level is a _POINTER_ to the input image, the other ones
      	// are created. Must be considered when destroying the pyramid
      	images[0] = X;     
	}
}

/****************************************************************************
 * Build the pyramid
 ****************************************************************************/

void FloatPyramid::build () 
{				
	switch (type) 
	{
		case PYRT_FIL33_RED22:
			gaussian33red22 (images, 0, noLevels-1);
			break;
		case PYRT_FIL35_RED22:
			gaussian35red22 (images, 0, noLevels-1);
			break;
		case PYRT_FIL53_RED22:
			gaussian53red22 (images, 0, noLevels-1);
			break;
		case PYRT_FIL73_RED22:
			gaussian73red22 (images, 0, noLevels-1);
			break;
		case PYRT_FIL37_RED22:
			gaussian37red22 (images, 0, noLevels-1);
			break;				
			
		// the mask is in each case different (see constuctor)
		case PYRT_FIL77LD_RED22:		
		case PYRT_FIL77RD_RED22:
			gaussianFMred22 (images, 0,noLevels-1, mask);
			break;
			
		case PYRT_MAX_RED22:
			maxRed22 (images, 0,noLevels-1);
			break;
			
		case PYRT_FIL77_RED22:
		default:
			ERR_THROW ("FloatPyramid::build(): unknown type or not yet implemented!\n");
	}				  
}


/****************************************************************************
 * Descend the pyramid
 ****************************************************************************/

void FloatPyramid::collapse () 
{
	
	// Descend the pyramid and rebuild the noLevels
	for (int l=noLevels-2; l>=0; --l) 
	{
		int xs=images[l]->xsize;
		int ys=images[l]->ysize;		
		int par_level=l+1;		
		
		FloatMatrix *Xp = images[par_level];
		FloatMatrix *Xl = images[l];
					
		for (int x=0; x<xs; ++x)
		for (int y=0; y<ys; ++y) 
		{
			float v;
			
			switch (type) 
			{
				case PYRT_FIL33_RED22:
					v = accessParents33red22 (Xp, x, y);
					break;
				case PYRT_FIL35_RED22:
					v = accessParents35red22 (Xp, x, y);
					break;
				case PYRT_FIL53_RED22:
					v = accessParents53red22 (Xp, x, y);
					break;
				case PYRT_FIL73_RED22:
					v = accessParents73red22 (Xp, x, y);
					break;
				case PYRT_FIL37_RED22:
					v = accessParents37red22 (Xp, x, y);
					break;				
				case PYRT_FIL77RD_RED22:
					v = accessParents77RDred22 (Xp, x, y);
					break;
				case PYRT_FIL77LD_RED22:
					v = accessParents77LDred22 (Xp, x, y);    		
					break;
					
				case PYRT_MAX_RED22:
					ERR_THROW ("Collapse not yet supported for that pyramid type!");
					
				case PYRT_FIL77_RED22:
				default:
					ERR_THROW ("FloatPyramid::collapse(): unknown type or not yet implemented!\n");
			}				  
				  
    		// Write it back
    		Xl->set(x,y, v);
    	}
	}		
}




