// C
#include <math.h>

// From this module
#include <CIL.h>

// From the main module
#include <CIL.h>

// From the IMAGE module
#include "FloatMatrix.h"
#include "FloatMatrix.h"
#include "Image.h"

#include "Pyramid.h"

// ***************************************************************
// Constructor
// ***************************************************************

Pyramid::Pyramid (Image *inputimage, int nLevels) 
{
	_alloc (inputimage, nLevels);	
}

// ***************************************************************
// Constructor
// The number of levels is only given indirectly by the desired
// smallest image size.
// ***************************************************************

Pyramid::Pyramid (Image *inputimage, int maxxsize, int maxysize, int minxsize, int minysize) 
{
    int levels;
	int xs,ys;

	// Calculate the number of levels
 	levels = 1;
  	xs = inputimage->xsize;
	ys = inputimage->ysize;

	/*
	cerr << "Pyramid base dimensions: " << inputimage->xsize << "x" << inputimage->ysize << " pixels.\n"
		 << "level: " << levels-1 << " xs=" << xs << " ys=" << ys << endl;
	*/
	
  	while (xs>maxxsize || ys>maxysize) 
  	{   		

   		if (xs/2>=minxsize && ys/2>=minysize) 
   		{
			xs /= 2;
  			ys /= 2;
    		++levels;
     	}
      	else
       		break;

       	// cerr << "level: " << levels-1 << " xs=" << xs << " ys=" << ys << endl;
	} 	

	_alloc (inputimage, levels);	
}

// ***************************************************************
// Constructor - help function
// ***************************************************************

void Pyramid::_alloc (Image *inputimage, int nLevels) 
{
	for (int i=0;i<nLevels; ++i)
		levelVector.push_back(NULL);
	
	reductionBuffer = NULL;
		
	changeImage (inputimage);
}

// ***************************************************************
// Return the image of the level and unlink it,
// i.e. do not destroy it when freeing the pyramid
// ***************************************************************

Image * Pyramid::getAndUnlinkLevel (int level)
{
	Image *rv = levelVector[level];
	levelVector[level] = NULL;
	return rv;
}

// ***************************************************************
// Destructor
// ***************************************************************

Pyramid::~Pyramid () {

	// The first level is a _POINTER_ to the input image, the other ones
	// are created.
	for (int i=1; i<levels(); ++i) 
	{
		if (levelVector[i]!=NULL)
			delete levelVector[i];
	}
	
	if (reductionBuffer!=NULL)
		delete reductionBuffer;
}

// ***************************************************************
// Change the baseimage of a pyramid and rebuild it.
// ***************************************************************

void Pyramid::changeImage (Image *inputimage) 
{

	if (inputimage==NULL)
		return;

	if ((levelVector[0]!=NULL) &&
		((levelVector[0]->xsize != inputimage->xsize) ||
		 (levelVector[0]->ysize != inputimage->ysize))) 
	{
		ERR_THROW ("Internal ERROR in Pyramid::changeImage()!!!!\n"
				"Changing size of a pyramid not yet implemented!!!\n");	
	}

	// Due to the inheritance properties of the pyramid,
	// this method can be called multiple times.	
	// (E.g. when a TextPyramid has been created)
	if (levelVector[0]!=inputimage) 
	{
    	
      	// The first level is a _POINTER_ to the input image, the other ones
      	// are created. Must be considered when destroying the pyramid
      	levelVector[0] = inputimage;
      	imagetype = inputimage->type;
      	
      	gaussian (0,levels()-1);	
	}
}

// ***************************************************************
// Write the pyramid to	files
// ***************************************************************

void Pyramid::write	(char *filename) 
{
	char buf[256];

	for	(int i=0; i<levels(); ++i) 
	{

		if ((*this)[i] == NULL)
			continue;

		sprintf	(buf, "%s%d.ppm", filename,	i);
		(*this)[i]->write(buf);
	}
}

/****************************************************************************
 * Build the Gaussian pyramid
 ****************************************************************************/

void Pyramid::gaussian (int begLevel, int endLevel)	
{
	int begx, begy;	
	FloatMatrix *BUF;
	
	// If not already done, allocate the reduction buffer
	if (reductionBuffer==NULL)
		reductionBuffer = new FloatMatrix (levelVector[0]->xsize, levelVector[0]->ysize);
	BUF = reductionBuffer;

	/* recentrage des niveau ? */
	if (begLevel<0)
		begLevel=0;
	if (endLevel>=levels())
		endLevel=levels()-1;
		
#ifdef CHECK_CODE
	if (levelVector[begLevel] == NULL) {
		cerr << "Internal error in Pyramid::gaussian()\n";
		exit (1);
	}		
#endif		
	
	// Travers the levels
	for(int l=begLevel+1;l<=endLevel;l++) 
	{
		float val;
		byte **CHI, **PAR;
				
		// The border
		int	topysize=levelVector[l-1]->ysize;
		int	topxsize=levelVector[l-1]->xsize;

		// Allocation of the level
		if (levelVector[l] == NULL)
			levelVector[l] = new Image (topxsize/2, topysize/2, imagetype);		
			
		// Travers the color planes
		for (int col=1; col<=baseImage().nbColorPlanes(); ++col) 
		{
		
			CHI = levelVector[l]->getPlane(col);
			PAR = levelVector[l-1]->getPlane(col);
			
     		// -------------------------------------------------------
     		// If the the line/column number is odd, then we just perform
     		// the reduction for this dimension. If the line/column number
     		// is even, then we begin the reduction with the third line/row
     		// (index 2), the first line/row (index 0) is a special case
     		// treated after the general case: the line/row must be
     		// duplicated.
     		// -------------------------------------------------------
     				
     		begx = (topxsize%2 ? 1 : 2);
     					
     		// --------------------------------------------
     		// Reduce the image - direction x
     		// --------------------------------------------			
     						
     		for(int y=0;y<topysize;++y)
     			for(int x=begx;x<topxsize;x+=2)
     				BUF->set(x,y, PAR[x-1][y] + 2*PAR[x][y] + PAR[x+1][y]);
     				
     		// duplicate first x column
     		if (begx==2) 
     		{
     			if (topxsize>1) 
     			{
     				for(int y=0;y<topysize;++y)
     					BUF->set(0, y, 3*PAR[0][y] +	PAR[1][y]);		
     			}
     			else 
     			{
     				for(int y=0;y<topysize;++y)
     					BUF->set(0, y, 4*PAR[0][y]);
     			}			
     		}
     				
     		// --------------------------------------------
     		// Reduce the image - direction y
     		// --------------------------------------------
     		
     		// If the column number is even, then begin with the
     		// first row (index 0) this time, since we already
     		// compensated the even columns above.
     		
     		begx = topxsize%2;			
     		begy = (topysize%2 ? 1 : 2);		
     		
     		for(int y=begy;y<topysize;y+=2) 
     		{
     			for(int x=begx;x<topxsize;x+=2)  
     			{
     				val = (int) rint ((BUF->get(x,y-1)+2*BUF->get(x,y)+BUF->get(x,y+1))/16.0);
                 		TRIMGRAY(val);
     				CHI[x/2][y/2] = (byte)val;
     			}
     		}
     		
     		// duplicate first y row
     		if (begy==2) 
     		{     		
     			if (topysize>1) 
     			{
         			for(int x=begx;x<topxsize;x+=2) 
         			{
         				val = (int) rint ((3*BUF->get(x,0)+BUF->get(x,1))/16.0);
                     		TRIMGRAY(val);
         				CHI[x/2][0] = (byte)val;
         			}
     			}
     			else 
     			{
     				for(int x=begx;x<topxsize;x+=2) 
     				{
         				val = (int) rint (BUF->get(x,0)/4.0);
                     		TRIMGRAY(val);
         				CHI[x/2][0] = (byte)val;
         			}
     			}
     		}
     		
     	} // plane
	}    // level
}



