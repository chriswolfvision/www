/***************************************************************************
                          FloatPyramid_Reduction.cc  -  description
                             -------------------                               	
    begin                : Thu May 30 2002
    copyright            : (C) 2002 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 *
 *	The different reduction functions for the AI Pyramids
 *  The filter kernel is not always matched to the reduction function.
 *
 ***************************************************************************/

// C
#include <math.h>

// From the IMAGE library
#include <FloatMatrix.h>
#include <Image.h>
#include <ImageFuncs.h>

// From the IMAGE PROCESSING library
#include <ImageProc.h>
#include <FilterMask.h>

// From this library
#include "FloatPyramid.h"

/****************************************************************************
 * The filtering function filtering a float Arr with a given filter Mask.
 * For reduction 2x2:
 * Calculate only every 2nd pixel of each row, and only every second row.
 * store the result in the reduction buffer
 *
 ****************************************************************************/

void FloatPyramid::filterForRed2x2 (FloatMatrix &parent, FloatMatrix &buf, FilterMask &fm) {
	int	ix, iy, bx, fx, fy, by;
	float newVal;
	int begx, begy,
		endx, endy,
		borderx, bordery;
	int fmx = fm.xsize,
		fmy = fm.ysize;
	int fmx2 = fmx/2,
		fmy2 = fmy/2;
	int xsize = parent.xsize,
		ysize = parent.ysize;
	
	// Treat only odd pixels
	begx = fmx2; if (begx%2==0) ++begx;
	begy = fmx2; if (begy%2==0) ++begy;
	endx = xsize-1; if (endx%2==0) --endx;
	endy = ysize-1; if (endy%2==0) --endy;
	borderx = xsize-fmx2;
	bordery = ysize-fmy2;
	
	// *************************************************************
	// The standard case	
	// *************************************************************	
	
	for	(int y=begy; y<(ysize-fmy2); y+=2)
	for	(int x=begx; x<(xsize-fmx2); x+=2) 
	{

		// Calculate the filter	boundaries
		bx = x-fmx2;
		by = y-fmy2;

		// convolve	
		newVal = 0;
		for	(iy=by,	fy=0; fy<fmy; ++iy, ++fy) 
		for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx) 				
			newVal += (parent.get(ix,iy)*fm.get(fx,fy));				
		
		buf.set(x,y,newVal);	
	}
	
	// *************************************************************
	// Border treatment.
	// Basically the same code as the standard
	// case, except that the borders of the image are checked, and
	// if the access exceeds the borders, the border values are
	// taken.
	// *************************************************************
	
	// horizontal borders
	// ===================
	
	for (int y=1; y<ysize; y+=2) 
	{
	
		// The left border
   	    // ---------------
   	
		for (int x=1; x<fmx2&&x<xsize; x+=2) 
		{
	    		    		
    		// Calculate the filter	boundaries
    		bx = x-fmx2;
    		by = y-fmy2;

    		// convolve	
    		newVal = 0;
    		for	(iy=by,	fy=0; fy<fmy; ++iy, ++fy) 
			{
    			int checkedy = (iy<0 ? 0 : (iy >= ysize ? (ysize-1) : iy));
    			for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx) 
				{
    				int checkedx = (ix<0 ? 0 : (ix >= xsize ? xsize-1 : ix));				
    				newVal += (parent.get(checkedx,checkedy)*fm.get(fx,fy));
    			}
    		}
    		buf.set(x,y,newVal);
    		
    	}
    	
    	// The right border
    	// ---------------

    	for (int x=endx; x>=borderx&&x>=0; x-=2) 
		{
    	
    		// Calculate the filter	boundaries
    		bx = x-fmx2;
    		by = y-fmy2;

    		// convolve	
    		newVal = 0;
    		for	(iy=by,	fy=0; fy<fmy; ++iy, ++fy) 
			{
    			int checkedy = (iy<0 ? 0 : (iy >= ysize ? (ysize-1) : iy));
    			for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx) 
				{
    				int checkedx = (ix<0 ? 0 : (ix >= xsize ? xsize-1 : ix));				
    				newVal += (parent.get(checkedx,checkedy)*fm.get(fx,fy));
    			}
    		}
    		buf.set(x,y,newVal);
    		
    	}
	}
	
	// vertical borders
	// ===================	
	
	for (int x=1; x<xsize; x+=2)
	{
		
		// The upper border
		// -----------------
		
		for (int y=1; y<fmx2&&y<xsize; y+=2) 
		{
    	
    		// Calculate the filter	boundaries
    		bx = x-fmx2;
    		by = y-fmy2;

    		// convolve	
    		newVal = 0;
    		for	(iy=by,	fy=0; fy<fmy; ++iy, ++fy) 
			{
    			int checkedy = (iy<0 ? 0 : (iy >= ysize ? (ysize-1) : iy));
    			for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx) 
				{
    				int checkedx = (ix<0 ? 0 : (ix >= xsize ? xsize-1 : ix));
    				newVal += (parent.get(checkedx,checkedy)*fm.get(fx,fy));
    			}
    		}
    		buf.set(x,y,newVal);
    	}
    	
    	// The lower border
		// -----------------

		for (int y=endy; y>=bordery&&y>=0; y-=2) 
		{    	
    	
    		// Calculate the filter	boundaries
    		bx = x-fmx2;
    		by = y-fmy2;

    		// convolve	
    		newVal = 0;
    		for	(iy=by,	fy=0; fy<fmy; ++iy, ++fy) 
			{
    			int checkedy = (iy<0 ? 0 : (iy >= ysize ? (ysize-1) : iy));
    			for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx) 
				{
    				int checkedx = (ix<0 ? 0 : (ix >= xsize ? xsize-1 : ix));		
    				newVal += (parent.get(checkedx,checkedy)*fm.get(fx,fy));
    			}
    		}
    		buf.set(x,y,newVal);
    	}
    }		
}

/****************************************************************************
 * Build the Gaussian pyramid
 * The general version accepting a filter mask as argument
 * reductions is always22
 * This version is slower than the other, highly specialized versions,
 * (it does not use seperable filters), but it must be used for the
 * diagonal smoothing filters
 ****************************************************************************/

void FloatPyramid::gaussianFMred22 (vector <FloatMatrix *> &vec, int begLevel, int endLevel,
	FilterMask &fm)	{	
	
	FloatMatrix *BUF;
	
	// If not already done, allocate the reduction buffer
	if (reductionBuffer==NULL)
		reductionBuffer = new FloatMatrix (vec[0]->xsize, vec[0]->ysize);
	BUF = reductionBuffer;

	/* recentrage des niveau ? */
	if (begLevel<0)
		begLevel=0;
	if (endLevel>=noLevels)
		endLevel=noLevels-1;
		
#ifdef CHECK_CODE
	if (vec[begLevel] == NULL) {
		cerr << "Internal error in FloatPyramid::gaussianFMred22()\n";
		exit (1);
	}		
#endif		
	
	// Travers the noLevels
	for(int l=begLevel+1;l<=endLevel;l++) {
		FloatMatrix *CHI, *PAR;
				
		// The border
		int	topysize=vec[l-1]->ysize;
		int	topxsize=vec[l-1]->xsize;

		// Allocation of the level
		if (vec[l]==NULL)
			vec[l] = new FloatMatrix (topxsize/2, topysize/2);
					
		CHI = vec[l];
		PAR = vec[l-1];
			
 		// ------------------------------------------------------- 		
 		// Filter the image, but calculate only every 2nd pixel
 		// of each row, and only every second row.
 		// --------------------------------------------			

 		filterForRed2x2 (*PAR, *BUF,fm);
 		
 		// ------------------------------------------------------- 		
 		// Reduce the image
 		// --------------------------------------------			
 		     	
 		for(int y=1;y<=topysize-1;y+=2)
 		for(int x=1;x<=topxsize-1;x+=2) {
 			CHI->set(x/2,y/2, BUF->get(x,y));     		
 		}
	}
}	

/****************************************************************************
 * Build the Gaussian pyramid
 * The an-isotropic 2x2 reduction version using a separable
 * 7x3 Gaussian filter [ 1 6 15 20 15 6 1 ] * [ 1 2 1 ]T
 ****************************************************************************/

void FloatPyramid::gaussian73red22 (vector <FloatMatrix *> &vec, int begLevel, int endLevel)	{	
	FloatMatrix *BUF;
	
	// If not already done, allocate the reduction buffer
	if (reductionBuffer==NULL)
		reductionBuffer = new FloatMatrix (vec[0]->xsize, vec[0]->ysize);
	BUF = reductionBuffer;

	/* recentrage des niveau ? */
	if (begLevel<0)
		begLevel=0;
	if (endLevel>=noLevels)
		endLevel=noLevels-1;
		
#ifdef CHECK_CODE
	if (vec[begLevel] == NULL) {
		cerr << "Internal error in FloatPyramid::gaussian32()\n";
		exit (1);
	}		
#endif		
	
	// Travers the levels
	for(int l=begLevel+1;l<=endLevel;l++) {
		float val;
		FloatMatrix *CHI, *PAR;
		int begx, endx, endy;
				
		// The border
		int	topysize=vec[l-1]->ysize;
		int	topxsize=vec[l-1]->xsize;

		// Allocation of the level
		if (vec[l]==NULL)
			vec[l] = new FloatMatrix (topxsize/2, topysize/2);
					
		CHI = vec[l];
		PAR = vec[l-1];
			
 		// ------------------------------------------------------- 		
 		// See research notebook 30.5.2002, page 125
		//
 		// Reduce the image - direction x
 		// (standard processing)
 		// --------------------------------------------			
 		
 		begx = 3;
 		endx = topxsize-5;
     						
 		for(int y=0;y<topysize;++y)
 			for(int x=begx;x<=endx;x+=2)
 				BUF->set(x,y, (
 					PAR->get(x-3,y) + 6.0*PAR->get(x-2,y) + 15.0*PAR->get(x-1,y) +
 					20.0*PAR->get(x,y) + 15.0*PAR->get(x+1,y) + 6.0*PAR->get(x+2,y)
 					+PAR->get(x+3,y)));
		
 		if (topxsize<5) {
 			cerr << "Internal error in FloatPyramid::gaussian73red32():\n"
 				"special case not yet implemented!\n";
 			exit (1);
 		}
 		else {
 		
     		// duplicate the first two x columns (always necessary)
 			for(int y=0;y<topysize;++y)
 				BUF->set(1, y, (22.0*PAR->get(0,y) + 20.0*PAR->get(1,y) +
 					15.0*PAR->get(2,y) + 6.0*PAR->get(3,y) + PAR->get(4,y)));	 			
	 					
 			// adjust the left border after the special treatment
 			begx = 1;
 			
 			// duplicate:
 			// position xsize-1 the last 3 columns
 			// position xsize-3 the last column
     		if (topxsize%2==0) {
     			int x=topxsize-1;
 				for(int y=0;y<topysize;++y)
 					BUF->set(x, y, (PAR->get(x-3,y) + 6.0*PAR->get(x-2,y) + 15.0*PAR->get(x-1,y) +
 						42.0*PAR->get(x,y)));
 						
 				x=topxsize-3;
 				for(int y=0;y<topysize;++y)
 					BUF->set(x, y, (PAR->get(x-3,y) + 6.0*PAR->get(x-2,y) + 15.0*PAR->get(x-1,y) +
 						20.0*PAR->get(x,y)+15.0*PAR->get(x+1,y)+ 7.0*PAR->get(x+2,y)));
	 					
	 			// adjust the right border after the special treatment
	 			endx = topxsize-1;
 			}
 			
 			// duplicate the last 2 columns
 			else {
 				int x=topxsize-2;
 				for(int y=0;y<topysize;++y)
 					BUF->set(x, y, (PAR->get(x-3,y) + 6.0*PAR->get(x-2,y) + 15.0*PAR->get(x-1,y) +
 						20.0*PAR->get(x,y)+22.0*PAR->get(x+1,y)));
	 					
	 			// adjust the right border after the special treatment
	 			endx = x;
 			}
 		}
 		
 		// --------------------------------------------
 		// Reduce the image - direction y
 		// --------------------------------------------	
 		
 		endy = (topysize%2==1 ? topysize-2 : topysize-3);
     		
 		// reduce, standard case
 		for(int y=1;y<=endy;y+=2) {
 			for(int x=begx;x<=endx;x+=2)  {
 				val = (BUF->get(x,y-1)+2.0*BUF->get(x,y)+BUF->get(x,y+1))/256.0;             	
 				CHI->set(x/2,y/2, val);
 			}
 		}
     		
 		// duplicate last y row
 		if (topysize%2==0) {     		 		 		     		
 			if (topysize>1) {
 				int y = topysize-1;
     			for(int x=begx;x<=endx;x+=2) {
     				val = (3.0*BUF->get(x,y)+BUF->get(x,y-1))/256.0;
     				CHI->set(x/2,y/2, val);
     			}
 			}
 			else {
 				for(int x=begx;x<=endx;x+=2) {
     				val = BUF->get(x,0)/64.0;
     				CHI->set(x/2, 0, val);
     			}
 			}
 		}     	
	}
}	

/****************************************************************************
 * Build the Gaussian pyramid
 * The an-isotropic 2x2 reduction version using a separable
 * 3x7 Gaussian filter [ 1 2 1 ] * [ 1 6 15 20 15 6  1 ]T
 * Similar to the case above (just rotated).
 ****************************************************************************/

void FloatPyramid::gaussian37red22 (vector <FloatMatrix *> &vec, int begLevel, int endLevel){
	FloatMatrix *BUF;
	
	// If not already done, allocate the reduction buffer
	if (reductionBuffer==NULL)
		reductionBuffer = new FloatMatrix (vec[0]->xsize, vec[0]->ysize);
	BUF = reductionBuffer;

	// recentrage des niveau
	if (begLevel<0)
		begLevel=0;
	if (endLevel>=noLevels)
		endLevel=noLevels-1;
		
#ifdef CHECK_CODE
	if (vec[begLevel] == NULL) {
		cerr << "Internal error in FloatPyramid::gaussian23()\n";
		exit (1);
	}		
#endif		
	
	// Travers the levels
	for(int l=begLevel+1;l<=endLevel;l++) {
		float val;
		FloatMatrix *CHI, *PAR;
		int begx=0, begy=0, endx, endy;
						
		// The border
		int	topysize=vec[l-1]->ysize;
		int	topxsize=vec[l-1]->xsize;

		// Allocation of the level
		if (vec[l]==NULL) {
			vec[l] = new FloatMatrix (topxsize/2, topysize/2);			
			cerr.flush();
		}
					
		CHI = vec[l];
		PAR = vec[l-1];
			
		// -------------------------------------------------------
 		// See research notebook 29.5.2002, p 129
 		// -------------------------------------------------------

 		endx = (topxsize%2==1 ? topxsize-2 : topxsize-3);
     					
 		// --------------------------------------------
 		// Reduce the image - direction x
 		// --------------------------------------------			
     						
 		for(int y=0;y<topysize;++y)
 			for(int x=1;x<=endx;x+=2)
 				BUF->set(x,y, (PAR->get(x-1,y) + 2.0*PAR->get(x,y) + PAR->get(x+1,y)));
     				
 		// duplicate last x column
 		if (topxsize%2==0) {
 			if (topxsize>1) {
 				int x = topxsize-1;
 				for(int y=0;y<topysize;++y)
 					BUF->set(x, y, PAR->get(x-1,y) + 3.0*PAR->get(x,y));		
 					
 				// adjust the right border after the special treatment
 				endx = x;
 			}
 			else {
 				for(int y=0;y<topysize;++y)
 					BUF->set(0, y, 4.0*PAR->get(0,y));
 				endx = 0;
 			}
 		}

		// --------------------------------------------
 		// Reduce the image - direction y
 		// Gaussian filter [ 1 6 15 20 15 6 1 ]T
 		// --------------------------------------------	
 		
 		begy = 3;
 		endy = topysize-5;
     		
 		// reduce, standard case
 		for(int y=begy;y<=endy;y+=2) {
 			for(int x=begx;x<=endx;x+=2)  {
 				int yh=y/2;
 				val = (BUF->get(x,y-3)+6.0*BUF->get(x,y-2)+15.0*BUF->get(x,y-1)+
 					 20.0*BUF->get(x,y)+15.0*BUF->get(x,y+1)+6.0*BUF->get(x,y+2)+
 					 BUF->get(x,y+3))/256.0;
 				CHI->set(x/2,yh, val);
 			}
 		}
 		
 		if (topysize<5) {
 			cerr << "Internal error in FloatPyramid::gaussian37red22():\n"
 				"special case not yet implemented!\n";
 			exit (1);
 		}
 		else {
 		
	 		// duplicate the first two rows (always necessary)
     		for(int x=begx;x<=endx;x+=2) {
     			val = (22.0*BUF->get(x,0)+20.0*BUF->get(x,1)+15.0*BUF->get(x,2)+
     				6.0*BUF->get(x,3)+BUF->get(x,4))/256.0;
     			CHI->set(x/2,0, val);
     		}
     		
     		// duplicate the last rows
     		if (topysize%2==0) {
 				int y = topysize-1;
 				int yh = y/2;
     			for(int x=begx;x<=endx;x+=2) {
     				val = (BUF->get(x,y-3)+6.0*BUF->get(x,y-2)+15.0*BUF->get(x,y-1)+
     					42.0*BUF->get(x,y))/256.0;
     				CHI->set(x/2,yh, val);
     			}
     			
     			y = topysize-3;
 				yh = y/2;
     			for(int x=begx;x<=endx;x+=2) {
     				val = (BUF->get(x,y-3)+6.0*BUF->get(x,y-2)+15.0*BUF->get(x,y-1)+
     					20.0*BUF->get(x,y)+15.0*BUF->get(x,y+1)+7.0*BUF->get(x,y+2))/256.0;
     				CHI->set(x/2,yh, val);
     			}
     		}
     		
     		else {
     			int y = topysize-2;
 				int yh = y/2;
     			for(int x=begx;x<=endx;x+=2) {
     				val = (BUF->get(x,y-3)+6.0*BUF->get(x,y-2)+15.0*BUF->get(x,y-1)+
     					20.0*BUF->get(x,y)+22.0*BUF->get(x,y+1))/256.0;
     				CHI->set(x/2,yh, val);
     			}
     		}
     	}
	}
}


/****************************************************************************
 * Build the Gaussian pyramid
 * The isotropic 2x2 reduction version using a separable
 * Gaussian filter [ 1 2 1 ] * [ 1 2 1 ]T
 ****************************************************************************/

void FloatPyramid::gaussian33red22 (vector <FloatMatrix *> &vec, int begLevel, int endLevel)	{
	int begx, endx, endy;	
	FloatMatrix *BUF;	
	
	// If not already done, allocate the reduction buffer
	if (reductionBuffer==NULL)
		reductionBuffer = new FloatMatrix (vec[0]->xsize, vec[0]->ysize);
	BUF = reductionBuffer;

	/* recentrage des niveau ? */
	if (begLevel<0)
		begLevel=0;
	if (endLevel>=noLevels)
		endLevel=noLevels-1;
		
#ifdef CHECK_CODE
	if (vec[begLevel] == NULL) {
		cerr << "Internal error in FloatPyramid::gaussian22()\n";
		exit (1);
	}		
#endif		
	
	// Travers the levels
	for(int l=begLevel+1;l<=endLevel;l++) {
		float val;
		FloatMatrix *CHI, *PAR;
				
		// The border
		int	topysize=vec[l-1]->ysize;
		int	topxsize=vec[l-1]->xsize;

		// Allocation of the level
		if (vec[l] == NULL)
			vec[l] = new FloatMatrix (topxsize/2, topysize/2);
			
		cerr << "level: " << vec[l]->xsize << ", " << vec[l]->ysize << endl;
					
		CHI = vec[l];
		PAR = vec[l-1];
			
 		// -------------------------------------------------------
 		// See research notebook 29.5.2002, p 129
 		// -------------------------------------------------------

 		begx = 1;
 		endx = (topxsize%2==1 ? topxsize-2 : topxsize-3);
     					
 		// --------------------------------------------
 		// Reduce the image - direction x
 		// --------------------------------------------			
     						
 		for(int y=0;y<topysize;++y)
 			for(int x=begx;x<=endx;x+=2)
 				BUF->set(x,y, (PAR->get(x-1,y) + 2.0*PAR->get(x,y) + PAR->get(x+1,y)));
     				
 		// duplicate last x column
 		if (topxsize%2==0) {
 			if (topxsize>1) {
 				int x = topxsize-1;
 				for(int y=0;y<topysize;++y)
 					BUF->set(x, y, 3.0*PAR->get(x,y) + PAR->get(x-1,y));		
 					
 				// adjust the right border after the special treatment
	 			endx = x;
 			}
 			else {
 				for(int y=0;y<topysize;++y)
 					BUF->set(0, y, 4.0*PAR->get(0,y));
 				endx = 0;
 			}
 		}
     				
 		// --------------------------------------------
 		// Reduce the image - direction y
 		// --------------------------------------------
     		
 		endy = (topysize%2==1 ? topysize-2 : topysize-3);
     		
 		// reduce, standard case
 		for(int y=1;y<=endy;y+=2) {
 			for(int x=1;x<=endx;x+=2)  {
 				val = (BUF->get(x,y-1)+2.0*BUF->get(x,y)+BUF->get(x,y+1))/16.0;             	
 				CHI->set(x/2,y/2, val);
 			}
 		}
     		
 		// duplicate last y row
 		if (topysize%2==0) {
     		
 			if (topysize>1) {
 				int y = topysize-1;
     			for(int x=begx;x<=endx;x+=2) {
     				val = (3.0*BUF->get(x,y)+BUF->get(x,y-1))/16.0;
     				CHI->set(x/2,y/2, val);
     			}
 			}
 			else {
 				for(int x=begx;x<topxsize;x+=2) {
     				val = BUF->get(x,0)/4.0;
     				CHI->set(x/2, 0, val);
     			}
 			}
 		}     	
	}
}

/****************************************************************************
 * Build the Gaussian pyramid
 * The an-isotropic 2x2 reduction version using a separable
 * Gaussian filter [ 1 4 6 4 1 ] * [ 1 2 1 ]T
 ****************************************************************************/

void FloatPyramid::gaussian53red22 (vector <FloatMatrix *> &vec, int begLevel, int endLevel)	{	
	FloatMatrix *BUF;
	
	// If not already done, allocate the reduction buffer
	if (reductionBuffer==NULL)
		reductionBuffer = new FloatMatrix (vec[0]->xsize, vec[0]->ysize);
	BUF = reductionBuffer;

	/* recentrage des niveau ? */
	if (begLevel<0)
		begLevel=0;
	if (endLevel>=noLevels)
		endLevel=noLevels-1;
		
#ifdef CHECK_CODE
	if (vec[begLevel] == NULL) {
		cerr << "Internal error in FloatPyramid::gaussian32()\n";
		exit (1);
	}		
#endif		
	
	// Travers the levels
	for(int l=begLevel+1;l<=endLevel;l++) {
		float val;
		FloatMatrix *CHI, *PAR;
		int begx, endx, endy;
				
		// The border
		int	topysize=vec[l-1]->ysize;
		int	topxsize=vec[l-1]->xsize;

		// Allocation of the level
		if (vec[l]==NULL)
			vec[l] = new FloatMatrix (topxsize/2, topysize/2);
					
		CHI = vec[l];
		PAR = vec[l-1];
			
 		// ------------------------------------------------------- 		
 		// See research notebook 29.5.2002, page 129
 		// The first line/column always needs to duplicated
 		// The last line/column as well, and if xsize%2==0 then
 		// also the second to last line/column
 		// -------------------------------------------------------
 		
 		begx = 3;
 		endx = (topxsize%2==1 ? topxsize-4 : topxsize-3);
     					
 		// --------------------------------------------
 		// Reduce the image - direction x
 		// (standard processing)
 		// --------------------------------------------			
     						
 		for(int y=0;y<topysize;++y)
 			for(int x=begx;x<=endx;x+=2)
 				BUF->set(x,y, (PAR->get(x-2,y) + 4.0*PAR->get(x-1,y) +
 					6.0*PAR->get(x,y) + 4.0*PAR->get(x+1,y) + PAR->get(x+2,y)));
     				
 		if (topxsize<4) {
 			cerr << "Internal error in FloatPyramid::gaussian53red32():\n"
 				"special case not yet implemented!\n";
 			exit (1);
 		}
 		else {
 		
     		// duplicate the first x column (always necessary)
 			for(int y=0;y<topysize;++y)
 				BUF->set(1, y, (5.0*PAR->get(0,y) +
 					6.0*PAR->get(1,y) + 4.0*PAR->get(2,y) + PAR->get(3,y)));	 			
	 					
 			// adjust the left border after the special treatment
 			begx = 1;
 			
 			// duplicate the last column
     		if (topxsize%2==1) {
     			int x=topxsize-2;
 				for(int y=0;y<topysize;++y)
 					BUF->set(x, y, (PAR->get(x-2,y) + 4.0*PAR->get(x-1,y) +
 						6.0*PAR->get(x,y) + 5.0*PAR->get(x+1,y)));
	 					
	 			// adjust the right border after the special treatment
	 			endx = x;
 			}
 			
 			// duplicate the last AND the second to last column			
 			else {
 				int x=topxsize-1;
 				for(int y=0;y<topysize;++y)
 					BUF->set(x, y, (PAR->get(x-2,y) + 4.0*PAR->get(x-1,y) +
 						11.0*PAR->get(x,y)));
	 					
	 			// adjust the right border after the special treatment
	 			endx = x;
 			}
 		}
     				
 		// --------------------------------------------
 		// Reduce the image - direction y
 		// --------------------------------------------	
 		
 		endy = (topysize%2==1 ? topysize-2 : topysize-3);
     		
 		// reduce, standard case
 		for(int y=1;y<=endy;y+=2) {
 			for(int x=begx;x<=endx;x+=2)  {
 				val = (BUF->get(x,y-1)+2.0*BUF->get(x,y)+BUF->get(x,y+1))/64.0;             	
 				CHI->set(x/2,y/2, val);
 			}
 		}
     		
 		// duplicate last y row
 		if (topysize%2==0) {     		 		 		     		
 			if (topysize>1) {
 				int y = topysize-1;
     			for(int x=begx;x<=endx;x+=2) {
     				val = (3.0*BUF->get(x,y)+BUF->get(x,y-1))/64.0;
     				CHI->set(x/2,y/2, val);
     			}
 			}
 			else {
 				for(int x=begx;x<=endx;x+=2) {
     				val = BUF->get(x,0)/16.0;
     				CHI->set(x/2, 0, val);
     			}
 			}
 		}     	
	}
}

/****************************************************************************
 * Build the Gaussian pyramid
 * The an-isotropic 2x2 reduction version using a separable
 * Gaussian filter [ 1 2 1 ] * [ 1 4 6 4 1 ]T
 * Similar to the case above (just rotated).
 ****************************************************************************/

void FloatPyramid::gaussian35red22 (vector <FloatMatrix *> &vec, int begLevel, int endLevel){
	FloatMatrix *BUF;
	
	// If not already done, allocate the reduction buffer
	if (reductionBuffer==NULL)
		reductionBuffer = new FloatMatrix (vec[0]->xsize, vec[0]->ysize);
	BUF = reductionBuffer;

	// recentrage des niveau
	if (begLevel<0)
		begLevel=0;
	if (endLevel>=noLevels)
		endLevel=noLevels-1;
		
#ifdef CHECK_CODE
	if (vec[begLevel] == NULL) {
		cerr << "Internal error in FloatPyramid::gaussian23()\n";
		exit (1);
	}		
#endif		
	
	// Travers the levels
	for(int l=begLevel+1;l<=endLevel;l++) {
		float val;
		FloatMatrix *CHI, *PAR;
		int begx=0, begy=0, endx, endy;
						
		// The border
		int	topysize=vec[l-1]->ysize;
		int	topxsize=vec[l-1]->xsize;

		// Allocation of the level
		if (vec[l]==NULL) {
			vec[l] = new FloatMatrix (topxsize/2, topysize/2);			
			cerr.flush();
		}
					
		CHI = vec[l];
		PAR = vec[l-1];
			
		// -------------------------------------------------------
 		// See research notebook 29.5.2002, p 129
 		// -------------------------------------------------------

 		endx = (topxsize%2==1 ? topxsize-2 : topxsize-3);
     					
 		// --------------------------------------------
 		// Reduce the image - direction x
 		// --------------------------------------------			
     						
 		for(int y=0;y<topysize;++y)
 			for(int x=1;x<=endx;x+=2)
 				BUF->set(x,y, (PAR->get(x-1,y) + 2.0*PAR->get(x,y) + PAR->get(x+1,y)));
     				
 		// duplicate last x column
 		if (topxsize%2==0) {
 			if (topxsize>1) {
 				int x = topxsize-1;
 				for(int y=0;y<topysize;++y)
 					BUF->set(x, y, PAR->get(x-1,y) + 3.0*PAR->get(x,y));		
 					
 				// adjust the right border after the special treatment
 				endx = x;
 			}
 			else {
 				for(int y=0;y<topysize;++y)
 					BUF->set(0, y, 4.0*PAR->get(0,y));
 				endx = 0;
 			}
 		}

		// --------------------------------------------
 		// Reduce the image - direction y
 		// Gaussian filter [ 1 4 6 4 1 ]T
 		// --------------------------------------------	
 		
 		begy = 3;
 		endy = (topysize%2==1 ? topysize-4 : topysize-3);
     		
 		// reduce, standard case
 		for(int y=begy;y<=endy;y+=2) {
 			for(int x=begx;x<=endx;x+=2)  {
 				int yh=y/2;
 				val = (BUF->get(x,y-2)+4.0*BUF->get(x,y-1)+6.0*BUF->get(x,y)+
 					 4.0*BUF->get(x,y+1)+BUF->get(x,y+2))/64.0;             	
 				CHI->set(x/2,yh, val);
 			}
 		}
 		
 		if (topysize<4) {
 			cerr << "Internal error in FloatPyramid::gaussian35red22():\n"
 				"special case not yet implemented!\n";
 			exit (1);
 		}
 		else {
 		
	 		// duplicate the first row (always necessary)
     		for(int x=begx;x<=endx;x+=2) {
     			val = (5.0*BUF->get(x,0)+6.0*BUF->get(x,1)+4.0*BUF->get(x,2)+BUF->get(x,3))/64.0;
     			CHI->set(x/2,0, val);
     		}
     		
     		// duplicate last y row
     		if (topysize%2==1) {     		 		 		     		
 				int y = topysize-2;
 				int yh = y/2;
     			for(int x=begx;x<=endx;x+=2) {
     				val = (BUF->get(x,y-2)+       
					   4.0*BUF->get(x,y-1)+
					   6.0*BUF->get(x,y)+
					   5.0*BUF->get(x,y+1))/64.0;
     				CHI->set(x/2,yh, val);
     			}
     		}
     		
     		// duplicate the last and the second to last y row
     		else {
     			int y = topysize-1;
 				int yh = y/2;
     			for(int x=begx;x<=endx;x+=2) {
     				val = (BUF->get(x,y-2)+4.0*BUF->get(x,y-1)+11.0*BUF->get(x,y))/4.0;
     				CHI->set(x/2,yh, val);
     			}
     		}
     	}
	}
}

/****************************************************************************
 * Build a maximum pyramid
 ****************************************************************************/

#define MAX3(x,y,z)	((x>y)?(x>z?x:z):(y>z?y:z))
#define MAX2(x,y)	((x>y)?x:y)
#define MAX1(x)		(x)

void FloatPyramid::maxRed22 (vector <FloatMatrix *> &vec, int begLevel, int endLevel)	
{
	int begx, endx, endy;	
	FloatMatrix *BUF;	
	
	// If not already done, allocate the reduction buffer
	if (reductionBuffer==NULL)
		reductionBuffer = new FloatMatrix (vec[0]->xsize, vec[0]->ysize);
	BUF = reductionBuffer;

	/* recentrage des niveau ? */
	if (begLevel<0)
		begLevel=0;
	if (endLevel>=noLevels)
		endLevel=noLevels-1;
		
#ifdef CHECK_CODE
	if (vec[begLevel] == NULL) {
		cerr << "Internal error in FloatPyramid::gaussian22()\n";
		exit (1);
	}		
#endif		
	
	// Travers the levels
	for(int l=begLevel+1;l<=endLevel;l++) {
		float val;
		FloatMatrix *CHI, *PAR;
				
		// The border
		int	topysize=vec[l-1]->ysize;
		int	topxsize=vec[l-1]->xsize;

		// Allocation of the level
		if (vec[l] == NULL)
			vec[l] = new FloatMatrix (topxsize/2, topysize/2);
					
		CHI = vec[l];
		PAR = vec[l-1];
			
 		// -------------------------------------------------------
 		// See research notebook 29.5.2002, p 129
 		// -------------------------------------------------------

 		begx = 1;
 		endx = (topxsize%2==1 ? topxsize-2 : topxsize-3);
     					
 		// --------------------------------------------
 		// Reduce the image - direction x
 		// --------------------------------------------			
     						
 		for(int y=0;y<topysize;++y)
 			for(int x=begx;x<=endx;x+=2)
 				BUF->set(x,y, (MAX3(PAR->get(x-1,y), PAR->get(x,y), PAR->get(x+1,y))));
     				
 		// duplicate last x column
 		if (topxsize%2==0) {
 			if (topxsize>1) {
 				int x = topxsize-1;
 				for(int y=0;y<topysize;++y)
 					BUF->set(x, y, MAX2(PAR->get(x,y), PAR->get(x-1,y)));
 					
 				// adjust the right border after the special treatment
	 			endx = x;
 			}
 			else {
 				for(int y=0;y<topysize;++y)
 					BUF->set(0, y, MAX1(PAR->get(0,y)));
 				endx = 0;
 			}
 		}
     				
 		// --------------------------------------------
 		// Reduce the image - direction y
 		// --------------------------------------------
     		
 		endy = (topysize%2==1 ? topysize-2 : topysize-3);
     		
 		// reduce, standard case
 		for(int y=1;y<=endy;y+=2) {
 			for(int x=1;x<=endx;x+=2)  {
 				val = MAX3(BUF->get(x,y-1),BUF->get(x,y),BUF->get(x,y+1));
 				CHI->set(x/2,y/2, val);
 			}
 		}
     		
 		// duplicate last y row
 		if (topysize%2==0) {
     		
 			if (topysize>1) {
 				int y = topysize-1;
     			for(int x=begx;x<=endx;x+=2) {
     				val = MAX2(BUF->get(x,y),BUF->get(x,y-1));
     				CHI->set(x/2,y/2, val);
     			}
 			}
 			else {
 				for(int x=begx;x<topxsize;x+=2) {
     				val = MAX1(BUF->get(x,0));
     				CHI->set(x/2, 0, val);
     			}
 			}
 		}     	
	}
}

