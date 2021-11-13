/*********************************************************************
 * Statistical routines
 *
 * Christian Wolf, chriswolf@gmx.at
 *********************************************************************/

// C
#include <errno.h>
#include <math.h>

// C++
#include <map>

// From the main module
#include <CIL.h>

// From the MATHEMATICS module
#include <ConfusionMatrix.h>

// From the IMAGE module
#include <FloatMatrix.h>

// From our own module
#include "ImageProc.h"

/***********************************************************************
 * Create a confusion matrix from a segmented image and its Groundtruth
 ***********************************************************************/
 
ConfusionMatrix * CreateConfusionMatrix (Image *seg, Image *gt)
{
	ConfusionMatrix *rv;
	map<unsigned char,int> labels_seg, labels_gt;
	unsigned int cnt_seg, cnt_gt;
	unsigned int i;
	
	rv = new ConfusionMatrix;
	
	if (seg->xsize != gt->xsize || 
	    seg->ysize != gt->ysize)
		ERR_THROW ("CreateConfusionMatrix(): image sizes do not match!");
	
	// Get the number of labels in each image. Just set the index
	// to zero at first, it will be set to the real index in a 
	// second pass
	for (int y=0; y<seg->ysize; ++y)
	for (int x=0; x<seg->xsize; ++x)
		labels_seg[seg->get(x,y)]=0;
		
	for (int y=0; y<gt->ysize; ++y)
	for (int x=0; x<gt->xsize; ++x)
		labels_gt[gt->get(x,y)]=0;
		
	// Allocate the label vectors and the matrix
	cnt_seg = labels_seg.size();
	cnt_gt = labels_gt.size();
	rv->gv_seg.resize(cnt_seg);
	rv->gv_gt.resize(cnt_gt);	
	rv->mat.resize(cnt_seg, cnt_gt);
	rv->mat.setZero();
	
	// Fill the label vectors
	i=0;
	for (map<unsigned char,int>::iterator iter=labels_seg.begin(); 
	     iter!=labels_seg.end(); ++iter, ++i)
	{
		labels_seg[iter->first]=i;
		rv->gv_seg[i]=iter->first;
	}
	i=0;
	for (map<unsigned char,int>::iterator iter=labels_gt.begin(); 
	     iter!=labels_gt.end(); ++iter, ++i)
	{
		labels_gt[iter->first]=i;
		rv->gv_gt[i]=iter->first;
	}
	
	// Calculate the matrix
	for (int y=0; y<seg->ysize; ++y)
	for (int x=0; x<seg->xsize; ++x)
	{
		unsigned char pix_seg=labels_seg[seg->get(x,y)];
		unsigned char pix_gt =labels_gt [gt->get(x,y)];
		rv->mat(pix_seg,pix_gt) = rv->mat(pix_seg,pix_gt)+1;		
	}
}	


/***************************************************************************
 * Generate a random image
 ***************************************************************************/

Image *uniformNoise (int xsize, int ysize, int planes) 
{
	RandomNumberGenerator randgen (true, 0, 255, 4096);
	Image *rv;
	int np;
	
	// Allocate
	rv = new Image (xsize, ysize, planes);
	
	// The different color planes
	np = rv->nbColorPlanes();
	for (int col=1; col<=np; ++col) 
	{		
		// Travers the pixels and create them
		for (int y=0; y<ysize; ++y)
		for (int x=0; x<xsize; ++x)
			rv->set(col, x, y, randgen.nextByte());
	}	
	
	return rv;
}

/***************************************************************************
 * Generate an image with Gaussian noise
 * Equation taken from:
 * The Scientist and Engineer's Guide to Digital Signal Processing
 * by Steven W. Smith, Ph.D.California Technical Publishing
 * http://www.dspguide.com
 * 
 * Exists in versions for integer and for float images
 ***************************************************************************/

Image *gaussianNoise (int xsize, int ysize, int planes, float mean, float stddev) 
{
	RandomNumberGenerator randgen (true, 0, 1, 4096);
	Image *rv;
	int np;
	
	// Allocate
	rv = new Image (xsize, ysize, planes);
	
	// The different color planes
	np = rv->nbColorPlanes();
	for (int col=1; col<=np; ++col) 
	{		
		// Travers the pixels and create them
		for (int y=0; y<ysize; ++y)
		for (int x=0; x<xsize; ++x)
		{
			float r1 = randgen.next();
			float r2 = randgen.next();
			r1 = sqrt(-2.0*log(r1))*cos(2*M_PI*r2);
			r1 = r1*stddev + mean;
			TRIMGRAY(r1);
			rv->set(col, x, y, (unsigned char) rint(r1));
		}
	}	
	
	return rv;
}

FloatMatrix *gaussianNoise (int xsize, int ysize, float mean, float stddev) 
{
	RandomNumberGenerator randgen (true, 0, 1, 4096);
	FloatMatrix *rv;
	
	// Allocate
	rv = new FloatMatrix (xsize, ysize);
	
	// Travers the pixels and create them
	for (int y=0; y<ysize; ++y)
	for (int x=0; x<xsize; ++x)
	{
		float r1 = randgen.next();
		float r2 = randgen.next();
		r1 = sqrt(-2.0*log(r1))*cos(2*M_PI*r2);
		r1 = r1*stddev + mean;
		rv->set(x, y, r1);
	}
	
	return rv;
}

// *************************************************************
// Returns the probability a pixel having the value val
// *************************************************************

double getProb (Image &im, byte val) 
{
	int count=0;

    // Count the times this combination takes place
    for	(int y=0; y<im.ysize; ++y)
	for	(int x=0; x<im.xsize; ++x) 
	{
		if (im.get(1,x,y)==val)
			++count;
	}

	return (double) count / ((double) im.xsize*im.ysize);
}

// *************************************************************
// Returns the probability of two pixels having the two values
// val1 and val2, where the distance vector of the pixels is
// x+d_x,y+d_y.
// *************************************************************

double getProb (Image &im, byte val1, byte val2, int d_x, int d_y) 
{
	int count=0;
	int beg_x, beg_y, end_x, end_y;

	// Calculate the window border
    if (d_x < 0) 
	{
		beg_x = abs(d_x);
		end_x = 0;
    }
    else 
	{
    	beg_x = 0;
    	end_x = d_x;
    }

	if (d_y < 0) 
	{
		beg_y = abs(d_y);
		end_y = 0;
    }
    else 
	{
    	beg_y = 0;
    	end_y = d_y;
    }

    // Count the times this combination takes place
    for	(int y=beg_y; y<im.ysize-end_y; ++y)
	for	(int x=beg_x; x<im.xsize-end_x; ++x) 
	{
		if ((im.get(1,x,y)==val1) && (im.get(1,x+d_x,y+d_y)==val2))
			++count;
	}

	return (double) count / ((double) (im.xsize-beg_x-end_x) * (im.ysize-beg_y-end_y));
}


