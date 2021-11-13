/**********************************************************
 * Binarization.cpp
 *
 * Thresholding algorithms for the Image class
 **********************************************************/

 // C
#include <string.h>
#include <math.h>
#include <errno.h>

 // C++
#include <iostream>
#include <set>

// From the IMAGE module
#include <Rect.h>
#include <FloatMatrix.h>
#include <color.h>

// From the MATH module
#include <Histogram.h>

// From the IMAGE PROCESSING module
#include <ImageProc.h>

// From this module
#include "Binarization.h"

// *************************************************************
// Threshold the image with a fixed value
// *************************************************************

void thresholdFixed (Image &im, int t) {
	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++)
		im.R[i][j]	= (im.R[i][j] > t ? 255 : 0) ;
}

// *************************************************************
// Threshold the image with a fixed value ->
// only threshold lower part, pixels above threshold stay.
// *************************************************************

void thresholdFixedLower (Image &im, int t) {
	byte r;
	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++) {
		r = im.R[i][j];
		im.R[i][j]	= (r > t ? r : 0) ;
	}
}

// *************************************************************
// Threshold the image with a fixed value ->
// only threshold upper part, pixels below threshold stay.
// *************************************************************

void thresholdFixedUpper (Image &im, int t) {
	byte r;
	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++) {
		r = im.R[i][j];
		im.R[i][j]	= (r > t ? 255 : r) ;
	}
}

// *************************************************************
// Threshold the image using a surface
// *************************************************************

void thresholdWithSurface (Image &im, FloatMatrix &surface) 
{
	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++)
		im.R[i][j]	= (im.R[i][j] >= surface.get(i,j) ? 255 : 0);
}

// *************************************************************
// Take a binarised image and calculate (in a windowed version)
// the mean of text and non-text and the variance of the
// noise process (compound variance of both distributions).
// *************************************************************

void getDistributions (Image &im, double khigh, double klow, int winx, int winy,
	FloatMatrix *&outim_t, FloatMatrix *& outim_b, FloatMatrix *& outim_var) {

	int wxh	= winx/2;
	int wyh	= winy/2;
	int x_firstth= wxh;
	int x_lastth = im.xsize-wxh-1;
	int y_firstth= wyh;
	int y_lastth = im.ysize-wyh-1;
	double t_sum, t_sum_sq, t_m=0, t_v,
		   b_sum, b_sum_sq, b_m=0, b_v,
		   var=0;
	int t_count, b_count;
	FloatMatrix *map_m, *map_s;

	outim_t   = new FloatMatrix (im.xsize, im.ysize);
	outim_b   = new FloatMatrix (im.xsize, im.ysize);
	outim_var = new FloatMatrix (im.xsize, im.ysize);
	outim_t->setZero();
	outim_b->setZero();
	outim_var->setZero();

	calcLocalStats (im, winx, winy, map_m, map_s);

	for	(int j = y_firstth ; j<=y_lastth; j++) {
		for	(int i=0 ; i <= im.xsize-winx; i++) {

			// NORMAL, NON-BORDER AREA IN THE MIDDLE OF THE WINDOW:
			t_sum = t_sum_sq = b_sum = b_sum_sq = t_count = b_count = 0;
			for	(int wy=0 ; wy<winy; wy++) {
    			for	(int wx=0 ; wx<winx; wx++) {
    				double m,s,thlow, thhigh;
    				double gv;

    				m  = map_m->get(i+wxh, j);
    				s  = map_s->get(i+wxh, j);
    				thlow  = m + klow * s;
    				thhigh = m + khigh* s;

    				gv = im.get(PLANE_RED,i+wx, j-wyh+wy);

    				// Pixel is text
    				if (gv < thlow) {
    					t_sum    += gv;
    					t_sum_sq += gv*gv;
    					++t_count;
    				}

    				// Pixel is background
    				if (gv > thhigh) {
    					b_sum    += gv;
    					b_sum_sq += gv*gv;
    					++b_count;
    				}

    			}
    		}
		    t_m = t_sum / (double) t_count;
		    b_m = b_sum / (double) b_count;
		    t_v = (t_sum_sq - (t_sum*t_sum)/(double)t_count)/(double) t_count;
		    b_v = (b_sum_sq - (b_sum*b_sum)/(double)b_count)/(double) b_count;
		    var = ((double) t_count * t_v + (double) b_count * b_v) /
		    	  	(double) (t_count + b_count);

		    outim_t->set  (i+wxh, j, t_m);
		    outim_b->set  (i+wxh, j, b_m);
		    outim_var->set(i+wxh, j, var);

		    if (i==0) {
        		// LEFT BORDER
        		for (int i=0; i<=x_firstth; ++i) {
                	outim_t->set  (i, j, t_m);
        		    outim_b->set  (i, j, b_m);
        		    outim_var->set(i, j, var);
        		}

        		// LEFT-UPPER CORNER
        		if (j==y_firstth)
        			for (int u=0; u<y_firstth; ++u)
        			for (int i=0; i<=x_firstth; ++i) {
        				outim_t->set  (i, u, t_m);
            		    outim_b->set  (i, u, b_m);
            		    outim_var->set(i, u, var);
            		 }

        		// LEFT-LOWER CORNER
        		if (j==y_lastth)
        			for (int u=y_lastth+1; u<im.ysize; ++u)
        			for (int i=0; i<=x_firstth; ++i) {
        				outim_t->set  (i, u, t_m);
            		    outim_b->set  (i, u, b_m);
            		    outim_var->set(i, u, var);
            		}
    		}

			// UPPER BORDER
			if (j==y_firstth)
				for (int u=0; u<y_firstth; ++u) {
					outim_t->set  (i+wxh, u, t_m);
        		    outim_b->set  (i+wxh, u, b_m);
        		    outim_var->set(i+wxh, u, var);
        		}

			// LOWER BORDER
			if (j==y_lastth)
				for (int u=y_lastth+1; u<im.ysize; ++u) {
					outim_t->set  (i+wxh, u, t_m);
        		    outim_b->set  (i+wxh, u, b_m);
        		    outim_var->set(i+wxh, u, var);
        		}

		}

		// RIGHT BORDER
		for (int i=x_lastth; i<im.xsize; ++i) {
        	outim_t->set  (i, j, t_m);
  		    outim_b->set  (i, j, b_m);
  		    outim_var->set(i, j, var);
  		}

  		// RIGHT-UPPER CORNER
		if (j==y_firstth)
			for (int u=0; u<y_firstth; ++u)
			for (int i=x_lastth; i<im.xsize; ++i) {
				outim_t->set  (i, u, t_m);
    		    outim_b->set  (i, u, b_m);
    		    outim_var->set(i, u, var);
    		}

		// RIGHT-LOWER CORNER
		if (j==y_lastth)
			for (int u=y_lastth+1; u<im.ysize; ++u)
			for (int i=x_lastth; i<im.xsize; ++i) {
				outim_t->set  (i, u, t_m);
    		    outim_b->set  (i, u, b_m);
    		    outim_var->set(i, u, var);
    		}
	}

	// Clean up
	delete map_m;
	delete map_s;
}

// *************************************************************
// Compute the optimal threshold with the fisher method
// and perform the thresholding
// *************************************************************

int thresholdFisher (Image &im, int kmin, int kmax, int deltaToMax) {
	Histogram<int> *hist;
	int kopt;

	hist = buildHistogram<int> (im);
	kopt = hist->getThresholdValueFisher (kmin, kmax);
	kopt +=	deltaToMax;
	thresholdFixed (im, kopt);

	delete hist;
	return kopt;
}

// *************************************************************
// Compute the optimal threshold with the fisher method
// and perform the thresholding.
// Used a windowed algorithm (pseudo-adaptiv).
// *************************************************************

FloatMatrix * surfaceFisherWindowed (Image &im, int kmin, int kmax,  int winx, int winy) {
	int wxh	= winx/2;
	int wyh	= winy/2;
	int x_firstth= wxh;
	int x_lastth = im.xsize-wxh-1;
	int y_lastth = im.ysize-wyh-1;
	int y_firstth= wyh;
	FloatMatrix *ret_im;
	Histogram<int> *h;
	Image *tb;
	int th=0;

	ret_im = new FloatMatrix (im.xsize, im.ysize);

	for	(int j = y_firstth ; j<=y_lastth; j++) {

		// NORMAL, NON-BORDER AREA IN THE MIDDLE OF THE WINDOW:
		for	(int i=0 ; i <= im.xsize-winx; i++) {

			// Cut the window out of the image and calcuate the fisher value
    		Rect r (j-wyh,j+wyh,i,i+winx-1);
    		tb  = im.cutSubImage (r,0,0);
    		h = buildHistogram<int> (*tb);
    		th = h->getThresholdValueFisher (kmin, kmax);

    		delete h;
    		delete tb;

    		ret_im->set(i+wxh,j,th);

    		if (i==0) {
        		// LEFT BORDER
        		for (int i=0; i<=x_firstth; ++i)
                	ret_im->set(i,j,th);

        		// LEFT-UPPER CORNER
        		if (j==y_firstth)
        			for (int u=0; u<y_firstth; ++u)
        			for (int i=0; i<=x_firstth; ++i)
        				ret_im->set(i,u,th);

        		// LEFT-LOWER CORNER
        		if (j==y_lastth)
        			for (int u=y_lastth+1; u<im.ysize-1; ++u)
        			for (int i=0; i<=x_firstth; ++i)
        				ret_im->set(i,u,th);
    		}

			// UPPER BORDER
			if (j==y_firstth)
				for (int u=0; u<y_firstth; ++u)
					ret_im->set(i+wxh,u,th);

			// LOWER BORDER
			if (j==y_lastth)
				for (int u=y_lastth+1; u<im.ysize-1; ++u)
					ret_im->set(i+wxh,u,th);
		}

		// RIGHT BORDER
		for (int i=x_lastth; i<im.xsize; ++i)
        	ret_im->set(i,j,th);

  		// RIGHT-UPPER CORNER
		if (j==y_firstth)
			for (int u=0; u<y_firstth; ++u)
			for (int i=x_lastth; i<im.xsize; ++i)
				ret_im->set(i,u,th);

		// RIGHT-LOWER CORNER
		if (j==y_lastth)
			for (int u=y_lastth+1; u<im.ysize-1; ++u)
			for (int i=x_lastth; i<im.xsize; ++i)
				ret_im->set(i,u,th);
	}

    return ret_im;
}



