/**************************************************************
 * Return a threshold surface for
 * Niblacks method
 * Also contains code which improves it:
 * (1) by Sauvola & Co.
 *     ICDAR 1997, pp 147-152
 * (2) by myself - Christian Wolf
 *     Research notebook 19.4.2001, page 129
 * (3) by myself - Christian Wolf
 *     20.4.2007
 *
 * See also:
 * Research notebook 24.4.2001, page 132 (Calculation of s)
 **************************************************************/
 
// From the IMAGE module
#include <Image.h> 
#include <FloatMatrix.h>

// From the PYRAMID module
#include <FloatPyramid.h>

// From the IMAGE PROCESSING module
#include <ImageProc.h>
 
// From this module
#include "Binarization.h" 

FloatMatrix * surfaceNiblackImproved (Image &im, NiblackVersion version,
	int winx, int winy, double k, double dR,
	FloatMatrix *& out_map_m, FloatMatrix *& out_map_s) {

	double m, s, max_s;
	double th=0;
	byte min_I;
	FloatMatrix *map_m, *map_s;
	FloatMatrix *ret_im;
	int wxh	= winx/2;
	int wyh	= winy/2;
	int x_firstth= wxh;
	int x_lastth = im.xsize-wxh-1;
	int y_lastth = im.ysize-wyh-1;
	int y_firstth= wyh;
	FloatPyramid *pyr=NULL;
	int pyr_access_level=0;
	float pyr_access_xfactor=0;
	float pyr_access_yfactor=0;
	FloatMatrix *pyr_access_top=NULL;
	int mx, my;

	ret_im = new FloatMatrix (im.xsize, im.ysize);

	// Create the local stats and store them in a map
	max_s = calcLocalStats (im, winx, winy, map_m, map_s);
	min_I = imageMin(im, PLANE_RED);
		
	// Create a pyramid which subsamples the standard dev
	if (version==NIBLACK_WOLF_2007)
	{
		int bigwin = (int) 2*(winx+winy);
		if (bigwin>im.xsize)
			bigwin=im.xsize;
		if (bigwin>im.ysize)
			bigwin=im.ysize;
		cerr << "Setting maximum(s) windows size to " << bigwin << endl;
		pyr = new FloatPyramid (map_s, (int) floor(log2(bigwin)), PYRT_MAX_RED22);
		pyr->build();
		
		pyr_access_level = pyr->levels()-1;
		pyr_access_xfactor = 1./pow(2.,pyr_access_level);
		pyr_access_yfactor = 1./pow(2.,pyr_access_level);
		pyr_access_top = (*pyr)[pyr_access_level];		
	}

	// In a second step, create the surface.
	// ----------------------------------------------------

	for	(int j = y_firstth ; j<=y_lastth; j++) {

		// NORMAL, NON-BORDER AREA IN THE MIDDLE OF THE WINDOW:
		for	(int i=0 ; i <= im.xsize-winx; i++) {

			m  = map_m->get(i+wxh, j);
    		s  = map_s->get(i+wxh, j);

    		// Calculate the threshold
    		switch (version) {

    			case NIBLACK_CLASSIC:
    				th = m + k*s;
    				break;

    			case NIBLACK_SAUVOLA:
	    			th = m * (1 + k*(s/dR-1));
	    			break;

	    		case NIBLACK_WOLF1:
    				th = m * (1 + k*(s/max_s-1));
    				break;

    			case NIBLACK_WOLF2:
    				th = m + k * (s/max_s-1) * (m-min_I);
    				break;
    				
    			case NIBLACK_WOLF_2007:
    				mx = (int) rint(i*pyr_access_xfactor);
    				my = (int) rint(j*pyr_access_yfactor);
    				if (mx>=pyr_access_top->xsize) mx=pyr_access_top->xsize-1;
    				if (my>=pyr_access_top->ysize) my=pyr_access_top->ysize-1;
    				max_s = pyr_access_top->get(mx,my);     					
    				th = m + k * (s/max_s-1) * (m-min_I);
    				break;

    			default:
    				cerr << "Unknown threshold type in ImageThresholder::surfaceNiblackImproved()\n";
    				exit (1);
    		}

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
        			for (int u=y_lastth+1; u<im.ysize; ++u)
        			for (int i=0; i<=x_firstth; ++i)
        				ret_im->set(i,u,th);
    		}

			// UPPER BORDER
			if (j==y_firstth)
				for (int u=0; u<y_firstth; ++u)
					ret_im->set(i+wxh,u,th);

			// LOWER BORDER
			if (j==y_lastth)
				for (int u=y_lastth+1; u<im.ysize; ++u)
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
			for (int u=y_lastth+1; u<im.ysize; ++u)
			for (int i=x_lastth; i<im.xsize; ++i)
				ret_im->set(i,u,th);
	}

	/*
	cerr << "Niblack: Writing maps ...\n";
	Image tmp1 (*map_m);
	Image tmp2 (*map_s);
	Image tmp3 (*ret_im);
	tmp1.write ("x_map_m");
	tmp2.write ("x_map_s");
	tmp3.write ("x_map_surface");
    */
    
    // Clean up;
    if (version==NIBLACK_WOLF_2007)
    	delete pyr;

	if (out_map_m!=NULL)
		out_map_m = map_m;
	if (out_map_s!=NULL)
		out_map_s = map_s;
	return ret_im;
}
