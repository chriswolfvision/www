// *************************************************************
// ImageDataOperator.cc
// Methods for conversion et. al for the image class
//
// author: Christian Wolf
// *************************************************************

// C
#include <errno.h>
#include <string.h>
#include <math.h>

// C++
#include <iostream>
#include <vector>

using namespace std;

// Port to Visual C++
#include <visualbug.h>

// From the main module
#include <CIL.h>

// From this module
#include "Image.h"
#include "FloatMatrixColor.h"
#include "BicubicPolynom.h"
#include "Rect.h"
#include "color.h"

#define MAXSIZE		4092

// *************************************************************
// Help function for the interpolation method
// *************************************************************

inline double weight_func (int x, int x0) {
	double d = ::sqrt(x*x + x0*x0);
	return d == 0 ? 1 : 1/d;
}

// *************************************************************
// Interpolate an image, i.e. increase its resolution,
// After Sato, Kanade et. al,
// Multimedia Systems 7(5) 1999, pp 385-395
// *************************************************************

Image *Image::interpolate (int factor) {
	Image *newim;
	int nys=factor*ysize, nxs=factor*xsize;
	int old_x_l, old_x_r,
		old_y_u, old_y_d;
	int new_x_l, new_x_r,
		new_y_u, new_y_d;

#ifdef CHECK_CODE
	if (factor<1||type!=3) {
		cerr <<	"Internal error in Image::interpolate()!\n";
		CORE_DUMP;
	}
#endif

	newim = new Image (nxs, nys, type);

	// Travers the pixels of the new interpolated image
	// and fill in

	for (int y=0; y<nys; ++y) {

		double ry = ((double) y / (double) factor);

		old_y_u = (int) floor(ry);
		old_y_d = (int) ceil(ry);
		new_y_u = 4*old_y_u;
		new_y_d = 4*old_y_d;

		for (int x=0; x<nxs; ++x) {
			double L=0;
			double w1, w2, w3, w4, weight_tot=0;

			double rx = ((double) x / (double) factor);

			// Fast special case if we hit the pixel exactly
			if (rx==rint(rx) && ry==rint(ry)) {
				newim->set(1,x,y, get(1,(int)rx,(int)ry));
				newim->set(2,x,y, get(2,(int)rx,(int)ry));
				newim->set(3,x,y, get(3,(int)rx,(int)ry));
				continue;
			}

			old_x_l = (int) floor(rx);
			old_x_r = (int) ceil(rx);
			new_x_l = 4*old_x_l;
			new_x_r = 4*old_x_r;

			if (old_x_r >= xsize) old_x_r = xsize-1;
			if (old_y_d >= ysize) old_y_d = ysize-1;

			// The weights
			w1 = weight_func (x-new_x_l, y-new_y_u);
			w2 = weight_func (x-new_x_r, y-new_y_u);
			w3 = weight_func (x-new_x_r, y-new_y_u);
			w4 = weight_func (x-new_x_l, y-new_y_d);
			weight_tot = w1 + w2 + w3 + w4;

			// For all 3 color planes
			DO_COLORS {
    			L = (w1 * get (col,old_x_l,old_y_u) +
    				 w2 * get (col,old_x_r,old_y_u) +
                     w3 * get (col,old_x_r,old_y_d) +
                     w4 * get (col,old_x_l,old_y_d)) / weight_tot;

    			if (L<0)   L = 0;
    			if (L>255) L = 255;
    			newim->set(col,x,y,(int) rint(L));
			}
	    }
	}

	return newim;
}

// *************************************************************
// Interpolate an image, i.e. increase its resolution,
// but just copy pixles (the "blocky" thing to do).
// *************************************************************

Image * Image::interpolateBlocks (int factor) {
	Image *newim;
	int nx, ny;

	nx = factor*xsize;
	ny = factor*ysize;
	newim = new Image (nx, ny, type);

	for (int y=0; y<ny; ++y) {
		for (int x=0; x<nx; ++x) {
			DO_AVAILABLE_COLORS {
				newim->set(col, x, y, get(col, x/factor, y/factor));
			}
		}
	}

	return newim;
}

// *************************************************************
// Interpolate an image, i.e. increase its resolution,
// Robust Method after Jean-Michel Jolion,
// Research notebook 23.1.2001, page 62-63.
//
// The mean and variance images are optional. If they are not
// specified (=NULL) then they are not used.
// *************************************************************

FloatMatrixColor *Image::interpolateJolion (FloatMatrixColor *avgImage,
									 	   FloatMatrixColor *varImage,
									       int xoffset, int yoffset, int factor) {

	FloatMatrixColor *newim;
	int nys=factor*ysize, nxs=factor*xsize;
	int old_x_l, old_x_r,
		old_y_u, old_y_d;
	int new_x_l, new_x_r,
		new_y_u, new_y_d;

#ifdef CHECK_CODE
	if (factor<1) {
		cerr <<	"Internal error in Image::interpolateJolion()!\n";
		CORE_DUMP;
	}
#endif

	newim = new FloatMatrixColor (nxs, nys);

	// Travers the pixels of the new interpolated image
	// and fill in

	for (int y=0; y<nys; ++y) {

		double ry = ((double) y / (double) factor);

		old_y_u = (int) floor(ry);
		old_y_d = (int) ceil(ry);
		new_y_u = 4*old_y_u;
		new_y_d = 4*old_y_d;

		for (int x=0; x<nxs; ++x) {

			double w1, w2, w3, w4, weight_tot=0;
			double g1, g2, g3, g4;
			double L;
			double v, rx;

			rx = ((double) x / (double) factor);
			// Fast special case if we hit the pixel exactly
			/*
			if (rx==rint(rx) && ry==rint(ry)) {
				newim->set(1,x,y, get(1,(int)rx,(int)ry));
				newim->set(2,x,y, get(2,(int)rx,(int)ry));
				newim->set(3,x,y, get(3,(int)rx,(int)ry));
				continue;
			}
			*/

			old_x_l = (int) floor(rx);
			old_x_r = (int) ceil(rx);
			new_x_l = 4*old_x_l;
			new_x_r = 4*old_x_r;

			if (old_x_r >= xsize) old_x_r = xsize-1;
			if (old_y_d >= ysize) old_y_d = ysize-1;

			// The weights depending on the distance to the corner points
			w1 = weight_func (x-new_x_l, y-new_y_u);
			w2 = weight_func (x-new_x_r, y-new_y_u);
			w3 = weight_func (x-new_x_r, y-new_y_d);
			w4 = weight_func (x-new_x_l, y-new_y_d);

			// Travers all 3 color planes
			DO_AVAILABLE_COLORS {

    			// The weights depending on the value of the pixel,
    			// i.e. if it is an outlier or not.
    			if (avgImage!=NULL && varImage!=NULL) {
        			g1 = get (col,old_x_l,old_y_u) - avgImage->get(col,old_x_l+xoffset, old_y_u+yoffset);
        			g2 = get (col,old_x_r,old_y_u) - avgImage->get(col,old_x_r+xoffset, old_y_u+yoffset);
        			g3 = get (col,old_x_r,old_y_d) - avgImage->get(col,old_x_r+xoffset, old_y_d+yoffset);
        			g4 = get (col,old_x_l,old_y_d) - avgImage->get(col,old_x_l+xoffset, old_y_d+yoffset);
        			if (g1<0) g1 *= -1.0;
        			if (g2<0) g2 *= -1.0;
        			if (g3<0) g3 *= -1.0;
        			if (g4<0) g4 *= -1.0;
        			v = varImage->get(col, old_x_l+xoffset, old_y_u+yoffset); if (v!=0) g1 /= v;
        			v = varImage->get(col, old_x_r+xoffset, old_y_u+yoffset); if (v!=0) g2 /= v;
        			v = varImage->get(col, old_x_r+xoffset, old_y_d+yoffset); if (v!=0) g3 /= v;
        			v = varImage->get(col, old_x_l+xoffset, old_y_d+yoffset); if (v!=0) g4 /= v;
        			g1 = 1.0 / ( 1.0 + g1);
        			g2 = 1.0 / ( 1.0 + g2);
        			g3 = 1.0 / ( 1.0 + g3);
        			g4 = 1.0 / ( 1.0 + g4);
        		}
        		else
        			g1 = g2 = g3 = g4 = 1.0;

    			weight_tot = w1*g1+w2*g2+w3*g3+w4*g4;
    			if (weight_tot == 0)
    				weight_tot = 1;

    			L = (w1*g1*get (col,old_x_l,old_y_u) +
    				 w2*g2*get (col,old_x_r,old_y_u) +
    				 w3*g3*get (col,old_x_r,old_y_d) +
    				 w4*g4*get (col,old_x_l,old_y_d)) / weight_tot;

    			if (L<0)   L=0;
    			if (L>255) L = 255;

    			newim->set(col,x,y,(int) rint(L));
			}
	    }
	}

	return newim;
}

// *************************************************************
// Interpolate an image, i.e. increase its resolution,
// The method uses a Bicubic interpolation scheme, combined with
// the robust approach Jean-Michel Jolion developped for the
// bilinear interpolation.
//
// See Research notebook 21.2.2001, p.71.
// The variables are those described in the notebook
//
// The mean and variance images are optional. If they are not
// specified (=NULL) then they are not used.
// *************************************************************

FloatMatrixColor *Image::interpolateBicubicRobust (FloatMatrixColor *avgImage,
									 	   		  FloatMatrixColor *varImage,
									              int xoffset, int yoffset, int factor) {



#ifdef CHECK_CODE
	if (factor<1) {
		cerr <<	"Internal error in Image::interpolateBicubicRobust()!\n";
		CORE_DUMP;
	}
#endif

	int nys=factor*ysize, nxs=factor*xsize;
	int oysize=ysize, oxsize=xsize;
	FloatMatrixColor *newim = new FloatMatrixColor (nxs, nys);
	BicubicPolynom *poly = BicubicPolynom::instance(factor);

	// Travers the pixels of the new interpolated image
	// and fill in
	for (int y=0; y<nys; ++y) {

		double My = ((double) y / (double) factor);
		int Ry = (int) floor(My);

		for (int x=0; x<nxs; ++x) {
			double a,b;

			double Mx = ((double) x / (double) factor);
			int Rx = (int) floor(Mx);

			if (Rx >= xsize) Rx = xsize-1;
			if (Ry >= ysize) Ry = ysize-1;

			a = Mx-(double)Rx;
			b = My-(double)Ry;

			// All colors
			DO_AVAILABLE_COLORS {
				double L, weight_tot;

    			// Travers the 16 neighbor pixels
    			L=0;
    			weight_tot = 0;

    			for (int m=-1; m<=2; ++m)
    			for (int n=-1; n<=2; ++n) {
    				int nx,ny;
    				double w,v;

    				// Treatment of the border
    				nx = Rx+m;
    				ny = Ry+n;
    				if (nx<0) nx = 0; if (nx>=oxsize) nx=oxsize-1;
    				if (ny<0) ny = 0; if (ny>=oysize) ny=oysize-1;

    				// The "robustness" weights, which depend on the value of the pixel,
    				// i.e. if it is an outlier or not.
    				if (avgImage != NULL) {
    					w = get (col,nx,ny) - avgImage->get(col,nx+xoffset, ny+yoffset);
            			if (w<0) w *= -1.0;
            			v = varImage->get(col, nx+xoffset, ny+yoffset);
            			if (v!=0) w /= v;
            			w = 1.0 / (1.0 + w);
            		}
            		else
            			w = 1.0;

  					L += get (col,nx,ny) * (*poly)(m-a) * (*poly)(b-n) * w;
   					weight_tot += w;
   				}

   				L *= 16.0;
   				if (weight_tot != 0)
	   				L /= weight_tot;

				if (L<0)   L=0;
    			if (L>255) L=255;
   				newim->set (col,x,y, L);
			}
	    }
	}

	return newim;
}

// *************************************************************
// Rescale the image by a factor 2.
// Code taken from the Pyramid library of Jean-Michel Jolion
// and rewritten for the Image class
// *************************************************************

Image * Image::subSample (int factor) {
	Image *rv;
    unsigned int nxs=xsize/factor,
    		     nys=ysize/factor;
	

	// Allocation of the new image
	rv = new Image (nxs, nys, type);

	DO_AVAILABLE_COLORS {

		for (unsigned int y=0; y<nys; ++y)
		for (unsigned int x=0; x<nxs; ++x) {

			rv->set(col,x,y, get(col,x*factor,y*factor));
		}
	}

	return rv;
}

// *************************************************************
// Rescale the image by a factor 2.
// Code taken from the Pyramid library of Jean-Michel Jolion
// and rewritten for the Image class
// *************************************************************

Image * Image::reScale2x2 ()	{
	FloatMatrix *C;
	Image *rv;
	int	i, j;
	int val;

#ifdef CHECK_CODE
	if ((xsize>MAXSIZE)||(ysize>MAXSIZE)) {
		cerr << "Internal error in Image::reScale2x2()!\n";
		exit (1);
	}
#endif

	// Allocation of the new image
	rv = new Image (xsize/2, ysize/2, type);

	// A temporary buffer
	C = new FloatMatrix(xsize+2, ysize+2);
	
	DO_AVAILABLE_COLORS {
	
		int xoffset, yoffset;
		int XS=xsize,
			YS=ysize;

    	// Where needs the image to be copied in the buffer
    	// depending if height and width are even or odd
    	if ((YS%2) && (YS!=1)) {
    		yoffset=0;
    		if ((XS%2) &&	(XS!=1))
    			xoffset=0;
    		else
    			xoffset=1;
    	}
    	else {
    		yoffset=1;
    		if ((XS%2) &&	(XS!=1))
    			xoffset=0;
    		else
    			xoffset=1;
    	}

    	// Copy the original image
    	for(i=0;i<YS;i++)
    		for(j=0;j<XS;j++)
    			C->set(j+xoffset, i+yoffset, get(col, j, i));


        // The height is even, copy the border line
		if (!(YS%2)) {
			for(j=0;j<XS;j++)
				C->set(j+xoffset, 0, C->get(j+xoffset, 1));

			if (!(XS%2))
				C->set(0, 0, C->get(1, 1));
		}

		// The width is even, copy the border line
		if (!(XS%2)) {
			for(i=0;i<YS;i++)
				C->set(0, i+yoffset, C->get(0, i+yoffset));
		}	
		
    	// readjustment	of the borders
    	if (!(YS%2))
    		YS++;
    	if (!(XS%2))
    		XS++;


    			
    	// OK, maintenant C	contient l'image a reduire par masque 3x3
    	for(i=1;i<YS;i+=2)
    		for(j=0;j<XS;j++)
    			C->set(j, i, C->get(j,i-1) +
    					   2*C->get(j,i) +
    						 C->get(j,i+1));
    								
    	/* reduction */
    	if (YS!=1) {
    		if (XS!=1) {
    			for(i=1;i<YS;i+=2)
    				for(j=1;j<XS;j+=2)  {
    					val = (int) rint ((C->get(j-1,i)+
    									 2*C->get(j,  i)+
    									   C->get(j+1,i))/16.0);
                    	TRIMGRAY(val);
    					rv->set(col, j/2, i/2, (byte)val);
    				}
    		}
    		else {
    			for(i=1;i<YS;i+=2) {
    				val = (int) rint ((C->get(1,i-1)+
    								 2*C->get(1,i)+
    								   C->get(1,i+1))/4.0);
    				TRIMGRAY(val);
    				rv->set(col, 0,i /2, (byte)val);
    			}
    		}
    	}
    	else  {
    		if (XS!=1) {
    			for(j=1;j<XS;j+=2) {
    				val =  (int) rint((C->get(j-1,1)+
    								 2*C->get(j,1)+
    								   C->get(j+1,1))/4.0);
    				TRIMGRAY(val);
    				rv->set(col, j/2, 0, (byte)val);
    			}
    					
    		}
    		else {
    			val = (int) rint (C->get(1,1)/4.0);
    			TRIMGRAY(val);
    			rv->set(col, 0, 0, (byte) val);				
    		}
    	}
	
	}
	
	return rv;
}

// *************************************************************
// Paste another image at a	specified position
// - The whole image
// *************************************************************

void Image::paste (Image &other, int xpos, int ypos) {	
	paste (other, xpos, ypos, 0, 0, other.xsize, other.ysize);
}

// *************************************************************
// Paste another image at a	specified position
// - Parts of the image
// *************************************************************

void Image::paste (Image &other,int dxpos,int dypos,int sxpos,int sypos,int xlen,int ylen) {
	int rxlen, rylen;
	
#ifdef CHECK_CODE
	if ((dxpos >= xsize) || (dypos >= ysize)) {
		cerr <<	"Error in Image::paste()!!!\n"
			 << "dxpos: " << dxpos << " xsize: " << xsize << endl
			 << "dypos: " << dypos << " ysize: " << ysize << endl;
		cerr << "Provoking core dump...\n";
		CORE_DUMP;
	}
#endif	

	// Check if we have enough space in the image.
	// If not, cut the source image
	rxlen = xsize-dxpos;
	rylen = ysize-dypos;
	rxlen = rxlen < xlen ? rxlen : xlen;
	rylen = rylen < ylen ? rylen : ylen;	

	for	(int y=0; y<rylen; ++y) {
		for	(int x=0; x<rxlen; ++x) {

			R[dxpos+x][dypos+y] = other.R[sxpos+x][sypos+y];

			if (type==3)	{
				
				// Both are color images, just copy
				if (other.type==3)	{
					G[dxpos+x][dypos+y] = other.G[sxpos+x][sypos+y];
					B[dxpos+x][dypos+y] = other.B[sxpos+x][sypos+y];
				}
				
				// The whole image is color, the other one is grayscale
				else {
					G[dxpos+x][dypos+y] = other.R[sxpos+x][sypos+y];
					B[dxpos+x][dypos+y] = other.R[sxpos+x][sypos+y];
				}
			}
		}
	}
}

// *************************************************************
// Copy	a color	plane
// *************************************************************

void Image::copyPlane (int dst,	int	src) {

	unsigned char **D, **S;

	D =	(dst==1	? R	: (dst == 2	? G	: B));
	S =	(src==1	? R	: (src == 2	? G	: B));

	int	x,y;
	for	(y=0; y<ysize; ++y) {
		for	(x=0; x<xsize; ++x) {
			D[x][y]	= S[x][y] ;
		}
	}
}

// *************************************************************
// Set a plane to a	constant value
// *************************************************************

void Image::setPlaneToValue	(int plane,	byte val) {
	unsigned char **D;

	D =	(plane==1 ?	R :	(plane == 2	? G	: B));

	for	(int y=0; y<ysize; ++y) {
		for	(int x=0; x<xsize; ++x) {
			D[x][y]	= val;
		}
	}
}
// *************************************************************
// clear the border
// *************************************************************

void Image::setBorderToValue (int plane, int bordersize, byte value) {

	for (int x=0; x<xsize; ++x) {
		for (int y=0; y<bordersize; ++y) {		
			set (plane, x,y, value);

			set (plane, x,ysize-y-1, value);
		}
	}
	for (int y=0; y<ysize; ++y) {
		for (int x=0; x<bordersize; ++x) {		
			set (plane, x, y , value);
			set (plane, xsize-x-1, y, value);
		}
	}
}

// *************************************************************
// Convert the image to	grayscale
// *************************************************************

void Image::convertRGB2GrayScale ()	{

	// No color	image
	if (type!=3)
		return;

	for	(int y=0; y<ysize; ++y) {
		for	(int x=0; x<xsize; ++x) {
			R[x][y]	= (int)	(0.299*R[x][y] + 0.587*G[x][y] + 0.114*B[x][y])/1;
		}
	}
	type = 2;

	FREE_IMAGE (G);
	FREE_IMAGE (B);
}

// *************************************************************
// Convert the image to	grayscale
// *************************************************************

Image * Image::convertRGB2GrayScaleReturn ()	{
    Image *ret;
	
	// No color	image
	if (type!=3)
		return NULL;
		
	ret = new Image (xsize, ysize, 2);

	for	(int y=0; y<ysize; ++y) {


		for	(int x=0; x<xsize; ++x) {
			ret->R[x][y]	= (int)	(0.299*R[x][y] + 0.587*G[x][y] + 0.114*B[x][y])/1;
		}
	}
	
	return ret;
}

// *************************************************************
// Delete the G and B color planes
// *************************************************************

void Image::convertR__2GrayScale ()	{

	// No color	image
	if (type!=3)
		return;

	type = 2;


	FREE_IMAGE (G);
	FREE_IMAGE (B);
}

// *************************************************************
// Convert the image to	a color	image
// But only allocate color planes, do not copy the pixels.
// *************************************************************

void Image::convertGrayScale2R__ ()	{

	// Image is	already	a color	image
	if (type==3)
		return;

	G =	CREATE_IMAGE (ysize,xsize);
	B =	CREATE_IMAGE (ysize,xsize);

	type = 3;
}

// *************************************************************
// Convert the image to	a color	image
// *************************************************************

void Image::convertGrayScale2RGB ()	{
	convertGrayScale2R__();		
	for	(int y=0; y<ysize; ++y) {
		for	(int x=0; x<xsize; ++x) {
			G[x][y] = B[x][y] = R[x][y];
		}
	}
}

// *************************************************************
// Convert into	LUV	and	scale afterwards
// *************************************************************

void Image::convertRGB2LUV () {
	int	x,y;
	double r,g,b,cl,cu,cv;

	// No color	image
	if (type!=3) {
		cerr <<	"Error in ImageDataOperator:convertRGB2LUV\n";
		exit (0);
	}

	init_color ();

	for	(y=0; y<ysize; ++y) {
		for	(x=0; x<xsize; ++x) {

			r =	get(1,x,y);
			g =	get(2,x,y);
			b =	get(3,x,y);
			r /= MAX_RED;
			g /= MAX_GREEN;
			b /= MAX_BLUE;

			rgb2luv	(r,g,b,&cl,&cu,&cv);

			set	(1,x,y,(int) scale_value(cl, MINVALUE_L, MAXVALUE_L,
				0, 255)/1);
			set	(2,x,y,(int) scale_value(cu, MINVALUE_U, MAXVALUE_U,
				0, 255)/1);
			set	(3,x,y,(int) scale_value(cv, MINVALUE_V, MAXVALUE_V,
				0, 255)/1);
		}
	}
}

// *************************************************************
// Convert into	RGB	and	scale the values
// *************************************************************

void Image::convertLUV2RGB () {
	int	x,y;
	double r,g,b,l,u,v;
	byte br, bb, bg;

	// No color	image
	if (type!=3) {
		cerr <<	"Error in Image:convertRGB2LUV\n";
		exit (1);
	}

	init_color ();


	for	(y=0; y<xsize; ++y) {
		for	(x=0; x<xsize; ++x) {

			l =	get(1,x,y);
			u =	get(2,x,y);
			v =	get(3,x,y);

			l =	scale_value	(l,	0, 255,	MINVALUE_L,	MAXVALUE_L);
			u =	scale_value	(u,	0, 255,	MINVALUE_U,	MAXVALUE_U);
			v =	scale_value	(v,	0, 255,	MINVALUE_V,	MAXVALUE_V);

			luv2rgb	(l,u,v,&r,&g,&b);

			r *= MAX_RED;
			g *= MAX_GREEN;
			b *= MAX_BLUE;

			if (r>255) br=255; else	br = (int) r/1;
			if (g>255) bg=255; else	bg = (int) g/1;
			if (b>255) bb=255; else	bb = (int) b/1;


			set	(1,x,y,br);
			set	(2,x,y,bg);
			set	(3,x,y,bb);
		}
	}
}


// *************************************************************
// Superimpose a second	image;
// *************************************************************

void Image::superImpose	(Image &other, int val)	{
	unsigned char **I, **J;

	J=R;
	I=other.R;

	for	(int i = 0 ; i < xsize ; i++)
		for	(int j = 0 ; j < ysize ; j++)
			if (I[i][j]	== 255)
				J[i][j]	= val ;
			else
				J[i][j]	= (I[i][j] > J[i][j] ? I[i][j] : J[i][j]) ;
}

// *************************************************************
// Change a	color
// *************************************************************

void Image::colorChange	(byte oldr,	byte oldg, byte	oldb,
						 byte newr,	byte newg, byte	newb) {

	// No color	image
	if (type!=3) {
		cerr <<	"Error in Image:colorChange()\n";
		exit (1);
	}

	for	(int y=0; y<ysize; ++y) {
		for	(int x=0; x<xsize; ++x) {

			if ((get(1,x,y)==oldr) && (get(2,x,y)==oldg) &&	(get(3,x,y)==oldb))	{
				set	(1,x,y,newr);
				set	(2,x,y,newg);
				set	(3,x,y,newb);
			}
		}
	}
}

// *************************************************************
// Change a	gray value
// *************************************************************

void Image::grayvalueChange	(byte oldr,	byte newr) {

	for	(int y=0; y<ysize; ++y) {
		for	(int x=0; x<xsize; ++x) {

			if (get(1,x,y)==oldr)
				set	(1,x,y,newr);		
		}
	}
}

// **********************************************************
// Take	the	image and cut a	subimages out
// **********************************************************

void Image::crop (Rect &r, int growX, int growY) {
	Image *im;
	Rect b;
	int width, height;

	b =	r;
	b.growAndClip (growX, growY, xsize, ysize);

	// Check the size
	width = b.width();
	height = b.height();
	if ((width==0) || (height==0))
		return;
		
	// Create the image		
	im = new Image (width, height, type);
	
	// Copy	the	part of	the	big	original image
	for	(int y=0; y<b.height();	++y) {
		for	(int x=0; x<b.width(); ++x)	{

			im->set	(1,x,y,	get(1,b.left+x,b.top+y));

			// Color images
			if (im->type ==	3) {
				im->set	(2,x,y,	get(2,b.left+x,b.top+y));
				im->set	(3,x,y,	get(3,b.left+x,b.top+y));
			}
		}
	}
	
	*this = *im;
}

// **********************************************************
// Copy the interior pixels into the border region
// **********************************************************

void Image::copyBorderPixels(int border_x, int border_y) {
	int ul, ur, ll, lr;

    // VERTICAL FLANKS
    for (int x=border_x; x<xsize-border_x; ++x) {
 		byte upper = get(PLANE_RED, x,border_y);
 		byte lower = get(PLANE_RED, x,ysize-border_y-1);
 		for (int y=0; y<border_y; ++y) {		
 			set (PLANE_RED, x,y, upper);
 			set (PLANE_RED, x,ysize-y-1, lower);
 		}
 	}
 	
 	// HORIZONTAL FLANKS
 	for (int y=border_y; y<ysize-border_y; ++y) {
 		byte left = get(PLANE_RED, border_x,y);
 		byte right = get(PLANE_RED, xsize-border_x-1,y);
 		for (int x=0; x<border_x; ++x) {		
 			set (PLANE_RED, x, y , left);
 			set (PLANE_RED, xsize-x-1, y, right);
 		}
 	}

 	ul = get (PLANE_RED, border_x, border_y);
 	ur = get (PLANE_RED, xsize-border_x-1, border_y);
 	ll = get (PLANE_RED, border_x, ysize-border_y-1);
 	lr = get (PLANE_RED, xsize-border_x-1, ysize-border_y-1);
 	for (int x=0; x<border_x; ++x) {
	 	for (int y=0; y<border_y; ++y) {
			set (PLANE_RED, x,y, ul);
			set (PLANE_RED, xsize-1-x,y, ur);
			set (PLANE_RED, x,ysize-1-y, ll);
			set (PLANE_RED, xsize-1-x,ysize-1-y, lr);
	 	}
	}
}

// *************************************************************
// Take an image which consists of one column of tiled images
// and produce an image with a maximum height and more than one
// column of tiled images if necessary
// *************************************************************

void Image::separateColumns (int maxheight, int border, byte background) {
	vector<int> colxpos;
	vector<int> colxsize;
	int max_width;
	Image *full_image;
	int full_x, full_y;
	int run_x, run_y;
	bool tiles_left;
	
	// Nothing to do
	if (ysize<=maxheight)
		return;
		
	// Determine the maximum width of a tile
	max_width = 0;
	for (int y=0; y<ysize; ++y) {
		for (int x=0; x<xsize; ++x) {
			if (R[x][y]!=background)
				if (x>max_width)
					max_width=x;
		}
	}		
	
	// Allocate an image containing the maximum possible image
	// We will crop it afterwards if we need less space.
	full_image = new Image ((ysize/maxheight+2)*(max_width+border),maxheight,type);		
	full_image->setPlaneToValue (PLANE_RED, background);
	if (type==3) {
        full_image->setPlaneToValue (PLANE_GREEN, background);
        full_image->setPlaneToValue (PLANE_BLUE, background);
	}
	full_x = 0;
	full_y = 0;	
	
	run_y = 0;
	run_x = 0;
	do {
	
		int y_max_section,
			y_empty_line,
			col_width;
				
		// Search the first empty line in this section.
		y_max_section = run_y+maxheight-1 < ysize-1 ? run_y+maxheight-1 : ysize-1;
		y_empty_line = y_max_section;
    	for (int y=y_max_section; y>=run_y; --y) {
    		bool isempty=true;
    		for (int x=0; x<xsize; ++x) {
    			if (R[x][y]!=background)
    				isempty=false;
    		}
    		if (isempty) {
    			y_empty_line = y;
    			break;
    		}
    	}
				
    	// Travers the lines until the first empty line and determine the
    	// column width;
    	col_width = 0;
    	for (int y=run_y; y<y_empty_line; ++y) {
    		for (int x=0; x<xsize; ++x) {
    			if (R[x][y]!=background)

    				if (x>col_width)
    					col_width=x;
    		}    		
    	}
    	col_width +=2;
    	
    	// Copy the column into the output image
    	for (int y=0; y<y_empty_line-run_y; ++y) {
    		for (int x=0; x<col_width; ++x) {
    			full_image->set(PLANE_RED,full_x+x,full_y+y, R[run_x+x][run_y+y]);
    			if (type==3) {
    				full_image->set(PLANE_GREEN,full_x+x,full_y+y, G[run_x+x][run_y+y]);
    				full_image->set(PLANE_BLUE, full_x+x,full_y+y, B[run_x+x][run_y+y]);
    			}    			
    		}
    	}
       	
    	run_y = y_empty_line;
    	full_x += col_width + border;
    	full_y = 0;    	    	
    	
    	// Check if there are still image tiles in the remaining part
    	tiles_left = false;
    	for (int y=run_y; y<=y_max_section; ++y) {    		
    		for (int x=0; x<xsize; ++x) {
    			if (R[x][y]!=background) {
    				tiles_left=true;
    				break;
    			}
    		}
    	}
	}	
	while (tiles_left);	
	
	// We copied all tiles, not crop the new image to its _real_ size
	Rect r(0,maxheight,0,full_x);
	full_image->crop (r,0,0);
	
	*this = *full_image;
}

// *************************************************************
// Raise the image to the next power of 2
// Set the new pixels to the mean value of the old image
// returns whether the size was changed or not.
// *************************************************************

bool Image::toNextPowerOf2 () 
{
	int xs = xsize,
		ys = ysize;
	int newxs=xs,newys=ys;
	bool ok=true;
	float c;
		
	// Check whether we already a power of 2
	c = logf((float)xsize)/logf(2.0);
	if (c!=truncf(c))
	{
		ok=false;
		c=ceilf(c);
		newxs = (int) rint(powf(2.0,ceilf(c)));
	}
	c = logf((float)ysize)/logf(2.0);
	if (c!=truncf(c))
	{
		ok=false;
		newys = (int) rint(powf(2.0,ceilf(c))); 
	}
	
	// We are not. Calculate the new size
	if (!ok)
	{
		// cerr << "New image size: " << newxs << "x" << newys << endl;
		Image newimage(newxs,newys,nbColorPlanes());
		float m[3];
		
		// Calculate the mean value
		m[0]=m[1]=m[2]=0;
		for (int y=0; y<ys; ++y)
		for (int x=0; x<xs; ++x)
		{
			DO_COLORS_IF(nbColorPlanes())
				m[col-1] += get(col,x,y);
		}
		DO_COLORS_IF(nbColorPlanes())
			m[col-1] /= (float) (xs*ys);
			
		// Set the new background value to the mean value
		for (int y=0; y<newimage.ysize; ++y)
		for (int x=0; x<newimage.xsize; ++x)
		{
			DO_COLORS_IF(nbColorPlanes())
				newimage.set(col,x,y,(byte) rint(m[col-1]));
		} 
		
		// Paste
		newimage.paste(*this,0,0);
		*this = newimage;
	}
	
	return !ok;
}
