/***************************************************************************
                          ImageFuncs.cc  -  description
                             -------------------
    begin                : Tue Mar 26 2002
    copyright            : (C) 2002 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 ***************************************************************************/

// C
#include <math.h>
#include <errno.h>
#include <string.h>

// C++
#include <iostream>

// If FITS PROCESSING is enabled ...
#ifdef HAVE_LIBCFITSIO 	  
#include <fitsio.h>
#endif

// From this module
#include "Image.h"
#include "FloatMatrix.h"
#include "MultiBandFloatMatrix.h"
#include "ImageFuncs.h"

// ***************************************************************
// Checks the image type of a given file and returns it
// ***************************************************************

ImageFileType checkFileType (const char *filename)
{
	char shortbuf[256];
	FILE * inpic;

	if ((inpic	= fopen(filename,"r")) == NULL)	
		ERR_THROW ("Cannot open file " << filename << ": " << strerror(errno) << endl);
	
	// Read the type flag
	fgets(shortbuf,250,inpic);
	fclose (inpic);

	// ----------------------------------------------------------------------------
	// Check if it's a NETPBM image
	
	if (shortbuf[0] == 'P') 
	{
		if (shortbuf[1]>='1' && shortbuf[1]<='6')
			return (ImageFileType) (IFT_PPM_P1 + shortbuf[1] - '1');
		else
			return IFT_UNKNOWN;
	}
	
	// ----------------------------------------------------------------------------
	// Check whether it is a float matrix		
	
	if (strncmp(shortbuf,MAGIC_NUMBER_FI,strlen(MAGIC_NUMBER_FI))==0) 
		return IFT_FLOATMATRIX;
		
	// ----------------------------------------------------------------------------
	// Check if it's a FITS image
#ifdef HAVE_LIBCFITSIO 	
	{	fitsfile *fp;
		int status=0;
		 	
		if (!fits_open_file(&fp, filename, READONLY, &status))
			return IFT_FITS;
	}
#endif	

	// ----------------------------------------------------------------------------
	// Check if it's a JPEG image
#ifdef HAVE_LIBJPEG				
	if (Image::check4JPEG (filename))
		return IFT_JPEG;
#endif			
	

	
	return IFT_UNKNOWN;
}

// ***************************************************************
// Create a 8bit gray scale image from a FLoatMatrix
// alwaysScale ... If true, then we always scale the image, even
// 				    if the maximum is smaller then 255.
// ***************************************************************

void float2image (FloatMatrix *f, Image *im, bool signedValues, bool alwaysScale) {
	int xs=f->xsize;
	int ys=f->ysize;
	float max=f->max();
	float min=f->min();
	float factor;
	byte **R=im->R;
	
	// ----------------------------------------------
	// Signed: 0=gray127, neg=black, positive=white
	// ----------------------------------------------
	if (signedValues) {
		float ma = fabs(max);
		float mi = fabs(min);
		factor = 255.0 / (ma > mi ? ma : mi);
		
		for (int y=0; y<ys; ++y) {
    		for (int x=0; x<xs; ++x) {
    			R[x][y] = (byte) rint(127.5 + factor*(f->get(x,y)));
    		}
    	}	
	}

	// ----------------------------------------------
	// Unsigned: 0=black, positive=white
	// ----------------------------------------------
	
	else {
    	factor=255.0/(max-min);
    	
    	// Do not scale, only truncate and copy the values
    	if (!alwaysScale && max<=255) {
        	for (int y=0; y<ys; ++y) {
        		for (int x=0; x<xs; ++x) {
        			R[x][y] = (byte) rint(f->get(x,y));		
        		}
        	}	
    	}
		
    	// Scale
    	else  {
        	for (int y=0; y<ys; ++y) {
        		for (int x=0; x<xs; ++x) {
        			float v = f->get(x,y);
        			v = rint((v-min)*factor);		
        			R[x][y] = (byte) v;
        		}
        	}
        }
	}
}

void float2floatScaled (FloatMatrix *f, FloatMatrix *R) {
	int xs=f->xsize;
	int ys=f->ysize;
	float max=f->max();
	float min=f->min();
	float factor=255.0/(max-min);
	
	// cerr << "Floatimage to 8bit image: maximum = " << max << endl;
	
	for (int y=0; y<ys; ++y) {
		for (int x=0; x<xs; ++x) {
			R->set(x,y, (f->get(x,y)-min)*factor);		
		}
	}
}

