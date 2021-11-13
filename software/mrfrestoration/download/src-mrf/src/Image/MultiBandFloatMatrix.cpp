/**********************************************************
 * MultiBandFloatMatrix.cpp
 *
 * Christian Wolf, chriswolf@gmx.at
 * Beginn: 22.3.2005
 **********************************************************/
 
// If FITS PROCESSING is enabled ...
#ifdef HAVE_LIBCFITSIO 	  
#include <fitsio.h>
#endif

// From the main module
#include <CIL.h>
 
// From this module
#include "MultiBandFloatMatrix.h"

/************************************************************************
 * Constructor - create an uninitialized image
 ************************************************************************/
 
MultiBandFloatMatrix::MultiBandFloatMatrix ()
{
	noBands = 0;
	xsize = -1;
	ysize = -1;
}

/************************************************************************
 * Constructor - create an uninitialized image
 ************************************************************************/
 
MultiBandFloatMatrix::MultiBandFloatMatrix (int xno, int xxs, int xys)
{
	noBands = xno;
	xsize = xxs;
	ysize = xys;
	
	matrices.resize(xno);
	for (int i=0; i<noBands; ++i)
		matrices[i] = new FloatMatrix (xsize, ysize);	
}

/************************************************************************
 * Constructor - read from a FITS file
 ************************************************************************/
 
#ifdef HAVE_LIBCFITSIO 	  
MultiBandFloatMatrix::MultiBandFloatMatrix (const char *filename)
{	
	FloatMatrix *f;
	fitsfile *fp;
    int status=0, hdupos, hdutype, foo;
	long dims[2];
	int index;
	
	if (fits_open_file(&fp, filename, READONLY, &status))
		throw EError ("Cannot open FITS file!\n");	
		
	// Traverse the images (HDUs)
	index=0;
	do {		
		float *scanline;
		long coords[2];
		
		// Get the current HDU position and type
		fits_get_hdu_num(fp, &hdupos);  
		fits_get_hdu_type(fp, &hdutype, &status);
		
		if (hdutype!=IMAGE_HDU)   
			throw EError ("FITS tables are not yet supported!");			
			
		// Check the value type
		fits_get_img_type(fp, &foo, &status);
		switch (foo)
		{
			case BYTE_IMG:
			case SHORT_IMG:
			case LONG_IMG:
			case FLOAT_IMG:
				break;
			default:
				ERR_THROW ("FITS: unsupported image type: " << foo << "!\n");	
		}
	
		// Check the number of dimensions
		fits_get_img_dim(fp, &foo, &status);
		if (foo!=2)
			throw EError ("FITS: only 2D images are supported!\n");
				
		// Get the image size
		fits_get_img_size(fp, 2, dims, &status);		
		
		// Allocate the Image and a scanline buffer
		f = new FloatMatrix (dims[0], dims[1]);	
		scanline = new float[dims[0]];
			
		coords[0]=1;
		for (coords[1]=1; coords[1]<=dims[1]; coords[1]++)
		{
			fits_read_pix(fp, TFLOAT, coords, 
				dims[0], NULL, scanline, NULL, &status);
			
			for (int x=0; x<dims[0]; ++x)	
				f->set(x, coords[1]-1, scanline[x]);
		}
	
		// Add the image plane to the stack	
		matrices.push_back(f);
		delete scanline;
		
		// Try to move to the next hdu
		fits_movrel_hdu(fp, 1, NULL, &status);  
		
	
	} while (status==0);
	
	// Close the file
	fits_close_file(fp, &status);
	 
	noBands=matrices.size();
	xsize=matrices[0]->xsize;
	ysize=matrices[0]->ysize;
	
	cerr << "FITS: " <<xsize << "," << ysize << "  " << noBands << " band(s)." << endl;
}
#endif // #ifdef HAVE_LIBCFITSIO

/************************************************************************
 * save into a FITS file
 ************************************************************************/
 
#ifdef HAVE_LIBCFITSIO 	  
void MultiBandFloatMatrix::saveFITS (const char *filename)
{
	fitsfile *fp;
    int status=0;

	if (fits_create_file(&fp, filename, &status))
		throw EError ("Cannot create FITS file!\n");	
	
	// Traverse all bands	
	ITERATE
	{
		long dims[2], coords[2];
		float *scanline;
		FloatMatrix *f=*iter;
		
		// Create the HDU
		dims[0] = f->xsize;
		dims[1] = f->xsize;		
		if (fits_create_img (fp, FLOAT_IMG, 2, dims, &status))
			throw EError ("Cannot create FITS HDU!\n");
					
		// Write the image band
		coords[0]=1;
		scanline = new float[dims[0]];
		for (coords[1]=1; coords[1]<=dims[1]; coords[1]++)
		{
			for (int x=0; x<dims[0]; ++x)	
				scanline[x] = f->get(x, coords[1]-1);
		
			fits_write_pix(fp, TFLOAT, coords, 
				dims[0], scanline, &status);							
		}
		delete scanline;
	}	
	
	// Close the file
	fits_close_file(fp, &status);
}

#endif // #ifdef HAVE_LIBCFITSIO

/************************************************************************
 * Destructor
 ************************************************************************/
 
MultiBandFloatMatrix::~MultiBandFloatMatrix ()
{
	for (int i=0; i<noBands; ++i)
		if (matrices[i]!=NULL)
			delete matrices[i];
}

