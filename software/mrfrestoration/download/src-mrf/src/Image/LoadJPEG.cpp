/***************************************************************************
                          LoadJPEG.cpp  -  description
                             -------------------
    begin                : Thu Apr 10 2003
    copyright            : (C) 2003 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 ***************************************************************************/

#ifdef HAVE_LIBJPEG

// C
#include <stdio.h>

// C++
#include <iostream>

// From the MAIN module of the library
#include <CIL.h>

// The JPEG library
#include <jpeglib.h>

// This module
#include "Image.h"

using namespace std;

/************************************************************************
 * Check if the file contains a jpeg image
 ************************************************************************/

bool Image::check4JPEG (const char *filename) 
{ 
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	FILE *infile;
	int i=-1;
	
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);

	if ((infile = fopen(filename, "rb")) == NULL) 
		ERR_THROW ("can't open '" << filename << "' for reading!\n");
		
	jpeg_stdio_src(&cinfo, infile);
	
	i = jpeg_read_header(&cinfo, TRUE);
	return (i==JPEG_HEADER_OK);
}

/************************************************************************
 * load a JPEG image
 ************************************************************************/

bool Image::readJPEG (const char *filename) 
{ 
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	FILE *infile;
	JSAMPARRAY buffer;		
	int row_stride, y;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);

	if ((infile = fopen(filename, "rb")) == NULL) 
		ERR_THROW ("can't open '" << filename << "' for reading!\n");
	
	jpeg_stdio_src(&cinfo, infile);
	if (jpeg_read_header(&cinfo, TRUE)!= JPEG_HEADER_OK ) 
		ERR_THROW ("invalid jpeg header\n");         		

	jpeg_start_decompress(&cinfo);

	// Allocate the image
	xsize = cinfo.output_width;
	ysize = cinfo.output_height;
	type =  cinfo.output_components;
	_alloc (xsize, ysize, type);        

	// Allocate the scanline buffer
	row_stride = cinfo.output_width * cinfo.output_components;
	buffer = (*cinfo.mem->alloc_sarray) ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

	y=0;
	while (cinfo.output_scanline < cinfo.output_height) 
	{
		jpeg_read_scanlines(&cinfo, buffer, 1);
		// Store the scanline in the image

		switch (type) 
		{

			// color image
			case 3:
				for (int x=0; x<xsize; ++x) 
				{
					set(PLANE_RED,x,y,   buffer[0][3*x+0]);
					set(PLANE_GREEN,x,y, buffer[0][3*x+1]);
					set(PLANE_BLUE,x,y,  buffer[0][3*x+2]);
				}
				break;

			case 1:
			case 2:
				for (int x=0; x<xsize; ++x) 
				{
					set(PLANE_RED,x,y, buffer[0][x]);
				}
				break;

			default:
				jpeg_destroy_decompress (&cinfo);
				ERR_THROW ("Internal error in Image::readJPEG()\n");
		}
		++y;
	}

	jpeg_finish_decompress(&cinfo);
	fclose (infile);
	jpeg_destroy_decompress(&cinfo);

	return true;
}

#endif
