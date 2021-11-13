/***************************************************************************
                          ImageFuncs.h  -  description
                             -------------------
    begin                : Tue Mar 26 2002
    copyright            : (C) 2002 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 ***************************************************************************/

#ifndef _WOLF_IMAGEFUNCS_H_
#define _WOLF_IMAGEFUNCS_H_

// From this module
#include "FloatMatrix.h"
#include "Image.h"

enum ImageFileType 
{
	IFT_PPM_P1=1,
	IFT_PPM_P2,
	IFT_PPM_P3,
	IFT_PPM_P4,
	IFT_PPM_P5,
	IFT_PPM_P6,
	IFT_FLOATMATRIX,
	IFT_JPEG,
	IFT_FITS,
	IFT_UNKNOWN
};

ImageFileType checkFileType (const char *filename);

void float2image (FloatMatrix *f, Image *im, bool signedValues, bool alwaysScale);
void float2floatScaled (FloatMatrix *f, FloatMatrix *R);

#endif
