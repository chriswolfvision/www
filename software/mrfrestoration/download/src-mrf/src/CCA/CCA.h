/***************************************************************************
                          CCA.h  -  description
                             -------------------
    begin                : Wed Oct 24 2001
    copyright            : (C) 2001 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 ***************************************************************************/

#ifndef _WOLF_CCA_H_
#define _WOLF_CCA_H_

// From the IMAGE module
#include <Image.h>

// From the MATHEMATICS module
#include <Matrix.h>

// From this module
#include "CPixel.h"
#include "CComponentList.h"

enum CCANeighborhoodType
{
	CCA_NBHD_4=0,
	CCA_NBHD_8
};

// Run a connected components analysis and return the list of components
CComponentList * connectedComponents (Image &im, CCANeighborhoodType nhType, bool ignoreBG);

void componentInfo (Image im, Matrix<int> *height, Matrix<int> *top, Matrix<int> *compid,
	Matrix<int> *cnt, CCANeighborhoodType nhType, bool ignoreBG);

// Get approximatvie information on components but do _NOT_ run a real CCA
void componentColumnHeight (Image &im, Matrix<int> &dst);
void componentLineWidth (Image &im, Matrix<int> &dst);

#endif
