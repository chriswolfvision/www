// C
#include <stdlib.h>
#include <stdio.h>

// From the main library
#include <CIL.h>

// From this library
#include "CComponentList.h"

// ********************************************************************
// Get a Matrix of the size of the original image (the size has to be
// specified), where each element (pixel) is either NULL if it does not
// belong to any cluster or a pointer to the respective component
// ********************************************************************

Matrix<CComponent *> * CComponentList::getPointerMatrix (int xsize, int ysize) 
{
	Matrix<CComponent *> *matrix;
	CPixel *pp;
	
	// Allocate the matrix
	matrix = new Matrix<CComponent *> (ysize, xsize);
	
	// Set all pixel-pointers to NULL
	for (int y=0; y<ysize; ++y)
	for (int x=0; x<xsize; ++x)
		matrix->set(x,y,NULL);
	
	// Travers all components
	ITERATE 
	{	
		pp = (*iter)->pixels;
		
		// Travers the pixels of the component and set them in the matrix
		while (pp!=NULL) 
		{
			matrix->set (pp->x, pp->y, *iter);
			pp = pp->next;
		}				
	}
	
	return matrix;
}

// **********************************************************
// Print the list to stdout
// **********************************************************

void CComponentList::print (FILE *fp) 
{
	fprintf	(stderr, "Connected	components list: %d	components\n", size());

	ITERATE
		(*iter)->print (fp);
}

