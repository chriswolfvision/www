/***********************************************************************
 * MRFLabelCounter2Labels2Labels.cpp
 * Count the clique labelings in a binary image
 *
 * Author: Christian Wolf
 * Begin: 1.6.2005
 ***********************************************************************/

// C
#include <math.h>
 
// From this Module
#include "MRFLabelCounter2Labels.h"

// **********************************************************************
// Constructor
// **********************************************************************

MRFLabelCounter2Labels::MRFLabelCounter2Labels (unsigned int cx, unsigned int cy)
	: MRFLabelCounter(cx,cy)
{
	if (cx*cy>32)
		throw EError ("MRFLabelCounter2Labels::MRFLabelCounter2Labels(): "
			"Only cliques of maximum 32 pixels are supported at this time!");
	cliques.resize((int) pow(2,cx*cy),0);
}

// **********************************************************************
// Destructor
// **********************************************************************

MRFLabelCounter2Labels::~MRFLabelCounter2Labels()
{	
}

void MRFLabelCounter2Labels::countInImage (Image &trainimage) 
{	
	unsigned int x_firstth= 0;
	unsigned int x_lastth = trainimage.xsize-cliquesizex;
	unsigned int y_firstth= 0;	
	unsigned int y_lastth = trainimage.ysize-cliquesizey;	
	
	// Shift a window across the trainimage and cut out the window
	for	(unsigned int y = y_firstth ; y<=y_lastth; y++) 
	for	(unsigned int x = x_firstth ; x<=x_lastth; x++) 
	{
		unsigned int c=0,bitadd = 1;
		
		// Calculate the integer value of the windows content
		for	(unsigned int wy=0 ; wy<cliquesizey; wy++) 
		for	(unsigned int wx=0 ; wx<cliquesizex; wx++) 
		{    			    			   				
			if (trainimage.get(x+wx, y+wy)>0)
				c |= bitadd;
			bitadd = bitadd << 1;
		}	
						
		++cliques[c];
	}			
}

// *********************************************************************
// Debug output
// *********************************************************************

ostream & operator << (ostream &os, MRFLabelCounter2Labels &c) 
{
	os << "MRFLabelCounter2Labels2Classes: \n";
	for (unsigned int i=0; i<c.cliques.size(); ++i) 
		os << i << ": " << c.cliques[i] << endl;
	return os;
}

