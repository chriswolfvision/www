/***********************************************************************
 * MRFLabelCounterMLabels.cpp
 * Count the clique labelings in a binary image
 *
 * Author: Christian Wolf
 * Begin: 1.6.2005
 ***********************************************************************/

// C
#include <math.h>
 
// From this Module
#include "MRFLabelCounterMLabels.h"

// **********************************************************************
// Constructor
// **********************************************************************

MRFLabelCounterMLabels::MRFLabelCounterMLabels (unsigned int cx, unsigned int cy)
	: MRFLabelCounter(cx,cy)
{			
}

// **********************************************************************
// Destructor
// **********************************************************************

MRFLabelCounterMLabels::~MRFLabelCounterMLabels()
{	
	for (map<unsigned char *,unsigned int>::iterator iter=cliques.begin();
		 iter!=cliques.end();
		 ++iter)
		delete iter->first;
}

// **********************************************************************
// Returns the count 
// **********************************************************************

unsigned int MRFLabelCounterMLabels::get (unsigned char *lab)
{
	map<unsigned char *,unsigned int>::iterator iter=cliques.find(lab);
	if (iter==cliques.end())		
	{
#ifdef PEDANTIC_CHECK_CODE
		ERR_THROW ("MRFLabelCounterMLabels::getMLabels(): labeling not found!");	
#endif		
		return 0;		
	}
	else
		return iter->second;
}

void MRFLabelCounterMLabels::countInImage (Image &trainimage) 
{	
	unsigned int x_firstth= 0;
	unsigned int x_lastth = trainimage.xsize-cliquesizex;
	unsigned int y_firstth= 0;	
	unsigned int y_lastth = trainimage.ysize-cliquesizey;
	unsigned char *p; 
	map<unsigned char *,unsigned int>::iterator iter;
	
	// Shift a window across the trainimage and cut out the window
	for	(unsigned int y = y_firstth ; y<=y_lastth; y++) 
	for	(unsigned int x = x_firstth ; x<=x_lastth; x++) 
	{
		unsigned int i=0;
		p = new unsigned char [cliquesizex*cliquesizey+1];
		
		// Create the vector containing the clique labeling
		for	(unsigned int wy=0 ; wy<cliquesizey; wy++) 
		for	(unsigned int wx=0 ; wx<cliquesizex; wx++) 
		{    			    			   				
			p[i]=trainimage.get(x+wx, y+wy);		
			++i;
		}		
		p[i]=CLIQUE_LABEL_TERMINATION;
		
		iter=cliques.find(p);
		if (iter==cliques.end())		
			cliques[p]=0;
		else
			++(iter->second);			
	}			
}

// *********************************************************************
// Debug output
// *********************************************************************

ostream & operator << (ostream &os, MRFLabelCounterMLabels &c) 
{
	os << "MRFLabelCounterMLabels:\n";
	for (map<unsigned char *,unsigned int>::iterator iter=c.cliques.begin();
			iter!=c.cliques.end();
			++iter)
	{
		for (unsigned int i=0; i<c.cliquesizex*c.cliquesizey; ++i)
			cerr << (int) (iter->first)[i];
		cerr << " : " << iter->second << endl;	
	}
	return os;
}

