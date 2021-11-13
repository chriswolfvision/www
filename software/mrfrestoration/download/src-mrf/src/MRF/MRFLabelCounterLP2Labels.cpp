/***********************************************************************
 * MRFLabelCounterLP2Labels.cpp
 * Count the clique labelings in a binary image
 *
 * Author: Christian Wolf
 * Begin: 7.2.2007
 ***********************************************************************/

// C
#include <math.h>

// From this Module
#include "MRFLabelCounterLP2Labels.h"
#include "MRFLineProcesses.h"

// **********************************************************************
// Constructor
// NsMask is the mask which extracts the neighbor Ns from a 
// complete clique labelling
// **********************************************************************

MRFLabelCounterLP2Labels::MRFLabelCounterLP2Labels (unsigned int cx, unsigned int cy,
	unsigned int xNsMask)
	:MRFLabelCounter(cx,cy)
{
	if (cx*cy>32)
		throw EError ("MRFLabelCounterLP2Labels::MRFLabelCounterLP2Labels(): "
			"Only cliques of maximum 32 pixels are supported at this time!");
	NsMask = xNsMask;
	NsMaskInv = ~NsMask;
}

// **********************************************************************
// Destructor
// **********************************************************************

MRFLabelCounterLP2Labels::~MRFLabelCounterLP2Labels()
{	
	map<unsigned int, NsMap *>::iterator iter=clique_sets.begin();
	while (iter!=clique_sets.end())
	{
		delete iter->second;
		++iter;
	}
}

// **********************************************************************
// Returns the count (the case of multiple labels per pixel)
// **********************************************************************

unsigned int MRFLabelCounterLP2Labels::get (unsigned int lab)
{
	unsigned int Ns=lab&NsMask;
	map<unsigned int, NsMap *>::iterator iter=clique_sets.find(Ns);
	if (iter==clique_sets.end())
		return 0;
		
	// There is a set with this neighbor hood. 
	// but is this labeling in the set?
	else
	{
		NsMap *nsmap=iter->second;
		map<unsigned int,unsigned int>::iterator iter2=nsmap->m.find(lab);
		if (iter2==nsmap->m.end())				
			return 0;		
		else
			return iter2->second;
	}
}

// **********************************************************************
// Increase the count for a given clique labeling
// **********************************************************************

void MRFLabelCounterLP2Labels::inc (unsigned int lab)
{
	unsigned int Ns=lab&NsMask;
	NsMap *nsmap;
	map<unsigned int, NsMap *>::iterator iter=clique_sets.find(Ns);
	
	// The clique set does not exist, create it
	if (iter==clique_sets.end())
	{
		nsmap = new NsMap;
		nsmap->m[lab]=1;
		clique_sets[Ns] = nsmap;
	}
		
	// The clique set already exists, search if the labeling exists
	else
	{	
		nsmap = iter->second;
		// Increment the count for this clique
		// or add it if it does not exist yet
		map<unsigned int, unsigned int>::iterator iter2=nsmap->m.find(lab);
		if (iter2==nsmap->m.end())		
			nsmap->m[lab]=1;
		else
			++(iter2->second);	
	}
}

// **********************************************************************
// Counts the clique labelings in a single training image and
// associated line process
// See research notebook 7.2.2007, p. 
// **********************************************************************

void MRFLabelCounterLP2Labels::countInImage(Image &ip, Image &lp) 
{	
	unsigned int x_firstth= 0;
	unsigned int x_lastth = ip.xsize-cliquesizex;
	unsigned int y_firstth= 0;	
	unsigned int y_lastth = ip.ysize-cliquesizey;
	unsigned int lpcliquesizex = 2*cliquesizex-1;
	unsigned int lpcliquesizey = 2*cliquesizey-1;
	map<unsigned int, unsigned int>::iterator iter;
	
	// Shift a window across the trainimage and cut out the window
	for	(unsigned int y = y_firstth ; y<=y_lastth; y++) 
	for	(unsigned int x = x_firstth ; x<=x_lastth; x++) 
	{
		unsigned int c=0,bitadd = 1;
		unsigned int lpx=2*x;
		unsigned int lpy=2*y;
		
		// Calculate the integer value of the window content
		// --> Intensity process
		for	(unsigned int wy=0 ; wy<cliquesizey; wy++) 
		for	(unsigned int wx=0 ; wx<cliquesizex; wx++) 
		{    			    			   				
			if (ip.get(x+wx, y+wy)>0)
				c |= bitadd;
			bitadd = bitadd << 1;
		}	
		
		// Calculate the integer value of the windows content
		// --> Line process
		for	(unsigned int wy=0 ; wy<lpcliquesizey; wy++) 
		for	(unsigned int wx=0 ; wx<lpcliquesizex; wx++) 
		{    
			// cerr << "{" << lpx+wx << "," << lpy+wy << "}";		
			if (IS_ELEMENT_OF_MRFLP(lpx+wx,lpy+wy))
			{				
				if (lp.get(lpx+wx,lpy+wy)>0)
					c |= bitadd;
				bitadd = bitadd << 1;
			}
		}
		
		inc(c);		
	}			
}

// *********************************************************************
// Debug output
// *********************************************************************

ostream & operator << (ostream &os, MRFLabelCounterLP2Labels &c) 
{
	os << "MRFLabelCounterLP2Labels: \n";
	for (map<unsigned int,NsMap *>::iterator iter=c.clique_sets.begin();
		 iter!=c.clique_sets.end();
		 ++iter)
	{
		NsMap *nsmap=iter->second;
		cerr << "Set for Neighborhood " << iter->first 
			 << "(" << nsmap->m.size() << " entries)" << endl;
		for (map<unsigned int,unsigned int>::iterator iter2=nsmap->m.begin();
				iter2!=nsmap->m.end();
				++iter2)
		{
			cerr << iter2->first << " : " << iter2->second << " times." << endl;	
		}		
	}
	return os;
}

