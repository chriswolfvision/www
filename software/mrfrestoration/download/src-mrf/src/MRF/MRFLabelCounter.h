/***********************************************************************
 * MRFLabelCounter.h
 * Count the clique labelings in a binary image
 *
 * Author: Christian Wolf
 * Begin: 1.6.2005
 ***********************************************************************/
 
#ifndef _WOLF_MRFLABELCOUNTER_H_
#define _WOLF_MRFLABELCOUNTER_H_

// C++
#include <vector>
#include <map>

// From the IMAGE module
#include <Image.h>

// From the main module
#include <CIL.h>

class MRFLabelCounter
{
	public:

		// Constructor
		MRFLabelCounter (unsigned int cliquesizex, unsigned int cliquesizey);
		virtual ~MRFLabelCounter() {}
		virtual void countInImage(Image &im)=0;
			
	protected: // DATA
	
		unsigned int cliquesizex, cliquesizey;
};

// **********************************************************************
// Constructor
// **********************************************************************

inline MRFLabelCounter::MRFLabelCounter (unsigned int cx, unsigned int cy)
{
	if (cx*cy>32)
		throw EError ((char *) "MRFLabelCounter::MRFLabelCounter(): "
			"Only cliques of maximum 32 pixels are supported at this time!");
			
	cliquesizex=cx;
	cliquesizey=cy;
}

#endif
