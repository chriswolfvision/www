/***********************************************************************
 * MRFLabelCounter2Labels.h
 * Count the clique labelings in a binary image
 *
 * Author: Christian Wolf
 * Begin: 1.6.2005
 ***********************************************************************/
 
#ifndef _WOLF_MRFLABELCOUNTER2LABELS_H_
#define _WOLF_MRFLABELCOUNTER2LABELS_H_

// C++
#include <vector>
#include <map>

// From the IMAGE module
#include <Image.h>

// From the main module
#include <CIL.h>

// From this module
#include "MRFLabelCounter.h"

class MRFLabelCounter2Labels : public MRFLabelCounter
{
	public:

		// Constructor
		MRFLabelCounter2Labels (unsigned int cliquesizex, unsigned int cliquesizey);
		~MRFLabelCounter2Labels();
		
		// Acessors
		unsigned int get (unsigned int i)		{ return cliques[i]; }
			
		// Count the cliques labelings in an image
		void countInImage(Image &im);		
		
		// Debug output
		friend ostream & operator << (ostream &os, MRFLabelCounter2Labels &c);
	
	private: // DATA
		
		vector<unsigned int> cliques;
};

#endif
