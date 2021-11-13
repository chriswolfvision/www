/***********************************************************************
 * MRFLabelCounterMLabelsMLabels.h
 * Count the clique labelings in a binary image
 *
 * Author: Christian Wolf
 * Begin: 1.6.2005
 ***********************************************************************/
 
#ifndef _WOLF_MRFLABELCOUNTERMLABELS_H_
#define _WOLF_MRFLABELCOUNTERMLABELS_H_

// C++
#include <vector>
#include <map>

// From the IMAGE module
#include <Image.h>

// From the main module
#include <CIL.h>

// From this module
#include "MRFLabelCounter.h"

#define CLIQUE_LABEL_TERMINATION	255

struct lesstthancl
{
	// Since the clique labeling vectors are terminated.
  	bool operator()(const unsigned char* c1, const unsigned char* c2) const
  	{
  		unsigned int i=0;
  		while (true)
  		{
  			if (c1[i]==CLIQUE_LABEL_TERMINATION)
  				return false;
  			else
  				if (c1[i]!=c2[i])
  					return c1[i]<c2[i];
  				else
  					++i;
  		}
  	}
};

class MRFLabelCounterMLabels : public MRFLabelCounter
{
	public:

		// Constructor
		MRFLabelCounterMLabels (unsigned int cliquesizex, unsigned int cliquesizey);
		~MRFLabelCounterMLabels();
		
		// Acessors
		unsigned int get (unsigned char *lab);
		
		// Count the cliques labelings in an image
		void countInImage(Image &im);
		// Debug output
		friend ostream & operator << (ostream &os, MRFLabelCounterMLabels &c);
	
	private: // DATA
	
		map<unsigned char *,unsigned int,lesstthancl> cliques;
};

#endif
