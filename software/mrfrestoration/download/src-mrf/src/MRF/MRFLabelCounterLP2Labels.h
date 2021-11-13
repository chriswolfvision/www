/***********************************************************************
 * MRFLabelCounterLP2Labels.h
 * Count the clique labelings in a binary image
 *
 * Author: Christian Wolf
 * Begin: 7.2.2007
 ***********************************************************************/
 
#ifndef _WOLF_MRFLABELCOUNTERLP2LABELS_H_
#define _WOLF_MRFLABELCOUNTERLP2LABELS_H_

// C++
#include <vector>
#include <map>

// From the IMAGE module
#include <Image.h>

// From the main module
#include <CIL.h>

// From this module
#include "MRFLabelCounter.h"

// Contains the clique counts for one given 
// constant Ns labeling
class NsMap
{
	public:
		map<unsigned int,unsigned int> m;
};

class MRFLabelCounterLP2Labels : public MRFLabelCounter
{
	public:		

		// Constructor
		MRFLabelCounterLP2Labels (unsigned int cx, unsigned int cy,	unsigned int xNsMask);
		~MRFLabelCounterLP2Labels();
		
		// Acessors
		unsigned int get (unsigned int i);
		
		// Iterator
		typedef map<unsigned int,NsMap *>::iterator iterator;
		iterator begin() { return clique_sets.begin(); }
		iterator end() { return clique_sets.end(); }
		
		// Increase the count for a clique labeling
		void inc (unsigned int lab);
		
		// Count the cliques labelings in an image
		void countInImage(Image &im, Image &lp);
		void countInImage(Image &im) {};
		
		// Debug output
		friend ostream & operator << (ostream &os, MRFLabelCounterLP2Labels &c);
	
	private: // DATA
			
		map<unsigned int, NsMap *> clique_sets;
		unsigned int NsMask,NsMaskInv;
};

#endif
