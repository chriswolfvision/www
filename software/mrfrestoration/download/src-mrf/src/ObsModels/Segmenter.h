/***********************************************************************
 * The generic segmenter class
 *
 * Author: Christian Wolf
 * Begin: 12.3.2007
 ***********************************************************************/
 
#ifndef _WOLF_SEGMENTER_H_
#define _WOLF_SEGMENTER_H_

// From the main module
#include <CIL.h>

// From the IMAGE module
#include <Image.h>


template <class TI>
class Segmenter 
{		
	public:
		TI *getObs() 				{ return obs; }
		Image *getLab()				{ return lab; }
		unsigned int getNoClasses() { return this->noClasses; }

	protected: 
	
		TI *obs;			// The observation field
		Image *lab; 		// The label field
		
		unsigned int noClasses;
};

#endif


