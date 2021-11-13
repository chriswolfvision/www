/***********************************************************************
 * Document Bleed-Through Restoration (separation recto/verso)
 * Using two Markov Random Fields.
 * Includes a separate line process for edges.
 *
 * Author: Christian Wolf
 * Begin: 30.1.2007
 ***********************************************************************/
 
#warning LINE PROCESSES CURRENTLY NOT SUPPORTED
#if 0
 
// From this module
#include "MRFSegmenterBleedThroughLP.h" 

// static class members. They need to be instantiated explicitely for the given
// type (e.g. <Image>) since the compiler does not know the instantiation
// which happens in another .cpp file.

// Different decisions taken on the state of different LP neighbors
// See research notebook 30.1.2007, p. 133

template <>
unsigned int MRFSegmenterBleedThroughLP<Image>::VirtualLabelKeepsBondRDConn[] = {0, 1, 0, 1};

template <>
unsigned int MRFSegmenterBleedThroughLP<Image>::VirtualLabelKeepsBondLDConn[] = {0, 0, 1, 1};

template <>
unsigned int MRFSegmenterBleedThroughLP<Image>::KeepsBondRDFromNeighbors[] = 
	{ 1, 1, 1, 1,
	  1, 0, 0, 0,
	  1, 0, 0, 0,
	  1, 0, 0, 0};
	  
template <>	 
unsigned int MRFSegmenterBleedThroughLP<Image>::KeepsBondLDFromNeighbors[] = 
	{ 1, 1, 1, 0,
	  1, 1, 0, 0,
	  1, 0, 1, 0,
	  0, 0, 0, 0};

template <>
unsigned int MRFSegmenterBleedThroughLP<Image>::BlockingCodeFromNeighbors[] = 
	{ 0, 0, 0, 2,
	  0, 1, 3, 3,
	  0, 3, 1, 3,
	  2, 3, 3, 3};
template <>	  
unsigned int MRFSegmenterBleedThroughLP<Image>::ParameterIndexFromNeighbors[] = 
	{ 0, 3, 3, 2,
	  3, 2, 1, 2,
	  3, 1, 2, 2,
	  2, 2, 2, 3};
	  

#endif
