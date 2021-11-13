/**********************************************************
 * BTKMeansClusterer.h
 * A special version of the K-Means Clusterer for
 * Bleed Through Removal
 *
 * Christian Wolf
 * Beginn: 22.12.2005
 **********************************************************/

#ifndef	_WOLF_BT_KMEANS_H_
#define	_WOLF_BT_KMEANS_H_

// from the MATHEMATICS module
#include "KMeansClusterer.h"

// The two possible labels for the bleed through mask
enum {
	BTMASK_BG=0,
	BTMASK_FG
};

// ******************************************************
// Class definition
// T... the datatype storing a scalar, not a vector!!!!
// ******************************************************

template <class T>
class BTKMeansClusterer : public KMeansClusterer<T>
{
	public:
		
		typedef Vector<T> VectorType;
		typedef Vector<float> VectorTypeFloat;	
	
		// Constructor
		BTKMeansClusterer (unsigned int d, unsigned int nC)	:
			KMeansClusterer<T>(d,nC) {}						
			
		void add (VectorTypeFloat v);
	
		// the [] operator (inherited) returns the normal labels
		// the () operator returns the mask labels!!!
		unsigned int operator() (unsigned int i) const 		{ return labels_mask[i]; }

		void  calcMaskLabels(unsigned int bgLabel, float threshold);
	
	private:
		
		// These are not the "normal" classified labels,
		// but the labels for the classification into 
		// "background" and "rest"
		vector<unsigned int> labels_mask;
};

// ******************************************************
// Add data
// ******************************************************

template <class T>
inline void BTKMeansClusterer<T>::add (VectorTypeFloat v)
{
	KMeansClusterer<T>::add(v);
	labels_mask.push_back(0);
}

// ******************************************************
// Create a matrix designed to mask the background pixels
// which are too clearly background pixels, i.e. we do 
// not want to mess around trying to relax them with 
// MRF algorithms.
// The criterion is ratio of the distances to the 2 nearest 
// cluster centers (the nearest is the background cluster,
// of course).
// ******************************************************

template <class T>
void BTKMeansClusterer<T>::calcMaskLabels(unsigned int bgLabel, float threshold)
{
	unsigned int no=this->data.size();
	
	// Traverse the data
	for (unsigned int i=0; i<no; ++i)
	{	
		float dist, minDist[] = {0.,0.};
		unsigned int minIndices[] = {0,0};
		unsigned int minsAvailable;
		
		
		// the label is background
		// -> let's see if we satisfy the additional constraint
		if (this->labels[i]==bgLabel)
		{
			
			// For each point determine the nearest cluster center
			// and the second nearest cluster center
			minsAvailable=0;
			for (unsigned int c=0; c<this->noClasses; ++c)
			{			
				dist = this->data[i].euclideanDistanceSquared(this->centers[c]);				
				if (minsAvailable==0 || (minsAvailable>0 && dist < minDist[0]))
				{
					minDist[1]=minDist[0];
					minDist[0]=dist;
					minIndices[1]=minIndices[0];
					minIndices[0]=c;
					++minsAvailable;
				} else
				{
					if (minsAvailable==1 || (minsAvailable>1 && dist < minDist[1]))
					{
						minDist[1]=dist;
						minIndices[1]=c;
						++minsAvailable;
					}
				}				
			}
#ifdef CHECK_CODE						
			if (minIndices[0]!=bgLabel || minsAvailable<2)
				ERR_THROW ("Internal Error in BTKMeansClusterer::calcMaskLabels()\n"
					<< "minsAvailable=" << minsAvailable << endl
					<< "minDist[0]=" << minDist[0] << endl
					<< "minDist[1]=" << minDist[1] << endl
					<< "minIndices[0]=" << minIndices[0] << endl
					<< "minIndices[1]=" << minIndices[1] << endl);
#endif	
			if (minDist[0]/minDist[1] < threshold)			
				labels_mask[i]=BTMASK_BG;
			else
				labels_mask[i]=BTMASK_FG;
				
			// cout << "RATIO " << minDist[0]/minDist[1] << endl;
		}
		
		// the label is NOT background
		// --> the mask pixel is NOT background
		else
			labels_mask[i]=BTMASK_FG;
	}
}

#endif
