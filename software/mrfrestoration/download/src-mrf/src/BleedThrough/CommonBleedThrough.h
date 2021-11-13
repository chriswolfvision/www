/***********************************************************************
 * Common stuff for all BleedThrough observation models
 *
 * Author: Christian Wolf
 * Begin: 24.6.2005
 ***********************************************************************/
 
#ifndef _WOLF_COMMONBLEEDTHROUGH_H_
#define _WOLF_COMMONBLEEDTHROUGH_H_

// C++
#include <set>

// From the IMAGE PROCESSING module
#include <ImageProc.h>

// From the CCA module
#include <CCA.h>

/***********************************************************************
 * PARAMETERS
 ***********************************************************************/

#define BLDTHR_MASK_THR		0.2

/***********************************************************************/
 
// There is an undef at the end of the file
#define NO_CLASSES			3

#define CL_IND_BG			0
#define CL_IND_RECTO		1
#define CL_IND_VERSO		2 
// EVALUATION ONLY: (!!!!!)
#define CL_IND_RECTO_VERSO	3
#define CL_IND_COUNT		4

#define NO_EROSIONS			0

#define KMEANS_INIT_BG		240
#define KMEANS_INIT_RECTO	10
#define KMEANS_INIT_VERSO	80

/***********************************************************************
 * PROPERTIES OF THE OBSERVATION MODEL
 ***********************************************************************/

enum PROPERTIES
{
	PROP_BG_MEAN=0
};	

 /***********************************************************************
 * Get the class index from the two field labels
 ***********************************************************************/
 
inline unsigned int getClass (unsigned char recto, unsigned char verso)
{
	if (recto)
		return CL_IND_RECTO;
	else
		if (verso)
			return CL_IND_VERSO;
		else
			return CL_IND_BG;
} 

 /***********************************************************************
 * Get the class index from the two field labels
 * - all different possibilities (recto does NOT cover verso!)
 ***********************************************************************/
 
inline unsigned int getClassAllComb (unsigned char recto, unsigned char verso)
{
	if (recto)
		if (verso)
			return CL_IND_RECTO_VERSO;
		else
			return CL_IND_RECTO;
	else
		if (verso)
			return CL_IND_VERSO;
		else
			return CL_IND_BG;
} 

 /***********************************************************************
 * Get the class index from a single label field
 * i.e. treat the two fields separately
 ***********************************************************************/
 
inline unsigned int getClassFromSingleField (unsigned char label, bool isRecto)
{
	if (isRecto)
		return (label? CL_IND_RECTO : CL_IND_BG);
	else
		return (label? CL_IND_VERSO : CL_IND_BG);
} 
 
 /***********************************************************************
 * Get the two field labels from the class index
 ***********************************************************************/
 
inline void getLabelsFromClass (unsigned int classLabel, 
	unsigned char &recto, unsigned char &verso)
{
	switch (classLabel)
	{
		case CL_IND_BG:
			recto=1;
			verso=1;
			return;
		case CL_IND_RECTO:
			recto=1;
			verso=0;
			return;
		case CL_IND_VERSO:
			recto=0;
			verso=1;
			return;
#ifdef PEDANTIC_CHECK_CODE
		default:
			ERR_THROW("Internal error in getLabelsFromClass");	
#endif
	}
} 

/***********************************************************************
 * Get the label of the background class
 ***********************************************************************/
 
inline unsigned int getBGLabel(Image *lab)
{
	Vector<unsigned int> count (NO_CLASSES);
	unsigned int max, argMax;
	
	HRULE; cerr << "Determining background label..." << endl;
	
	// Count the number of pixels for each label
	count.setZero();	
	for (int y=0; y<lab->ysize; ++y) 
	for (int x=0; x<lab->xsize; ++x) 	
		++count(lab->get(x,y));
		
	// Determine the label with the most pixels
	argMax=0;
	max=count(0);
	for (int l=1; l<NO_CLASSES; ++l) 
	{
		if (count(l)>max)
		{
			max=count(l);
			argMax=l;
		}
	}
	return argMax;
}	
	
/***********************************************************************
 * Reorder the labels such that the given label is 0
 ***********************************************************************/
 
inline void reOrderLabels(Image *lab, unsigned int labelToBeFirst)
{		
	// nothing to do
	if (labelToBeFirst==0)
		return;
	
	// Redefine the labels such that this label is at position 0
	for (int y=0; y<lab->ysize; ++y) 
	for (int x=0; x<lab->xsize; ++x) 	
	{
		if (lab->get(x,y)==0)
			lab->set(x,y,labelToBeFirst);
		else 
			if (lab->get(x,y)==labelToBeFirst)
				lab->set(x,y,0);
	}
}

/***********************************************************************
 * Determine the recto label
 * Method: use a scaneline based algorithm
 * lineNr ... if != -1, then use only this line
 ***********************************************************************/

#define THR_SIZE				20

#define LAB2INDEX(x)			(x)-1
#define INDEX2LAB(x)			(x)+1
 
inline unsigned char determineRectoLabel (Image *input, int lineNr)
{
	int xs = input->xsize;
	int ys = input->ysize;
	Matrix<int> mId(ys,xs);
	Matrix<int> mCnt(ys,xs);
	unsigned int lastLabel, label;
	set<int> idsL1, idsL2;
	Image *I;
	
	cerr << "Determining recto label..." << endl;
	
	// First heavily filter the image, in order to remove the borders around
	// each character.
	I = new Image(*input);
	filterMedian(*I,PLANE_RED,3);
			
	// Perform a connected components analysis and 
	// back project the results into different matrices
	componentInfo (*I, NULL, NULL, &mId, &mCnt, CCA_NBHD_8, false);
	
	// Traverse the lines
	for (int y=0; y<ys; ++y)
	{
		if (lineNr!=-1 && y!=lineNr)
			continue;
	
		// Traverse the pixels of a line
		lastLabel = I->get(0,y);
		for (int x=1; x<xs; ++x)
		{			
			label = I->get(x,y);		
			
			// the label changes between recto and verso
			if (label!=lastLabel && label!=0 && lastLabel!=0)
			{
				// We only count this label change if both connected components
				// exceed a certain size
				if ((mCnt.get(x,y)>THR_SIZE) && (mCnt.get(x-1,y)>THR_SIZE))
				{
					if (label==1) idsL1.insert (mId(y,x));
					if (label==2) idsL2.insert (mId(y,x));
					if (lastLabel==1) idsL1.insert (mId(y,x-1));
					if (lastLabel==2) idsL2.insert (mId(y,x-1));
				}
			}
			
			lastLabel=label;
		}		
	}	
	
	cerr << "L1: " << idsL1.size() << " components" << endl
		 << "L2: " << idsL2.size() << " components" << endl;
	
	// Clean up
	delete I;
	
	return idsL1.size() < idsL2.size() ? 1 : 2;
}

/***********************************************************************
 * Determine the recto label
 * Using the gray value only
 * We already know that 0 is back ground and that 
 * the two candidate labels are 1 and 2.
 ***********************************************************************/
  
inline unsigned char determineRectoLabelWithGV (Image *lab, Image *obs)
{
	float means[] = {0.,0.};
	unsigned cnt[] = {0,0};
		
	cerr << "Determining recto label (with grayvalue) ..." << endl;
	
	for (int y=0; y<lab->ysize; ++y)
	for (int x=0; x<lab->xsize; ++x)
	{
		if (lab->get(x,y)==1)
		{
			means[0]+=obs->get(x,y);
			++cnt[0];
		}
		if (lab->get(x,y)==2)
		{	
			means[1]+=obs->get(x,y);
			++cnt[1];
 		}
	}	
	means[0]/=(float)cnt[0];
	means[1]/=(float)cnt[1];
	
	cerr << "Recognition of recto label with gray value: " << endl
	     << "L1 has mean " << means[0] << endl
	     << "L2 has mean " << means[1] << endl;
	
	return (means[0]<=means[1]) ? 1 : 2;
}

/***********************************************************************
 * Classify the labels and separate the label field into 2 fields
 * and then update the recto and the verso image such that each of these
 * two label images contains only 2 different labels.
 * 
 * if lineNr!=-1 then determine the recto label using only this line
 ***********************************************************************/
 
inline void recto2RectoVerso(Image *lab, Image *labv, Image *obs, 
	bool useGV4RectoRecognition, int lineNr)
{
	// use the label with the maximum vote count
	unsigned char labelRecto,labelVerso;	  	  
	if (useGV4RectoRecognition)
		labelRecto=determineRectoLabelWithGV(lab, obs);
	else
		labelRecto=determineRectoLabel(lab, lineNr);
		
	labelVerso = (labelRecto == 1 ? 2 : 1);
	
	cerr << "recto label: " << (int) labelRecto << endl
		 << "verso label: " << (int) labelVerso << endl;
	
	// Using the label information, create the two label fields:
	for (int y=0; y<lab->ysize; ++y)
	for (int x=0; x<lab->xsize; ++x)
	{
		unsigned char l=lab->get(x,y);		
		
		if (l==0)
		{
			lab->set (x,y, 0); 	
			labv->set(x,y, 0);
		}
				
		if (l==labelRecto)
		{
			lab->set (x,y, 1); 	
			labv->set(x,y, 0);
		}
		
		if (l==labelVerso)
		{
			lab->set (x,y, 0); 	
			labv->set(x,y, 1);
		}
	}
}

#endif
