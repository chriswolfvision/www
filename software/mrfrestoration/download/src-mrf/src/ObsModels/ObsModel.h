/***********************************************************************
 * An observation model
 *
 * Author: Christian Wolf
 * Begin: 25.5.2005
 ***********************************************************************/
 
#ifndef _WOLF_MRFOBSMODEL_H_
#define _WOLF_MRFOBSMODEL_H_

// From the main module
#include <CIL.h>		
		 
template <class TI> 
class Segmenter;

/***********************************************************************
 * The abstract base class
 ***********************************************************************/

template <class TI>
class ObsModel
{	
	public:
		typedef typename TI::PixelType PixelType;
		
		// Constructor & Destructor
		ObsModel (bool doBlur); 
		virtual ~ObsModel();
		
		// This method is called by the segmenter object 
		// in order to make itself known
		virtual void init(Segmenter<TI> *xseg);
	
		// The _negative_ term in the exponential
		virtual float condEnergy (int x, int y)=0;
		
		// the complete log likelihood
		virtual float logLikelihood (unsigned char obs, unsigned char lab)=0;
		virtual float logLikelihood (unsigned char obs, unsigned char lab,
			unsigned int x, unsigned int y)=0;
		
		
		virtual void estimateParams()=0;
		virtual void updateParams()=0;	
		virtual float getProperty(unsigned int whichProperty)=0;
		
		virtual Image *visualizeSegmentation()=0;
		
		virtual void printParams()=0;
		
		unsigned int getNoClasses() { return segmenter->getNoClasses(); }

	protected:		
	
		Segmenter<TI> *segmenter;	
		
		bool doBlur;
		FloatMatrix *blurred;
};

/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>  
ObsModel<TI>::ObsModel (bool xDoBlur)
{
	segmenter = NULL;
	doBlur = xDoBlur;	
	blurred=NULL;
}

/***********************************************************************
 * Destructor
 ***********************************************************************/

template <class TI>  
ObsModel<TI>::~ObsModel ()
{
	if (doBlur && blurred!=NULL)
		delete blurred;
}

/***********************************************************************
 * This method is called by the segmenter object 
 * in order to make itself known
 ***********************************************************************/

template <class TI>  
void ObsModel<TI>::init(Segmenter<TI> *xsegmenter) 
{ 
	segmenter=xsegmenter; 	
	if (doBlur)
		blurred = new FloatMatrix (segmenter->getLab()->xsize, segmenter->getLab()->ysize);
}



#endif

