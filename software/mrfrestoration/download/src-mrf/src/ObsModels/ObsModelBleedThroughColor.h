/***********************************************************************
 * An observation model:
 * the model for bleed through removal
 * COLOR PROCESSING!!!!!
 *
 * Author: Christian Wolf
 * Begin: 20.6.2005
 ***********************************************************************/
 
#ifndef _WOLF_MRFOBSMODELBLEEDTHROUGHCOLOR_H_
#define _WOLF_MRFOBSMODELBLEEDTHROUGHCOLOR_H_

// From the main module
#include <CIL.h>		

// From the MRF module
#include <ObsModel.h>

// From this module
#include "CommonBleedThrough.h"
#include "MRFSegmenterBleedThrough.h"

#define NO_COLOR_PLANES			3
		
template <class TI> 
class ObsModelBleedThroughColor : public ObsModel<TI>
{		
	public:	
		typedef typename TI::PixelType PixelType;
		
		// Constructor and destructor
		ObsModelBleedThroughColor(bool doBlur);
		~ObsModelBleedThroughColor();
		
		void init(MRFSegmenter<TI> *xmrf);
		
		float condEnergy (int x, int y);
		float logLikelihood (unsigned char obs, unsigned char lab);
		float logLikelihood (unsigned char obs, unsigned char lab, 
			unsigned x, unsigned y);
		void estimateParams();
		void updateParams();
		void beforeIteration();
		float getProperty(unsigned int whichProperty) { ERR_THROW ("getProperty() not implemented"); }
		void printParams();
		
		TI * visualizeSegmentation();

	private:		

		// put a copy of the pointers here for speed reasons			
		TI *obs;
		Image *lab, *labv;
	
		Vector<float> *means[NO_COLOR_PLANES];	
		Matrix<float> *covariances[NO_COLOR_PLANES];
		
		// Help vectors; Kept here instead of in the methods
		// for speed reasons only.
		Vector<float> *foo, *bar;
};

/***********************************************************************
 * Constructor
 *  
 ***********************************************************************/

template <class TI>  
ObsModelBleedThroughColor<TI>::ObsModelBleedThroughColor (bool doBlur)
	:ObsModel<TI>(doBlur)
{ 
	if (doBlur)
		ERR_THROW("Blurring not yet supported for ObsModelBleedThroughColor!");
		
 	for (unsigned int i=0; i<NO_CLASSES; ++i)
	{
		means[i]        = new Vector<float> (NO_COLOR_PLANES);
		covariances[i] = new Matrix<float> (NO_COLOR_PLANES, NO_COLOR_PLANES);
	}
	foo = new Vector<float> (NO_COLOR_PLANES);
	bar = new Vector<float> (NO_COLOR_PLANES);
}

/***********************************************************************
 * Destructor
 ***********************************************************************/

template <class TI>  
ObsModelBleedThroughColor<TI>::~ObsModelBleedThroughColor()
{
	for (unsigned int i=0; i<NO_CLASSES; ++i)
	{
		delete means[i];
		delete covariances[i];
	}
	delete foo;
	delete bar;
}

/***********************************************************************
 * Initialize
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThroughColor<TI>::init (MRFSegmenter<TI> *xmrf)
{
	cerr << "Observation model: bleedthrough Gaussian color.\n";
	
	// In our heart we know, that the segmenter actually is 
	// of type MRFSegmenterBleedThrough, since only this one
	// is compatible to this model.
	// So lets help out the compiler a little bit ...
	MRFSegmenterBleedThrough<TI> *s = (MRFSegmenterBleedThrough<TI> *) xmrf;

	ObsModel<TI>::init(xmrf);
	obs =s->getObs();
	lab =s->getLab();
	labv=s->getLabv();	
	
	if (obs->nbColorPlanes()<3)
		throw EError ("Can't do color processing on a grayscale image!\n");
}

/***********************************************************************
 * This method is called before each iteration during
 * the estimation process
 ***********************************************************************/

template <class TI>  
void ObsModelBleedThroughColor<TI>::beforeIteration () 
{
	cerr << "\nNo blurring included!!!\n";
}

/***********************************************************************
 * Calculate the likelihood
 ***********************************************************************/
 
template <class TI>  
float ObsModelBleedThroughColor<TI>::condEnergy (int x, int y)
{
	unsigned char lr = this->lab->get (x,y);
	unsigned char lv = this->labv->get(x,y);
	int l=getClass(lr,lv);

	// The difference vector
	for (unsigned int i=0; i<NO_COLOR_PLANES; ++i)
		(*foo)(i) = this->obs->get(i+1,x,y); 
	*foo -= (*means[l]);
	
	// multiply with the inverse of the covariance matrix
	for (unsigned int i=0; i<NO_COLOR_PLANES; ++i)
	{
		float sum=0;
		for (unsigned int j=0; j<NO_COLOR_PLANES; ++j)
			sum += (*foo)(j) * ((*covariances[l])(j,i));
		(*bar)(i)=sum;
	}
		
	// And then dot-multiply it with the original difference
	return bar->dotProduct(*foo);
}

/***********************************************************************
 * Calculate the log likelihood
 * directly from a label and an observed value
 ***********************************************************************/
 
template <class TI>  
float ObsModelBleedThroughColor<TI>::logLikelihood (unsigned char obs, unsigned char lab)
{		
	ERR_THROW ("logLikelihood(): not implemented");
	return 0;		
}

template <class TI>  
float ObsModelBleedThroughColor<TI>::logLikelihood (unsigned char obs, unsigned char lab,
	unsigned int x, unsigned int y)
{		
	ERR_THROW ("logLikelihood(): not implemented");
	return 0;		
}

/***********************************************************************
 * Estimate the parameters from the current labeling 
 * and the observation:
 * Use the standard ML estimators
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThroughColor<TI>::estimateParams()
{
	unsigned int xs=this->obs->xsize;
	unsigned int ys=this->obs->ysize;
	unsigned int noPixels[NO_CLASSES];
	
	for (unsigned int i=0; i<NO_CLASSES; ++i)
	{
		means[i]->setZero();
		covariances[i]->setZero();
		noPixels[i]=0;
	}
	
	// Calculate the mean vectors
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		unsigned char lr = this->lab->get (x,y);
		unsigned char lv = this->labv->get(x,y);
		int l=getClass(lr,lv);		
		for (unsigned int i=0; i<NO_COLOR_PLANES; ++i)
			(*means[l])(i) += this->obs->get(i+1,x,y); 
		++noPixels[l];
	}		
	for (unsigned int i=0; i<NO_CLASSES; ++i)
		(*means[i]) /= (float) noPixels[i];
	
	// Calculate the covariance matrices
	// and invert them
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{		
		unsigned char lr = this->lab->get (x,y);
		unsigned char lv = this->labv->get(x,y);
		int l=getClass(lr,lv);
		
		// The difference vector
		for (unsigned int i=0; i<NO_COLOR_PLANES; ++i)
			(*foo)(i) = (*means[l])(i) - this->obs->get(i+1,x,y); 
		
		for (unsigned int i=0; i<NO_COLOR_PLANES; ++i)
		for (unsigned int j=0; j<NO_COLOR_PLANES; ++j)
			(*covariances[l])(i,j) += (*foo)(i)*(*foo)(j);		
	}
	for (unsigned int i=0; i<NO_CLASSES; ++i)
	{
		(*covariances[i]) /= (float) noPixels[i];
		covariances[i]->invert();		
	}
	
	// Debug output	
	for (unsigned int i=0; i<NO_CLASSES; ++i)
		cerr << "class " << i << ": " << *(means[i]);
	for (unsigned int i=0; i<NO_CLASSES; ++i)
		cerr << "class " << i << ": " << *covariances[i];
}

/***********************************************************************
 * Update the parameters, 
 * i.e. calculate the precomputed values which accelerate the 
 * calculations
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThroughColor<TI>::updateParams()
{
}

/***********************************************************************
 * Visualize the result by creating an image where each pixel
 * contains the color of its label
 ***********************************************************************/
 
template <class TI>  
TI * ObsModelBleedThroughColor<TI>::visualizeSegmentation()
{ 
	ERR_THROW ("ObsModelBleedThroughColor<TI>::visualizeSegmentation:"
		"the general template has not yet been implemented.\n");
}

/***********************************************************************
 * Print the parameters to stderr
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThroughColor<TI>::printParams()
{
	ERR_THROW ("Not yet implemented.");	
}

#endif
