/***********************************************************************
 * An observation model:
 * the model for bleed through removal
 *
 * Author: Christian Wolf
 * Begin: 7.7.2006
 ***********************************************************************/
 
#ifndef _WOLF_MRFOBSMODELBLEEDTHROUGHNONPARAMETRIC_H_
#define _WOLF_MRFOBSMODELBLEEDTHROUGHNONPARAMETRIC_H_


/***********************************************************************
 * PARAMETERS */
 
#define NO_FILTER_ITS		15

/***********************************************************************/

// From the main module
#include <CIL.h>		

// From the IMAGE PROCESSING module
#include <ImageProc.h>

// From the VISUALIZATION module
#include <Visualize.h>

// From the MRF module
#include <ObsModel.h>

// From this module
#include "CommonBleedThrough.h"
#include "MRFSegmenterBleedThrough.h"
		
template <class TI> 
class ObsModelBleedThroughNonParametric : public ObsModel<TI>
{		
	public:	
		typedef typename TI::PixelType PixelType;
		
		// Constructor & destructor
		ObsModelBleedThroughNonParametric (bool doBlur);
		~ObsModelBleedThroughNonParametric();
		
		void init(MRFSegmenter<TI> *xmrf);
		
		float condEnergy (int x, int y);
		float logLikelihood (unsigned char obs, unsigned char lab);
		float logLikelihood (unsigned char obs, unsigned char lab, 
			unsigned x, unsigned y);
			
		void estimateParams();		
		void updateParams()						{}
		float getProperty(unsigned int whichProperty)	{ ERR_THROW ("getProperty() not implemented"); }
		void printParams();
		TI * visualizeSegmentation();

	private:		

		// put a copy of the pointers here for speed reasons			
		TI *obs;
		Image *lab, *labv;
	
		Histogram<float> *hist[NO_CLASSES];
};

/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>  
ObsModelBleedThroughNonParametric<TI>::ObsModelBleedThroughNonParametric (bool doBlur)
	:ObsModel<TI>(doBlur)
{
	if (doBlur)
		ERR_THROW("Blurring not yet supported for ObsModelBleedThroughColor!"); 
}

/***********************************************************************
 * Destructor
 ***********************************************************************/

template <class TI>  
ObsModelBleedThroughNonParametric<TI>::~ObsModelBleedThroughNonParametric ()
{	
	for (int i=0; i<NO_CLASSES; ++i)
		delete hist[i];
}

/***********************************************************************
 * Initialize
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThroughNonParametric<TI>::init (MRFSegmenter<TI> *xmrf)
{
	cerr << "Observation model: bleedthrough non parametric gray.\n";

	// In our heart we know, that the segmenter actually is 
	// of type MRFSegmenterBleedThrough, since only this one
	// is compatible with this model.
	// So lets help out the compiler a little bit ...
	MRFSegmenterBleedThrough<TI> *s = (MRFSegmenterBleedThrough<TI> *) xmrf;

	ObsModel<TI>::init(xmrf);
	obs =s->getObs();
	lab =s->getLab();
	labv=s->getLabv();
	
	for (int i=0; i<NO_CLASSES; ++i)
		hist[i] = new Histogram<float> (256, 0, 255, true);
}

/***********************************************************************
 * Calculate the likelihood
 ***********************************************************************/
 
template <class TI>  
float ObsModelBleedThroughNonParametric<TI>::condEnergy (int x, int y)
{
	unsigned char lr = this->lab->get (x,y),
				  lv = this->labv->get(x,y);
				  
//#warning CONSTANT COND ENERGY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!				  
	//return hist[0]->fromValue(this->obs->get(x,y)); 
				 
	return hist[getClass(lr,lv)]->fromValue(this->obs->get(x,y));
}

/***********************************************************************
 * Calculate the log likelihood
 * directly from a label and an observed value
 ***********************************************************************/
 
template <class TI>  
float ObsModelBleedThroughNonParametric<TI>::logLikelihood (unsigned char obs, unsigned char lab)
{		
	ERR_THROW ("logLikelihood(): not implemented");
	return 0;		
}

template <class TI>  
float ObsModelBleedThroughNonParametric<TI>::logLikelihood (unsigned char obs, unsigned char lab,
	unsigned int x, unsigned int y)
{		
	ERR_THROW ("logLikelihood(): not implemented");
	return 0;		
}

/***********************************************************************
 * Estimate the parameters from the current labeling 
 * and the observation:
 ***********************************************************************/

template <class TI>  
void ObsModelBleedThroughNonParametric<TI>::estimateParams()
{	
#ifdef HAVE_VISUALIZE 
#ifdef WRITE_DEBUG_IMAGES
	static int funcCallCount=0;	
	++funcCallCount;
#endif
#endif
	DataSerie gaussFilter;
	gaussFilter.add (1);
	gaussFilter.add (4);
	gaussFilter.add (6);
	gaussFilter.add (4);
	gaussFilter.add (1);
	gaussFilter.normalize();

	for (int i=0; i<NO_CLASSES; ++i)
		hist[i]->clear();
		
	for (int y=0; y<this->lab->ysize; ++y)
	for (int x=0; x<this->lab->xsize; ++x)
	{
		unsigned char lr = this->lab->get (x,y),
					  lv = this->labv->get(x,y);					
		int label=getClass(lr,lv);				
		if (label<0 || label>=NO_CLASSES)
			ERR_THROW ("Internal error in ObsModelBleedThroughNonParametric<>::estimateParams(1)");
		hist[label]->add(this->obs->get(x,y));
	}

	for (int i=0; i<NO_CLASSES; ++i)
	{
		for (int f=0; f<NO_FILTER_ITS; ++f)
			hist[i]->filter(gaussFilter);		
		hist[i]->normalize();
	}
	
#ifdef HAVE_VISUALIZE 
#ifdef WRITE_DEBUG_IMAGES
	{
		ostringstream s; 
		s << "likelihoodplot_nonparam_" << setfill ('0') << setw(3) << funcCallCount << ".txt"; 
		ofstream st (s.str().c_str(), ios::out);
		if (!st.good())
			ERR_THROW ("Cannot open file " << s.str() << "for writing!\n");
	
		// Debug: plot the class histograms as well as the fitted distributions
		for (unsigned int c=0; c<NO_CLASSES; ++c)
		{				
			st << *hist[c];
			st << endl;		
		}
		st.close();	
	}
#endif	
#endif
}

/***********************************************************************
 * Visualize the result by creating an image where each pixel
 * contains the gray value of its label
 * Flip the right (verso) part of the image
 ***********************************************************************/
 
template <class TI>  
TI * ObsModelBleedThroughNonParametric<TI>::visualizeSegmentation()
{ 
	ERR_THROW ("ObsModelBleedThroughNonParametric<TI>::visualizeSegmentation:"
		"the general template has not yet been implemented.\n");
}

template <>
Image * ObsModelBleedThroughNonParametric<Image>::visualizeSegmentation()
{	
	Image *rv = new Image (2*this->obs->xsize, this->obs->ysize, 
		this->obs->nbColorPlanes());
	unsigned int label;
	float means[NO_CLASSES];
	
	for (int i=0; i<NO_CLASSES; ++i)
		means[i]=hist[i]->getMean();
		
	for (int y=0; y<this->obs->ysize; ++y)
	for (int x=0; x<this->obs->xsize; ++x)
	{	
		label = (this->lab->get(x,y)  ? CL_IND_RECTO : CL_IND_BG);
		rv->set(PLANE_RED, x, y, (byte) TRIMMEDGRAY(rint(means[label])));
		label = (this->labv->get(x,y) ? CL_IND_RECTO : CL_IND_BG);
		rv->set(PLANE_RED, rv->xsize-x-1, y, (byte) TRIMMEDGRAY(rint(means[label])));
	}
	return rv;
}

/***********************************************************************
 * Print the parameters to stderr
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThroughNonParametric<TI>::printParams()
{
	ERR_THROW ("Not yet implemented.");	
}

#endif
