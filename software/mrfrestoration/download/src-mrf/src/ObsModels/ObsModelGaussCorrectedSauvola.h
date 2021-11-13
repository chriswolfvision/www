/***********************************************************************
 * ObsModelGaussCorrectedSauvola.h
 * An observation model:
 * Gaussian noise added to a threshold image obtained with Sauvola et al.'s
 * adaptaive method 
 *
 * Author: Christian Wolf
 * Begin: 31.3.2009
 ***********************************************************************/
 
#ifndef _WOLF_MRFOBSMODELGAUSSCORRECTEDSAUVOLA_H_
#define _WOLF_MRFOBSMODELGAUSSCORRECTEDSAUVOLA_H_

// C++
#include <vector>

// From the main module
#include <CIL.h>		

// From the IMAGE module
#include <Image.h>

// From the BINARIZATION module
#include <Binarization.h>

// From the MRF module
#include <MRFPriorModels/MRFSegmenterPotts.h>

// From this module
#include <Segmenter.h>
#include <ObsModel.h>
		
template <class TI> 
class ObsModelGaussCorrectedSauvola : public ObsModel<TI>
{		
	public:	
		typedef typename TI::PixelType PixelType;
		
		ObsModelGaussCorrectedSauvola ();
		~ObsModelGaussCorrectedSauvola();
		
		void init(Segmenter<TI> *xmrf);
		
		float condEnergy (int x, int y);  
		float logLikelihood (unsigned char obs, unsigned char lab);
		float logLikelihood (unsigned char obs, unsigned char lab, 
			unsigned x, unsigned y);

		void estimateParams();
		void updateParams();	
		void printParams();
		float getProperty(unsigned int whichProperty)	{ ERR_THROW ("getProperty() not implemented"); }
		TI * visualizeSegmentation();

	private:		

		// put a copy of the pointers here for speed reasons			
		TI *obs;
		Image *lab;
		
		// The Sauvola threshold surface
		FloatMatrix *tsurf;
		
		// The only compatible prior model is the Potts model
		MRFSegmenterPotts<TI> *pottsmodel;
	
		vector<float> means, 
					  sigma_squares;
		
		vector<float> add_for_log,muls_in_exp;
};

/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>  
ObsModelGaussCorrectedSauvola<TI>::ObsModelGaussCorrectedSauvola ()
	: ObsModel<TI>(false)
{
	cerr << "Observation model: Gauss-corrected Sauvola" << endl;
}

/***********************************************************************
 * Destructor
 ***********************************************************************/

template <class TI>  
ObsModelGaussCorrectedSauvola<TI>::~ObsModelGaussCorrectedSauvola()
{
}

/***********************************************************************
 * Initialize
 ***********************************************************************/
 
template <class TI>  
void ObsModelGaussCorrectedSauvola<TI>::init (Segmenter<TI> *xsegmenter)
{
	pottsmodel = (MRFSegmenterPotts<TI> *) xsegmenter;		
	ObsModel<TI>::init(xsegmenter);
	obs=xsegmenter->getObs();
	lab=xsegmenter->getLab();
	
	if (this->getNoClasses()>2)
		ERR_THROW ("Observation model 'ObsModelGaussCorrectedSauvola' requires 2 classes");	
	
	means.resize(this->getNoClasses());
	sigma_squares.resize(this->getNoClasses());
	muls_in_exp.resize(this->getNoClasses());
	add_for_log.resize(this->getNoClasses());
}

/***********************************************************************
 * Calculate the energy
 * The version for standard MRFs, i.e. a flat image
 *
 * The energy is related to the log likelihood, 
 * however, ATTENTION! the energy is proportional to the 
 * NEGATIVE (!!!!) of the likelihood!
 ***********************************************************************/
 
template <class TI>  
float ObsModelGaussCorrectedSauvola<TI>::condEnergy (int x, int y)
{	
	ERR_THROW ("ObsModelGaussCorrectedSauvola:condEnergy not yet implemented");
}

/***********************************************************************
 * Calculate the log likelihood
 * directly from a label and an observed value
 ***********************************************************************/
 
template <class TI>  
float ObsModelGaussCorrectedSauvola<TI>::logLikelihood (unsigned char obs, unsigned char lab)
{		
	ERR_THROW ("ObsModelGaussCorrectedSauvola<TI>::logLikelihood() cannot be used: non-stationary model!");
}

template <class TI>  
float ObsModelGaussCorrectedSauvola<TI>::logLikelihood (unsigned char obs, unsigned char lab, 
	unsigned x, unsigned y)
{	
	float rv;	
		
	/*
	if (obs>=tsurf->get(x,y))
		return (lab ? 1 : 1000);
	else
		return (lab ? 1000 : 1);
	*/
	
	if (lab==0)
		rv = (0-obs+tsurf->get(x,y)-127.5);
	else
		rv = (255-obs+tsurf->get(x,y)-127.5);
	
	rv=rv*rv*muls_in_exp[lab];	
	return add_for_log[lab] - rv;	
}

/***********************************************************************
 * estimate the parameters from the current labeling 
 * and the observation:
 * Use the standard ML estimators
 ***********************************************************************/
 
template <class TI>  
void ObsModelGaussCorrectedSauvola<TI>::estimateParams()
{
	unsigned int xs=this->obs->xsize;
	unsigned int ys=this->obs->ysize;
	vector<unsigned int> noPixels(this->getNoClasses());
	unsigned int l;
	float d;
// 	FloatMatrix *dummy1=NULL, *dummy2=NULL;
	
	cerr << "Estimating observation model parameters (Gauss-corrected Sauvola)" << endl;
// 		 << "Calculating threshold surface." << endl;
	
// 	tsurf = surfaceNiblackImproved (*(this->obs), NIBLACK_WOLF2, 40, 40, 0.5, 128, dummy1, dummy2);
	
// 	{
// 		Image tmp (*(this->obs));
// 		thresholdNiblackImproved (tmp, NIBLACK_WOLF2, 40, 40, 0.5, 128);
// 		tmp.write ("x_sauvola.pgm");
// 	}

	tsurf = pottsmodel->getThresholdSurface();
	
// 	for (unsigned int y=0; y<tsurf->ysize; ++y)
// 	for (unsigned int x=0; x<tsurf->xsize; ++x)
// 	{
// 		cerr << "[" << tsurf->get(x,y) << "]"; cerr.flush();
// 	}

	
	// init
	for (unsigned int i=0; i<this->getNoClasses(); ++i)
	{
		means[i]=0;
		sigma_squares[i]=0;
		noPixels[i]=0;
	}
	
	// Calculate the means
	// Used for the calculation of the variance only,
	// the means will be replace with 0 and 255 
	// AFTERWARDS!!!!!!!
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
// 		l=((this->obs->get(x,y) >= tsurf->get(x,y))?1:0);
// 		this->lab->set(x,y,l);
		l=this->lab->get(x,y);
		d = this->obs->get(x,y) - means[l];
		if (l>=this->getNoClasses())
			ERR_THROW ("Internal error in ObsModelGauss1VCorr<TI>::estimateParams():"
				" invalid label: " << l);
		means[l] += this->obs->get(x,y);
		++noPixels[l];
	}		
	
	for (unsigned int i=0; i<this->getNoClasses(); ++i)
		means[i] /= (float) noPixels[i];
	
	// Calculate a unique variance but over multiple
	// means
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{		
		l=this->lab->get(x,y);
		d = this->obs->get(x,y) - means[l];
		
		// we used the first variance (index 0)
		// for the calculation of the unique
		// variance
		sigma_squares[0] += d*d;
	}
	sigma_squares[0] /= (xs*ys);
	sigma_squares[1] = sigma_squares[0];
		
	// The means are fixed to 0 and 255 here!
  	means[0] = 0;
  	means[1] = 255;
		
	// DEBUG
	this->printParams();
				 	
	updateParams();
}

/***********************************************************************
 * Update the parameters, 
 * i.e. calculate the precomputed values which accelerate the 
 * calculations
 ***********************************************************************/
 
template <class TI>  
void ObsModelGaussCorrectedSauvola<TI>::updateParams()
{
	for (unsigned int i=0; i<this->getNoClasses(); ++i)
	{
		if (sigma_squares[i]==0)
			ERR_THROW ("ObsModelGaussCorrectedSauvola: noise variance is zero!");
						
		muls_in_exp[i]=1.0/(2.0*sigma_squares[i]);
		add_for_log[i] = 1.0/(sqrtf(2.*M_PI*sigma_squares[i]));		
		if (add_for_log[i]<=0)
			ERR_THROW ("Internal error in "
			"ObsModelGaussCorrectedSauvola::updateParams: "
			"invalid model parameters");
		add_for_log[i] = logf(add_for_log[i]);
	}
}

/***********************************************************************
 * Visualize the result by creating an image where each pixel
 * contains the gray value of its label
 ***********************************************************************/
 
template <class TI>  
TI * ObsModelGaussCorrectedSauvola<TI>::visualizeSegmentation()
{ 
	ERR_THROW ("ObsModelGaussCorrectedSauvola<TI>::visualizeSegmentation():"
		"the general template has not yet been implemented.\n");
}
 
template <>
Image * ObsModelGaussCorrectedSauvola<Image>::visualizeSegmentation()
{	
	Image *rv = new Image (this->obs->xsize, this->obs->ysize, 
		this->obs->nbColorPlanes());
		
	for (int y=0; y<this->obs->ysize; ++y)
	for (int x=0; x<this->obs->xsize; ++x)
	{
		unsigned int l = this->segmenter->getLab()->get(x,y);
		rv->set(PLANE_RED, x, y, (byte) TRIMMEDGRAY(rint(means[l])));
	}
	return rv;
}

/***********************************************************************
 * Print the parameters to stderr
 ***********************************************************************/
 
template <class TI>  
void ObsModelGaussCorrectedSauvola<TI>::printParams()
{
	for (unsigned int i=0; i<this->getNoClasses(); ++i)
		cerr << "class " << i << ": mean=" << means[i]
			 << " sigma=" << sqrt(sigma_squares[i]) << endl;	
}

#endif
