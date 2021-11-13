/***********************************************************************
 * An observation model:
 * the model for bleed through removal
 *
 * Author: Christian Wolf
 * Begin: 10.6.2005
 ***********************************************************************/
 
#ifndef _WOLF_MRFOBSMODELBLEEDTHROUGH_H_
#define _WOLF_MRFOBSMODELBLEEDTHROUGH_H_

// From the main module
#include <CIL.h>		

// From the MATHEMATICS module
//#include <EMClustererGauss1V.h>

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
class ObsModelBleedThrough : public ObsModel<TI>
{		
	public:	
		typedef typename TI::PixelType PixelType;
		
		ObsModelBleedThrough (bool doBlur);
		~ObsModelBleedThrough();
		
		void init(Segmenter<TI> *xmrf);
		
		float condEnergy (int x, int y);
		float logLikelihood (unsigned char obs, unsigned char lab);
		float logLikelihood (unsigned char obs, unsigned char lab, 
			unsigned x, unsigned y);
		void estimateParams();
		void estimateParamsMLorEM(bool doEM);
		void updateParams();
		
		float getProperty(unsigned int whichProperty);
		void printParams();
		
		TI * visualizeSegmentation();

	private:		

		// put a copy of the pointers here for speed reasons			
		TI *obs;
		Image *lab, *labv;
		Image *mask;
			
		float means[NO_CLASSES];	
		float sigma_squares[NO_CLASSES];
		
		// For speedup
		float muls_in_exp[NO_CLASSES],
			  add_for_log[NO_CLASSES];	
};


/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>  
ObsModelBleedThrough<TI>::ObsModelBleedThrough (bool doBlur)
	:ObsModel<TI>(doBlur)
{
}

/***********************************************************************
 * Destructor
 ***********************************************************************/

template <class TI>  
ObsModelBleedThrough<TI>::~ObsModelBleedThrough ()
{
}

/***********************************************************************
 * Initialize
 ***********************************************************************/

template <class TI>  
float ObsModelBleedThrough<TI>::getProperty(unsigned int whichProperty)
{
	if (whichProperty==PROP_BG_MEAN)
	 	return means[CL_IND_BG];
	
	ERR_THROW ("ObsModelBleedThrough::getProperty(): unknown property"
		<< whichProperty);
}

/***********************************************************************
 * Initialize
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThrough<TI>::init (Segmenter<TI> *xmrf)
{
	cerr << "Observation model: bleedthrough Gaussian gray.\n";

	// In our heart we know, that the segmenter actually is 
	// of type MRFSegmenterBleedThrough, since only this one
	// is compatible to this model.
	// So lets help out the compiler a little bit ...
	MRFSegmenterBleedThrough<TI> *s = (MRFSegmenterBleedThrough<TI> *) xmrf;

	ObsModel<TI>::init(xmrf);
	obs =s->getObs();
	lab =s->getLab();
	labv=s->getLabv();
	mask=s->getMask();
}

/***********************************************************************
 * Calculate the log likelihood
 * directly from a label and an observed value
 ***********************************************************************/
 
template <class TI>  
float ObsModelBleedThrough<TI>::logLikelihood (unsigned char o, unsigned char l,
	unsigned int x, unsigned int y)
{
	return logLikelihood(o,l);		
}
 
template <class TI>  
float ObsModelBleedThrough<TI>::logLikelihood (unsigned char o, unsigned char l)
{		
#ifdef PEDANTIC_CHECK_CODE
	if (l>=this->getNoClasses)
		ERR_THROW ("Internal error in ObsModelBleedThrough<TI>::logLikelihood()!");
#endif
	float rv=means[l] - (float) o;
	rv*=(rv*muls_in_exp[l]);	
	return add_for_log[l] - rv;		
}

/***********************************************************************
 * Calculate the energy potential
 * *
 * The energy is related to the log likelihood, 
 * however, ATTENTION! the energy is proportional to the 
 * NEGATIVE (!!!!) of the likelihood!
 ***********************************************************************/
 
#define M(i,j)	means[getClass(this->lab->get (i,j),this->labv->get(i,j))]

template <class TI>  
float ObsModelBleedThrough<TI>::condEnergy (int x, int y)
{		
	int l=getClass(this->lab->get (x,y),this->labv->get(x,y));
	float lab2obs, rv;
	
	if (this->doBlur)
	{
		if (x>0 && y>0 && x<this->lab->xsize-1 && y<this->lab->ysize-1)
		{					
			lab2obs = 0.0625*M(x-1,y-1) + 0.1250*M(x,y-1) + 0.0625*M(x+1,y-1) +
					  0.1250*M(x-1,y)   + 0.2500*M(x,y)   + 0.1250*M(x+1,y) + 
					  0.0625*M(x-1,y+1) + 0.1250*M(x,y+1) + 0.0625*M(x+1,y+1);
					  		
			// lab2obs = 0.5*M(x,y) + 0.5*M(x-1,y);			
			// lab2obs = 1.0*M(x,y);		
			DEBUGIFPX(x,y)
			{
				cerr << DEBUGSTRPX << " blurredlab=" << lab2obs 
					 << " ; means[" << l << "]=" << M(x,y);
			}
		}
		else		
			lab2obs = means[l];
	}
	else	
		lab2obs = means[l];
	
	rv=lab2obs - (float) this->obs->get(x,y),
	rv=rv*rv*muls_in_exp[l];
			
	/*
	DEBUGIFPX(x,y)
	{
		cerr << "cond(" <<  (int) lr << "-" << (int) lv << ":" << l 
			<< "; v=" << (int) this->obs->get(x,y) 
			<< "; mean=" << means[l]
			<< "; result=" << rv << ")\n";
	}	      
	*/
	      
	//The energy is the NEGATIVE (!!!!) of the likelihood!
	return rv - add_for_log[l];			
}
#undef M

/***********************************************************************
 * Estimate the parameters from the current labeling 
 * and the observation:
 * Use the standard ML estimators.
 * if doEM is true, then
 * we interpret the ensemble of the three classes as a 
 * Gaussian Mixture Model (GMM) and estimate its parameters
 * using the Expectation-Maximization Algorithm (EM).
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThrough<TI>::estimateParamsMLorEM(bool doEM)
{
	unsigned int xs=this->obs->xsize;
	unsigned int ys=this->obs->ysize;	
	vector<Image *> erIms(NO_CLASSES);
#ifdef HAVE_VISUALIZE 
#ifdef WRITE_DEBUG_IMAGES
	static int funcCallCount=0;	
	++funcCallCount;
#endif
#endif

	// Create a segmented+eroded image for each class 
	for (unsigned int c=0; c<NO_CLASSES; ++c)
	{					
		// Create the image
		erIms[c] = new Image (xs, ys, 1);
		for (unsigned int y=0; y<ys; ++y)
		for (unsigned int x=0; x<xs; ++x)
			erIms[c]->set(x,y, (
				getClass(this->lab->get (x,y),this->labv->get(x,y))==c ? 1 : 0));
				
		// Perform the erosions if requested
		if (NO_EROSIONS>0)
			erode(*erIms[c], NO_EROSIONS, false, 1);
		
		// Take the mask into account if it exists
		if (mask!=NULL)
		{
			for (unsigned int y=0; y<ys; ++y)
			for (unsigned int x=0; x<xs; ++x)
				if (this->mask->get(x,y)==0)
					erIms[c]->set(x,y,0);
		}	
	}
	
	// Then fit a gaussian to each class
	for (unsigned int c=0; c<NO_CLASSES; ++c)
	{									
		fitGaussian (*(this->obs), *erIms[c], 1, means[c], sigma_squares[c]);
		cerr << "class " << c << " (G): mu=" << means[c]
			 << "\ts2=" << sigma_squares[c] << endl;
	}		
					
#ifdef HAVE_VISUALIZE 
#ifdef WRITE_DEBUG_IMAGES
	{
		ostringstream s; 
		s << "likelihoodplot_ml_" << setfill ('0') << setw(3) << funcCallCount << ".txt"; 
		ofstream st (s.str().c_str(), ios::out);
		if (!st.good())
			ERR_THROW ("Cannot open file " << s.str() << "for writing!\n");
	
		// Debug: plot the class histograms as well as the fitted distributions
		for (unsigned int c=0; c<NO_CLASSES; ++c)
		{							 								
			plotClassDistributions (this->obs, erIms[c], st, 0, 0, 1.);
			st << endl;		
			plotGaussian (st, means[c], sigma_squares[c], 0, 255, 1, 1, false);
			st << endl;		
		}
		st.close();	
	}
	{
		ostringstream s; 
		s << "likelihoodplot_ml_all_" << setfill ('0') << setw(3) << funcCallCount << ".txt"; 
		ofstream st (s.str().c_str(), ios::out);
		if (!st.good())
			ERR_THROW ("Cannot open file " << s.str() << "for writing!\n");
	
		plotClassDistributions (this->obs, NULL, st, 0, 0, 1.);		
		st.close();	
	}
#endif	
#endif	

    
    {
    	/*
    	#warning MU AND SIGMA SET MANUALLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	means[0] = 196;
    	means[1] = 51;
    	means[2] = 97;    	
    	sigma_squares[0] = sigma_squares[1] = sigma_squares[2] = 100;
    	*/
    }

 
	// If requested, then use the ML estimators as an 
	// initialization to the EM algorithm
	if (doEM)
	{
		ERR_THROW ("EMClusterer deactivated for the moment (compilation problems only)\n");
/*		
		vector<int> noPixels(NO_CLASSES);
		int noPixTot=0;
		EMClustererGauss1V<float> clusterer(*(this->obs), NO_CLASSES);
		
		// Count the pixels for each class
		for (unsigned int y=0; y<ys; ++y)
		for (unsigned int x=0; x<xs; ++x)
			++noPixels[getClass(this->lab->get (x,y),this->labv->get(x,y))];	
		for (unsigned int j=0; j<NO_CLASSES; ++j)
			noPixTot+=noPixels[j];
		
		// Initialize the EM components with the results from
		// the K-Means clustering and the ML estimation	
		for (unsigned int j=0; j<NO_CLASSES; ++j)
			clusterer.initComponent(j, (float)noPixels[j]/(float)noPixTot, 
				means[j], sigma_squares[j]);
			
		// Perform the clustering. 
		// The result will be applied to the observation image
		clusterer.doCluster(30);	
								
#ifdef HAVE_VISUALIZE 
#ifdef WRITE_DEBUG_IMAGES
		{
			ostringstream s; 
			s << "likelihoodplot_em_end_" << setfill ('0') << setw(3) << funcCallCount << ".txt"; 
			ofstream st (s.str().c_str(), ios::out);
			if (!st.good())
				ERR_THROW ("Cannot open file " << s.str() << "for writing!\n");
		
			// Debug: plot the class histograms as well as the fitted distributions
			for (unsigned int c=0; c<NO_CLASSES; ++c)
			{							 								
				plotClassDistributions (this->obs, erIms[c], st, 0, 0,
					noPixels[c]/(float)noPixTot);
				st << endl;		
				plotGaussian (st, means[c], sigma_squares[c], 0, 255, 1, 
					noPixels[c]/(float)noPixTot, false);
				st << endl;		
			}
			st.close();	
		}		
#endif	
#endif	
*/
	}
						
	// Clean up
	for (unsigned int c=0; c<NO_CLASSES; ++c)
		delete erIms[c];		
}

/***********************************************************************
 * Estimate the parameters from the current labeling 
 * and the observation:
 ***********************************************************************/

template <class TI>  
void ObsModelBleedThrough<TI>::estimateParams()
{	
	estimateParamsMLorEM(false);
	updateParams();	
}

/***********************************************************************
 * Update the parameters, 
 * i.e. calculate the precomputed values which accelerate the 
 * calculations
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThrough<TI>::updateParams()
{
	for (unsigned int i=0; i<NO_CLASSES; ++i)
	{
		if (sigma_squares[i]==0)
			ERR_THROW ("ObsModelGauss1V: noise variance is zero!");
			
		muls_in_exp[i]=1/(2.0*sigma_squares[i]);
		
		add_for_log[i] = 1.0/(sqrtf(2.*M_PI*sigma_squares[i]));		
		if (add_for_log[i]<=0)
			ERR_THROW ("Internal error in "
			"ObsModelBleedThrough::updateParams: "
			"invalid model parameters");
		add_for_log[i] = logf(add_for_log[i]);
					
		// DEBUG	
		cerr << "class " << i << ": mean=" << means[i]
			 << " sigma=" << sqrt(sigma_squares[i])
			 << " muls_in_exp=" << muls_in_exp[i]
			 << " add_for_log=" << add_for_log[i] << endl;
	}
}

/***********************************************************************
 * Visualize the result by creating an image where each pixel
 * contains the gray value of its label
 * Flip the right (verso) part of the image
 ***********************************************************************/
 
template <class TI>  
TI * ObsModelBleedThrough<TI>::visualizeSegmentation()
{ 
	ERR_THROW ("ObsModelBleedThrough<TI>::visualizeSegmentation:"
		"the general template has not yet been implemented.\n");
}

template <>
Image * ObsModelBleedThrough<Image>::visualizeSegmentation()
{	
	Image *rv = new Image (2*this->obs->xsize, this->obs->ysize, 
		this->obs->nbColorPlanes());
	unsigned int label;
		
	for (int y=0; y<this->obs->ysize; ++y)
	for (int x=0; x<this->obs->xsize; ++x)
	{	
		label = (this->lab->get(x,y)  ? CL_IND_RECTO : CL_IND_BG);
		rv->set(PLANE_RED, x, y, (byte) TRIMMEDGRAY(rint(means[label])));
		label = (this->labv->get(x,y) ? CL_IND_RECTO : CL_IND_BG);
		
		// Flip the verso image
		rv->set(PLANE_RED, rv->xsize-x-1, y, (byte) TRIMMEDGRAY(rint(means[label])));
		// rv->set(PLANE_RED, this->obs->xsize+x, y, (byte) TRIMMEDGRAY(rint(means[label])));
	}
	return rv;
}

/***********************************************************************
 * Print the parameters to stderr
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThrough<TI>::printParams()
{
	for (unsigned int i=0; i<NO_CLASSES; ++i)
		cerr << "class " << i << ": mean=" << means[i]
			 << " sigma=" << sqrt(sigma_squares[i]) << endl;	
}

#endif
