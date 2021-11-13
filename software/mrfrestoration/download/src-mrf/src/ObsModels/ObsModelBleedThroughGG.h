/***********************************************************************
 * An observation model:
 * the model for bleed through removal
 * Generalized Gaussian, uni-variate
 *
 * Author: Christian Wolf
 * Begin: 24.6.2005
 ***********************************************************************/
 
#ifndef _WOLF_MRFOBSMODELBLEEDTHROUGHGG_H_
#define _WOLF_MRFOBSMODELBLEEDTHROUGHGG_H_

// From the main module
#include <CIL.h>		

// From the NUMERICAL RECIPES module (includes THIRD PARTY CODE!!!)
#ifdef HAVE_NUM_RECIPES
#include <NumericalRecipes.h>
#endif

// From the MATHEMATICS module
#include <MathFuncs.h>

// From the IMAGE PROCESSING module
#include <ImageProc.h>

// From the VISUALIZATION module
#ifdef HAVE_VISUALIZE
#include <Visualize.h>
#endif

// From the MRF model
#include <ObsModel.h>

// From this module
#include "CommonBleedThrough.h"
#include "MRFSegmenterBleedThrough.h"
		
template <class TI> 
class ObsModelBleedThroughGG : public ObsModel<TI>
{		
	public:	
		typedef typename TI::PixelType PixelType;
		
		ObsModelBleedThroughGG (bool doBlur);
		void init(MRFSegmenter<TI> *xmrf);
		
		float condEnergy (int x, int y);
		float logLikelihood (unsigned char obs, unsigned char lab);
		float logLikelihood (unsigned char obs, unsigned char lab, 
			unsigned x, unsigned y);
		void estimateParams();
		void updateParams();
		float getProperty(unsigned int whichProperty)	{ ERR_THROW ("getProperty() not implemented");
		 }
		void printParams();
		
		TI *visualizeSegmentation();

	private:		

		// put a copy of the pointers here for speed reasons			
		TI *obs;
		Image *lab, *labv;
	
		float shapes[NO_CLASSES];
		float means[NO_CLASSES];	
		float sigma_squares[NO_CLASSES];
		
		// For speedup
		float muls_in_exp[NO_CLASSES];	
};

/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>  
ObsModelBleedThroughGG<TI>::ObsModelBleedThroughGG (bool doBlur) 
	: ObsModel<TI>(doBlur) 
{
	if (doBlur)
		ERR_THROW ("Blurring not yet supported by ObsModelBleedThroughGG!");
}

/***********************************************************************
 * Initialize
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThroughGG<TI>::init (MRFSegmenter<TI> *xmrf)
{
	cerr << "Observation model: bleed through generalized Gaussian.\n";
	// In our heart we know, that the segmenter actually is 
	// of type MRFSegmenterBleedThrough, since only this one
	// is compatible to this model.
	// So lets help out the compiler a little bit ...
	MRFSegmenterBleedThrough<TI> *s = (MRFSegmenterBleedThrough<TI> *) xmrf;

	ObsModel<TI>::init(xmrf);
	obs =s->getObs();
	lab =s->getLab();
	labv=s->getLabv();
}

/***********************************************************************
 * Calculate the conditional energy
 * the formula is 
 * en = ((gamma(3/p)/(sigma^2*gamma(1/p)))^1/2 * (z-mu))^p
 * but since the numerical recipes library provides us with
 * ln(gamma(x)) instead of gamma(x), its better to rewrite the above 
 * equation as
 * en = ((exp(lngamma(3/p)-2*lngamma(sigma)-lngamma(1/p)))^1/2 * (z-mu))^p
 ***********************************************************************/
 
template <class TI>  
float ObsModelBleedThroughGG<TI>::condEnergy (int x, int y)
{
	unsigned char lr = this->lab->get (x,y);
	unsigned char lv = this->labv->get(x,y);
	int l=getClass(lr,lv);
	float f,p=shapes[l];		
	f = expf(numrec::gammln(3./p)-numrec::gammln(sigma_squares[l])-numrec::gammln(1./p));
	f = (sqrt(f) * fabs ((float) this->obs->get(x,y) - means[l]));
	f = pow(f,p);
    ERR_THROW("ObsModelBleedThroughGG<TI>::condEnergy(): No lognormal implemented!");
	return f;	
}

/***********************************************************************
 * Calculate the log likelihood
 * directly from a label and an observed value
 ***********************************************************************/
 
template <class TI>  
float ObsModelBleedThroughGG<TI>::logLikelihood (unsigned char obs, unsigned char lab)
{		
	ERR_THROW ("logLikelihood(): not implemented");
	return 0;		
}

template <class TI>  
float ObsModelBleedThroughGG<TI>::logLikelihood (unsigned char obs, unsigned char lab,
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
 
 /*
class ShapeFunctional
{
	public:
		ShapeFunctional (float a)				{ alpha=a; }
		float operator() (float x) 
		{
			return exp(gammln(5./x)+gammln(1./x)-2.*gammln(5./x))-alpha;	
		}
	private:
		float alpha;
};  
 
template <class T>
void fitGeneralizedGaussianForMean (T &im, Image &segI, unsigned char label, 
	float inmean, float &outvar, float &outshape)
{
	DataSerie s;
	float m4;
	dataSerieFromImage (im,segI,label,s);	
	
	// Get the necessary moments
	outvar  =     s.getVariance();	
	
	// Calculate the shape parameters
	m4=0.;
	for (DataSerie::iterator iter = s.begin(); iter!=s.end(); ++iter)
		m4 += pow(inmean-*iter,4);
	m4 /= s.size();			
	
	cerr << "fitGeneralizedGaussianWithMean(): m=" << inmean << " v=" << outvar << " m4=" << m4 << endl;
	
	try 
	{
		outshape = rootWithBisection (ShapeFunctional(m4/(outvar*outvar)), 0.1, 10, 0.001);
	}
	catch (EError &e)
	{
		ERR_THROW ("fitGeneralizedGaussian(): fit is numerically unstable!!");		
	} 						
}
*/
 
 
template <class TI>  
void ObsModelBleedThroughGG<TI>::estimateParams()
{
	unsigned int xs=this->obs->xsize;
	unsigned int ys=this->obs->ysize;
	vector<Image *> erIms(3);
	float minMean=0;
	unsigned int minLabel=0;
#ifdef HAVE_VISUALIZE 
#ifdef WRITE_DEBUG_IMAGES	
	static int funcCallCount=0;	
#endif
#endif
	
	cerr << "OBS MODEL: fitting:" << endl;
	
	/* Decomment this to
	 * LOAD MANUAL GROUND TRUTH 	 
	{
		char *fname="image4_200x200_gt_strict.pgm";
		cerr << "MANUAL LOAD OF GROUND TRUTH FROM FILE: " << fname << endl;
		Image mgt (fname);
		if (mgt.xsize!=(int) xs || mgt.ysize!=(int)ys)
			ERR_THROW ("Manually loaded ground truth doesn't fit image size!");
				
		erIms[0] = new Image (xs, ys, 1);
		for (unsigned int y=0; y<ys; ++y)
		for (unsigned int x=0; x<xs; ++x)
			erIms[0]->set(x,y, (mgt.get(x,y)?1:0) );
				
		plotClassDistributions (*(this->obs), *erIms[0], cout, 1, 0, 1.);
		cout << endl;		
				
		fitLogNormal (*(this->obs), *erIms[0], 1, 1,
			means[0], sigma_squares[0]);			
		plotLogNormal (cout, means[0], sigma_squares[0], 0, 255, 1);
		
		cerr << "log-normal: m=" << means[0] << " s=" << sqrt(sigma_squares[0]) << endl;
						
		cout << endl;		
		exit(0);
	}
	*/
	
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
	}
	
	// Determine the label with the minimum mean gray level
	for (unsigned int c=0; c<NO_CLASSES; ++c)
	{
		DataSerie s;
		dataSerieFromImage (*(this->obs), *erIms[c], c, s);		
		means[c] = s.getMean();
		if (c==0 || means[c]<minMean)
		{
			minMean=means[c];
			minLabel=c;
		}
	}
		
	// Travers all classes and create a binary image for each class
	// Then fit a generalized gaussian to each class
	for (unsigned int c=0; c<NO_CLASSES; ++c)
	{									
		// Fit the distribution, but the staticical law
		// depends on the class
		
		
		ERR_THROW ("ObsModelBleedThroughGG<TI>::estimateParams(): Lognormal disabled!");
		/*
		if (c==minLabel)
		{
			fitLogNormal (*(this->obs), *erIms[c], 1, 1, 
				means[c], sigma_squares[c]);
			cerr << "class " << c << " (LogN): mu=" << means[c]
					<< "\ts2=" << sigma_squares[c] << endl;				
		}
		
		else 
		*/
		
		
		{
			fitGeneralizedGaussian (*(this->obs), *erIms[c], 1, 
				means[c], sigma_squares[c], shapes[c]);
			cerr << "class " << c << " (GG): mu=" << means[c]
					<< "\ts2=" << sigma_squares[c]
					<< "\tp=" << shapes[c] << endl;
		}							
		/*
		case 2:
			#warning "mean value manually set"
			means[c] = 132.;
			fitGeneralizedGaussianForMean (*(this->obs), *erIms[c], 1, 
				means[c], sigma_squares[c], shapes[c]);	
			cerr << "class " << c << " (GG): mu=" << means[c]
					<< "\ts2=" << sigma_squares[c]
					<< "\tp=" << shapes[c] << endl;
			break;
		*/
	}		
		
#ifdef HAVE_VISUALIZE 
#ifdef WRITE_DEBUG_IMAGES
	// Create a new file
	ostringstream s; 
	s << "likelihoodplot_" << funcCallCount++ << ".txt"; 
	ofstream st (s.str().c_str(), ios::out);
	if (!st.good())
		ERR_THROW ("Cannot open file " << s.str() << "for writing!\n");

	// Debug: plot the class histograms as well as the fitted distributions
	for (unsigned int c=0; c<NO_CLASSES; ++c)
	{							 								
	 	plotClassDistributions (this->obs, erIms[c], st, 1, 0, 1.);
		st << endl;
						 
		 #warning LOG NORMAL DISABLED
		/*
		if (c==minLabel)		
			plotLogNormal (st, means[c], sigma_squares[c], 0, 255, 1);
		else
		*/
		
		
			plotGeneralizedGaussian (st, means[c], sigma_squares[c], shapes[c], 0, 255, 1, false);
		st << endl;		
	}
	st.close();	
#endif	
#endif		
	// Clean up
	for (unsigned int c=0; c<NO_CLASSES; ++c)
		delete erIms[c];		
}

/***********************************************************************
 * Update the parameters, 
 * i.e. calculate the precomputed values which accelerate the 
 * calculations
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThroughGG<TI>::updateParams()
{
}

/***********************************************************************
 * Visualize the result by creating an image where each pixel
 * contains the color of its label
 ***********************************************************************/
 
template <class TI>  
TI * ObsModelBleedThroughGG<TI>::visualizeSegmentation()
{ 
	ERR_THROW ("ObsModelBleedThroughGG<TI>::visualizeSegmentation:"
		"the general template has not yet been implemented.\n");
}

/***********************************************************************
 * Print the parameters to stderr
 ***********************************************************************/
 
template <class TI>  
void ObsModelBleedThroughGG<TI>::printParams()
{
	ERR_THROW ("Not yet implemented.");	
}

#endif
