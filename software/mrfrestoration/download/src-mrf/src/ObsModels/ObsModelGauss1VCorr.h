/***********************************************************************
 * ObsModelGauss1VCorr.h
 * An observation model:
 * _CORRELATED_ single variate Gaussian Noise 
 *
 * Author: Christian Wolf
 * Begin: 25.5.2005
 ***********************************************************************/
 
#ifndef _WOLF_MRFOBSMODELGAUSS1VCORR_H_
#define _WOLF_MRFOBSMODELGAUSS1VCORR_H_

// C++
#include <vector>

// From the main module
#include <CIL.h>		

// From the IMAGE module
#include <Image.h>

// From this module
#include <Segmenter.h>
#include <ObsModel.h>
		
template <class TI> 
class ObsModelGauss1VCorr : public ObsModel<TI>
{		
	public:	
		typedef typename TI::PixelType PixelType;
		
		ObsModelGauss1VCorr (bool doBlur);
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
	
		vector<float> means, 
					  sigma_squares;
		
		vector<float> add_for_log,muls_in_exp;
};

/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>  
ObsModelGauss1VCorr<TI>::ObsModelGauss1VCorr (bool doBlur)
	: ObsModel<TI>(doBlur)
{
}

/***********************************************************************
 * Initialize
 ***********************************************************************/
 
template <class TI>  
void ObsModelGauss1VCorr<TI>::init (Segmenter<TI> *xsegmenter)
{
	ObsModel<TI>::init(xsegmenter);
	obs=xsegmenter->getObs();
	lab=xsegmenter->getLab();
	
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
 
#define M(i,j)	(means[this->lab->get(i,j)])
 
template <class TI>  
float ObsModelGauss1VCorr<TI>::condEnergy (int x, int y)
{		
	int l=this->lab->get(x,y);
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
	
	//The energy is the NEGATIVE (!!!!) of the likelihood!
	return rv - add_for_log[l];		
}
#undef M

/***********************************************************************
 * Calculate the log likelihood
 * directly from a label and an observed value
 ***********************************************************************/
 
template <class TI>  
float ObsModelGauss1VCorr<TI>::logLikelihood (unsigned char obs, unsigned char lab,
	unsigned int x, unsigned int y)
{		
	return logLikelihood(obs,lab,x,y);
} 
 
template <class TI>  
float ObsModelGauss1VCorr<TI>::logLikelihood (unsigned char obs, unsigned char lab)
{		
	float rv=means[lab] - (float) obs;
	rv=rv*rv*muls_in_exp[lab];	
	return add_for_log[lab] - rv;		
}

/***********************************************************************
 * estimate the parameters from the current labeling 
 * and the observation:
 * Use the standard ML estimators
 ***********************************************************************/
 
class M 
{ 
	public:
		float m; 
		unsigned int ind; 
		M (float xm, unsigned int xind) { m=xm; ind=xind; }
		bool operator < (const M & o) const { return m<o.m;} 
}; 
 
template <class TI>  
void ObsModelGauss1VCorr<TI>::estimateParams()
{
	unsigned int xs=this->obs->xsize;
	unsigned int ys=this->obs->ysize;
	vector<unsigned int> noPixels(this->getNoClasses());
	unsigned int l;
	float d;
	
	cerr << "Estimating observation model parameters (Gauss1VCorr)" << endl;
	
	// init
	for (unsigned int i=0; i<this->getNoClasses(); ++i)
	{
		means[i]=0;
		sigma_squares[i]=0;
		noPixels[i]=0;
	}
	
	// Calculate the means
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		l=this->lab->get(x,y);
		if (l>=this->getNoClasses())
			ERR_THROW ("Internal error in ObsModelGauss1VCorr<TI>::estimateParams():"
				" invalid label: " << l);
		means[l] += this->obs->get(x,y);
		++noPixels[l];
	}		
	
	for (unsigned int i=0; i<this->getNoClasses(); ++i)
		means[i] /= (float) noPixels[i];
	
	// Calculate the variances
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{		
		l = this->lab->get(x,y);
		d = this->obs->get(x,y) - means[l];
		sigma_squares[l] += d*d;
	}
	for (unsigned int i=0; i<this->getNoClasses(); ++i)
		sigma_squares[i] /= (float) noPixels[i];
		
	// DEBUG
	this->printParams();
			 	
	/*
	{
		#warning "CHANGING MEANS AND SIGMA SQUARES!!!!!"
		cerr << "Using fixed params from the ground truth!!!!!!!!\n";
		
		// Sort the means		
		set<M> s;
		s.insert(M(means[0],0));
		s.insert(M(means[1],1));
		s.insert(M(means[2],2));
		set<M>::iterator iter=s.begin();
		
		means[iter->ind] = 51; 	++iter;
		means[iter->ind] = 97; 	++iter;
		means[iter->ind] =196;
		
		for (unsigned int i=0; i<this->getNoClasses(); ++i)
			// sigma_squares[i] /= 10;			
			sigma_squares[i] = 100;
	}
	*/
			 	
	updateParams();
}

/***********************************************************************
 * Update the parameters, 
 * i.e. calculate the precomputed values which accelerate the 
 * calculations
 ***********************************************************************/
 
template <class TI>  
void ObsModelGauss1VCorr<TI>::updateParams()
{
	for (unsigned int i=0; i<this->getNoClasses(); ++i)
	{
		if (sigma_squares[i]==0)
			ERR_THROW ("ObsModelGauss1V: noise variance is zero!");
						
		muls_in_exp[i]=1.0/(2.0*sigma_squares[i]);
		add_for_log[i] = 1.0/(sqrtf(2.*M_PI*sigma_squares[i]));		
		if (add_for_log[i]<=0)
			ERR_THROW ("Internal error in "
			"ObsModelGuass1VCorr::updateParams: "
			"invalid model parameters");
		add_for_log[i] = logf(add_for_log[i]);
	}
}

/***********************************************************************
 * Visualize the result by creating an image where each pixel
 * contains the gray value of its label
 ***********************************************************************/
 
template <class TI>  
TI * ObsModelGauss1VCorr<TI>::visualizeSegmentation()
{ 
	ERR_THROW ("ObsModelGauss1VCorr<TI>::visualizeSegmentation():"
		"the general template has not yet been implemented.\n");
}
 
template <>
Image * ObsModelGauss1VCorr<Image>::visualizeSegmentation()
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
void ObsModelGauss1VCorr<TI>::printParams()
{
	for (unsigned int i=0; i<this->getNoClasses(); ++i)
		cerr << "class " << i << ": mean=" << means[i]
			 << " sigma=" << sqrt(sigma_squares[i]) << endl;	
}

#endif
