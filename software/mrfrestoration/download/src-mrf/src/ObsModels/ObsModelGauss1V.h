/***********************************************************************
 * An observation model:
 * Single variate Gaussian model, 
 * where the labels represent the following original gray values:
 * label 0 = gv 0
 * label 1 = gv 255
 *
 * The only parameter is the variance sigma^2
 *
 * Author: Christian Wolf
 * Begin: 25.5.2005
 ***********************************************************************/
 
#ifndef _WOLF_MRFOBSMODELGAUSS1V_H_
#define _WOLF_MRFOBSMODELGAUSS1V_H_

// From the main module
#include <CIL.h>		

// From this module
#include <ObsModel.h>
		
template <class TI> 
class ObsModelGauss1V : public ObsModel<TI>
{		
	public:	
		typedef typename TI::PixelType PixelType;
		
		ObsModelGauss1V (bool doBlur);
		
		float condEnergy (int x, int y);
		float logLikelihood (unsigned char obs, unsigned char lab);
		
		void estimateParams();
		void updateParams();

	private:	// METHODS
			
	private: 	// DATA

		float sigma_square;	
		float mul, exp_div;
};

/***********************************************************************
 * Calculate the likelihood
 ***********************************************************************/
 
template <class TI>  
float ObsModelGauss1V<TI>::condEnergy (int x, int y)
{
	float f=255.0*(float) this->lab->get(x,y) - (float) this->obs->get(x,y); 
	f*=f;
	return mul*exp(f/exp_div);		
}

/***********************************************************************
 * Calculate the log likelihood
 * directly from a label and an observed value
 ***********************************************************************/
 
template <class TI>  
float ObsModelGauss1V<TI>::logLikelihood (unsigned char obs, unsigned char lab)
{		
	ERR_THROW ("logLikelihood: not implemented");
	return 0;		
}

/***********************************************************************
 * estimate the parameters from the current labeling 
 * and the observation:
 * Use the standard ML estimators
 ***********************************************************************/
 
template <class TI>  
void ObsModelGauss1V<TI>::estimateParams()
{
	unsigned int xs=this->obs->xsize;
	unsigned int ys=this->obs->ysize;
	float d;
	
	// Calculate the variance
	sigma_square=0;
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		// Class 0 = gv 0
		if (this->lab->get(x,y)==0)
			d = this->obs->get(x,y);
			
		// Class 1 = gv 255
		else
			d = (255 - this->obs->get(x,y));
		
		sigma_square += d*d;
	}
	sigma_square /= (float) (xs*ys);

	updateParams();
}

/***********************************************************************
 * Update the parameters, 
 * i.e. calculate the precomputed values which accelerate the 
 * calculations
 ***********************************************************************/
 
template <class TI>  
void ObsModelGauss1V<TI>::updateParams()
{
	if (sigma_square==0)
		ERR_THROW ("ObsModelGauss1V: noise variance is zero!");
		
	mul=1.0/sqrt(2.0*M_PI*sigma_square);
	exp_div=2.0*sigma_square;
}

#endif
