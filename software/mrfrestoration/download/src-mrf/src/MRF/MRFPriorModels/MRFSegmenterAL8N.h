/***********************************************************************
 * Segmentation using Markov Random Field Models:
 * The Auto Logistic Model 
 * with 8-connected neighborhood
 *
 * AL = 2 classes !!!!
 * if you need more, use the MLL8N model
 *
 * Author: Christian Wolf
 * Begin: 24.5.2005
 ***********************************************************************/
 
#ifndef _WOLF_MRFSEGMENTERALOG8N_H_
#define _WOLF_MRFSEGMENTERALOG8N_H_

// from the IMAGE PROCESSING module
#include <ImageProc.h>

// From this module
#include "MRFSegmenter.h"
#include "Neighborhood3x3.h"

/***********************************************************************
 * Needed for the clique potential function
 ***********************************************************************/

#define GAMMA(x,y)			((x)==(y)?1:-1)
#define GAMMAMASK(x,y)		((((x)>0 && (y)>0) || ((x)==0 && (y)==0)) ? 1 : -1)

enum {
	MRFP_AL8N_SINGLE=0,
	
	MRFP_FIRST_PAIRWISE,
	MRFP_AL8N_1_H=MRFP_FIRST_PAIRWISE,
	MRFP_AL8N_1_V,
	MRFP_AL8N_1_RD,
	MRFP_AL8N_1_LD,
	
	MRFP_AL8N_AFTER_LAST
};


	
template <class TI>
class MRFSegmenterAL8N : public MRFSegmenter<TI>
{

	public:
	
		// Constructor and Destructor
		MRFSegmenterAL8N(TI *o, Image *l, ObsModel<TI> *m, char *xForceMRFParams);
		~MRFSegmenterAL8N();
		
		// The properties of this model
		virtual void _init();
			
	protected:	// methods		
	
		// The properties of this model
		unsigned int getCliqueSizeX()				{ return 3; }
		unsigned int getCliqueSizeY()				{ return 3; }
		unsigned int getCliqueMaskCenterPixel()		{ return CL_MASK_S; }
		unsigned int getCliqueArrIndexCenterPixel()	{ return 4; }
		unsigned int firstSingleCliqueParameter()	{ return MRFP_AL8N_SINGLE; }
		unsigned int firstPairwiseCliqueParameter()	{ return MRFP_FIRST_PAIRWISE; }
		float priorEnergy (int x, int y);
		void parameterLessTermOfEnergy2Labels (unsigned int clab, Vector<float> &rv);		
		void parameterLessTermOfEnergyMLabels (unsigned char *lab, Vector<float> &rv);
		// help functions
		float _priorEnergy (Image *i, Vector<float> &xMRFParams, int x, int y);		
		void  _mirrorParamsX(Vector<float> &xMRFParams);
		float _sumOfPairwiseParams(Vector<float> &xMRFParams);
};

/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>
MRFSegmenterAL8N<TI>::MRFSegmenterAL8N(TI *o, Image *l, ObsModel<TI> *m, char *xForceMRFParams)
	: MRFSegmenter<TI> (o,l,m,xForceMRFParams,2)
{
	cerr << "New model: MRFSegmenterAL8N" << endl;
	
	// The clique potential parameters
	this->noMRFParams=MRFP_AL8N_AFTER_LAST;
	this->MRFParams.resize(this->noMRFParams);
	if (this->noClasses!=2)
		ERR_THROW ("MRFSegmenterAL8N<>: n.o.classes!=2 not supported in this model!");
}

/***********************************************************************
 * Destructor
 ***********************************************************************/

template <class TI>
MRFSegmenterAL8N<TI>::~MRFSegmenterAL8N()
{
}

/***********************************************************************
 * The initial segmentation. We use a K-means
 ***********************************************************************/

template <class TI>
void MRFSegmenterAL8N<TI>::_init()
{
	typedef typename TI::PixelType PixelType;	
	unsigned int dim=this->obs->nbColorPlanes();
	unsigned int xs=this->obs->xsize;
	unsigned int ys=this->obs->ysize;	
	unsigned int pixnum;	
	Vector<float> v(dim);		
	KMeansClusterer<PixelType> kmc(dim, this->noClasses);	
		
	// perform K-means clustering
	HRULE; cerr << "K-Means Initialisation (MRFSegmenterAL8N):" << endl;
		
	// Collect the pixels and put them into the
	// K-Means clusterer
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		for (unsigned int d=0; d<dim; ++d)
			v.set(d, this->obs->get(d+1,x,y));
		kmc.add (v);
	}
	
	kmc.initRandom();	
	kmc.doCluster();
	
	// Collect the results ->
	// initialise the label field. Only the recto label field
	// is initialised to 3 classes at this moment.
	pixnum=0;	
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		// Initialise the labeling
		this->lab->set(x,y, kmc[pixnum]);
		++pixnum;
	}	
}

/***********************************************************************
 * Change the parameters so that they reflect the model for a page 
 * which has been mirrored in X direction.
 * In this case we just need to swap the two diagonal parameters
 ***********************************************************************/

template <class TI>
void MRFSegmenterAL8N<TI>::_mirrorParamsX(Vector<float> &xMRFParams)
{
	float foo=xMRFParams[MRFP_AL8N_1_RD];
	xMRFParams[MRFP_AL8N_1_RD] = xMRFParams[MRFP_AL8N_1_LD];
	xMRFParams[MRFP_AL8N_1_LD] = foo;
}

/***********************************************************************
 * 
 ***********************************************************************/

template <class TI>
float MRFSegmenterAL8N<TI>::_sumOfPairwiseParams(Vector<float> &xMRFParams)
{
	float rv=0;
	for (unsigned int i=MRFP_FIRST_PAIRWISE; i<MRFP_AL8N_AFTER_LAST; ++i)
		rv += xMRFParams[i];
	
	return rv;
}

/***********************************************************************
 * Calculate the prior probability 
 ***********************************************************************/
 
template <class TI>  
float MRFSegmenterAL8N<TI>::priorEnergy (int x, int y) 
{
	return _priorEnergy (this->lab, this->MRFParams, x, y);
}

template <class TI>  
float MRFSegmenterAL8N<TI>::_priorEnergy (Image *labim, Vector<float> &xMRFParams, int x, int y) 
{
	float sum;
	int f=labim->get(x,y);
	
	/*
	#warning PRIOR IS ZERO!!!!!!
	return 0;
	*/
		
	// The contribution of the site wise clique
	sum = (f==0 ? 0 : xMRFParams[MRFP_AL8N_SINGLE]);
	
	// The contribution of the pair wise cliques: 
	// the different directions
	sum += xMRFParams[MRFP_AL8N_1_H] * 
		(GAMMA(f,labim->get(x-1,y)) + GAMMA(f,labim->get(x+1,y)));
	sum += xMRFParams[MRFP_AL8N_1_V] * 
		(GAMMA(f,labim->get(x,y-1)) + GAMMA(f,labim->get(x,y+1)));
	sum += xMRFParams[MRFP_AL8N_1_RD] * 
		(GAMMA(f,labim->get(x-1,y+1)) + GAMMA(f,labim->get(x+1,y-1)));
	sum += xMRFParams[MRFP_AL8N_1_LD] * 
		(GAMMA(f,labim->get(x-1,y-1)) + GAMMA(f,labim->get(x+1,y+1)));
		
#ifdef PEDANTIC_CHECK_CODE		
	if (isnan(sum)) 
	{
		ERR_THROW ("MRFSegmenterAL8N::priorEnergy(): energy is NaN; details: " <<
		GAMMA(f,labim->get(x-1,y)) << " " << GAMMA(f,labim->get(x+1,y)) << " " << 
		GAMMA(f,labim->get(x,y-1)) << " " << GAMMA(f,labim->get(x,y+1)) << " " << 
		GAMMA(f,labim->get(x-1,y+1)) << " " << GAMMA(f,labim->get(x+1,y-1)) << " " << 
		GAMMA(f,labim->get(x-1,y-1)) << " " << GAMMA(f,labim->get(x+1,y+1)) );
	}
#endif
	
	return sum;
}

/***********************************************************************
 * Calculate Ni()
 * i.e. the function returning the different contributions to each 
 * parameter
 * Needed for the least squares estimation of the parameters
 ***********************************************************************/
 
template <class TI>  
void MRFSegmenterAL8N<TI>::parameterLessTermOfEnergyMLabels (unsigned char *lab, Vector<float> &rv)
{
	ERR_THROW ("MRFSegmenterAL8N<TI>::parameterLessTermOfEnergyMLabels(): the model allows only 2 labels!");
}
 
template <class TI>  
void MRFSegmenterAL8N<TI>::parameterLessTermOfEnergy2Labels (unsigned int cLab, Vector<float> &rv)
{
	unsigned int f=cLab&CL_MASK_S;
	
	// The contribution of the single site clique
	rv[MRFP_AL8N_SINGLE] = (f==0?0:1);
	
	// The contribution of the pair wise cliques: 
	// the different directions
	rv[MRFP_AL8N_1_H] = (GAMMAMASK(f,cLab&CL_MASK_WE) + GAMMAMASK(f,cLab&CL_MASK_EA));
	rv[MRFP_AL8N_1_V] =	(GAMMAMASK(f,cLab&CL_MASK_NO) + GAMMAMASK(f,cLab&CL_MASK_SO));
	rv[MRFP_AL8N_1_RD]= (GAMMAMASK(f,cLab&CL_MASK_SW) + GAMMAMASK(f,cLab&CL_MASK_NE));
	rv[MRFP_AL8N_1_LD]=	(GAMMAMASK(f,cLab&CL_MASK_NW) + GAMMAMASK(f,cLab&CL_MASK_SE));
	
#ifdef PEDANTIC_CHECK_CODE
	for (int i=0; i<MRFP_AL8N_AFTER_LAST; ++i)
	{
		if (isnan(rv[i]))
			ERR_THROW ("MRFSegmenterAL8N::parameterLessTermOfEnergy(): the value with index "
				<< i << " is NaN");
	}	
#endif	
}

#endif
