/***********************************************************************
 * Segmentation using Markov Random Field Models:
 * The Multi Level Logistic Model 
 * with 8-connected neighborhood
 *
 * multiple classes. If you need only 2, 
 * use the AL8N model
 *
 * Author: Christian Wolf
 * Begin: 19.6.2006
 ***********************************************************************/
 
#ifndef _WOLF_MRFSEGMENTERMLL8N_H_
#define _WOLF_MRFSEGMENTERMLL8N_H_

// From this module
#include "MRFSegmenterAL8N.h"
#include "Neighborhood3x3.h"

/***********************************************************************
 * Needed for the clique potential function
 ***********************************************************************/

#define GAMMA(x,y)			((x)==(y)?1:-1)

enum {		
	MRFP_MLL8N_FIRST_PAIRWISE=0,// The pairwise clique parameters
	MRFP_MLL8N_1_H=MRFP_MLL8N_FIRST_PAIRWISE,			
	MRFP_MLL8N_1_V,
	MRFP_MLL8N_1_RD,
	MRFP_MLL8N_1_LD,
		
	MRFP_MLL8N_FIRST_SINGLE		// The single clique parameters
								// ATTENTION: the first clique label
								// (label 0) does _NOT_ have a parameter
};

template <class TI>
class MRFSegmenterMLL8N : public MRFSegmenter<TI>
{

	public:
	
		// Constructor and Destructor
		MRFSegmenterMLL8N(TI *o, Image *l, ObsModel<TI> *m, 
			char *xForceMRFParams, int noClasses, 
			bool xOptDoGaussFilterBeforeKmeans);
		~MRFSegmenterMLL8N();
		
		// The properties of this model
		void _init();
			
	protected:	// methods		
	
		// The properties of this model
		unsigned int getCliqueSizeX()				{ return 3; }
		unsigned int getCliqueSizeY()				{ return 3; }
		unsigned int getCliqueMaskCenterPixel()		{ return CL_MASK_S; }
		unsigned int firstSingleCliqueParameter()	{ return MRFP_MLL8N_FIRST_SINGLE; }
		unsigned int firstPairwiseCliqueParameter()	{ return MRFP_MLL8N_FIRST_PAIRWISE; }
		float priorEnergy (int x, int y);
		void parameterLessTermOfEnergy2Labels (unsigned int clab, Vector<float> &rv);		
		void parameterLessTermOfEnergyMLabels (unsigned char *clab, Vector<float> &rv);		
		unsigned int getCliqueArrIndexCenterPixel()	{ return 4; }
		
		// help functions
		float _priorEnergy (Image *i, Vector<float> &xMRFParams, int x, int y);		
		void  _mirrorParamsX(Vector<float> &xMRFParams);
		
	protected:	// DATA
	
		bool optDoGaussFilterBeforeKmeans;
	
};

/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>
MRFSegmenterMLL8N<TI>::MRFSegmenterMLL8N(TI *o, Image *l, ObsModel<TI> *m,
	char *xForceMRFParams, int xNoClasses, bool xOptDoGaussFilterBeforeKmeans)
	: MRFSegmenter<TI> (o,l,m, xForceMRFParams, xNoClasses)
{
	cerr << "New model: MRFSegmenterMLL8N" << endl;
	
	// The clique potential parameters
	this->noMRFParams=MRFP_MLL8N_FIRST_SINGLE+this->noClasses-1;
	this->MRFParams.resize(this->noMRFParams);
	optDoGaussFilterBeforeKmeans = xOptDoGaussFilterBeforeKmeans;
}

/***********************************************************************
 * Destructor
 ***********************************************************************/

template <class TI>
MRFSegmenterMLL8N<TI>::~MRFSegmenterMLL8N()
{
}

/***********************************************************************
 * The initial segmentation. We use a K-means
 ***********************************************************************/

template <class TI>
void MRFSegmenterMLL8N<TI>::_init()
{
	typedef typename TI::PixelType PixelType;	
	unsigned int dim=this->obs->nbColorPlanes();
	unsigned int xs=this->obs->xsize;
	unsigned int ys=this->obs->ysize;	
	unsigned int pixnum;	
	Vector<float> v(dim);		
	KMeansClusterer<PixelType> kmc(dim, this->noClasses);
	Image *filtered;
		
	// perform K-means clustering
	HRULE; cerr << "K-Means Initialisation (MRFSegmenterMLL8N):" << endl;
	
	if (optDoGaussFilterBeforeKmeans)
	{
		cerr << "Low pass filtering image for k-means input.\n";
		filtered = new Image (*this->obs);	
		FilterMask fm;
		fm  << 1.0 <<  4.0  << 1.0 <<  4.0 <<  1.0 << '\n'
			<< 4.0 << 16.0 << 26.0 << 16.0 <<  4.0 << '\n'
			<< 7.0 << 26.0 << 41.0 << 26.0 <<  7.0 << '\n'
			<< 4.0 << 16.0 << 26.0 << 16.0 <<  4.0 << '\n'
			<< 1.0 <<  4.0  << 7.0 <<  4.0 <<  1.0;	
 		fm.normalize();
		filterMask (*filtered, fm);
	}
	else
	{
		cerr << "_No_ low pass filtering.\n";
		filtered = this->obs;
	}
		
	// Collect the pixels and put them into the
	// K-Means clusterer
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		for (unsigned int d=0; d<dim; ++d)
			v.set(d, filtered->get(d+1,x,y));
		kmc.add (v);
	}
	
	kmc.initRandom();	
	kmc.doCluster();
	
	// Collect the results
	pixnum=0;	
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		// Initialise the labeling
		this->lab->set(x,y, kmc[pixnum]);
		++pixnum;
	}	
	
	// Debug output
	{	Image tmp (*(this->lab));
		brightenImage (tmp);
		reverseVideo (tmp);		
		tmp.write ("x_005_kmeans_reordered_rev.pgm");
	}
	
	if (optDoGaussFilterBeforeKmeans)
		delete filtered;
}

/***********************************************************************
 * Change the parameters so that they reflect the model for a page 
 * which has been mirrored in X direction.
 * In this case we just need to swap the two diagonal parameters
 ***********************************************************************/

template <class TI>
void MRFSegmenterMLL8N<TI>::_mirrorParamsX(Vector<float> &xMRFParams)
{
	float foo=xMRFParams[MRFP_MLL8N_1_RD];
	xMRFParams[MRFP_MLL8N_1_RD] = xMRFParams[MRFP_MLL8N_1_LD];
	xMRFParams[MRFP_MLL8N_1_LD] = foo;
}

/***********************************************************************
 * Calculate the prior probability 
 ***********************************************************************/
 
template <class TI>  
float MRFSegmenterMLL8N<TI>::priorEnergy (int x, int y) 
{
	return _priorEnergy (this->lab, this->MRFParams, x, y);
}

template <class TI>  
float MRFSegmenterMLL8N<TI>::_priorEnergy (Image *labim, Vector<float> &xMRFParams, int x, int y) 
{
	float sum=0;
	int f=labim->get(x,y);
	
	/*
	#warning PRIOR IS ZERO in MLL8N!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	return 0;
	*/
	
	// The contribution of the site wise clique	
	// Label 0 does _NOT_ have a parameter!!!
	if (f>0)
		xMRFParams[MRFP_MLL8N_FIRST_SINGLE+f-1];
	
	// The contribution of the pair wise cliques: 
	// the different directions
	sum += xMRFParams[MRFP_MLL8N_1_H] * 
		(GAMMA(f,labim->get(x-1,y)) + GAMMA(f,labim->get(x+1,y)));
	sum += xMRFParams[MRFP_MLL8N_1_V] * 
		(GAMMA(f,labim->get(x,y-1)) + GAMMA(f,labim->get(x,y+1)));
	sum += xMRFParams[MRFP_MLL8N_1_RD] * 
		(GAMMA(f,labim->get(x-1,y+1)) + GAMMA(f,labim->get(x+1,y-1)));
	sum += xMRFParams[MRFP_MLL8N_1_LD] * 
		(GAMMA(f,labim->get(x-1,y-1)) + GAMMA(f,labim->get(x+1,y+1)));
		
#ifdef PEDANTIC_CHECK_CODE		
	if (isnan(sum)) 
	{
		ERR_THROW ("MRFSegmenterMLL8N::priorEnergy(): energy is NaN; details: " <<
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
void MRFSegmenterMLL8N<TI>::parameterLessTermOfEnergyMLabels (unsigned char *cLab, Vector<float> &rv)
{
	unsigned char f=cLab[getCliqueArrIndexCenterPixel()];
	
	// The contribution of the single site clique
	// Label 0 does _NOT_ have a parameter!!!
	for (unsigned int i=1; i<this->noClasses; ++i)
		rv[MRFP_MLL8N_FIRST_SINGLE+i-1] = (f==i?1:0);
	
	// The contribution of the pair wise cliques: 
	// the different directions
	rv[MRFP_MLL8N_1_H] = GAMMA(f,cLab[CL_INDEX_WE]) + GAMMA(f,cLab[CL_INDEX_EA]);
	rv[MRFP_MLL8N_1_V] = GAMMA(f,cLab[CL_INDEX_NO]) + GAMMA(f,cLab[CL_INDEX_SO]);
	rv[MRFP_MLL8N_1_RD]= GAMMA(f,cLab[CL_INDEX_SW]) + GAMMA(f,cLab[CL_INDEX_NE]);
	rv[MRFP_MLL8N_1_LD]= GAMMA(f,cLab[CL_INDEX_NW]) + GAMMA(f,cLab[CL_INDEX_SE]);
	
#ifdef PEDANTIC_CHECK_CODE
	for (unsigned int i=0; i<this->noMRFParams; ++i)
	{
		if (isnan(rv[i]))
			ERR_THROW ("MRFSegmenterMLL8N::parameterLessTermOfEnergy(): the value with index "
				<< i << " is NaN");
	}	
#endif	
}
  
template <class TI>  
void MRFSegmenterMLL8N<TI>::parameterLessTermOfEnergy2Labels (unsigned int cLab, Vector<float> &rv)
{
	ERR_THROW ("MRFSegmenterMLL8N<TI>::parameterLessTermOfEnergy2Labels: "
		"does not make sense with this model!");

	/*
	unsigned int f=cLab&CL_MASK_S;
	
	// The contribution of the single site clique
	for (unsigned int i=0; i<this->noClasses; ++i)
		rv[MRFP_MLL8N_FIRST_SINGLE+i] = (f==i?0:1);
	
	// The contribution of the pair wise cliques: 
	// the different directions
	rv[MRFP_MLL8N_1_H] = (GAMMA(f,cLab&CL_MASK_WE) + GAMMA(f,cLab&CL_MASK_EA));
	rv[MRFP_MLL8N_1_V] = (GAMMA(f,cLab&CL_MASK_NO) + GAMMA(f,cLab&CL_MASK_SO));
	rv[MRFP_MLL8N_1_RD]= (GAMMA(f,cLab&CL_MASK_SW) + GAMMA(f,cLab&CL_MASK_NE));
	rv[MRFP_MLL8N_1_LD]= (GAMMA(f,cLab&CL_MASK_NW) + GAMMA(f,cLab&CL_MASK_SE));
	
#ifdef PEDANTIC_CHECK_CODE
	for (int i=0; i<this->noMRFParams; ++i)
	{
		if (isnan(rv[i]))
			ERR_THROW ("MRFSegmenterMLL8N::parameterLessTermOfEnergy(): the value with index "
				<< i << " is NaN");
	}	
#endif	
	*/
}
	
#endif
