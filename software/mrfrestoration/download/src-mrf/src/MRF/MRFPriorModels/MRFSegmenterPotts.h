/***********************************************************************
 * Segmentation using Markov Random Field Models:
 * The potts model.
 * Very similar to the Multi Level Logistic Model 
 *
 * Author: Christian Wolf
 * Begin: 12.3.2008
 ***********************************************************************/
 
#ifndef _WOLF_MRFSEGMENTERPOTTS_H_
#define _WOLF_MRFSEGMENTERPOTTS_H_

// C++
#include <sstream>

// From Kolmogorov's graph cut code
#include "graphcut.h"

// From the BINARIZATION module
#include "Binarization.h"

// From this module
#include "MRFSegmenter.h"
#include "Neighborhood3x3.h"

using namespace std;

/***********************************************************************
 * MISC
 ***********************************************************************/

#define KR_DELTA(x,y)				((x)==(y)?1:0)
#define REV_KR_DELTA(x,y)			((x)==(y)?0:1)
#define WEIGHT_INFINITY 			(1.0e38)
#define NODENR(x,y)					((y)*xs+(x))
#define PARAM2WEIGHT_GRCT(x)		(-1.0*this->betaAdjust*(x))
#define PARAM2WEIGHT_EXPMOVE(x)		(-1.0*this->betaAdjust*(x))
#define LOGLH2WEIGHT(x)				(-1.0*(x))
#define WEIGHT_C2_H(l1,l2) 			\
	(PARAM2WEIGHT_GRCT(this->MRFParams[MRFP_POTTS_1_H]*REV_KR_DELTA((l1),(l2))))
#define WEIGHT_C2_V(l1,l2) 			\
	(PARAM2WEIGHT_GRCT(this->MRFParams[MRFP_POTTS_1_V]*REV_KR_DELTA((l1),(l2))))

/***********************************************************************/

enum {		
	MRFP_POTTS_FIRST_PAIRWISE=0,// The pairwise clique parameters
	MRFP_POTTS_1_H=MRFP_POTTS_FIRST_PAIRWISE,			
	MRFP_POTTS_1_V,	
		
	MRFP_POTTS_FIRST_SINGLE		// The single clique parameters
								// ATTENTION: the first clique label
								// (label 0) does _NOT_ have a parameter
};

// For the graph cut 
typedef Graph<float,float,float> GraphType;

template <class TI>
class MRFSegmenterPotts : public MRFSegmenter<TI>
{

	public:
	
		// Constructor and Destructor
		MRFSegmenterPotts(TI *o, Image *l, ObsModel<TI> *m, 
			TI *xEventualColorObs,
			char *xForceMRFParams, int noClasses, 
			bool xOptDoGaussFilterBeforeKmeans, float xBetaAdjust,
			ColorSpaceCode xColorSpaceCode);
		~MRFSegmenterPotts();
		
		// The properties of this model
		void _init();
		void _initKMeans();
		void _initSauvola();
		
		FloatMatrix *getThresholdSurface() { return thresholdSurface; }
		
		void graphCutAdaptedToProblem(unsigned int noIterations);
			
	protected:	// methods		
	
		// The properties of this model
		unsigned int getCliqueSizeX()				{ return 3; }
		unsigned int getCliqueSizeY()				{ return 3; }
		unsigned int getCliqueMaskCenterPixel()		{ return CL_MASK_S; }
		unsigned int firstSingleCliqueParameter()	{ return MRFP_POTTS_FIRST_SINGLE; }
		unsigned int firstPairwiseCliqueParameter()	{ return MRFP_POTTS_FIRST_PAIRWISE; }
		float priorEnergy (int x, int y);
		void parameterLessTermOfEnergy2Labels (unsigned int clab, Vector<float> &rv);		
		void parameterLessTermOfEnergyMLabels (unsigned char *clab, Vector<float> &rv);		
		unsigned int getCliqueArrIndexCenterPixel()	{ return 4; }
		
		void graphCut ();
		void expansionMove (unsigned int noCycles, bool fastKolmogorovVersion);
		float clique2GraphCutH (unsigned char lab1, unsigned char lab2);
		float clique2GraphCutV (unsigned char lab1, unsigned char lab2);
		float totalGraphCutEnergy();
		
		// The (fast) Kolmogorov version of the expansion move algorithm
		float singleExpansionMoveKolmogorov (unsigned int alpha);
		void addEdgeEMKolmogorov (GraphType *g, 
			unsigned int xs, unsigned int ys, 
			unsigned int x, unsigned int y, 
			unsigned int dx, unsigned int dy,
			unsigned char alpha, vector<float> &tweights);
			
		// The (slow) Boykov version of the expansion move algorithm
		float singleExpansionMoveBoykov (unsigned int alpha);
		void addEdgeEMBoykov(GraphType *g, 
			unsigned int nodenr1, unsigned int nodenr2,
			unsigned int lab1, unsigned int lab2,
			unsigned int alpha, unsigned int curIndex, bool isHoriz);
			
			// help functions
		float _priorEnergy (Image *i, Vector<float> &xMRFParams, int x, int y);		
		void  _mirrorParamsX(Vector<float> &xMRFParams);
		float _sumOfPairwiseParams(Vector<float> &xMRFParams);

		
	public: 	// DATA
		
		bool doSauvolaInsteadOfKmeans;
		
	protected:	// DATA
	
		bool optDoGaussFilterBeforeKmeans;
		float betaAdjust;
		
		// In case of Sauvola et al initialization:
		FloatMatrix *thresholdSurface;
		
		Image *eventualColorObs;
		ColorSpaceCode colspc;
};

/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>
MRFSegmenterPotts<TI>::MRFSegmenterPotts(TI *o, Image *l, ObsModel<TI> *m,
	TI *xEventualColorObs, char *xForceMRFParams, int xNoClasses, 
	bool xOptDoGaussFilterBeforeKmeans, float xBetaAdjust,
	ColorSpaceCode xColorSpaceCode)
	: MRFSegmenter<TI> (o,l,m, xForceMRFParams, xNoClasses)
{		
	// The clique potential parameters
	this->noMRFParams=MRFP_POTTS_FIRST_SINGLE+this->noClasses-1;
	this->MRFParams.resize(this->noMRFParams);
	optDoGaussFilterBeforeKmeans = xOptDoGaussFilterBeforeKmeans;
	doSauvolaInsteadOfKmeans = false;
	betaAdjust = xBetaAdjust;
	eventualColorObs = xEventualColorObs;
	colspc = xColorSpaceCode;
	
	cerr << "New model: MRFSegmenterPotts (" 
		 << this->noMRFParams << " params)" << endl;
}

/***********************************************************************
 * Destructor
 ***********************************************************************/

template <class TI>
MRFSegmenterPotts<TI>::~MRFSegmenterPotts()
{
	delete thresholdSurface;
}

/***********************************************************************
 * The initial segmentation. 
 ***********************************************************************/

template <class TI>
void MRFSegmenterPotts<TI>::_init()
{
	if (doSauvolaInsteadOfKmeans)
		_initSauvola();
	else
		_initKMeans();
}

/***********************************************************************
 * The initial segmentation, 
 * Special case: we use the improved Sauvola et al.
 ***********************************************************************/

template <class TI>
void MRFSegmenterPotts<TI>::_initSauvola()
{
	HRULE; cerr << "Wolf-improved Sauvola et al. initialisation (MRFSegmenterPotts):" << endl
		   		<< "Calculating preliminary first step threshold surface ...";
		   cerr.flush();
	
	*(this->lab) = *(this->obs);
	thresholdSurface = surfaceNiblackImproved (*(this->lab), NIBLACK_WOLF2, 40, 40, 0.5, 128);
	thresholdWithSurface (*(this->lab), *thresholdSurface);
	
	// Label 255 -> 1
	for (unsigned int y=0; y<this->lab->ysize; ++y)
	for (unsigned int x=0; x<this->lab->xsize; ++x)
		this->lab->set(x,y, (this->lab->get(x,y) ? 1 : 0));
		
	cerr << "Done." << endl;
}

/***********************************************************************
 * The initial segmentation,
 * Special case: we use a K-means
 ***********************************************************************/

template <class TI>
void MRFSegmenterPotts<TI>::_initKMeans()
{
	typedef typename TI::PixelType PixelType;
	Image *chosenObs = ((this->eventualColorObs==NULL) ? this->obs : this->eventualColorObs);
	unsigned int dim=chosenObs->nbColorPlanes();
	unsigned int xs=chosenObs->xsize;
	unsigned int ys=chosenObs->ysize;		
	unsigned int pixnum;	
	Vector<float> v(dim);		
	KMeansClusterer<PixelType> kmc(dim, this->noClasses);
	Image *filtered;
	LabCie94Distance *customDistance;
	double pl, pu, pv;
		
	// perform K-means clustering
	HRULE; cerr << "K-Means Initialisation (MRFSegmenterPotts):" << endl;
	
	if (optDoGaussFilterBeforeKmeans)
	{
		cerr << "Low pass filtering image for k-means input.\n";
		filtered = new Image (*chosenObs);	
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
		filtered = chosenObs;
	}
	
	if ((filtered->nbColorPlanes()!=3) && (this->colspc!=COLSPC_RGB))
		ERR_THROW ("Incompatible colorspace and dimension in MRFSegmenterPotts<TI>::_init()\n");

	// Collect the pixels and put them into the
	// K-Means clusterer
	switch (this->colspc)
	{
		case COLSPC_RGB:
			cerr << "Kmeans using RGB with " << chosenObs->nbColorPlanes() << "color planes\n";
			for (unsigned int y=0; y<ys; ++y)
			for (unsigned int x=0; x<xs; ++x)
			{
				for (unsigned int d=0; d<dim; ++d)
					v.set(d, filtered->get(d+1,x,y));
				kmc.add (v);
			}
			break;
			
		case COLSPC_LUV:
			cerr << "Kmeans using L*u*v* color space.\n";
			for (unsigned int y=0; y<ys; ++y)
			for (unsigned int x=0; x<xs; ++x)
			{
				rgb2luv (filtered->get(PLANE_RED,x,y)/255., 
						 filtered->get(PLANE_GREEN,x,y)/255., 
						 filtered->get(PLANE_BLUE,x,y)/255.,
						 &pl, &pu, &pv);
				v.set(0, pl);
				v.set(1, pu);
				v.set(2, pv);
				kmc.add (v);
				
			}
			break;	
			
		case COLSPC_HSV:
			for (unsigned int y=0; y<ys; ++y)
			for (unsigned int x=0; x<xs; ++x)
			{
				rgb2hsv (filtered->get(PLANE_RED,x,y)/255., 
						 filtered->get(PLANE_GREEN,x,y)/255., 
						 filtered->get(PLANE_BLUE,x,y)/255.,
						 &pl, &pu, &pv);
				pl = pl*255./360.;
				pu = pu*255.;
				pv = pv*255.;
				v.set(0, pl);
				v.set(1, pu);
				v.set(2, pv);		
				kmc.add (v);
			}
			break;
			
		case COLSPC_LAB:
		case COLSPC_LAB_CIE94:
			cerr << "Kmeans using L*a*b* color space.\n";
			for (unsigned int y=0; y<ys; ++y)
			for (unsigned int x=0; x<xs; ++x)
			{
				rgb2lab (filtered->get(PLANE_RED,x,y)/255., 
						 filtered->get(PLANE_GREEN,x,y)/255., 
						 filtered->get(PLANE_BLUE,x,y)/255.,
						 &pl, &pu, &pv);
				v.set(0, pl);
				v.set(1, pu);
				v.set(2, pv);
				kmc.add (v);				
			}
			break;
			
		default:
			ERR_THROW ("Unknown color space in segmentKMeans()\n");
	}
	
	// Do we use a color space with non euclidean distance?
	if (this->colspc==COLSPC_LAB_CIE94)
	{
		customDistance = new LabCie94Distance;
		kmc.setDistance(customDistance);	
		cerr << "Setting CIE94 distance function for kmeans.\n";
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
void MRFSegmenterPotts<TI>::_mirrorParamsX(Vector<float> &xMRFParams)
{
	// No diagonal neighbors, nothing to do ...
}

/***********************************************************************
 * 
 ***********************************************************************/

template <class TI>
float MRFSegmenterPotts<TI>::_sumOfPairwiseParams(Vector<float> &xMRFParams)
{
	float rv=0;
	for (unsigned int i=MRFP_POTTS_FIRST_PAIRWISE; i<MRFP_POTTS_FIRST_SINGLE; ++i)
		rv += xMRFParams[i];
		
	return rv;
}

/***********************************************************************
 * Calculate the prior probability 
 ***********************************************************************/

template <class TI>  
float MRFSegmenterPotts<TI>::priorEnergy (int x, int y) 
{
	return _priorEnergy (this->lab, this->MRFParams, x, y);
}

template <class TI>  
float MRFSegmenterPotts<TI>::_priorEnergy (Image *labim, Vector<float> &xMRFParams, int x, int y) 
{
	float sum=0;
	int f=labim->get(x,y);
	
	// The contribution of the site wise clique	
	// Label 0 does _NOT_ have a parameter!!!
	if (f>0)
		sum += xMRFParams[MRFP_POTTS_FIRST_SINGLE+f-1];
	
	// The contribution of the pair wise cliques: 
	// the different directions
	sum += xMRFParams[MRFP_POTTS_1_H] * 
		(KR_DELTA(f,labim->get(x-1,y)) + KR_DELTA(f,labim->get(x+1,y)));
	sum += xMRFParams[MRFP_POTTS_1_V] * 
		(KR_DELTA(f,labim->get(x,y-1)) + KR_DELTA(f,labim->get(x,y+1)));
	
		
#ifdef PEDANTIC_CHECK_CODE		
	if (isnan(sum)) 
	{
		ERR_THROW ("MRFSegmenterPotts::priorEnergy(): energy is NaN; details: " <<
		KR_DELTA(f,labim->get(x-1,y)) << " " << KR_DELTA(f,labim->get(x+1,y)) << " " << 
		KR_DELTA(f,labim->get(x,y-1)) << " " << KR_DELTA(f,labim->get(x,y+1)) << " ");
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
void MRFSegmenterPotts<TI>::parameterLessTermOfEnergyMLabels (unsigned char *cLab, Vector<float> &rv)
{
	unsigned char f=cLab[getCliqueArrIndexCenterPixel()];
	
	// The contribution of the single site clique
	// Label 0 does _NOT_ have a parameter!!!
	for (unsigned int i=1; i<this->noClasses; ++i)
		rv[MRFP_POTTS_FIRST_SINGLE+i-1] = (f==i?1:0);
	
	// The contribution of the pair wise cliques: 
	// the different directions
	rv[MRFP_POTTS_1_H] = KR_DELTA(f,cLab[CL_INDEX_WE]) + KR_DELTA(f,cLab[CL_INDEX_EA]);
	rv[MRFP_POTTS_1_V] = KR_DELTA(f,cLab[CL_INDEX_NO]) + KR_DELTA(f,cLab[CL_INDEX_SO]);	
	
#ifdef PEDANTIC_CHECK_CODE
	for (unsigned int i=0; i<this->noMRFParams; ++i)
	{
		if (isnan(rv[i]))
			ERR_THROW ("MRFSegmenterPotts::parameterLessTermOfEnergyMLabels(): the value with index "
				<< i << " is NaN");
	}	
#endif	
}
  
template <class TI>  
void MRFSegmenterPotts<TI>::parameterLessTermOfEnergy2Labels (unsigned int cLab, Vector<float> &rv)
{	
	unsigned int f=cLab&CL_MASK_S;
	
	// The contribution of the single site clique	
	rv[MRFP_POTTS_FIRST_SINGLE+0] = (f>0?1:0);
	
	// The contribution of the pair wise cliques: 
	// the different directions
	rv[MRFP_POTTS_1_H] = (KR_DELTA(f,cLab&CL_MASK_WE) + KR_DELTA(f,cLab&CL_MASK_EA));
	rv[MRFP_POTTS_1_V] = (KR_DELTA(f,cLab&CL_MASK_NO) + KR_DELTA(f,cLab&CL_MASK_SO));	
	
#ifdef PEDANTIC_CHECK_CODE
	for (int i=0; i<this->noMRFParams; ++i)
	{
		if (isnan(rv[i]))
			ERR_THROW ("MRFSegmenterPotts::parameterLessTermOfEnergy2Labels(): the value with index "
				<< i << " is NaN");
	}	
#endif	
}

/***********************************************************************
 * Graph cut optimization
 ***********************************************************************/

template <class TI>  
void MRFSegmenterPotts<TI>::graphCutAdaptedToProblem (unsigned int noIterations)
{		
	if (this->noClasses==2)
		graphCut();		
	else
	{
	/*
#warning Running 2 different versions of expansion move!!	
		Image initlab (*this->lab);
		expansionMove(noIterations, false);
		this->lab->write("x_expmove_boykov.pgm");
		*this->lab = initlab;
		*/
		expansionMove(noIterations, true);
		/*
		this->lab->write("x_expmove_kolmogorov.pgm");
		*/
	}
}

/***********************************************************************
 * Graph cut optimization
 ***********************************************************************/

template <class TI>  
void MRFSegmenterPotts<TI>::graphCut ()
{		
	unsigned int xs=this->lab->xsize;
	unsigned int ys=this->lab->ysize;
	float flow;
	GraphType *g;
	
	if (this->noClasses!=2)
		ERR_THROW ("Graph cut only works with 2 classes!");
	
	cerr << "Creating graph ..."; cerr.flush();
	
	g = new GraphType(xs*ys, 4*xs*ys);
		
	cerr << "\nDone. Adding nodes and edges ..."; cerr.flush();

	// Add the nodes and their t-edges
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		double l1=exp(this->obsModel->logLikelihood(this->obs->get(x,y),1,x,y)),
			   l0=exp(this->obsModel->logLikelihood(this->obs->get(x,y),0,x,y));
		double quot=l1/l0;
		double dweight=log(quot); 
			//+logf(expf(clique1GraphCut(1))/expf(clique1GraphCut(0)));
		float weight = (float) dweight;
			
		if (!finitef(weight))
			ERR_THROW ("Weight is not a finite number: l0=" << l0 << " l1="
				<< l1 << " quot=" << quot << " dweight=" << dweight 
				<< " weight=" << weight);
		
		g -> add_node(); 
			
/*		if (weight==0)
			ERR_THROW ("graphCut(): zero weight! obs=" << this->obs->get(x,y));		*/	
			
		if (weight>0)
			g -> add_tweights(NODENR(x,y), weight, 0);
		else
			g -> add_tweights(NODENR(x,y), 0, -1.*weight);
			
		// cerr << "\nweight: " << weight;
	}
		
	
	// Add the horizontal n-edges
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs-1; ++x)
		g -> add_edge(NODENR(x,y), NODENR(x+1,y), 
			PARAM2WEIGHT_GRCT(this->MRFParams[MRFP_POTTS_1_H]),
			PARAM2WEIGHT_GRCT(this->MRFParams[MRFP_POTTS_1_H]));
			
	// Add the vertical n-edges
	for (unsigned int y=0; y<ys-1; ++y)
	for (unsigned int x=0; x<xs; ++x)
		g -> add_edge(NODENR(x,y), NODENR(x,y+1), 
			PARAM2WEIGHT_GRCT(this->MRFParams[MRFP_POTTS_1_V]),
			PARAM2WEIGHT_GRCT(this->MRFParams[MRFP_POTTS_1_V]));
	
	cerr << "\nDone. Calculating max flow ..."; cerr.flush();

	flow = g -> maxflow();
	
	cerr << "\nDone. Flow = " << flow << endl;
	
	cerr << "Collecting results..."; cerr.flush();
		
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		this->lab->set(x,y,
			(g->what_segment(NODENR(x,y)) == GraphType::SOURCE)	? 1 : 0);
	}
	cerr << "Done." << endl;
	
	delete g;
}

/***********************************************************************
 * Expansion move optimization
 ***********************************************************************/

template <class TI>  
void MRFSegmenterPotts<TI>::expansionMove (unsigned int noCycles, bool fastKolmogorovVersion)
{
	float oldEnergy, newEnergy;
	bool didAccept;
	
	HRULE;
	cerr << "Starting alpha-expansion move algorithm.\n";
	HRULE;

	cerr << "Parameters:\n"
		 << " V_H(1=1) = " << clique2GraphCutH(1,1) 
		 << " V_V(1=1) = " << clique2GraphCutV(1,1) << endl
		 << " V_H(1=0) = " << clique2GraphCutH(1,0) 
		 << " V_V(1=0) = " << clique2GraphCutV(1,0) << endl
		 << " Adjusting factor (already included): " << betaAdjust << endl;
	//cerr << "SingleNode weights for label 0:" << endl;
	// for (unsigned int i=0; i<255; ++i)
	//	cerr << i << ": " << LOGLH2WEIGHT(this->obsModel->logLikelihood(i,0)) << endl;
	//this->lab->setZero();
		
	oldEnergy=totalGraphCutEnergy();
	for (unsigned int cycle=0; cycle<noCycles; ++cycle)
	{
		didAccept=false;
		for (unsigned int alpha=0; alpha<this->noClasses; ++alpha)
		{
			cerr << "expansion move: cycle=" << cycle << " alpha=" << alpha << endl;
			
			Image *labCopy = new Image (*this->lab);
			if (fastKolmogorovVersion)
				singleExpansionMoveKolmogorov(alpha);			
			else
				singleExpansionMoveBoykov(alpha);
			newEnergy=totalGraphCutEnergy();
			
			/*
			if (alpha==1 && cycle==0)
				return;
			#warning EXPANSION MOVE: ONLY ONE INTERATION
			*/
			
			if (newEnergy>=oldEnergy)
			{
				cerr << "revert.\n";
				*(this->lab) = *labCopy;
			}
			else
			{
				cerr << "accept.\n";
				oldEnergy=newEnergy;
				didAccept=true;
			}
			delete labCopy;
// 			{
// 				ostringstream name;
// 				Image tmp(*this->lab);
// 				brightenImage(tmp);
//     			name << "x_expansionmove_" << cycle << "_" << alpha << ".pgm";
//     			tmp.write (name.str().c_str());
// 			}
		}
		if (!didAccept)
		{
			cerr << "\nNo changes anymore, we are finished.\n";
			break;
		}
	}
	cerr << endl;
}


template <class TI>  
float MRFSegmenterPotts<TI>::clique2GraphCutH (unsigned char lab1, unsigned char lab2) 
{
	return PARAM2WEIGHT_EXPMOVE(
		this->MRFParams[MRFP_POTTS_1_H]*REV_KR_DELTA(lab1,lab2));
}

template <class TI>  
float MRFSegmenterPotts<TI>::clique2GraphCutV (unsigned char lab1, unsigned char lab2) 
{
	return PARAM2WEIGHT_EXPMOVE(
		this->MRFParams[MRFP_POTTS_1_V]*REV_KR_DELTA(lab1,lab2));
}

/***********************************************************************
 * Calculate the total energy of the label image
 ***********************************************************************/

template <class TI>  
float MRFSegmenterPotts<TI>::totalGraphCutEnergy()
{
	unsigned int xs=this->lab->xsize;
	unsigned int ys=this->lab->ysize;
	float e=0;
	
	// 1-cliques
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
		e += LOGLH2WEIGHT(this->obsModel->logLikelihood(
						this->obs->get(x,y),
						this->lab->get(x,y),
						x,y));

	// Horizontal 2-cliques
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs-1; ++x)		
		e += clique2GraphCutH(this->lab->get(x,y),this->lab->get(x+1,y)); 
	
	// Vertical 2-cliques
	for (unsigned int y=0; y<ys-1; ++y)
	for (unsigned int x=0; x<xs; ++x)
		e += clique2GraphCutV(this->lab->get(x,y),this->lab->get(x,y+1));
		
	cerr << "GraphCutEnergy = " << e << endl;	
		
	return e;
}

/*
 *
 * THE KOLMOGOROV VERSION OF THE EXPANSION MOVE ALGORITHM
 * FASTER THAN THE BOYKOV VERSION!
 *
 */

/***********************************************************************
 * Add an n-edge to the expansion move cut graph 
 ***********************************************************************/

template <class TI>
inline void MRFSegmenterPotts<TI>::addEdgeEMKolmogorov (GraphType *g, 
	unsigned int xs, unsigned int ys, 
	unsigned int x, unsigned int y, 
	unsigned int dx, unsigned int dy,
	unsigned char alpha, vector<float> &tweights) 
{	
	unsigned int nodenr1 = NODENR(x,y),
				 nodenr2 = NODENR(dx,dy);
	unsigned char lab1 = this->lab->get(x,y),
				  lab2 = this->lab->get(dx,dy);
	float weight;
	
#define E00	WEIGHT_C2_H(lab1,lab2)		  		  
#define E01	WEIGHT_C2_H(lab1,alpha)		  
#define E10	WEIGHT_C2_H(alpha,lab2)		  
#define E11	WEIGHT_C2_H(alpha,alpha)		  
		  
	// A horizontal clique
	if (x==xs)
	{
		weight = E01+E10-E00-E11;
		tweights[nodenr1] += E10-E00; 
		tweights[nodenr2] += E11-E10;
	}			 

#define E00	WEIGHT_C2_H(lab1,lab2)		  		  
#define E01	WEIGHT_C2_H(lab1,alpha)		  
#define E10	WEIGHT_C2_H(alpha,lab2)		  
#define E11	WEIGHT_C2_H(alpha,alpha)	
	
	// A vertical clique
	else
	{
		weight = E01+E10-E00-E11;
		tweights[nodenr1] += E10-E00; 
		tweights[nodenr2] += E11-E10;		  
	}	
		
	g -> add_edge(nodenr1, nodenr2, weight, weight);
	
	/*
	// A horizontal clique
	if (x==xs)
	{
		weight = WEIGHT_C2_H(lab1,alpha)+WEIGHT_C2_H(alpha,lab2)-
				 WEIGHT_C2_H(lab1,lab2) -WEIGHT_C2_H(alpha,alpha);
		tweights[nodenr1] += WEIGHT_C2_H(alpha,lab2)-WEIGHT_C2_H(lab1,lab2);
		tweights[nodenr2] += WEIGHT_C2_H(alpha,lab2)-WEIGHT_C2_H(alpha,alpha);			  
	}			 
	// A vertical clique
	else
	{
		weight = WEIGHT_C2_V(lab1,alpha)+WEIGHT_C2_V(alpha,lab2)-
				 WEIGHT_C2_V(lab1,lab2) -WEIGHT_C2_V(alpha,alpha);
		tweights[nodenr1] += WEIGHT_C2_V(alpha,lab2)-WEIGHT_C2_V(lab1,lab2);
		tweights[nodenr2] += WEIGHT_C2_V(alpha,lab2)-WEIGHT_C2_V(alpha,alpha);			  
	}	
	
	g -> add_edge(nodenr1, nodenr2, weight, weight);
	*/
}


template <class TI>  
float MRFSegmenterPotts<TI>::singleExpansionMoveKolmogorov (unsigned int alpha)
{
	typedef Graph<float,float,float> GraphType;
	unsigned int xs=this->lab->xsize;
	unsigned int ys=this->lab->ysize;
	float flow;
	GraphType *g;
	unsigned int curIndex;	
	vector<float> tweights(xs*ys);
	
	cerr << "Creating graph ..."; cerr.flush();
	
	g = new GraphType(2*xs*ys, 4*xs*ys);
		
	cerr << "Done.\nAdding nodes and edges ..."; cerr.flush();

	// Add the pixel nodes and the t-edges resulting from the
	// observation model
	// The t-edges are not yet added to the graph,
	// only added to a helper-structure. They will 
	// be combined with the t-edges resulting from the 
	// prior model.
	curIndex=0;
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		tweights[NODENR(x,y)] = LOGLH2WEIGHT(this->obsModel->logLikelihood(
					this->obs->get(x,y),
					alpha,x,y)) -
				  LOGLH2WEIGHT(this->obsModel->logLikelihood(
					this->obs->get(x,y),
					this->lab->get(x,y),x,y));
				  ;
				
		g -> add_node();
		++curIndex;
	}
	
	
	// Add the horizontal n-edges 
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs-1; ++x)
		addEdgeEMKolmogorov(g,xs,ys,x,y,x+1,y,alpha,tweights);
	
	// Add the vertical n-edges 
	for (unsigned int y=0; y<ys-1; ++y)
	for (unsigned int x=0; x<xs; ++x)
		addEdgeEMKolmogorov(g,xs,ys,x,y,x,y+1,alpha,tweights);
	
	// Add the t-edges
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		float weight = tweights[NODENR(x,y)];
		if (weight>0)
			g -> add_tweights(NODENR(x,y), weight, 0);
		else if (weight<0)
			g -> add_tweights(NODENR(x,y), 0, -1.*weight);		
	}
	
	cerr << "Done.\nCalculating max flow ..."; cerr.flush();

	flow = g -> maxflow();
	
	cerr << "Done.\nFlow = " << flow << endl;
	cerr << "Collecting results..."; cerr.flush();
		
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
		if (g->what_segment(NODENR(x,y)) == GraphType::SINK)	
			this->lab->set(x,y,alpha);
			
	cerr << "Done." << endl;
	delete g;	
	return flow;
}

 /*
  *
  *
  *
  * THE BOYKOV VERSION OF THE EXPANSION MOVE ALGORITHM
  * SLOW!!!!!!!!!
  *
  *
  *
  *
  */

/***********************************************************************
 * Add an n-edge to the expansion move cut graph
 ***********************************************************************/

#define PRIOR_FUNC(x,y)	(isHoriz ? (clique2GraphCutH((x),(y))) : (clique2GraphCutV((x),(y))))

template <class TI>  
void MRFSegmenterPotts<TI>::addEdgeEMBoykov(GraphType *g, 
	unsigned int nodenr1, unsigned int nodenr2,
	unsigned int lab1, unsigned int lab2,
	unsigned int alpha, unsigned int curIndex, bool isHoriz)
{
	float x; 
	if (lab1==lab2) 
	{ 
		x = PRIOR_FUNC(lab1,alpha); 
		if (x!=0) 
		{
			g -> add_edge(nodenr1, nodenr2, x, x); 
			// cerr << "n-edge: " << nodenr1 << "\t" << nodenr2
			// 	 << "\tcap: " << x << "," << x << endl;
		}
	} 
	else 
	{ 
		g->add_node(); 
		x= PRIOR_FUNC(lab1,lab2); 
		if (x==0)
		{
			ERR_THROW ("t-edge cannot be zero:"
				<< "\nlab1=" << lab1 
				<< "\nlab2=" << lab2 
				<< "\nprior=" << PRIOR_FUNC(lab1,lab2)
				<< "\nkrdelta=" << KR_DELTA(lab1,lab2)
				<< "\nweight=" << x
				<< "\nbeta-hori=" << this->MRFParams[MRFP_POTTS_1_H] 
		 		<< "\nbeta-vert=" << this->MRFParams[MRFP_POTTS_1_V]);
		}
		//cerr << "begin\n";
		if (x!=0) 
		{
			g->add_tweights (curIndex, 0, x); 
			// cerr << "t-edge: " << curIndex << "\tcap: " << 0 
			// 	 << "," << x << endl;;
		}
		x = PRIOR_FUNC(lab1,alpha);
		if (x!=0) 
		{
			g->add_edge (nodenr1,curIndex,x,x);
			// cerr << "n-edge: " << nodenr1 << "\t" << curIndex
			// 	 << "\tcap: " << x << "," << x << endl;
		}
		x = PRIOR_FUNC(alpha,lab2);
		if (x!=0) 
		{
			g->add_edge (curIndex, nodenr2, x, x);
			// cerr << "n-edge: " << curIndex << "\t" << nodenr2
			// 	 << "\tcap: " << x << "," << x << endl;
		}
		//cerr << "end\n";
		++curIndex;
	}
} 

template <class TI>  
float MRFSegmenterPotts<TI>::singleExpansionMoveBoykov (unsigned int alpha)
{
	typedef Graph<float,float,float> GraphType;
	unsigned int xs=this->lab->xsize;
	unsigned int ys=this->lab->ysize;
	float flow;
	GraphType *g;
	unsigned int curIndex;	
	
	cerr << "Creating graph ..."; cerr.flush();
	
	g = new GraphType(2*xs*ys, 4*xs*ys);
		
	cerr << "Done.\nAdding nodes and edges ..."; cerr.flush();

	// Add the pixel nodes and their t-edges
	curIndex=0;
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		float wsource, wsink;
		
		wsource = LOGLH2WEIGHT(this->obsModel->logLikelihood(
					this->obs->get(x,y),
					alpha,x,y))
				  ;//+clique1GraphCut(alpha);
				
		if (this->lab->get(x,y)==alpha)
			wsink = WEIGHT_INFINITY;
		else
			wsink = LOGLH2WEIGHT(this->obsModel->logLikelihood(
						this->obs->get(x,y),
						this->lab->get(x,y),x,y))
					;//+clique1GraphCut(this->lab->get(x,y));
		
		g -> add_node();
		g -> add_tweights(curIndex, wsource, wsink);
		//cerr << "T-edge: " << curIndex << "\tcap: " << wsource 
		//	 << "," << wsink << endl;;
		++curIndex;
	}
	
	// Add the horizontal n-edges 
	// and the corresponding auxiliary nodes
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs-1; ++x)
	{
		addEdgeEMBoykov(g, NODENR(x,y), NODENR(x+1,y), 
			this->lab->get(x,y), this->lab->get(x+1,y), 
			alpha, curIndex, true);
	}
	
	// Add the vertical n-edges 
	// and the corresponding auxiliary nodes
	for (unsigned int y=0; y<ys-1; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		addEdgeEMBoykov(g, NODENR(x,y), NODENR(x,y+1),
			this->lab->get(x,y), this->lab->get(x,y+1), 
			alpha, curIndex, false);
	}
	
	cerr << "Done.\nCalculating max flow ..."; cerr.flush();

	flow = g -> maxflow();
	
	cerr << "Done.\nFlow = " << flow << endl;
	cerr << "Collecting results..."; cerr.flush();
		
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
		if (g->what_segment(NODENR(x,y)) == GraphType::SINK)	
			this->lab->set(x,y,alpha);
			
	cerr << "Done." << endl;
	delete g;	
	return flow;
}


	
#endif
