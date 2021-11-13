
/***********************************************************************
 * Document Bleed-Through Restoration (separation recto/verso)
 * Using two Markov Random Fields.
 *
 * Author: Christian Wolf
 * Begin: 6.6.2005
 ***********************************************************************/
 
#ifndef _WOLF_MRFSEGMENTERBLEEDTHROUGH_H_
#define _WOLF_MRFSEGMENTERBLEEDTHROUGH_H_

// C++
#include <set>

// From the CONNECTED COMPONENTS ANALYSIS modules
#include <CCA.h>

// From the MARKOV RANDOM FIELDS module

#ifdef DBLMRF_WITH_POTTS_AND_GRAPHCUT
	#include <MRFPriorModels/MRFSegmenterPotts.h>
#else	
	#include <MRFPriorModels/MRFSegmenterAL8N.h>
#endif	

// From this module
#include "CommonBleedThrough.h"
#include "VersoReplacementPyramid.h"
#include "BTKMeansClusterer.h"

#ifdef DBLMRF_WITH_POTTS_AND_GRAPHCUT	
	#define U2pr(f1,f2,o)	((-1.0*this->obsModel->logLikelihood((o), getClass((f1),(f2)))) \
						 -1.0*(f1*f2)*versoConstraintParam)	
	#define BETAA(x)	((x)*this->betaAdjust)
#endif

	
template <class TI>
#ifdef DBLMRF_WITH_POTTS_AND_GRAPHCUT
class MRFSegmenterBleedThrough : public MRFSegmenterPotts<TI>
#else
class MRFSegmenterBleedThrough : public MRFSegmenterAL8N<TI>
#endif
{

	public:
	
		// Constructor and Destructor
		MRFSegmenterBleedThrough(TI *o, Image *lr, Image *lv, ObsModel<TI> *model,
			TI *xEventualColorObs,
			char *xForceMRFParams, bool xUseMask, bool xUseGV4RectoRecognition, 
			bool xOptDoGaussFilterBeforeKmeans,
			float xBetaAdjust, float xVersoConstraintParamFactor,
			ColorSpaceCode xColorSpaceCode);
		~MRFSegmenterBleedThrough();
		
		// Redefinition of virtual methods
		virtual void _init();	
		
		float priorEnergy(int x, int y);
		bool estimateMRFParamsLS(bool doMedianFilter, bool singleCliquesManually);
		unsigned int getCliqueArrIndexCenterPixel();
				
		void graphCutAdaptedToProblem(unsigned int noIts) { expansionMove(noIts); }
		void expansionMove (unsigned int noCycles);
		unsigned int gibbsSampler (double T, bool greedy);
					
		void restoreFromSegmentation();
		
		// Accessors
		Image *getLabv()	{ return labv; }
		
		/*
		void simulatedAnnealing (double T1, double C, int noIts, int noItsGreedy,
			char *debugString);
		*/
		
		void writeLabelField(char *filename, DebugWriteImageOptions options);
		
	protected:	// methods	
	
		void singleExpansionMove (bool rectoIsFixed);
		inline float calcTEdgeFromNeighbors(unsigned int x, unsigned int y,
			bool rectoIsFixed, Matrix<unsigned char> &H);
		float totalGraphCutEnergy();
			
		
	protected:	// data
	
		
		Image *labv,				// The label field for the verso side
			  *labvModif;			// The label field for the verso side, currently modified		
		
		// The MRF Parameters for the verso side
		unsigned int noMRFParamsv;
		Vector<float> MRFParamsv;	
		
 		bool useMask, useGV4RectoRecognition, optDoGaussFilterBeforeKmeans;
		unsigned char maskLabelRecto, maskLabelVerso;
		
		float versoConstraintParam, versoConstraintParamFactor;
};

/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>
MRFSegmenterBleedThrough<TI>::MRFSegmenterBleedThrough(TI *o, Image *lr, Image *lv, 
	ObsModel<TI> *m, TI *xEventualColorObs,
	char *xForceMRFParams, bool xUseMask, bool xUseGV4RectoRecognition,
	bool xOptDoGaussFilterBeforeKmeans, float xBetaAdjust, 
	float xVersoConstraintParamFactor, ColorSpaceCode xColorSpaceCode)
#ifdef DBLMRF_WITH_POTTS_AND_GRAPHCUT	
	: MRFSegmenterPotts<TI> (o,lr,m,xEventualColorObs, 
	xForceMRFParams,2,xOptDoGaussFilterBeforeKmeans,xBetaAdjust,xColorSpaceCode)
#else
	: MRFSegmenterAL8N<TI> (o,lr,m,xForceMRFParams)
#endif		
{
	cerr << "New model: MRFSegmenterBleedThrough" << endl;
	// Init the verso side
	labv = lv;
	MRFParamsv.resize(this->noMRFParams);
	labvModif = new Image (lv->xsize, lv->ysize, 1);
	useMask = xUseMask;	
	useGV4RectoRecognition = xUseGV4RectoRecognition;
	optDoGaussFilterBeforeKmeans = xOptDoGaussFilterBeforeKmeans;
	versoConstraintParamFactor = xVersoConstraintParamFactor;
}

/***********************************************************************
 * Destructor
 ***********************************************************************/

template <class TI>
MRFSegmenterBleedThrough<TI>::~MRFSegmenterBleedThrough()
{
	delete labvModif;
}

/***********************************************************************
 * A couple of virtual methods we need to redefine
 ***********************************************************************/

template <class TI>
unsigned int MRFSegmenterBleedThrough<TI>::getCliqueArrIndexCenterPixel()
{
	ERR_THROW ("MRFSegmenterBleedThrough<TI>::getCliqueArrIndexCenterPixel():"
		"does not make sense for 2 labels!");
}

/***********************************************************************
 * Calculate the prior probability 
 ***********************************************************************/
 
template <class TI>  
inline float MRFSegmenterBleedThrough<TI>::priorEnergy (int x, int y) 
{
	float rv;
	
	// The priors coming from the two independent MRFs
	rv = this->_priorEnergy (this->lab,  this->MRFParams,  x, y)+
		 this->_priorEnergy (this->labv, this->MRFParamsv, x, y);	  
		 
	// The prior from the joined information
	// it is necessary to constrain the verso pixels which are 
	// hidden from the recto text pixels.
	if (this->lab->get(x,y) && this->labv->get(x,y))
		rv += versoConstraintParam;
	
	return rv;
}

/***********************************************************************
 * Estimate the Parameters
 * Estimate them on the recto field only and then adapt them
 * for the verso side
 * returns true if the estimation was successful
 ***********************************************************************/

template <class T>  
bool MRFSegmenterBleedThrough<T>::estimateMRFParamsLS(bool doMedianFilter, bool singleCliquesManually)
{
	if (!_estimateMRFParamsLS (this->lab,  this->MRFParams, doMedianFilter, singleCliquesManually))
		return false;
	for (unsigned int i=0; i<this->noMRFParams; ++i)
		this->MRFParamsv[i] = this->MRFParams[i];	
	_mirrorParamsX(this->MRFParamsv);
	
	versoConstraintParam=versoConstraintParamFactor*this->_sumOfPairwiseParams(this->MRFParamsv);	
	cerr << "versoConstraint=x*" << versoConstraintParamFactor 
	 	 << "=" << versoConstraintParam << endl;
	return true;
}

/***********************************************************************
 * Calculate a specific contribution to certain t-edge weights
 * for the graph cut optimization
 ***********************************************************************/
 
#define ADD_WEIGHT(nx,ny,beta) \
{	unsigned int f1s=(rectoIsFixed?this->lab->get(nx,ny):this->labv->get(nx,ny)); \
	if  (H(ny,nx)==0) \
		weight += (f1s==0?-1.:1.)*beta; \
}		

template <typename TI>
inline float MRFSegmenterBleedThrough<TI>::calcTEdgeFromNeighbors(
	unsigned int x, unsigned int y, bool rectoIsFixed, Matrix<unsigned char> &H)
{
	float weight;
	unsigned int xs=this->lab->xsize;
	unsigned int ys=this->lab->ysize;		
	Vector<float> *p = (rectoIsFixed?(&this->MRFParams):(&this->MRFParamsv));
	
#ifdef DBLMRF_WITH_POTTS_AND_GRAPHCUT	
	weight=0;
	if (x>0) 	ADD_WEIGHT(x-1,y,BETAA((*p)[MRFP_POTTS_1_H]));
	if (x<xs-1)	ADD_WEIGHT(x+1,y,BETAA((*p)[MRFP_POTTS_1_H]));
	if (y>0) 	ADD_WEIGHT(x,y-1,BETAA((*p)[MRFP_POTTS_1_V]));	
	if (y<ys-1)	ADD_WEIGHT(x,y+1,BETAA((*p)[MRFP_POTTS_1_V]));
	return weight;
#else
	ERR_THROW ("Internal error in MRFSegmenterBleedThrough<TI>::calcTEdgeFromNeighbors!");
#endif
}

/***********************************************************************
 * Graph cut optimization
 *
 * One of the fields is "fixed", i.e. the variables corresponding to 
 * nonregular sites are fixed and not estimated. The other field
 * is estimated.
 ***********************************************************************/
 
template <class T>  
void MRFSegmenterBleedThrough<T>::singleExpansionMove (bool rectoIsFixed)
{		
	unsigned int xs=this->lab->xsize;
	unsigned int ys=this->lab->ysize;
	float flow, alphaFix, alphaEst;
	GraphType *g;
	Matrix<unsigned char> H (this->lab->ysize, this->lab->xsize);
	Matrix<unsigned int> FixedNodeIndex (this->lab->ysize, this->lab->xsize);
	unsigned int regularCount, curIndex;
	Vector<float> *paramFix, *paramEst;
		
#ifdef DBLMRF_WITH_POTTS_AND_GRAPHCUT				
	if (rectoIsFixed)
	{		
		paramFix = &this->MRFParams;
		paramEst = &this->MRFParamsv;		
	}
	else
	{		
		paramFix = &this->MRFParamsv;
		paramEst = &this->MRFParams;		
	}
	alphaEst = (*paramEst)[MRFP_POTTS_FIRST_SINGLE];
	alphaFix = (*paramFix)[MRFP_POTTS_FIRST_SINGLE];

	cerr << "Searching regular sites ... "; cerr.flush();
	regularCount=0;
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		unsigned char o = this->obs->get(x,y);		
		/*
		if (o>=this->obsModel->getProperty(PROP_BG_MEAN))
		{
			cerr <<  "U2pr(0,0,o)=" << U2pr(0,0,o) << endl
				<<  "U2pr(1,1,o)=" << U2pr(1,1,o) << endl
				<<  "U2pr(0,1,o)=" << U2pr(0,1,o) << endl
				<<  "U2pr(1,0,o)=" << U2pr(1,0,o) << endl
				<<  "obs=" << (int) o << endl;
			exit(1);
		}
		*/
		if ((U2pr(0,0,o)+U2pr(1,1,o))<=(U2pr(1,0,o)+U2pr(0,1,o)))
		{
			H(y,x)=1;
			++regularCount;
		}
		else
			H(y,x)=0;
	}

	cerr << "Done.\nNumber of regular sites: " << regularCount
		 << " (" << 100.0*(float)regularCount/(float)(xs*ys) 
		 << "%).\nCreating graph ... "; cerr.flush();


#warning ONLY NORMAL ALPHA EXPANSION MOVE!
	cerr << "Done...\nONLY NORMAL ALPHA EXPANSION MOVE!";
	H.setZero();	
	
	g = new GraphType((int) rint(1.5*xs*ys), 6*xs*ys);
		
	cerr << "Done.\nAdding nodes ... "; cerr.flush();
	
	// Add the nodes of the estimated field 
	curIndex=0;
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{				
		g->add_node();
		++curIndex; 
	}
	
	// Add the nodes of the fixed field 
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{				
		if (H(y,x)==1)
		{
			g->add_node(); 
			FixedNodeIndex(y,x)=curIndex;
			++curIndex;
		}
	}
	
	cerr << "Done.\nAdding t-edges ... "; cerr.flush();	

	// Add the t-edges
	curIndex=0;
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		float weightEstNode,weightFixNode;
		unsigned char o = this->obs->get(x,y);		
		
		weightEstNode = -1.0*alphaEst;
		weightFixNode = -1.0*alphaFix;
#warning NO PIXEL PRIORS!!!		
		weightEstNode = 0;
		weightFixNode = 0;

		if (H(y,x)==0)
		{			
			if (rectoIsFixed)
				weightEstNode +=
					(U2pr(this->lab->get(x,y),1,o) - U2pr(this->lab->get(x,y),0,o));
			else			
				weightEstNode += 
					(U2pr(1,this->labv->get(x,y),o) - U2pr(0,this->labv->get(x,y),o));
		}
		
		// H(y,x)==1
		else
		{	
			weightEstNode += (U2pr(1,1,o) - U2pr(1,0,o));			
			weightFixNode += (U2pr(1,0,o) - U2pr(0,0,o));
			weightFixNode += calcTEdgeFromNeighbors(x,y,rectoIsFixed,H);
		}		
			
		if (weightEstNode>0)
			g -> add_tweights(curIndex, weightEstNode, 0);
		if (weightEstNode<0)
			g -> add_tweights(curIndex, 0, -1.*weightEstNode);
			
		if (H(y,x)==1)
		{
			if (weightFixNode>0)
				g -> add_tweights(FixedNodeIndex(y,x), weightFixNode, 0);
			if (weightFixNode<0)
				g -> add_tweights(FixedNodeIndex(y,x), 0, -1.*weightFixNode);
		}
			
		++curIndex;			
	}	
			
	// Add the horizontal n-edges
	cerr << "Done.\nAdding horizontal n-edges (w=" 
	     << -1.0*BETAA((*paramEst)[MRFP_POTTS_1_H]) << ")... "; cerr.flush();	
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs-1; ++x)
	{
		g -> add_edge(NODENR(x,y), NODENR(x+1,y), 
			-1.0*BETAA((*paramEst)[MRFP_POTTS_1_H]),
			-1.0*BETAA((*paramEst)[MRFP_POTTS_1_H]));
			
		if ((H(y,x)==1) && (H(y,x+1)==1))
			g -> add_edge(FixedNodeIndex(y,x), FixedNodeIndex(y,x+1), 
				-1.0*BETAA((*paramFix)[MRFP_POTTS_1_H]),
				-1.0*BETAA((*paramFix)[MRFP_POTTS_1_H]));
	}
			
	// Add the vertical n-edges
	cerr << "Done.\nAdding vertical n-edges (w=" 
	     << -1.0*BETAA((*paramEst)[MRFP_POTTS_1_V]) << ")... "; cerr.flush();		     
	for (unsigned int y=0; y<ys-1; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		g -> add_edge(NODENR(x,y), NODENR(x,y+1), 
			-1.0*BETAA((*paramEst)[MRFP_POTTS_1_V]),
			-1.0*BETAA((*paramEst)[MRFP_POTTS_1_V]));
			
		if ((H(y,x)==1) && (H(y+1,x)==1))
			g -> add_edge(FixedNodeIndex(y,x), FixedNodeIndex(y+1,x), 
				-1.0*BETAA((*paramFix)[MRFP_POTTS_1_V]),
				-1.0*BETAA((*paramFix)[MRFP_POTTS_1_V]));	
	}
	
	// Add the n-edges between the two fields
	cerr << "Done.\nAdding inter-field n-edges ... "; cerr.flush();	
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		if (H(y,x)==1)
		{
			unsigned char o = this->obs->get(x,y);
			float weight = U2pr(0,1,o)+U2pr(1,0,o)-U2pr(0,0,o)-U2pr(1,1,o);
			if (weight<0)
			{
				ERR_THROW ("Non-regular weight for pixel marked as regular!\n"
					 <<  "U2pr(0,0,o)=" << U2pr(0,0,o) << endl
					 <<  "U2pr(1,1,o)=" << U2pr(1,1,o) << endl
					 <<  "U2pr(0,1,o)=" << U2pr(0,1,o) << endl
					 <<  "U2pr(1,0,o)=" << U2pr(1,0,o) << endl
					 <<  "obs=" << (int) o);
			}
			g -> add_edge(NODENR(x,y), FixedNodeIndex(y,x), weight, weight);
		}		
	}
	
	cerr << "Done.\nCalculating max flow ... "; cerr.flush();
	flow = g -> maxflow();	
	cerr << "\nDone. Flow = " << flow << endl;
	
	cerr << "Collecting results..."; cerr.flush();
		
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		if (rectoIsFixed)
			this->labv->set(x,y,
				(g->what_segment(NODENR(x,y)) == GraphType::SOURCE)	? 0 : 1);
		else
			this->lab->set(x,y,
				(g->what_segment(NODENR(x,y)) == GraphType::SOURCE)	? 0 : 1);
			
		if (H(y,x)==1)
		{
			if (rectoIsFixed)
				this->lab->set(x,y,
					(g->what_segment(FixedNodeIndex(y,x)) == GraphType::SOURCE)	? 0 : 1);
			else
				this->labv->set(x,y,
					(g->what_segment(FixedNodeIndex(y,x)) == GraphType::SOURCE)	? 0 : 1);			
		}
	}
	cerr << "Done.\n";
	
	delete g;
#else
	ERR_THROW("Graphcut optimization is not supported if the bleed-through class inherits "
		"from the autologistic model! define DBLMRF_WITH_POTTS_AND_GRAPHCUT to change this."); 	
#endif	
}


/***********************************************************************
 * The improved expansion move optimization:
 * We alternate, always keep one field fixed and estimate the other field.
 * However, not all variables of the "fixed" field are really fixed,
 * according to the observation, some variables of the fixed field are
 * estimated as well.
 ***********************************************************************/

template <class TI>  
void MRFSegmenterBleedThrough<TI>::expansionMove (unsigned int noCycles)
{
	float oldEnergy, newEnergy;
	unsigned int revertsInARow;
	
	HRULE;
	cerr << "Starting improved alpha-expansion move algorithm.\n";
	HRULE;
	
	{
		ostringstream name;				
		name << "x_expansionmove_ICE" << this->curICEIteration
				<< "_cy0_init.pgm";
		writeLabelField((char *) name.str().c_str(), DBG_WIMAGES_SEG);
	}
	
	cerr << "Parameters (not adjusted, not yet taken negative!!):\n"
		 << " alpha1  = " << this->MRFParams[MRFP_POTTS_FIRST_SINGLE] << endl
		 << " beta1_h = " << this->MRFParams[MRFP_POTTS_1_H] << endl
		 << " beta1_v = " << this->MRFParams[MRFP_POTTS_1_V] << endl
		 << " alpha2  = " << this->MRFParamsv[MRFP_POTTS_FIRST_SINGLE] << endl
		 << " beta2_h = " << this->MRFParamsv[MRFP_POTTS_1_H] << endl
		 << " beta2_v = " << this->MRFParamsv[MRFP_POTTS_1_V] << endl
		 << " Adjusting factor: " << this->betaAdjust << endl;
		
	oldEnergy=totalGraphCutEnergy();
	revertsInARow=0;
	for (unsigned int cycle=0; cycle<2*noCycles; ++cycle)
	{		
		bool rectoIsFixed=((cycle%2)==1);
	
		cerr << "expansion move: cycle=" << cycle;
		if (rectoIsFixed)
			cerr << " recto=fixed, verso=estimated.\n";
		else
			cerr << " recto=estimated, verso=fixed.\n";
		
		Image *labCopy = new Image (*this->lab);
		singleExpansionMove(rectoIsFixed);			
		newEnergy=totalGraphCutEnergy();
		
		if (newEnergy>=oldEnergy)
		{
			cerr << "revert.\n";
			*(this->lab) = *labCopy;
			++revertsInARow;
		}
		else
		{
			cerr << "accept.\n";
		/*
		{
			cerr << "ALWAYS ACCEPTING!\n";
		*/	
			oldEnergy=newEnergy;
			revertsInARow=0;
		}
		delete labCopy;
		{
			ostringstream name;				
			name << "x_expansionmove_ICE" << this->curICEIteration
				 << "_cy" << cycle << ".pgm";
			writeLabelField((char *) name.str().c_str(), DBG_WIMAGES_SEG);
		}
		if (revertsInARow>1)
		{
			cerr << "No changes anymore, we are finished.\n";
			break;
		}
	}	
}


/***********************************************************************
 * Calculate the total energy of the label image
 ***********************************************************************/

template <class TI>  
float MRFSegmenterBleedThrough<TI>::totalGraphCutEnergy()
{
	unsigned int xs=this->lab->xsize;
	unsigned int ys=this->lab->ysize;
	float e=0;
	
	// 1-cliques
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		e += U2pr(this->lab->get(x,y),this->labv->get(x,y),this->obs->get(x,y));
#warning NO PIXEL PRIORS!!!
		/*		
		e += this->lab->get(x,y) *this->MRFParams[MRFP_POTTS_FIRST_SINGLE];
		e += this->labv->get(x,y)*this->MRFParamsv[MRFP_POTTS_FIRST_SINGLE];
		*/
	}

	// Horizontal 2-cliques
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs-1; ++x)
	{
		e += -1.0*BETAA(this->MRFParams[MRFP_POTTS_1_H])*
			REV_KR_DELTA(this->lab->get(x,y),this->lab->get(x+1,y));
		e += -1.0*BETAA(this->MRFParamsv[MRFP_POTTS_1_H])*
			REV_KR_DELTA(this->labv->get(x,y),this->labv->get(x+1,y)); 		
	}
	
	// Vertical 2-cliques
	for (unsigned int y=0; y<ys-1; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		e += -1.0*BETAA(this->MRFParams[MRFP_POTTS_1_V])*
			REV_KR_DELTA(this->lab->get(x,y),this->lab->get(x,y+1));
		e += -1.0*BETAA(this->MRFParamsv[MRFP_POTTS_1_V])*
			REV_KR_DELTA(this->labv->get(x,y),this->labv->get(x,y+1));
	}	
		
	cerr << "GraphCutEnergy = " << e << endl;	
		
	return e;
}

// **********************************************************************
// One iteration of the gibbs sampler
// Needs to be redefined since we have now two label images
// **********************************************************************

#define	PRINT_ENERGY(l,k)		{ this->lab->set (x, y, (l)); this->labv->set(x, y, (k)); \
								  float p=this->priorEnergy(x,y), c=this->condEnergy(x,y); \
								  cerr << "  |  " << (int) l << "-" << (int) k << ":  prior=" \
								  	   <<  p << " cond=" <<  c << " sum=" << p+c << endl; }

enum {
	STATE_UNCHANGED=0,
	STATE_CHANGE1,
	STATE_CHANGE2,
	STATE_CHANGE12,
	
	STATE_COUNT	
};

template <typename TI>
unsigned int MRFSegmenterBleedThrough<TI>::gibbsSampler (double T, bool greedy) 
{
	unsigned int nochangedpixels=0;	
	
	if (this->noClasses>2)
		ERR_THROW ("MRFSegmenterBleedThrough<TI>::gibbsSampler(): only works for 2 labels!");
	
	// ------------------------------------------------------------------------
	// We are in "greedy" mode, 
	// This corresponds to the ICM=Iterated Conditional Modes Algorithm 
	// proposed by Besag
	// ------------------------------------------------------------------------
	
	if (greedy)
	{
		// A full sweep of the image
		for (int y=1; y<this->lab->ysize-1; ++y) 
		for (int x=1; x<this->lab->xsize-1; ++x) 
		{
			double energy[STATE_COUNT], min;
			byte curvalue1, curvalue2,
				 newvalue1, newvalue2;
			int argMin;
			
			if (this->mask!=NULL && this->mask->get(x,y)==0)
			{
				this->labModif->set (x,y,maskLabelRecto);
				this->labvModif->set(x,y,maskLabelVerso);
				continue;
			}
			
			// ------------------------------------------------------------------------
			// DETAILED DEBUG OUTPUT: all energies for the 4 different possible states
			// ------------------------------------------------------------------------
			
			DEBUGIFPX(x,y)
			{
				curvalue1 = this->lab->get(x, y);		
				curvalue2 = this->labv->get(x, y);
				
				HRULE;
				cerr << DEBUGSTRPX 
					 << "observed gray value: " << (int) this->obs->get(x,y) << endl
					 << DEBUGSTRPX 
					 << "unch. state=" << (int) this->lab->get(x, y) 
					 << "-" << (int) this->labv->get(x, y) 
					 << endl;
					 
				PRINT_ENERGY(0,0);
				PRINT_ENERGY(0,1);
				PRINT_ENERGY(1,0);
				PRINT_ENERGY(1,1);
				HRULE;
																     
				this->lab->set(x, y, curvalue1);		
				this->labv->set(x, y, curvalue2);				     
			}			
			
			curvalue1 = this->lab->get(x, y);		
			curvalue2 = this->labv->get(x, y);					
			newvalue1 = (curvalue1 > 0 ? 0 : 1);
			newvalue2 = (curvalue2 > 0 ? 0 : 1);									

			// Calculate the current energy
			// as well as the energy for all three possible state changes
			energy[STATE_UNCHANGED] = this->priorEnergy(x,y) + this->condEnergy(x,y);	
			this->lab->set (x, y,  newvalue1);
			energy[STATE_CHANGE1]   = this->priorEnergy(x,y) + this->condEnergy(x,y);			
			this->labv->set (x, y, newvalue2);
			energy[STATE_CHANGE12]  = this->priorEnergy(x,y) + this->condEnergy(x,y);			
			this->lab->set (x, y,  curvalue1);
			energy[STATE_CHANGE2]   = this->priorEnergy(x,y) + this->condEnergy(x,y);
						
			// chose the best one
			min=energy[0];
			argMin=0;
			for (unsigned int i=1; i<STATE_COUNT; ++i)
			{
				if (energy[i]<min)
				{
					min=energy[i];
					argMin=i;
				}
			}			
			
			// Choose the best state
			switch (argMin)
			{		
				case STATE_UNCHANGED:
					this->labModif->set (x, y, curvalue1);
					this->labvModif->set(x, y, curvalue2);
					break;
				case STATE_CHANGE1:				
					this->labModif->set (x, y, newvalue1);
					this->labvModif->set(x, y, curvalue2);
					++nochangedpixels;
					break;
				case STATE_CHANGE2:
					this->labModif->set (x, y, curvalue1);
					this->labvModif->set(x, y, newvalue2);
					++nochangedpixels;
					break;
				case STATE_CHANGE12:
					this->labModif->set (x, y, newvalue1);
					this->labvModif->set(x, y, newvalue2);
					++nochangedpixels;
					break;
#ifdef PEDANTIC_CHECK_CODE					
				default:
					throw EError ("MRFSegmenterBleedThrough::gibbsSampler(): internal error(1)!\n");
#endif							
			}
				
			// Restore the old values. Necessary since we 
			// perform the modifications on a copy only!!
			// lab has already been done.	
			this->labv->set (x, y, curvalue2);
			
			DEBUGIFPX(x,y)
				cerr << endl;
		}		
	}
	
	// ------------------------------------------------------------------------
	// We are in the gibbs sampler mode
	// proposed by Geman & Geman
	// ------------------------------------------------------------------------
	
	else
	{
		// A full sweep of the image
		for (int y=1; y<this->lab->ysize-1; ++y) 
		for (int x=1; x<this->lab->xsize-1; ++x) 
		{
			double en_unchanged, en_changed, q;
			byte curvalue1, curvalue2,
				 newvalue1, newvalue2;
			bool takeNewValues;
			byte rb;
			
			if (this->mask!=NULL && this->mask->get(x,y)==0)
			{
				this->labModif->set (x,y,maskLabelRecto);
				this->labvModif->set(x,y,maskLabelVerso);
				continue;
			}
			
			// ------------------------------------------------------------------------
			// DETAILED DEBUG OUTPUT: all energies for the 4 different possible states
			// ------------------------------------------------------------------------
			
			DEBUGIFPX(x,y)
			{
				curvalue1 = this->lab->get(x, y);		
				curvalue2 = this->labv->get(x, y);
				
				HRULE;
				cerr << DEBUGSTRPX 
					 << "observed gray value: " << (int) this->obs->get(x,y) << endl
					 << DEBUGSTRPX 
					 << "unch. state=" << (int) this->lab->get(x, y) 
					 << "-" << (int) this->labv->get(x, y) 
					 << endl;
				HRULE;
					 
				PRINT_ENERGY(0,0);
				PRINT_ENERGY(0,1);
				PRINT_ENERGY(1,0);
				PRINT_ENERGY(1,1);
																     
				this->lab->set(x, y, curvalue1);		
				this->labv->set(x, y, curvalue2);				     
			}
			
			// ------------------------------------------------------------------------
	
			// Calculate the current energy
			en_unchanged = this->priorEnergy(x,y) + this->condEnergy(x,y);						
			
			// Chose randomly one of three new possible configurations
			// and set the label fields to it
			rb = (this->randgen->nextByte())%3;	
			curvalue1 = this->lab->get(x, y);		
			curvalue2 = this->labv->get(x, y);
			newvalue1 = curvalue1;
			newvalue2 = curvalue2;					
			switch (rb)
			{
				case 0:				
					newvalue1 = (curvalue1 > 0 ? 0 : 1);
					this->lab->set (x, y, newvalue1);						
					break;
				case 1:		
					newvalue2 = (curvalue2 > 0 ? 0 : 1);
					this->labv->set(x, y, newvalue2);
					break;
				case 2:			
					newvalue1 = (curvalue1 > 0 ? 0 : 1);
					newvalue2 = (curvalue2 > 0 ? 0 : 1);					
					this->lab->set (x, y, newvalue1);							
					this->labv->set(x, y, newvalue2);
					break;
#ifdef PEDANTIC_CHECK_CODE					
				default:
					throw EError ("MRFSegmenterBleedThrough::gibbsSampler(): internal error(2)!\n");
#endif					
			}
			
			// calculate the energy of the changed fields
			en_changed = this->priorEnergy(x,y) + this->condEnergy(x,y);			
	
			// The energy change
			q = exp(-1.0*(en_changed-en_unchanged)/T);
			
			// Energy changes to worse, so we should not change it,
			// but with a certain probability we do it anyway
			if (q<=1) 
				takeNewValues = (this->randgen->next() <= q);
			else
				takeNewValues=true;
			
			// accept the configuration
			if (takeNewValues)
			{
				this->labModif->set (x,y, newvalue1);
				this->labvModif->set(x,y, newvalue2);					
			}
			
			// revert the configuration
			else
			{
				this->labModif->set (x,y, curvalue1);
				this->labvModif->set(x,y, curvalue2);					
			}
			
			// Restore the old values. Necessary since we 
			// perform the modifications on a copy only!!
			this->lab->set  (x, y, curvalue1);
			this->labv->set (x, y, curvalue2);
			
			DEBUGIFPX(x,y)
				cerr << endl;
		}	
	}
	
	// Apply the changes
	for (int y=1; y<this->lab->ysize-1; ++y) 
	for (int x=1; x<this->lab->xsize-1; ++x)
	{ 
		this->lab->set(x,y,  this->labModif->get(x,y));
		this->labv->set(x,y, this->labvModif->get(x,y));
	}
		
	return nochangedpixels;	
}
#undef PRINT_ENERGY

/***********************************************************************
 * Initialise the two label images
 ***********************************************************************/
 
template <class TI>
void MRFSegmenterBleedThrough<TI>::_init()
{
	typedef typename TI::PixelType PixelType;	
	Image *chosenObs = ((this->eventualColorObs==NULL) ? this->obs : this->eventualColorObs);
	unsigned int dim=chosenObs->nbColorPlanes();
	unsigned int xs=chosenObs->xsize;
	unsigned int ys=chosenObs->ysize;	
	unsigned int pixnum;	
	Vector<float> v(dim);		
	BTKMeansClusterer<PixelType> kmc(dim, NO_CLASSES);
	unsigned int bgLabel;
	Image *filtered;
	LabCie94Distance *customDistance;
	double pl, pu, pv;
	
	// perform K-means clustering
	HRULE; cerr << "K-Means Initialisation (MRFSegmenterBleedThrough):" << endl;
		
	if (optDoGaussFilterBeforeKmeans)
	{
		cerr << "Low pass filtering image for k-means input.\n";
		filtered = new Image (*chosenObs);	
		FilterMask fm;
		fm  << 1.0 <<  4.0  << 7.0 <<  4.0 <<  1.0 << '\n'
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
		ERR_THROW ("Incompatible colorspace and dimension in MRFSegmenterBleedThrough<TI>::_init()\n");

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
	
	// Initialise the cluster centers	
	kmc.initZero();
	kmc.init((unsigned int) 0, (PixelType) KMEANS_INIT_BG);
	kmc.init((unsigned int) 1, (PixelType) KMEANS_INIT_RECTO);
	kmc.init((unsigned int) 2, (PixelType) KMEANS_INIT_VERSO);		
	
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
		
		
	/*	
	#warning Fisher initisalisation
	{
		HRULE; cerr << "Fisher/Otsu initialisation (2 thresholds):" << endl;
		Histogram<unsigned int> *H;
		int t1,t2;
	
		H = buildHistogram<unsigned int> (*(chosenObs), false, 0);
		H->getThresholdValueFisher2th (0, 256, t1,t2);
		delete H;
		
		for (unsigned int y=0; y<ys; ++y)
		for (unsigned int x=0; x<xs; ++x)
		{
			unsigned char v=chosenObs->get(x,y);
			if (v<=t1)
				this->lab->set(x,y,0);
			else
				if (v<=t2)
					this->lab->set(x,y,1);
				else
					this->lab->set(x,y,2);
		}	
	}	
	*/	
			
	// Debug: print the number of pixels per class
	{	Vector<unsigned int> count (NO_CLASSES);
		count.setZero();	
		for (int y=0; y<this->lab->ysize; ++y) 
		for (int x=0; x<this->lab->xsize; ++x) 	
			++count[this->lab->get(x,y)];
		for (unsigned int l=0; l<NO_CLASSES; ++l) 
			cerr << "class " << l << ": " << count[l] << " pixels." << endl;
	}
	
	// Get the background label
	bgLabel = getBGLabel(this->lab);
	cerr << "Background is label " << bgLabel << endl;
		
	if (useMask)
	{
		// Create the mask
		if (this->mask==NULL)
			this->mask = new Image (xs, ys, 1);
		maskLabelRecto=0;
		maskLabelVerso=0;
							
		kmc.calcMaskLabels(bgLabel, BLDTHR_MASK_THR);
		
		// Collect the results ->
		// initialise the mask
		pixnum=0;	
		for (unsigned int y=0; y<ys; ++y)
		for (unsigned int x=0; x<xs; ++x)
		{
			// Initialise the labeling
			this->mask->set(x,y, kmc(pixnum));
			++pixnum;
		}	
		dilate (*(this->mask), 1, false, BTMASK_FG); 
	}
	
	reOrderLabels(this->lab, bgLabel);

 	// Debug output
	{	Image tmp (2*this->lab->xsize,(useMask?2:1)*this->lab->ysize);
		tmp.setPlaneToValue(PLANE_RED,0);
				
		// UPPER RIGHT - labels
		// The recto and verso labels are NOT correctly ordered!!!!
		{
			Image m (*(this->lab));
			brightenImage(m);
			tmp.paste (m,this->lab->xsize,0);
		}
		// LOWER RIGHT - mask
		if (useMask)
		{
			Image m (*(this->mask));
			brightenImage(m);
			tmp.paste (m,0,this->lab->ysize);
		}
		brightenImage (tmp);
		reverseVideo (tmp);
		// UPPER LEFT		
		tmp.paste(*(chosenObs),0,0);		
		tmp.write ("x_005_kmeans_reordered_rev.pgm");
	}
	
	// Debug output
	/*
	if (this->mask!=NULL)
	{
		Histogram<float> *h = new Histogram<float> (256, 0, 255, true);	
		h->clear();	
		for (unsigned int y=0; y<ys; ++y)
		for (unsigned int x=0; x<xs; ++x)
		{
			if (this->mask->get(x,y)!=0)
				++(h->bins[chosenObs->get(x,y)]);
		}		
		h->normalize();
		cerr << "Masked histogram:" << endl;
		cerr << *h;
		delete h;
   	}
	*/
	recto2RectoVerso(this->lab, this->labv, chosenObs, useGV4RectoRecognition, -1);
		
	// Clean up
	if (optDoGaussFilterBeforeKmeans)
		delete filtered;	
	if (this->colspc==COLSPC_LAB_CIE94)
		delete customDistance;
}

/***********************************************************************
 * Create the restored image
 ***********************************************************************/

template <class TI>
void MRFSegmenterBleedThrough<TI>::restoreFromSegmentation()
{
	Image tmp (this->lab->xsize, this->lab->ysize, 1);
	for (int y=0; y<this->lab->ysize; ++y)
	for (int x=0; x<this->lab->xsize; ++x)
		tmp.set(x,y,getClass(this->lab->get(x,y),this->labv->get(x,y)));
		
	VersoReplacementPyramid<TI> p (this->obs, &tmp, CL_IND_BG, CL_IND_VERSO,
		6, 4, (unsigned char) this->obsModel->getProperty(PROP_BG_MEAN));	
	p.restore();
	
	/* Debug
	{ 	Image tmp (this->lab->xsize,this->lab->ysize);	
		for (int y=0; y<tmp.ysize; ++y)
		for (int x=0; x<tmp.xsize; ++x)		
		{
			switch (getClass(this->lab->get(x,y),this->labv->get(x,y)))
			{
				case CL_IND_BG: 
					tmp.set(x,y,197);
					break;		
				case CL_IND_RECTO: 
					tmp.set(x,y,57);
					break;
				case CL_IND_VERSO: 
					tmp.set(x,y,96);
					break;
			}
		}	
		tmp.write("x_rest_tmp_a.pgm");
		p.debugOutput("x_rest_pyr");
	} */
}

/***********************************************************************
 * Write the label field into a file
 * This version is used to compare the labels with a ground truth field
 * See also the previous method ...
 ***********************************************************************/
 
template <class TI>
void MRFSegmenterBleedThrough<TI>::writeLabelField(char *filename, DebugWriteImageOptions options)
{
	Image *tmp;	
	switch (options)
	{
		case DBG_WIMAGES_NO:
			break;
		
		case DBG_WIMAGES_SEG:	
			tmp = this->obsModel->visualizeSegmentation();		
			tmp->write(filename);
			delete tmp;
			break;
	
		case DBG_WIMAGES_LABELS:	
			tmp = new Image (this->lab->xsize,this->lab->ysize);
			for (int y=0; y<this->lab->ysize; ++y)
			for (int x=0; x<this->lab->xsize; ++x)
				tmp->set(x,y,getClassAllComb (this->lab->get(x,y),this->labv->get(x,y)));	
			tmp->write(filename);
			delete tmp;
			break;
			
		case DBG_WIMAGES_RESTORE:
			{ 	Image rest (*this->obs);
				Image tmp (this->lab->xsize, this->lab->ysize, 1);
				for (int y=0; y<this->lab->ysize; ++y)
				for (int x=0; x<this->lab->xsize; ++x)
					tmp.set(x,y,getClass(this->lab->get(x,y),this->labv->get(x,y)));
				VersoReplacementPyramid<TI> p (&rest, &tmp, CL_IND_BG, CL_IND_VERSO,
					6, 4, (unsigned char) this->obsModel->getProperty(PROP_BG_MEAN));	
				p.restore();
				rest.write(filename);			
			}
			break;
		
		default:
			break;
	}	
}

#endif

