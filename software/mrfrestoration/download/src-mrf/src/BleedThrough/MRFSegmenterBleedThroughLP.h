/***********************************************************************
 * Document Bleed-Through Restoration (separation recto/verso)
 * Using two Markov Random Fields.
 * Includes a separate line process for edges.
 *
 * Author: Christian Wolf
 * Begin: 30.1.2007
 ***********************************************************************/
 
#ifndef _WOLF_MRFSEGMENTERBLEEDTHROUGHLP_H_
#define _WOLF_MRFSEGMENTERBLEEDTHROUGHLP_H_

// From the MARKOV RANDOM FIELDS module
#include <MRFLineProcesses.h>
#include <MRFLabelCounterLP2Labels.h>

// From this module
#include "MRFSegmenterBleedThrough.h"

#define LSTSQ_MAX_CNT_COMBINATIONS		100
#define LSTSQ_MIN_CNT_COMBINATIONS		40
#define LSTSQ_MIN_IMPORTANCE			10

// Needed for the calculation of the prior energy
#define PHI_HORIZ_VERT(x,y)			(1-lineim->get((x),(y)))
#define PHI_HORIZ_VERT_MASK(i)		(1-(cLab&(i)))
#define PHI_RD(x,y)					(VirtualLabelKeepsBondRDConn[lineim->get((x),(y))])
#define PHI_LD(x,y)					(VirtualLabelKeepsBondLDConn[lineim->get((x),(y))])
#define PHI_RD_MASK(n,w,e,s)		KeepsBondRDFromNeighbors[(cLab&n)+2*(cLab&w)+4*(cLab&e)+8*(cLab&s)]
#define PHI_LD_MASK(n,w,e,s)		KeepsBondLDFromNeighbors[(cLab&n)+2*(cLab&w)+4*(cLab&e)+8*(cLab&s)]


// The Line process parameters
enum {
	MRFPLP_ZETA=0,
	MRFPLP_ETA,
	MRFPLP_KAPPA,
	
	MRFPLP_AFTER_LAST
};

template <class TI>
class MRFSegmenterBleedThroughLP : public MRFSegmenterBleedThrough<TI>
{
	public:
	
		// Constructor and Destructor
		MRFSegmenterBleedThroughLP(TI *o, Image *lr, Image *lv, ObsModel<TI> *model,
			bool xUseMask, bool xUseGV4RectoRecognition, 
			float xVersoConstraintParamFactor);
		~MRFSegmenterBleedThroughLP();
		
		// Redefinition of virtual methods
		void _init();	
		float priorEnergy(int x, int y);
		bool estimateMRFParamsLS(bool doMedianFilter, bool singleCliquesManually);
		unsigned int getCliqueArrIndexCenterPixel();
		unsigned int gibbsSampler (double T, bool greedy);							
		unsigned int gibbsSamplerLP (double T, bool greedy) ;
		
	protected:	// methods	
						
		float _priorEnergy (Image *labim, Image *lineim, Vector<float> &xMRFParams, int x, int y);
		float _priorEnergyLP (Image *lineim, int x, int y);
		float  priorEnergyLP (int x, int y);
		
		void parameterLessTermOfEnergy2Labels (unsigned int clab, Vector<float> &rv);		
		void parameterLessTermOfEnergyMLabels (unsigned char *lab, Vector<float> &rv);
	
		void updateVirtualLabel(Image *l, int x, int y);
		void updateVirtualLabels4LPLabel(Image *l, int x, int y);
		
		void _initLP(Image &ip, Image &lp);
		bool _estimateMRFParamsLS(Image *xip, Image *xlp, 
			Vector<float> &xMRFParams, bool doMedianFilter, bool singleCliquesManually);
		unsigned int prepareLSEquation2Labels(Image *ip, Image *lp, 
			Matrix<float> &lhs, Vector<float> &rhs);
			
	protected:	// data
			
		Image *line,				// The recto line process
			  *linev,				// The verso line process
			  *lineModif,
			  *linevModif;
			  
		// The Line process parameters. They are identical for 
		// both sides (recto and verso).
		unsigned int noMRFParamsLP;
		Vector<float> MRFParamsLP;	
			  				
	private:	// static data
	
		static unsigned int VirtualLabelKeepsBondRDConn[];
		static unsigned int VirtualLabelKeepsBondLDConn[];
		static unsigned int BlockingCodeFromNeighbors[];
		static unsigned int KeepsBondRDFromNeighbors[];
		static unsigned int KeepsBondLDFromNeighbors[];		
		static unsigned int ParameterIndexFromNeighbors[];
};

/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>
MRFSegmenterBleedThroughLP<TI>::MRFSegmenterBleedThroughLP(TI *o, Image *lr, Image *lv, 
	ObsModel<TI> *m, bool xUseMask, bool xUseGV4RectoRecognition,
	float xVersoConstraintParamFactor)
	: MRFSegmenterBleedThrough<TI>
	(o,lr,lv,m,xUseMask,xUseGV4RectoRecognition,false,xVersoConstraintParamFactor)
{
	// Init the line process	
	line       = new Image (2*lv->xsize, 2*lv->ysize, 1);
	linev      = new Image (2*lv->xsize, 2*lv->ysize, 1);
	lineModif  = new Image (2*lv->xsize, 2*lv->ysize, 1);
	linevModif = new Image (2*lv->xsize, 2*lv->ysize, 1);	
	noMRFParamsLP=MRFPLP_AFTER_LAST;
	MRFParamsLP.resize(noMRFParamsLP);
#warning LP parameters manually set!!
	MRFParamsLP[MRFPLP_ZETA] = 0.9;
	MRFParamsLP[MRFPLP_ETA] = 1.8;
	MRFParamsLP[MRFPLP_KAPPA] = 2.7;
}

/***********************************************************************
 * Destructor
 ***********************************************************************/

template <class TI>
MRFSegmenterBleedThroughLP<TI>::~MRFSegmenterBleedThroughLP()
{
	delete line;
	delete linev;
	delete lineModif;
	delete linevModif;
}

/***********************************************************************
 * A couple of virtual methods we need to redefine
 ***********************************************************************/

template <class TI>
unsigned int MRFSegmenterBleedThroughLP<TI>::getCliqueArrIndexCenterPixel()
{
	ERR_THROW ("MRFSegmenterBleedThroughLP<TI>::getCliqueArrIndexCenterPixel():"
		"does not make sense for 2 labels!");
}

/***********************************************************************
 * Calculate the prior probability for a single combination 
 * of intensity/line process
 * U(f|l)
 ***********************************************************************/

template <class TI>  
float MRFSegmenterBleedThroughLP<TI>::_priorEnergy (Image *labim, Image *lineim,
	Vector<float> &xMRFParams, int x, int y) 
{
	float sum;
	int f=labim->get(x,y);	
		
	// The contribution of the site wise clique
	sum = (f==0 ? 0 : xMRFParams[MRFP_AL8N_SINGLE]);
	
	// The contribution of the pair wise cliques: 
	// Horizontal & Vertical
	sum += xMRFParams[MRFP_AL8N_1_H]* 
		(GAMMA(f,labim->get(x-1,y)) * PHI_HORIZ_VERT(2*x-1,2*y)
		+GAMMA(f,labim->get(x+1,y)) * PHI_HORIZ_VERT(2*x+1,2*y));
	sum += xMRFParams[MRFP_AL8N_1_V]* 
		(GAMMA(f,labim->get(x,y-1)) * PHI_HORIZ_VERT(2*x,2*y-1) 
		+GAMMA(f,labim->get(x,y+1)) * PHI_HORIZ_VERT(2*x,2*y+1));
		
	// Right-diagonal
	sum += xMRFParams[MRFP_AL8N_1_RD] * 
		(GAMMA(f,labim->get(x-1,y+1)) * PHI_RD(2*x-1,2*y+1)
		+GAMMA(f,labim->get(x+1,y-1)) * PHI_RD(2*x+1,2*y-1));
	
	// Left-diagonal
	sum += xMRFParams[MRFP_AL8N_1_LD] * 
		(GAMMA(f,labim->get(x-1,y-1)) * PHI_LD(2*x-1,2*y-1) 
		+GAMMA(f,labim->get(x+1,y+1)) * PHI_LD(2*x+1,2*y+1));
		
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
 * Calculate the prior probability 
 * U(f|l)
 ***********************************************************************/
 
template <class TI>  
inline float MRFSegmenterBleedThroughLP<TI>::priorEnergy (int x, int y) 
{
	float rv;
	
	// The priors coming from the two independent MRFs
	rv = this->_priorEnergy (this->lab,  line,  this->MRFParams,  x, y)+
		 this->_priorEnergy (this->labv, linev, this->MRFParamsv, x, y);	  
		 
	// The prior from the joined information
	// it is necessary to constrain the verso pixels which are 
	// hidden from the recto text pixels.
	if (this->lab->get(x,y) && this->labv->get(x,y))
		rv += this->versoConstraintParam;
	
	return rv;
}

/***********************************************************************
 * Calculate the prior probability for a single combination 
 * of intensity/line process
 * U(l)
 ***********************************************************************/

template <class TI>  
float MRFSegmenterBleedThroughLP<TI>::_priorEnergyLP (Image *lineim, int x, int y) 
{
	unsigned int n,s,w,e;
	float sum=0;
	
	// Two horizontally neighboring cliques
	if (x%2==0)
	{		
		if (x>1) 
		{
			n = lineim->get(x-1,y-1),			
			w = lineim->get(x-2,y),
			e = lineim->get(x,y);	
			s = lineim->get(x-1,y+1),
			sum += MRFParamsLP[ParameterIndexFromNeighbors[n+2*w+4*e+8*s]];
		}
		if (x<lineim->xsize-2) 
		{
			n = lineim->get(x+1,y-1),			
			w = lineim->get(x,y),
			e = lineim->get(x+2,y);	
			s = lineim->get(x+1,y+1),
			sum += MRFParamsLP[ParameterIndexFromNeighbors[n+2*w+4*e+8*s]];
		}
	}
	
	// Two vertically neighboring cliques	
	else
	{
		if (y>1) 
		{
			n = lineim->get(x,y-2),			
			w = lineim->get(x-1,y-1),
			e = lineim->get(x+1,y-1);	
			s = lineim->get(x,y),
			sum += MRFParamsLP[ParameterIndexFromNeighbors[n+2*w+4*e+8*s]];
		}
		if (y<lineim->ysize-2) 
		{
			n = lineim->get(x,y),			
			w = lineim->get(x-1,y+1),
			e = lineim->get(x+1,y+1);	
			s = lineim->get(x,y+2),
			sum += MRFParamsLP[ParameterIndexFromNeighbors[n+2*w+4*e+8*s]];
		}
	}
	return sum;
}

/***********************************************************************
 * Calculate the prior probability 
 * U(l)
 ***********************************************************************/

template <class TI>  
float MRFSegmenterBleedThroughLP<TI>::priorEnergyLP (int x, int y) 
{
	return 
		this->_priorEnergyLP (line,  x, y)+
		this->_priorEnergyLP (linev, x, y);	
}

/***********************************************************************
 * Calculate Ni() i.e. the function returning the different 
 * contributions to each parameter
 * Needed for the least squares estimation of the parameters
 * The lineimage ist always the recto one.
 * For the positions of the clique labels inside the vector,
 * See research notebook 7.2.2007, p.145
 ***********************************************************************/

template <class TI>  
void MRFSegmenterBleedThroughLP<TI>::parameterLessTermOfEnergy2Labels (
	unsigned int cLab, Vector<float> &rv)
{	
	unsigned int f=cLab&CL_MASK_S;
		
	// The contribution of the single site clique
	rv[MRFP_AL8N_SINGLE] = (f==0?0:1);
	
	// The contribution of the pair wise cliques: 
	// the different directions
	rv[MRFP_AL8N_1_H] = 
		(GAMMAMASK(f,cLab&CL_MASK_WE) * PHI_HORIZ_VERT_MASK(14)
		+GAMMAMASK(f,cLab&CL_MASK_EA) * PHI_HORIZ_VERT_MASK(15));
	rv[MRFP_AL8N_1_V] =	
		(GAMMAMASK(f,cLab&CL_MASK_NO) * PHI_HORIZ_VERT_MASK(12)
		+GAMMAMASK(f,cLab&CL_MASK_SO) * PHI_HORIZ_VERT_MASK(17));
	rv[MRFP_AL8N_1_RD]= 
		(GAMMAMASK(f,cLab&CL_MASK_SW) * PHI_RD_MASK(14,16,17,19)
		+GAMMAMASK(f,cLab&CL_MASK_NE) * PHI_RD_MASK(10,12,13,15));
	rv[MRFP_AL8N_1_LD]=	
		(GAMMAMASK(f,cLab&CL_MASK_NW) * PHI_LD_MASK(9,11,12,14) 
		+GAMMAMASK(f,cLab&CL_MASK_SE) * PHI_LD_MASK(15,17,18,20));

#ifdef CHECK_CODE
	for (int i=0; i<MRFP_AL8N_AFTER_LAST; ++i)
	{
		if (isnan(rv[i]))
			ERR_THROW ("MRFSegmenterBleedThroughLP::parameterLessTermOfEnergy(): the value with index "
				<< i << " is NaN. Clab=" << cLab);
	}	
#endif	
}	

template <class TI> 
void MRFSegmenterBleedThroughLP<TI>::parameterLessTermOfEnergyMLabels (
	unsigned char *lab, Vector<float> &rv)
{
	ERR_THROW ("Internal error in MRFSegmenterBleedThroughLP<TI>::parameterLessTermOfEnergyMLabels");
}

/***********************************************************************
 * Some of the entries in the line process matrix are not actually part
 * of the line process but just "virtual" labels which depend on their
 * neighbors. Calculate them
 ***********************************************************************/

template <class TI>  
void MRFSegmenterBleedThroughLP<TI>::updateVirtualLabel(Image *l, int x, int y)
{
	unsigned int n = l->get(x,y-1),
				 s = l->get(x,y+1),
				 w = l->get(x-1,y),
				 e = l->get(x+1,y);
	l->set(x,y,BlockingCodeFromNeighbors[n+2*w+4*e+8*s]);	
}

template <class TI>  
void MRFSegmenterBleedThroughLP<TI>::updateVirtualLabels4LPLabel(Image *l, int x, int y)
{
	if (x>1) updateVirtualLabel(l,x-1,y);	
	if (y>1) updateVirtualLabel(l,x,y-1);
	if (x<l->xsize-2) updateVirtualLabel(l,x+1,y);
	if (y<l->ysize-2) updateVirtualLabel(l,x,y+1);
}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 ***********************************************************************
 ***********************************************************************
 * Estimate the Parameters of the Prior Model
 * Use the least squares method proposed by Derin and Elliott (1985, 1987)
 *
 * returns true if the estimation was successful
 ***********************************************************************
 ***********************************************************************
 *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
 
// A combination of two clique labelings and 
// their importance (=the product of the number 
// of times they occur in the training image).
class LabComb
{
	public:
		unsigned int lab1, lab2, cnt1, cnt2, importance;
		LabComb (unsigned int l1, unsigned int l2, 
			     unsigned int c1, unsigned int c2)
		{
			lab1=l1; lab2=l2;
			cnt1=c1; cnt2=c2;
			importance = cnt1*cnt2;
		}
		bool operator < (const LabComb &o) const
		{
			return importance < o.importance;
		}	
};

/***********************************************************************
 * Prepare the matrix N and the vector p for the equation solved
 * for the least squares estimation,
 * i.e. search for clique pairs useable for the least squares estimation
 * For each clique pair, prepare the RHS and the LHS of the
 * corresponding equation and store them in two data structures
 * The case of 2 labels
 ***********************************************************************/

template <class T>  
unsigned int MRFSegmenterBleedThroughLP<T>::prepareLSEquation2Labels(
	Image *ip, Image *lp, Matrix<float> &lhs, Vector<float> &rhs)
{
	unsigned int cliquesizex=this->getCliqueSizeX();
	unsigned int cliquesizey=this->getCliqueSizeY();	
	FloatVector plTermA(this->noMRFParams), plTermB(this->noMRFParams);
	unsigned int dataCount,dataline;	
	
	// Estimate the clique probabilities
	MRFLabelCounterLP2Labels lc (cliquesizex, cliquesizey, CL_MASK_NS);
	lc.countInImage(*ip,*lp);
	
	cerr << lc;

	// --------------------------------------------------
	// Create a sorted list of labeling combinations
	// --------------------------------------------------
	set<LabComb> combinations;

	// Travers the clique sets. Each set contains several 
	// labelings with identical fNs neighborhoods
	for (MRFLabelCounterLP2Labels::iterator iter=lc.begin(); iter!=lc.end(); ++iter)
	{
		NsMap *nsmap = iter->second;
		map<unsigned int, unsigned int>::iterator iter1,iter2;
				
		// Get all possible labeling combinations 
		// for each clique set
		for (iter1=nsmap->m.begin(); iter1!=nsmap->m.end(); ++iter1)
		for (iter2=nsmap->m.begin(); iter2!=nsmap->m.end(); ++iter2)
		{			
			if (iter1->first != iter2->first)
				combinations.insert(LabComb(iter1->first,iter2->first,
					iter1->second, iter2->second));
		}
	}

	cerr << "Found " << combinations.size() << " combinations.\n";

	// --------------------------------------------------
	// Travers the list of combinations a first time to
	// find out which equations are usefull.
	// --------------------------------------------------
	
	dataCount=0;
	set<LabComb>::iterator iter=combinations.end();	
	do 
	{
		--iter;
		cerr << "combination nr. " << dataCount << " : " << iter->cnt1 << "x" << iter->cnt2 
			 << " = " << iter->importance << endl;
			 
		if (iter->importance<LSTSQ_MIN_IMPORTANCE)
			break;
		if (dataCount>=LSTSQ_MAX_CNT_COMBINATIONS)
			break;
		++dataCount;
	} while (iter!=combinations.begin());
	
	if (dataCount < LSTSQ_MIN_CNT_COMBINATIONS)
		ERR_THROW ("We don't have enough data. Some special treatment should be here ...");
	
	// Allocate the matrix and the vector
	lhs.resize(dataCount,MRFP_AL8N_AFTER_LAST);
	rhs.resize(dataCount);
	
	// --------------------------------------------------
	// Travers the list a second time and build the
	// matrix and the vector
	// --------------------------------------------------
	
	dataline=0;
	iter=combinations.end();	
	do
	{
		--iter;
		if (iter->importance < LSTSQ_MIN_IMPORTANCE)
			break;
			
		unsigned int cLabA   = iter->lab1,
		             cLabB   = iter->lab2, 
		             cntLabA = iter->cnt1, 
		             cntLabB = iter->cnt2;
		             					
		parameterLessTermOfEnergy2Labels(cLabA, plTermA);
		parameterLessTermOfEnergy2Labels(cLabB,  plTermB);
		
		// unusable, since the equation would be 0=constant
		if (plTermA == 0 || plTermB == 0)
			continue;
		
		// Calculate the LHS
		for (unsigned int i=0; i<MRFP_AL8N_AFTER_LAST; ++i)
			lhs(dataline,i)=(plTermA[i]-plTermB[i]);
		
		// Calculate the RHS 
		rhs[dataline] = logf((float) cntLabB/(float) cntLabA);
		
		++dataline;
		if (dataline>=dataCount)
			break;
	} while (iter!=combinations.begin());
	
	// We were too optimistic, some data was not valid.
	// We need to crop the matrix and the vector.
	if (dataline<dataCount)
	{		
		cerr << "Matrix and vector cropped. " << endl
			 << "Old dimensions: lhs=" << lhs.rows() << "x" << lhs.columns() 
			 << ", rhs=" << rhs.size() << endl;
		lhs.crop(0, MRFP_AL8N_AFTER_LAST-1, 0, dataline-1);
		rhs = rhs.cropR(0, dataline-1);		
		cerr << "New dimensions: lhs=" << lhs.rows() << "x" << lhs.columns() 
			 << ", rhs=" << rhs.size() << endl;
			 
		
	}	
	return dataline;
}

template <class T>  
bool MRFSegmenterBleedThroughLP<T>::_estimateMRFParamsLS(
	Image *xip, Image *lp, Vector<float> &xMRFParams, bool doMedianFilter, bool singleCliquesManually)
{		
	unsigned int dataCount;		
	Matrix<float> *Npsinv;
	Vector<float> solution;	
	bool ok=true;
	unsigned int derinBeg;
	Image *ip;
	
	// Eventually filter the image if requested
	if (doMedianFilter)	
	{
		ip = new Image (*xip);	
		filterMedian (*ip, 1, 3);		
	}	
	else
		ip = xip;
	
	HRULE; cerr << "PRIOR: least squares estimation:" << endl;
	
	// stores the LHS and the RHS of the equation to solve
	Matrix<float> N;
	Vector<float> P;			
	
	dataCount = prepareLSEquation2Labels(ip, lp, N, P);
	
	/*
	cerr << "_BEFORE_ LEAST SQUARES ESTIMATION - DERIN & CO (1987):\n"
		 << "=== N: " << N
		 << "=== P: " << P
		 << endl;
	*/
	
	cerr << "Estimate all clique parameters with Derin et al." << endl;
	
	// Get the pseudo inverse of N
	Npsinv = N.pseudoInverse();
		
	// Get the solution
	solution = (*Npsinv)*P;	
	
	// Where part in the parameter vector is used by this solution?
	derinBeg = 0;
	
	// Copy the solution values into the parameter vector
	for (unsigned int p=0; p<solution.size(); ++p)
	{
		xMRFParams[derinBeg+p] = solution[p];
		if (isnan(xMRFParams[p]))
		{
			HRULE;
			cerr << "=== N: " << N
				 << "=== P: " << P
				 << "=== Npsinv: " << *Npsinv;
			HRULE;
			
			/* If we didn't have enough data, then the segmentation 
			 * is almost finished (not enough different labelings in
			 * the image), we therefore just end.
			 */
			if (dataCount < LSTSQ_MIN_DATA_WHEN_INSTABLE)
			{
				if (doMedianFilter)
				{				
					cerr << "The matrix N is unstable, trying to determine the parameters\n"
						 << "without median filtering!" << endl;
					_estimateMRFParamsLS(xip, lp, xMRFParams, false, singleCliquesManually);
				}
				else
					ERR_THROW ("MRFSegmenter::estimateMRFParamsLS(): parameter #" 
						<< p << " is NaN");
			}
		}
	}	
	
	// Debug output
	{		
		cerr << this->noMRFParams << " params: ";	
		for (unsigned int p=0; p<this->noMRFParams; ++p)
		{
			if (p>0) cerr << ", ";
			cerr << xMRFParams[p];
		}
		cerr << endl;
	}
					
	// Clean up
	delete Npsinv;	
	if (doMedianFilter)
		delete ip;		
		
	return ok;
}

/***********************************************************************
 * Estimate the Parameters
 * Estimate them on the recto field only and then adapt them
 * for the verso side
 * returns true if the estimation was successful
 ***********************************************************************/

template <class T>  
bool MRFSegmenterBleedThroughLP<T>::estimateMRFParamsLS(bool doMedianFilter, bool singleCliquesManually)
{		
	if (!_estimateMRFParamsLS (this->lab,  this->line, this->MRFParams, doMedianFilter, singleCliquesManually))
		return false;
	for (unsigned int i=0; i<this->noMRFParams; ++i)
		this->MRFParamsv[i] = this->MRFParams[i];	
	_mirrorParamsX(this->MRFParamsv);
	
	this->versoConstraintParam=0;
	for (unsigned int i=MRFP_FIRST_PAIRWISE; i<MRFP_AL8N_AFTER_LAST; ++i)
			this->versoConstraintParam += this->MRFParamsv[i];
	this->versoConstraintParam*=this->versoConstraintParamFactor;
	
	cerr << "versoConstraint=x*" << this->versoConstraintParamFactor 
	 	 << "=" << this->versoConstraintParam << endl;
	return true;	
}

// **********************************************************************
// One iteration of the gibbs sampler
// Needs to be redefined since we have now 4 label fields
// **********************************************************************

template <typename TI>
unsigned int MRFSegmenterBleedThroughLP<TI>::gibbsSampler (double T, bool greedy) 
{
	return
		// The gibbs sampler of the intensity process is equal to the inherited one
		// since the priorEnergy() method is virtual, thus it calls the new one.
 		MRFSegmenterBleedThrough<TI>::gibbsSampler   (T, greedy)
		                            + gibbsSamplerLP (T, greedy);
}

// **********************************************************************
// The gibbs sampler for the line process is different
// **********************************************************************

template <typename TI>
unsigned int MRFSegmenterBleedThroughLP<TI>::gibbsSamplerLP (double T, bool greedy) 
{
	unsigned int nochangedpixels=0;	
	
	if (this->noClasses>2)
		ERR_THROW ("MRFSegmenterBleedThroughLP<TI>::gibbsSampler(): only works for 2 labels!");
	
	// ------------------------------------------------------------------------
	// We are in "greedy" mode, 
	// This corresponds to the ICM=Iterated Conditional Modes Algorithm 
	// proposed by Besag
	// ------------------------------------------------------------------------
	
	if (greedy)
	{
		// A full sweep of the image
		for (int y=1; y<line->ysize-1; ++y) 
		for (int x=1; x<line->xsize-1; ++x) 
		{
			if (!IS_ELEMENT_OF_MRFLP(x,y))
				continue;
				
			double energy[STATE_COUNT], min;
			byte curvalue1, curvalue2,
				 newvalue1, newvalue2;
			int argMin;			
			
			curvalue1 = line->get(x, y);		
			curvalue2 = linev->get(x, y);					
			newvalue1 = (curvalue1 > 0 ? 0 : 1);
			newvalue2 = (curvalue2 > 0 ? 0 : 1);									

			// Calculate the current energy
			// as well as the energy for all three possible state changes
			energy[STATE_UNCHANGED] = this->priorEnergyLP(x,y);
			line->set (x, y,  newvalue1);
			updateVirtualLabels4LPLabel(line,x,y);
			energy[STATE_CHANGE1]   = this->priorEnergyLP(x,y);
			linev->set (x, y, newvalue2);
			updateVirtualLabels4LPLabel(linev,x,y);
			energy[STATE_CHANGE12]  = this->priorEnergyLP(x,y);
			line->set (x, y,  curvalue1);
			updateVirtualLabels4LPLabel(line,x,y);
			energy[STATE_CHANGE2]   = this->priorEnergyLP(x,y);
						
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
					lineModif->set (x, y, curvalue1);
					linevModif->set(x, y, curvalue2);
					break;
				case STATE_CHANGE1:				
					lineModif->set (x, y, newvalue1);
					linevModif->set(x, y, curvalue2);
					++nochangedpixels;
					break;
				case STATE_CHANGE2:
					lineModif->set (x, y, curvalue1);
					linevModif->set(x, y, newvalue2);
					++nochangedpixels;
					break;
				case STATE_CHANGE12:
					lineModif->set (x, y, newvalue1);
					linevModif->set(x, y, newvalue2);
					++nochangedpixels;
					break;
#ifdef PEDANTIC_CHECK_CODE					
				default:
					throw EError ("MRFSegmenterBleedThroughLP::gibbsSamplerLP(): internal error(1)!\n");
#endif							
			}
				
			// Restore the old values. Necessary since we 
			// perform the modifications on a copy only!!
			// lab has already been done.	
			linev->set (x, y, curvalue2);
			
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
		for (int y=1; y<line->ysize-1; ++y) 
		for (int x=1; x<line->xsize-1; ++x) 
		{
			if (!IS_ELEMENT_OF_MRFLP(x,y))
				continue;
				
			double en_unchanged, en_changed, q;
			byte curvalue1, curvalue2,
				 newvalue1, newvalue2;
			bool takeNewValues;
			byte rb;									
			
			// ------------------------------------------------------------------------
	
			// Calculate the current energy
			en_unchanged = this->priorEnergyLP(x,y);
			
			// Chose randomly one of three new possible configurations
			// and set the label fields to it
			rb = (this->randgen->nextByte())%3;	
			curvalue1 = line->get(x, y);		
			curvalue2 = linev->get(x, y);
			newvalue1 = curvalue1;
			newvalue2 = curvalue2;					
			switch (rb)
			{
				case 0:				
					newvalue1 = (curvalue1 > 0 ? 0 : 1);
					line->set (x, y, newvalue1);
					updateVirtualLabels4LPLabel(line,x,y);						
					break;
				case 1:		
					newvalue2 = (curvalue2 > 0 ? 0 : 1);
					linev->set(x, y, newvalue2);
					updateVirtualLabels4LPLabel(linev,x,y);
					break;
				case 2:			
					newvalue1 = (curvalue1 > 0 ? 0 : 1);
					newvalue2 = (curvalue2 > 0 ? 0 : 1);					
					line->set (x, y, newvalue1);							
					linev->set(x, y, newvalue2);
					updateVirtualLabels4LPLabel(line,x,y);
					updateVirtualLabels4LPLabel(linev,x,y);
					break;
#ifdef PEDANTIC_CHECK_CODE					
				default:
					throw EError ("MRFSegmenterBleedThroughLP::gibbsSamplerLP(): internal error(2)!\n");
#endif					
			}
			
			// calculate the energy of the changed fields
			en_changed = this->priorEnergyLP(x,y);
	
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
				lineModif->set (x,y, newvalue1);
				linevModif->set(x,y, newvalue2);					
			}
			
			// revert the configuration
			else
			{
				lineModif->set (x,y, curvalue1);
				linevModif->set(x,y, curvalue2);					
			}
			
			// Restore the old values. Necessary since we 
			// perform the modifications on a copy only!!
			line->set  (x, y, curvalue1);
			linev->set (x, y, curvalue2);
			updateVirtualLabels4LPLabel(line,x,y);
			updateVirtualLabels4LPLabel(linev,x,y);
			
			DEBUGIFPX(x,y)
				cerr << endl;
		}	
	}
	
	// Apply the changes
	for (int y=1; y<line->ysize-1; ++y) 
	for (int x=1; x<line->xsize-1; ++x)
	{ 
		if (IS_ELEMENT_OF_MRFLP(x,y))
		{				
			line->set(x,y,  lineModif->get(x,y));
			linev->set(x,y, linevModif->get(x,y));
		}
		if (IS_VIRTUAL_LABEL_MRFLP(x,y))
		{
			updateVirtualLabel(line,x,y);
			updateVirtualLabel(linev,x,y);
		}
	}
		
	return nochangedpixels;	
}
#undef PRINT_ENERGY

/***********************************************************************
 * Initialise a line process
 ***********************************************************************/

template <class TI>
void MRFSegmenterBleedThroughLP<TI>::_initLP(Image &ip, Image &lp)
{
	// Sobel filter the image:
	// a special version with a 3x2 mask. 
	// The results are 
	// directly stored in the line process (since it is gray value)
	// See research notebook 16.2.2007, p. 147
	
	// Border treatment. Could be done more properly, but we would 
	// need to implement a lot of special cases. 
	// So we just set the border edges to zero.
	for (int y=0; y<lp.ysize; ++y)	
	{
		lp.set(0,y,0);
		lp.set(lp.xsize-1,y,0);
	}			
	for (int x=0; x<lp.xsize; ++x)
	{
		lp.set(x,0,0);
		lp.set(x,lp.ysize-1,0);
	}	
		
	// The standard non border case
	for (int y=1; y<lp.ysize-1; ++y)
	for (int x=1; x<lp.xsize-1; ++x)
	{
		if (IS_HORIZ_EDGE_OF_MRFLP(x,y))
		{
			lp.set(x,y, (ip.get(x/2,y/2) != ip.get(x/2,y/2+1))?1:0);
		}				
		if (IS_VERT_EDGE_OF_MRFLP(x,y))
		{
			lp.set(x,y, (ip.get(x/2,y/2) != ip.get(x/2+1,y/2))?1:0);
		}
	}		
}

/***********************************************************************
 * Initialise everything
 ***********************************************************************/

template <class TI>
void MRFSegmenterBleedThroughLP<TI>::_init()
{
	// Init the intensity process
	MRFSegmenterBleedThrough<TI>::_init();
	
	// Init the line process
	cerr << "Initializing the two line processes ...";
	_initLP(*this->lab, *line);
	_initLP(*this->labv, *linev);
	cerr << "Done." << endl;
		
	lineModif->setZero();
	linevModif->setZero();
	
	// Debug output
	{	Image tmp (*line), tmp2 (*linev);
		brightenImage (tmp);
		brightenImage (tmp2);
		reverseVideo (tmp);
		reverseVideo (tmp2);
		tmp.write  ("x_007_init_lp_recto.pgm");
		tmp2.write ("x_007_init_lp_verso.pgm");
	}
}

#endif

