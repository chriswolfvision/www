/***********************************************************************
 * Segmentation using Markov Random Field Models
 *
 * Author: Christian Wolf
 * Begin: 24.5.2005
 *
 * WHAT SUBCLASSES OF MRFSEGMENTER NEED TO DO:
 * - All virtual methods need to implemented (you might have guessed...)
 * - noMRFParams must be set in the constructor
 ***********************************************************************/
 
#ifndef _WOLF_MRFSEGMENTER_H_
#define _WOLF_MRFSEGMENTER_H_

/***********************************************************************
 * PARAMETERS
 ***********************************************************************/
 
#define LSTSQ_DATA_PERCENTAGE			50
#define LSTSQ_MIN_DATA_WHEN_INSTABLE	30

/***********************************************************************
 * MISC
 ***********************************************************************/

#define OPTALG_SIMA		"sima"
#define OPTALG_GRCT		"grct"


// Options for the output of debug images
enum DebugWriteImageOptions 
{
	DBG_WIMAGES_NO=0,
	DBG_WIMAGES_SEG,
	DBG_WIMAGES_LABELS,
	DBG_WIMAGES_RESTORE
};
 
/***********************************************************************/
 
// C++
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>

// From the main module
#include <CIL.h>

// From the MATH module
#include <Vector.h>
#include <Matrix.h>
#include <RandomNumberGenerator.h>
#include <Histogram.h>

// From the IMAGE PROCESSING module
#include <ImageProc.h>

// From the MRF OBSERVATION MODELS module
#include <Segmenter.h>

// From this module
#include <MRFLabelCounter.h>
#include <MRFLabelCounter2Labels.h>
#include <MRFLabelCounterMLabels.h>
#include <ObsModel.h>

typedef Vector<float> FloatVector;

enum ObsModelType
{
	MRF_OMT_GAUSS1V=0,
	MRF_OMT_GAUSS1VCORR,
	MRF_OMT_GAUSSMV
};

// Used for the least squares estimation of the MRF parameters
struct LabelPair
{
	unsigned int dataIndex, count;
	float PDiv;	
	const bool operator < (const LabelPair &o) const { 	return count < o.count;	}
};

/***********************************************************************
 * Template parameters:
 * TI .... the type of observation image 
 *         (Image, FloatMatrix or MultiBandFloatMatrix)
 *
 * Constructor:
 * o ..... the observation image (input)
 * l ..... the label image (output)
 * obsModel .... The observation model (input)
 * mask ........ The background mask. The pixels which shall not be 
 *               treated need to be 0. (may be NULL -> ignored)
 * maskLabel ... the label of the pixels which are masked out by "mask"
 ***********************************************************************/

template <class TI>
class MRFSegmenter : public Segmenter<TI>
{
	public:
	
		// Constructor & Destructor
		MRFSegmenter (TI *o, Image *l, ObsModel<TI> *obsModel, 
			char *xForceMRFParams, unsigned int noClasses);
		virtual ~MRFSegmenter();
		
		// Initialisation techniques
		void initKMeans(vector<typename TI::PixelType> *initValues);
		
		// Parameter estimation
		virtual bool estimateMRFParamsLS(bool doMedianFilter, bool singleCliquesManually);
		virtual void estimateObsParams();
		
		// Optimization & sampling
		virtual void iteratedConditionalEstimation (int noItsICE, 
			double T1, double C, int noItsSA, int noItsSAGreedy, 
			bool doMedianFilter, bool singleCliquesManually, 
			char *optAlg,
			DebugWriteImageOptions debugOutputOptions, char *debugString);
		virtual void simulatedAnnealing (double T1, double C, int noIts, int noItsGreedy, 
 			DebugWriteImageOptions debugOutputOptions, char *debugString=NULL); 
		virtual unsigned int gibbsSampler (double T, bool isGreedy);
				
		// MUST BE CALLED!!!
		void init();
		
		// Accessors
		Image *getMask(){ return mask; }
		
		// Access to the results
		virtual void writeLabelField (char *filename, DebugWriteImageOptions options);
		
		// Sub classes supporting graph cut inference need to implement this method
		virtual void graphCutAdaptedToProblem(unsigned int noIterations)
		{
			ERR_THROW ("MRFSegmenter: Graph cut optimization not implemented!");
		}		
		
	protected:	// Methods
	
		float condEnergy  (int x, int y);					// likelihood
				
		// ----------------------------------------------------
		// These functions need to be defined by the subclasses
		// ----------------------------------------------------
		
		// The prior energy
		
		virtual float priorEnergy (int x, int y)=0;			// prior		
		
		// For the least squares estimation (the case of 2 labels)
		virtual unsigned int getCliqueSizeX()=0;
		virtual unsigned int getCliqueSizeY()=0;
		virtual unsigned int getCliqueMaskCenterPixel()=0;		
		virtual void  parameterLessTermOfEnergy2Labels (unsigned int lab, Vector<float> &rv)=0;
		
		// For the least squares estimation (the case of multiple labels)
		virtual unsigned int getCliqueArrIndexCenterPixel()=0;
		virtual unsigned int firstSingleCliqueParameter()=0;
		virtual unsigned int firstPairwiseCliqueParameter()=0;
		virtual void  parameterLessTermOfEnergyMLabels (unsigned char *lab, Vector<float> &rv)=0;
		
		virtual void _init()=0;		
		
		// ----------------------------------------------------
		
	protected:	// Help methods		
		
		
		void _estimateSingleCliqueParameters(Image *lab, Vector<float> &p);
		bool _estimateMRFParamsLS(Image *l, Vector<float> &xMRFParams,
			bool doMedianFilter, bool singleCliquesManually);		
		unsigned int  prepareLSEquation2Labels(Image *l, vector<FloatVector *> &lhs,
			multiset<LabelPair> &rhs);
		unsigned int prepareLSEquationMLabels(Image *l, vector<FloatVector *> &lhs,
			multiset<LabelPair> &rhs);
		unsigned int gibbsSampler2Labels (double T, bool greedy);
		unsigned int gibbsSamplerMLabels (double T, bool greedy);
			
	protected: 	// Data
	
		Image *labModif;	// The label field currently modified
				
		unsigned int noMRFParams;
		Vector<float> MRFParams;
		char *forceMRFParams;		
		
		// The observation model
		ObsModel<TI> *obsModel;
		
		// The mask which determines the pixel not 
		// to process
		Image *mask;
		unsigned int maskLabel;
		
		// misc stuff
		RandomNumberGenerator *randgen;
		
		// Counter of the current iteration number 
		// of iterated conditional estimation.
		// Will be used by the children.
		unsigned int curICEIteration;
};

/***********************************************************************
 * Constructor
 ***********************************************************************/

template <class TI>  
MRFSegmenter<TI>::MRFSegmenter (TI *o, Image *l, ObsModel<TI> *xObsModel,
	char *xForceMRFParams, unsigned int xNoClasses)
{
	cerr << "New model: MRFSegmenter" << endl;
	
	// the arrays containing clique labelings are terminated with the byte 255
	if (xNoClasses>=CLIQUE_LABEL_TERMINATION)
		ERR_THROW ("MRFSegmenter<TI>::MRFSegmenter: number of classes "
		"must not be greater than 254!");
	this->obs = o;
	this->lab = l;
	obsModel = xObsModel;	
	mask=NULL;
	maskLabel=0;
	this->noClasses = xNoClasses;
	randgen = new RandomNumberGenerator (true, 0, 1, 2048);	
	labModif = new Image (l->xsize, l->ysize, 1);
	forceMRFParams = xForceMRFParams;
}

/***********************************************************************
 * Destructor
 ***********************************************************************/

template <class TI>  
MRFSegmenter<TI>::~MRFSegmenter()
{
	delete randgen;
}

/***********************************************************************
 * Post constructor initialisation
 * here is some code which needs to be exectuted after ALL
 * constructors (including the subclass constructors) have been called
 ***********************************************************************/

template <class T>  
void MRFSegmenter<T>::init () 
{
	cerr << "Init of MRFSegmenter.\nInit of observation model:\n";
	obsModel->init(this);	
	cerr << "Init of prior model:\n";
	this->_init();	
}

/***********************************************************************
 * Calculate the likelihood of the data given the parameters
 ***********************************************************************/

template <class T>  
inline float MRFSegmenter<T>::condEnergy (int x, int y) 
{
	return obsModel->condEnergy (x,y);
}

/***********************************************************************
 * Initialise the labeling using the KMeans clustering algorithm
 ***********************************************************************/
 
template <class T>  
void MRFSegmenter<T>::initKMeans(vector<typename T::PixelType> *initValues)
{
	typedef typename T::PixelType PixelType;
	
	HRULE; cerr << "K-Means Initialisation (MRFSegmenter): " << endl;
	
	segmentKMeans (*this->obs, *this->lab, this->noClasses, initValues, NULL);
	
	// Debug: print the number of pixels per class
	Vector<unsigned int> count (this->noClasses);
	count.setZero();	
	for (int y=0; y<this->lab->ysize; ++y) 
	for (int x=0; x<this->lab->xsize; ++x) 	
		++count(this->lab->get(x,y));
	for (unsigned int l=0; l<this->noClasses; ++l) 
		cerr << "class " << l << ": " << count[l] << " pixels." << endl;
		
	// Debug output
	{	Image tmp (*(this->lab));
		brightenImage (tmp);
		reverseVideo (tmp);		
		tmp.write ("x_005_kmeans_reordered_rev.pgm");
	}
}

/***********************************************************************
 * Estimate the Parameters of the Observation Model
 ***********************************************************************/
 
template <class T>  
void MRFSegmenter<T>::estimateObsParams()
{
	obsModel->estimateParams();
}

// **********************************************************************
// Iterated Conditional Estimation
// **********************************************************************

template <typename TI>
void MRFSegmenter<TI>::iteratedConditionalEstimation (int noItsICE, 
	double T1, double C, int noItsSA, int noItsSAGreedy, 
	bool doMedianFilter, bool singleCliquesManually,
	char *optAlg, 
	DebugWriteImageOptions debugOutputOptions, char *debugString) 
{    
	char buf[100];
	char *debugStringTreated=NULL;	
	bool foundOptAlg;
	
	for (int j=0; j<noItsICE; ++j)
	{
		HRULE; cerr << "ICE ITERATION NR. " << j+1 << endl; cerr.flush();
		this->curICEIteration=j;
	
		if (debugString!=NULL)
		{						
			sprintf (buf, "%s_%2.2d",debugString,j);
			debugStringTreated=buf;
		}	
		else
			debugStringTreated=NULL;
				
		// --------------------------------------------- 
		// Parameter estimation
		// --------------------------------------------- 
				
		if (!estimateMRFParamsLS(doMedianFilter, singleCliquesManually))
			break;
			
		estimateObsParams();				
		if (debugString!=NULL)
			cerr << "debugString(1)=[" << debugString << "]" << endl;
							
		// --------------------------------------------- 
		// Maximization of the posterior probability
		// --------------------------------------------- 
							
		foundOptAlg=false;
		if (strcmp(optAlg,OPTALG_SIMA)==0)
		{					
			simulatedAnnealing (T1, C, noItsSA, noItsSAGreedy, debugOutputOptions, debugStringTreated);
			foundOptAlg=true;
		}
		if (strcmp(optAlg,OPTALG_GRCT)==0)
		{		
			graphCutAdaptedToProblem(noItsSA);				
			foundOptAlg=true;
		}
		if (!foundOptAlg)
		{
			ERR_THROW ("Unknown optimization algorithm: " << optAlg);
		}
	}	
}

// **********************************************************************
// Simulated annealing
// The image should already be initialized.
// **********************************************************************

template <typename TI>
void MRFSegmenter<TI>::simulatedAnnealing (double T1, double C, int noIts, int noItsGreedy,
	DebugWriteImageOptions debugOutputOptions, char *debugString) 
{    
    unsigned int nochangedpixels;
    double T;
    int k;
	
	HRULE;
	cerr << "T1=" << T1
		 << " C=" << C
		 << " noIt=" << noIts
		 << " noItsGreedy=" << noItsGreedy 
		 << " outputoptions=" << debugOutputOptions << endl;
	
	if ((debugString!=NULL) && (debugOutputOptions!=DBG_WIMAGES_NO))
	{			
		char buf[100];
		sprintf (buf, "%s_0000.pgm",debugString);		
		writeLabelField(buf, debugOutputOptions);			
	}

    // Iterate: One iteration contains a full sweep of the image
   	nochangedpixels = 1;
	for (k=1; k<=(noIts+noItsGreedy) ; ++k) 
	{		
    	// If we are in greedy mode and no pixels have been changed
    	// anymore, then exit
    	if (k>noIts && nochangedpixels==0)
    		break;    	

		// Set the temperature
		// Equation proposed by Daniel Dementhon, taken from Duda-Hart-Stork
		T = T1 * pow(C,k-1);		
		nochangedpixels = gibbsSampler(T, k>=noIts);
		
		// DEBUG Output of the number of iterations
		if ((k==1)||(k%1==0)||k>noIts) 
		{
			cerr << "[" << k; //<< ",T=" << T;
				
			// We are in greedy mode
			if (k>noIts)
				cerr << ":" << nochangedpixels;

			cerr << "]"; cerr.flush();
    	}
		
		if (debugString!=NULL)
		{			
			char buf[100];
			sprintf (buf, "%s_%4.4d.pgm",debugString,k);
			writeLabelField(buf, debugOutputOptions);				
		}
	}

	cerr << endl << k-noIts << " iterations in greedy mode.\n";
}

// **********************************************************************
// One iteration of the gibbs sampler
// There are two versions, one for 2 labels, one for multiple labels
// **********************************************************************

template <typename TI>
unsigned int MRFSegmenter<TI>::gibbsSampler (double T, bool greedy) 
{
	if (this->noClasses<=2)
		return gibbsSampler2Labels (T, greedy);
	else
		return gibbsSamplerMLabels (T, greedy);
}

// **********************************************************************
// One iteration of the gibbs sampler
// The version for 2 labels
// **********************************************************************

template <typename TI>
unsigned int MRFSegmenter<TI>::gibbsSampler2Labels (double T, bool greedy) 
{
	unsigned int nochangedpixels=0;
	
	// A full sweep of the image
	for (int y=1; y<this->lab->ysize-1; ++y) 
	for (int x=1; x<this->lab->xsize-1; ++x) 
	{
		double en_unchanged, en_changed, q;
		byte curvalue, newvalue;
		bool takeNewValue;
		
		if (mask!=NULL && mask->get(x,y)==0)
		{
			labModif->set(x,y,maskLabel);
			continue;
		}
	
		// Calculate the current energy
		en_unchanged = priorEnergy(x,y) + condEnergy(x,y);

		// Change the pixel
		curvalue = this->lab->get(x, y);
		newvalue = (curvalue > 0 ? 0 : 1);
		takeNewValue = true;
		this->lab->set (x, y, newvalue);

		// calc the energy of the changed image
		en_changed = priorEnergy(x,y) + condEnergy(x,y);

		// The energy change
		q = exp(-1.0*(en_changed-en_unchanged)/T);
		
		// Energy changes to worse, so we should not change it,
		// but with a certain probability we do it anyway
		if (q<1) 
		{
			double randomnumber = randgen->next();

			if (randomnumber > q)
				takeNewValue = false;
		}

		// In the last iteration, we force the "good"
		// decisions only, ---> switch into GREEDY mode.
		if (greedy) 
		{
			if (q<1)
				takeNewValue = false;
			else
				++nochangedpixels;
		}		
		
		labModif->set(x,y, (takeNewValue ? newvalue : curvalue));
		
		// Restore the old value. Necessary since we 
		// perform the modifications on a copy only!!
		this->lab->set(x, y, curvalue);		
	}	
	
	// Apply the changes
	for (int y=1; y<this->lab->ysize-1; ++y) 
	for (int x=1; x<this->lab->xsize-1; ++x) 
		this->lab->set(x,y, labModif->get(x,y));
	
	return nochangedpixels;	
}

// **********************************************************************
// One iteration of the gibbs sampler
// The version for multiple labels
// **********************************************************************

#define	PRINT_ENERGY(l)			{ this->lab->set (x, y, (l)); \
								  float p=this->priorEnergy(x,y), c=this->condEnergy(x,y); \
								  cerr << "  |  " << (int) l << ": prior=" \
								  	   <<  p << " cond=" <<  c << " sum=" << p+c << endl; }

template <typename TI>
unsigned int MRFSegmenter<TI>::gibbsSamplerMLabels (double T, bool greedy) 
{
	unsigned int nochangedpixels=0;	
	bool takeNewValue;
	
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
			double energy, en_min;			
			unsigned char st_unchanged, st_changed;
			
			// ------------------------------------------------------------------------
			// DETAILED DEBUG OUTPUT: all energies for the 4 different possible states
			// ------------------------------------------------------------------------
			
			DEBUGIFPX(x,y)
			{
				unsigned int curvalue1 = this->lab->get(x, y);						
				
				HRULE;
				cerr << DEBUGSTRPX 
					 << "observed gray value: " << (int) this->obs->get(x,y) << endl
					 << DEBUGSTRPX 
					 << "unch. state=" << (int) this->lab->get(x, y) 
					 << endl;
					 
				for (unsigned i=0; i<this->noClasses; ++i)					 
					PRINT_ENERGY(i);
				HRULE;
																     
				this->lab->set(x, y, curvalue1);					
			}
			
			// Calculate the current energy
			st_changed = st_unchanged = this->lab->get(x, y);		
			en_min = this->priorEnergy(x,y) + this->condEnergy(x,y);				
						
			// Calculate the energy for all possible states (changes and current)
			// choose the minimum
			for (unsigned int st=0; st<this->noClasses; ++st)
			{				
				if (st==st_changed) // do not recalculate the current state
					continue;
				this->lab->set (x, y, st);
				energy = this->priorEnergy(x,y) + this->condEnergy(x,y);				
				if (energy<en_min)
				{
					en_min=energy;
					st_changed=st;
				}
			}
				
			labModif->set(x,y, st_changed);
			if (st_unchanged!=st_changed)
				++nochangedpixels;			
			
			// Restore the old value. Necessary since we 
			// perform the modifications on a copy only!!
			this->lab->set(x, y, st_unchanged);	
				
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
			double en_changed, en_unchanged, q;			
			unsigned char st_unchanged,st_changed;
			
			// ------------------------------------------------------------------------
			// DETAILED DEBUG OUTPUT: all energies for the 4 different possible states
			// ------------------------------------------------------------------------
			
			DEBUGIFPX(x,y)
			{
				unsigned int curvalue1 = this->lab->get(x, y);						
				
				HRULE;
				cerr << DEBUGSTRPX 
					 << "observed gray value: " << (int) this->obs->get(x,y) << endl
					 << DEBUGSTRPX 
					 << "unch. state=" << (int) this->lab->get(x, y) 
					 << endl;
					 
				for (unsigned i=0; i<this->noClasses; ++i)					 
					PRINT_ENERGY(i);				
				HRULE;
																     
				this->lab->set(x, y, curvalue1);					
			}
			
			// Calculate the current energy
			st_unchanged = this->lab->get(x, y);		
			en_unchanged = this->priorEnergy(x,y) + this->condEnergy(x,y);							
			
			// Chose randomly one of three new possible configurations
			// and set the label fields to it. 
			// There are noClasses-1 possible state changes since the
			// current one does not count.
			st_changed = (this->randgen->nextByte())%(this->noClasses-1);
			if (st_changed>=st_unchanged)
				++st_changed;
				
			this->lab->set (x, y, st_changed);						
								
			// calculate the energy of the changed field
			en_changed = this->priorEnergy(x,y) + this->condEnergy(x,y);
	
			// The energy change
			q = exp(-1.0*(en_changed-en_unchanged)/T);
			
			// Energy changes to worse, so we should not change it,
			// but with a certain probability we do it anyway
			if (q<1) 
				takeNewValue = (this->randgen->next() <= q);
			else
			{
				takeNewValue=true;
				++nochangedpixels;
			}
			
			labModif->set(x,y, (takeNewValue ? st_changed : st_unchanged));	
			
			// Restore the old values. Necessary since we 
			// perform the modifications on a copy only!!
			this->lab->set  (x, y, st_unchanged);
			
			DEBUGIFPX(x,y)
				cerr << endl;
		}	
	}
	
	// Apply the changes
	for (int y=1; y<this->lab->ysize-1; ++y) 
	for (int x=1; x<this->lab->xsize-1; ++x)
		this->lab->set(x,y,  this->labModif->get(x,y));
		
	return nochangedpixels;	
}
#undef PRINT_ENERGY

/***********************************************************************
 * Write the label field into a file
 ***********************************************************************/
 
template <class TI>
void MRFSegmenter<TI>::writeLabelField(char *s, DebugWriteImageOptions options)
{	
	Image *ptmp;
	switch (options)
	{
		case DBG_WIMAGES_NO:			
			break;
		
		case DBG_WIMAGES_SEG:
			ptmp = this->obsModel->visualizeSegmentation();
			ptmp->write(s);
			delete ptmp;
			break;
	
		case DBG_WIMAGES_LABELS:
			this->lab->write(s);
			break;
		
		default:
			break;
	}		
}

/***********************************************************************
 * Estimate the single clique parameters with a simple histogramming
 * technique
 ***********************************************************************/

template <class T>  
void MRFSegmenter<T>::_estimateSingleCliqueParameters(Image *lab, Vector<float> &p)
{
	Histogram<unsigned int> *H;
	unsigned int sum;
	
	H = buildHistogram<unsigned int> (*lab);
	sum = H->sum();
		
	// Label 0 does not have a parameter, so we start with label 1
	// The parameter for label 1 is stored at array index firstSingleCliqueParameter()
	for (unsigned int l=1; l<this->getNoClasses(); ++l)	
		p[firstSingleCliqueParameter()+l-1] = -1.0*log((float)(*H)[l]/(float)sum);			
		
	/* DEBUG OUTPUT */
	{
		cerr << "Manually estimated single clique parameters:";
		for (unsigned int l=0; l<this->getNoClasses()-1; ++l)	
			cerr << "p[" << firstSingleCliqueParameter()+l << "]=" 
			     << p[firstSingleCliqueParameter()+l] << " ";
		cerr << endl;
	}
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
  
template <class T>  
bool MRFSegmenter<T>::estimateMRFParamsLS(bool doMedianFilter, bool singleCliquesManually)
{
	return _estimateMRFParamsLS (this->lab, this->MRFParams, doMedianFilter, 
		singleCliquesManually);
}

/***********************************************************************
 * Prepare the matrix N and the vector p for the equation solved
 * for the least squares estimation,
 * i.e. search for clique pairs useable for the least squares estimation
 * For each clique pair, prepare the RHS and the LHS of the
 * corresponding equation and store them in two data structures
 * The case of 2 labels
 ***********************************************************************/

template <class T>  
unsigned int MRFSegmenter<T>::prepareLSEquation2Labels(Image *l,
	vector<FloatVector *> &lhs,
	multiset<LabelPair> &rhs)
{
	unsigned int cliquesizex=getCliqueSizeX();
	unsigned int cliquesizey=getCliqueSizeY();
	unsigned int cliquecount=(int) pow(2,cliquesizex*cliquesizey);
	unsigned int centerMask=getCliqueMaskCenterPixel();
	FloatVector plTermA(noMRFParams), plTermB(noMRFParams);
	FloatVector *diff; 
	unsigned int dataCount;	
	LabelPair labPair;
	
	// Estimate the clique probabilities
	MRFLabelCounter2Labels lc (cliquesizex, cliquesizey);
	lc.countInImage(*l);

	dataCount=0;
	for (unsigned int cLabA=0; cLabA<cliquecount; ++cLabA)
	{
		unsigned int cLabB, cntLabA, cntLabB;
		
		// We only consider pairs, so a clique with center pixel 1 
		// is automatically considered with its sister clique
		if ((cLabA&centerMask)>0)
			continue;
						
		cLabB = cLabA|centerMask;
		cntLabA=lc.get(cLabA);
		cntLabB=lc.get(cLabB);
			
 		// unusable, there is no log(0)
		if (cntLabA==0 || cntLabB==0)
			continue;
			
		parameterLessTermOfEnergy2Labels(cLabA, plTermA);
		parameterLessTermOfEnergy2Labels(cLabB,  plTermB);
		
		// unusable, since the equation would be 0=constant
		if (plTermA == 0 || plTermB == 0)
			continue;
		
		// Calculate the LHS
		diff = new FloatVector (plTermA);
		*diff -= plTermB;
		lhs.push_back(diff);
		
		// Calculate the RHS and store it in a sorted container
		labPair.PDiv = logf((float) cntLabB/(float) cntLabA);
		labPair.dataIndex = dataCount; 
		labPair.count=cntLabA+cntLabB;
		rhs.insert(labPair);
		
		++dataCount;
	}	
	return dataCount;
}

/***********************************************************************
 * Prepare the matrix N and the vector p for the equation solved
 * for the least squares estimation
 * i.e. search for clique pairs useable for the least squares estimation
 * For each clique pair, prepare the RHS and the LHS of the
 * corresponding equation and store them in two data structures
 * The case of multiple labels 
 ***********************************************************************/

template <class T>  
unsigned int MRFSegmenter<T>::prepareLSEquationMLabels(Image *l,
	vector<FloatVector *> &lhs,
	multiset<LabelPair> &rhs)
{
	unsigned int cliquesizex=getCliqueSizeX();
	unsigned int cliquesizey=getCliqueSizeY();
	unsigned int cliquesize=cliquesizex*cliquesizey;
	unsigned int centerIndex=getCliqueArrIndexCenterPixel();
	FloatVector plTermA(noMRFParams), plTermB(noMRFParams);
	FloatVector *diff; 
	unsigned char *cLabA, *cLabB;
	unsigned int dataCount, pos;
	LabelPair labPair;
	
	cLabA = new unsigned char [cliquesize+1];
	cLabB = new unsigned char [cliquesize+1];
	
	// Estimate the clique probabilities
	MRFLabelCounterMLabels lc (cliquesizex, cliquesizey);
	lc.countInImage(*l);
	
	for (unsigned int i=0; i<cliquesize; ++i)
		cLabA[i]=0;
				
	// Travers all the possible clique labelings	
	dataCount=0;
	while (true)
	{					
		// Change to the next clique labeling
		pos=0;
		while (pos<cliquesize)
		{			
			if (cLabA[pos]<this->noClasses-1)
			{
				++cLabA[pos];				
				break;
			}
			else
			{	
				cLabA[pos]=0;
				++pos;				
			}			
		} 
		
		cLabA[cliquesize]=CLIQUE_LABEL_TERMINATION;		
		// cerr << "["; prClique(cLabA); cerr << "]" << endl;
			
		// We did all the different clique labelings
		if (pos>=cliquesize && cLabA[pos-1]==0)
			break;
		
		// We consider all labeling pairs
		for (unsigned char BCenter=cLabA[centerIndex]+1; BCenter<this->noClasses; ++BCenter)
		{	
			memcpy(cLabB,cLabA,cliquesize+1);
			cLabB[centerIndex]=BCenter;								
			unsigned int cntLabA=lc.get(cLabA);
			unsigned int cntLabB=lc.get(cLabB);
				
			// unusable, there is no log(0)
			if (cntLabA==0 || cntLabB==0)
				continue;
				
			parameterLessTermOfEnergyMLabels(cLabA, plTermA);
			parameterLessTermOfEnergyMLabels(cLabB, plTermB);
			
			// unusable, since the equation would be 0=constant
			if (plTermA == 0 || plTermB == 0)
				continue;
			
			// Calculate the LHS
			diff = new FloatVector (plTermA);
			*diff -= plTermB;
			lhs.push_back(diff);
			
			/*
			cerr << "[";
			prClique(cLabA);
			cerr << "/";
			prClique(cLabB);
			cerr << ";" << cntLabB << "," << cntLabA << "]" << endl; 
			*/
			
			// Calculate the RHS and store it in a sorted container
			labPair.PDiv = logf((float) cntLabB/(float) cntLabA);			
			
			labPair.dataIndex = dataCount; 
			labPair.count=cntLabA+cntLabB;		
			rhs.insert(labPair);
			
			++dataCount;
		}
	}	
	
	// Clean up
	delete cLabA;
	delete cLabB;
	
	return dataCount;
}


template <class T>  
bool MRFSegmenter<T>::_estimateMRFParamsLS(Image *xl, Vector<float> &xMRFParams,
	bool doMedianFilter, bool singleCliquesManually)
{		
	unsigned int dataCount, dataCountReduced;		
	Matrix<float> *Npsinv;
	Vector<float> solution;	
	bool ok=true;
	unsigned int derinBeg;
	Image *l;
	
	if (doMedianFilter)	
	{
		cerr << "Median filtering before least-squares estimation.\n";
		l = new Image (*xl);	
		filterMedian (*l, 1, 3);		
	}	
	else
		l = xl;
		
	// We don't estimate, we got the parameters by command line argument
	if (forceMRFParams!=NULL)
	{
		double *tmp = new double [noMRFParams];
		parseParameterString (forceMRFParams, tmp, noMRFParams, (char *) "-h");
		for (unsigned int p=0; p<noMRFParams; ++p)
			xMRFParams[p] = tmp[p];
		delete [] tmp;
			
		// Debug 
		cerr << noMRFParams << " params (forced from command line): ";	
		for (unsigned int p=0; p<noMRFParams; ++p)
		{
			xMRFParams[p] = tmp[p];
			if (p>0) cerr << ", ";
			cerr << xMRFParams[p];
		}
		cerr << endl;
		return true;
	}
	
	HRULE; cerr << "PRIOR: least squares estimation:" << endl;
	
	// stores the LHS of the equation to solve
	vector<FloatVector *> lhs;								
	
	// stores the RHS of the equation to solve
	// the index is the n.o. labels
	multiset<LabelPair> rhs;			
	
	if (this->noClasses<=2)
		dataCount = prepareLSEquation2Labels(l, lhs, rhs);
	else
		dataCount = prepareLSEquationMLabels(l, lhs, rhs);					
			
	// -------------------------------------------------------------------
	// Use only the label pairs with the largest number of found labelings
	// Take the RHS and build a vector
	// Take the LHS and build a matrix
	// -------------------------------------------------------------------
	
	dataCountReduced = (dataCount*LSTSQ_DATA_PERCENTAGE)/100;
	if (dataCountReduced<noMRFParams)
		dataCountReduced=dataCount;
	if (dataCountReduced>rhs.size())
		dataCountReduced=rhs.size();
		
	if (dataCount<noMRFParams)
		ERR_THROW ("MRFSegmenter::estimateMRFParamsLS(): not enough training data!!");
		
	cerr << "DataCount (reduced/total): " << dataCountReduced << "/" << dataCount << endl;
	
	multiset<LabelPair>::iterator iter=rhs.end();
	--iter;
	Matrix<float> N(dataCountReduced, noMRFParams);
	Vector<float> P(dataCountReduced);	
	
	for (unsigned int d=0; d<dataCountReduced; ++d,--iter)
	{
		// cerr << iter->count << ": " << iter->dataIndex << "=" << iter->PDiv << endl;
				
		P[d] = iter->PDiv;
		for (unsigned int p=0; p<noMRFParams; ++p)
		{
			unsigned int ind=iter->dataIndex;			
			N(d,p) = (*lhs[ind])[p];
		}
		
		if (iter==rhs.begin())
			ERR_THROW("MRFSegmenter<T>::_estimateMRFParamsLS(): internal error (1)\n"
				<< "d=" << d << " rhs.size()=" << rhs.size());
	}
	
	/*
	cerr << "_BEFORE_ LEAST SQUARES ESTIMATION - DERIN & CO (1987):\n"
		 << "=== N: " << N
		 << "=== P: " << P
		 << endl;
	*/
	
	// -------------------------------------------------------------------
	// Estimate the single clique parameters manually and estimate the 
	// rest with Derin et al.
	// See Research notebook of 19.10.2006, page 96
	// -------------------------------------------------------------------
	
	if (singleCliquesManually)
	{
		unsigned int noManParams=this->getNoClasses()-1;
		cerr << "Estimate single clique parameters manually (no Derin et al.)" << endl;
				
		_estimateSingleCliqueParameters(l, xMRFParams);
			
		Matrix<float> *Nk=N.cropR(0, noManParams-1, 0, N.rows()-1);
		Matrix<float> *Nu=N.cropR(noManParams, N.columns()-1, 0, N.rows()-1);
		Vector<float> Thetak = xMRFParams.cropR(0,noManParams-1);
		Vector<float> Pprime = ((*Nk)*Thetak);
		Pprime = P - Pprime;
	
		// Get the pseudo inverse of Nu
		Npsinv = Nu->pseudoInverse();
			
		// Get the solution
		solution = (*Npsinv)*Pprime;	
		
		// Where part in the parameter vector is used by this solution?
		derinBeg = firstPairwiseCliqueParameter();
				
		delete Nk;
		delete Nu;
	}	
	
	// -------------------------------------------------------------------
	// Estimate _ALL_ parameters using Derin et al.
	// -------------------------------------------------------------------
	
	else
	{
		cerr << "Estimate all clique parameters with Derin et al." << endl;
		
		// Get the pseudo inverse of N
		Npsinv = N.pseudoInverse();
			
		// Get the solution
		solution = (*Npsinv)*P;	
		
		// Where part in the parameter vector is used by this solution?
		derinBeg = 0;
	}	
	
	cerr << "First parameter estimated with Derin et al. at index " << derinBeg << endl;
	
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
			if (dataCountReduced < LSTSQ_MIN_DATA_WHEN_INSTABLE)
			{
				if (doMedianFilter)
				{				
					cerr << "The matrix N is unstable, trying to determine the parameters\n"
						 << "without median filtering!" << endl;
					_estimateMRFParamsLS(xl, xMRFParams, false, singleCliquesManually);
				}
				else
					ERR_THROW ("MRFSegmenter::estimateMRFParamsLS(): parameter #" 
						<< p << " is NaN");
			}
		}
	}	
	
	// Debug output
	{		
		cerr << noMRFParams << " params: ";	
		for (unsigned int p=0; p<noMRFParams; ++p)
		{
			if (p>0) cerr << ", ";
			cerr << xMRFParams[p];
		}		
	}
					
	// Clean up
	delete Npsinv;
	for (vector<FloatVector *>::iterator iter=lhs.begin(); iter!=lhs.end(); ++iter)
		delete *iter;		
	if (doMedianFilter)
		delete l;		
	
	cerr << " ;" << endl;
		
	return ok;
}

#endif


