/***************************************************************************
 * Segmentation using Markov Random Field Models:
 * The main program
 *
 * Author: Christian Wolf
 * Begin: 24.5.2005
 ***************************************************************************/
 
// Prior model numbers
enum {
	PRIORM_SINGLE_LOGISTIC=0,
	PRIORM_DOUBLE,
	PRIORM_DOUBLE_LP,
	PRIORM_SINGLE_POTTS,
	
	PRIORM_COUNT
};	 
 
 
// Observation model numbers
enum {
	OBSM_DOUBLE_GAUSS_GV=0,
	OBSM_DOUBLE_GAUSS_COLOR,
	OBSM_DOUBLE_GG_GV,
	OBSM_SINGLE_GAUSS_GV,
	OBSM_DOUBLE_NONPAR_GV,
	OBSM_SINGLE_GAUSS_CORR_SAUVOLA,
	
	OBSM_COUNT
};	

// Which prior model is compatible with which observation model?
// Rows = Priors
// Columns = Observation models
bool ModelCompatibility[][OBSM_COUNT] = 
{
	{false, false, false , true, false, true},
	{true,  true,  true,  false, true, false},
	{true,  true,  true,  false, true, false},
	{false, false, false , true, false, true}	
};
 
 /***************************************************************************
 * DEFAULT VALUES
 ***************************************************************************/

#define DEF_I_ICE 				1 
#define DEF_I_SA				300
#define DEF_I_EXPMOVE			10
#define DEF_BETA_ADJUST			1.0
#define DEF_I_ICM				10
#define DEF_C					0.97
#define DEF_T1					2.0
#define DEF_PRIORM				PRIORM_DOUBLE
#define DEF_OBSM				OBSM_DOUBLE_GAUSS_GV
#define DEF_CONSTR_PARAM		-0.5
#define DEF_COLSPC				0

// #define DEF_BG_SUB_ALPHA		0.80
#define DEF_NO_SINGLE_CLASSES	3				
#define DEF_OPTALG				OPTALG_SIMA

#define DEF(x)				"(def.: " << (x) << ")"

/***************************************************************************/

#define LABEL_FILE_NAME		"labels.pgm"

#define CHECK_NUMERIC		c=*optarg; if ((c=='-')||(c=='.')) c=optarg[1]; \
							if (!isdigit(c)) { usage(*argv); \
							cerr << "Numerical argument expected - got '" << optarg << "'\n"; exit (1); }

/***************************************************************************/

// C++
#include <new>
#include <iostream>

// From the IMAGE module
#include <Image.h>
#include <FloatMatrix.h>

// From the MRF (Markov Random Fields) module
#include "MRFSegmenter.h"
#include "MRFPriorModels/MRFSegmenterMLL8N.h"
#include "MRFPriorModels/MRFSegmenterPotts.h"

// From The ObsModels module
#include "ObsModelGauss1V.h"
#include "ObsModelGauss1VCorr.h"

// From the BLEEDTHROUGH module
#include "CommonBleedThrough.h"
#include "MRFSegmenterBleedThrough.h"
// #include "MRFSegmenterBleedThroughLP.h"
#include "ObsModelBleedThrough.h"
#include "ObsModelBleedThroughColor.h"
#ifdef HAVE_NUM_REC
#include "ObsModelBleedThroughGG.h"
#endif
#include "ObsModelBleedThroughNonParametric.h"
#include "ObsModelGaussCorrectedSauvola.h"

/**********************************************************
 * Global new handler for memory management
 **********************************************************/

static void globalNewHandler () 
{
	cerr << "***** ERROR: OUT OF MEMORY.\n"
 	     << "      PROVOKING CORE DUMP TO EASE DEBUGGING!\n";
 	int *crash_pointer=0;
  	*crash_pointer = 1;
   	exit (1); // exit should not be reached :)
}

/**********************************************************
 * The error function
 **********************************************************/

void Error (char *text)
{
	cerr << text;
}

/**********************************************************
 * The usage message
 **********************************************************/

static void usage (char *com)
{
	cerr << "usage: " << com << " [ options ] <in> <seg-out> <rest-out> <label-out>\n\n"
		 << "Options:\n"		 
		 << "[-P <model-number>]  Choose the prior MRF model " << DEF(DEF_PRIORM) << endl
		 << "                     " << PRIORM_SINGLE_LOGISTIC << " single MRF (autologistic or MLL) " << endl
		 << "                     " << PRIORM_DOUBLE << " double MRF (2x MLL / 2x Potts)" << endl
		 << "                     " << PRIORM_DOUBLE_LP << " double MRF w line process" << endl
		 << "                     " << PRIORM_SINGLE_POTTS << " single MRF (Potts)" << endl
		 << "[-O <model-numer>]   Choose the obervation model " << DEF(DEF_OBSM) << endl
		 << "                     " << OBSM_DOUBLE_GAUSS_GV << " double field - Gaussian grayvalue " << endl
		 << "                     " << OBSM_DOUBLE_GAUSS_COLOR << " double field - Gaussian color " << endl
		 << "                     " << OBSM_DOUBLE_GG_GV << " double field - Generalized Gaussian grayvalue " << endl
		 << "                     (only if numerical recipes are available)" 
		 << "                     " << OBSM_SINGLE_GAUSS_GV << " single field - Gaussian grayvalue" << endl
		 << "                     " << OBSM_DOUBLE_NONPAR_GV << " double field - non-parametric" << endl
		 << "                     " << OBSM_SINGLE_GAUSS_CORR_SAUVOLA << " single field - Gauss corrected Sauvola" << endl
		 << "[-I <iterations> ]   iterations ICE " << DEF(DEF_I_ICE) << endl
		 << "[-i <iterations> ]   expansion move:      n.o. cycles " << DEF(DEF_I_EXPMOVE) << endl
		 << "                     simulated annealing: n.o. iter. gibbs sampler " << DEF(DEF_I_SA) << endl
		 << "[-h x,x,x,...]       force values for the MRF parameters\n"
		 << "[-a <adjust>]        adj. factor for beta parameter (graph cut only) " << DEF(DEF_BETA_ADJUST) << endl
		 << "[-g <iterations> ]   iterations ICM (greedy mode) " << DEF(DEF_I_ICM) << endl
		 << "[-t <t1> ]           start temperature simulated annealing " << DEF(DEF_T1) << endl
		 << "[-c <c> ]            parameter c simulated annealing " << DEF(DEF_C) << endl
		 << "[-Z <o.algo>]        Chose the optimization algorithm " << DEF(DEF_OPTALG) << endl
		 << "                     " << OPTALG_SIMA << " simulated annealing and/or ICM" << endl
		 << "                     " << OPTALG_GRCT << " graph cut" << endl
		 << "[-B]                 Bleedthrough: separate result images into recto and verso " << endl
		 << "                     part, even if a single MRF model is chosen." << endl
		 << "[-b]                 do use blurring in the model (not all models)" << endl
		 << "[-D]                 use Derin et al. for a subset of parameters only" << endl
		 << "[-M]                 do median filter before least squares estimation" << endl
		 << "[-G]                 do Gauss filter before k-means" << endl
		 << "[-s]                 do Sauvola et al.'s binarization instead of k-means" << endl
		 << "                     (2 class segmentation, Potts model only)." << endl
		 << "[-m]                 use the background substraction mask" << endl
		 << "[-o color-code]      color space and distance for kmeans init. " << DEF(DEF_COLSPC) << endl
		 << "                     (0=rgb, 1=hsv, 2=L*u*v*, 3=L*a*b*, 4=L*a*b*+CIE94)\n"
		 << "[-C <n.o.cl>]        number of classes (single MRF only) " << DEF(DEF_NO_SINGLE_CLASSES) 
		 << endl
		 << "[-V <factor>]        constraint parm. for verso px (double MRF) " << DEF(DEF_CONSTR_PARAM) 
		 << endl
		 << "[-R]                 use gray.val. for the recognition of the recto class " << endl
		 << "[-l]                 write the label field into the file '" << LABEL_FILE_NAME << "'\n"
		 << "[-d <option>]        write intermediate debug images to disk or not: \n"
		 << "                     " << DBG_WIMAGES_NO << ": do not write\n"
		 << "                     " << DBG_WIMAGES_SEG << ": write visible segmentation output\n"
		 << "                     " << DBG_WIMAGES_LABELS << ": write label images\n"		 
		 << "                     " << DBG_WIMAGES_RESTORE << ": write restored images\n"
		 ;		 
}

/**********************************************************
 * The main program
 **********************************************************/

int	main (int argc,	char **argv) 
{
	char *ifname, *ofname_seg, *ofname_res, *ofname_lab, buf[512];
	int c;
	int   optIICE 		= DEF_I_ICE, 
		  optISA 		= DEF_I_SA,
		  optIExpMove	= DEF_I_EXPMOVE,
		  optIICM		= DEF_I_ICM,
		  optPRIORM	 	= DEF_PRIORM,
		  optOBSM		= DEF_OBSM,
		  optNOSingleClasses = DEF_NO_SINGLE_CLASSES,
		  optColSpc		= DEF_COLSPC;
	float optT1 		= DEF_T1, 
		  optC			= DEF_C,
		  optConstraintParamFactor = DEF_CONSTR_PARAM,
		  optBetaAdjust = DEF_BETA_ADJUST;
	// float optBGSubAlpha	= DEF_BG_SUB_ALPHA;
	bool  optDoUseMask	= false,
		  optDoBlur = false,
		  optDoWriteLabelField = false,		  
		  optDoDerinForAll = true,
		  optDoMedianFilter = false,
		  optDoGaussFilterBeforeKmeans = false,
		  optUseGV4RectoRecognition = false,
		  optForceBleedThroughVisualization = false,
		  optDoSauvolaInsteadOfKmeans = false;
	char *optOptAlg = (char *) DEF_OPTALG,
		 *optForceMRFParams = NULL;
	DebugWriteImageOptions optDebugImageOptions = DBG_WIMAGES_NO;
	ObsModel<Image> *obsModel;
	MRFSegmenter<Image> *ms;
	Image *eventualColorObs;
	unsigned int noIterations;
	char *cp;
	//Image *mask;
	
	// Install new handler for memory management
	set_new_handler	(globalNewHandler);
	
	//	***************************************************
	//	Arguments
	//	***************************************************

    while ((c =	getopt (argc, argv,	"mMbo:O:P:a:h:Z:I:i:t:c:g:lC:V:d:DRGBs")) != EOF) {

		switch (c) {

			case 'i':
				CHECK_NUMERIC;
				optISA = atol(optarg);
				optIExpMove = optISA;
				break;
				
			case 'I':
				CHECK_NUMERIC;
				optIICE = atol(optarg);
				break;
				
			case 'g':
				CHECK_NUMERIC;
				optIICM = atol(optarg);
				break;
				
			case 't':
				CHECK_NUMERIC;
				optT1 = atol(optarg);
				break;
				
			case 'c':
				CHECK_NUMERIC;
				optC = atof(optarg);
				break;
				
			case 'a':
				CHECK_NUMERIC;
				optBetaAdjust = atof(optarg);
				break;
			
			case 'O':
				CHECK_NUMERIC;
				optOBSM = atol(optarg);
				break;	
				
			case 'o':
				CHECK_NUMERIC;
				optColSpc = atof(optarg);
				break;
				
			case 'P':	
				CHECK_NUMERIC;
				optPRIORM = atol(optarg);
				break;	
				
			case 'Z':
				optOptAlg = optarg;
				break;
				
			case 'h':
				optForceMRFParams = optarg;
				break;
				
			case 'm':
				optDoUseMask=true;
				break;
				
			case 'B':
				optForceBleedThroughVisualization=true;
				break;
				
			case 'b':
				optDoBlur=true;
				break;
				
			case 'D':
				optDoDerinForAll=false;
				break;
				
			case 'M':
				optDoMedianFilter=true;
				break;
				
			case 's':
				optDoSauvolaInsteadOfKmeans=true;
				break;
				
			case 'G':
				optDoGaussFilterBeforeKmeans=true;
				break;
				
			case 'l':
				optDoWriteLabelField=true;
				break;
			
			case 'd':
				CHECK_NUMERIC;
				optDebugImageOptions=(DebugWriteImageOptions) atol(optarg);				
				break;
				
			case 'C':
				CHECK_NUMERIC;
				optNOSingleClasses = atol(optarg);
				break;	
				
			case 'V':
				CHECK_NUMERIC;
				optConstraintParamFactor = atof(optarg);
				break;
			
			case 'R':				
				optUseGV4RectoRecognition = true;
				break;
				
			case '?':
				usage (*argv);				
				cerr << "Syntax error!!\n";
				exit (1);
				
			default:
				usage (*argv);
				exit (1);
		}
	}

	if (argc-optind != 4)
	{		
		usage (*argv);
 		cerr << "\nPlease specify an input and 3 output filenames!\n";
		exit(1);
	}
	ifname  = argv[optind];
	ofname_seg = argv[optind+1];
	ofname_res = argv[optind+2];
	ofname_lab = argv[optind+3];
	
	if (!ModelCompatibility[optPRIORM][optOBSM])
	{	
		usage (*argv);
		cerr << "The prior model is not compatible with the observation model!\n";
		exit(1);
	}
	
	if ((optColSpc!=0) && (optPRIORM==PRIORM_SINGLE_LOGISTIC))
	{
		usage (*argv);
		cerr << "Non RGB color spaces are currently not supported for the logistic model!\n";
		exit (1);
	}
	
	if (optDoSauvolaInsteadOfKmeans && (optNOSingleClasses!=2))
	{
		usage (*argv);
		cerr << "Number of classes must be 2 for Sauvola et al.'s binarization\n";
		exit(1);
	}
	
	if (optDoSauvolaInsteadOfKmeans && (optPRIORM!=PRIORM_SINGLE_POTTS))
	{
		usage (*argv);
		cerr << "The single Potts model must be chosen for Sauvola et al.'s binarization\n";
		exit(1);
	}
	
	cerr << "############################################################"  << endl
	     << "# Image segmentation with Markov random fields" << endl
		 << "#  <- " << ifname << endl
		 << "#  -> seg: " << ofname_seg << endl
		 << "#  -> res: " << ofname_res << endl
		 << "#  -> lab: " << ofname_lab << endl
		 << "############################################################" << endl;
	
	// Main routine
	try
	{
		Image observation (ifname);
		Image labRecto (observation.xsize,observation.ysize,1);
		Image labVerso (observation.xsize,observation.ysize,1);				
		
		// -----------------------------------------------------------------------
		// Choice of the model and the method
		// -----------------------------------------------------------------------
				
		// A single MRF 
		/*
		ObsModelGauss1VCorr<Image> obsmodel;
		MRFSegmenterAL8N<Image> ms (&observation, &labRecto, &obsmodel);				
		ms.initKMeans(2, (vector<Image::PixelType> *) NULL);
		*/
		
		// In case the observation model is grayscale one, we will need
		// to convert an eventual color image into grayscale. In order 
		// to be able to use the color information nevertheless for
		// the kmeans initialisation, we keep the original color image
		if (observation.nbColorPlanes()==3)
			eventualColorObs = new Image (observation);
		else
			eventualColorObs = NULL;
			
		// Choose the observation model	
		switch (optOBSM)
		{
			case OBSM_DOUBLE_GAUSS_GV:
				obsModel = new ObsModelBleedThrough<Image> (optDoBlur);
				observation.convertRGB2GrayScale();	
				break;
			case OBSM_DOUBLE_GAUSS_COLOR:
				obsModel = new ObsModelBleedThroughColor<Image> (optDoBlur);
				break;
#ifdef HAVE_NUM_REC
			case OBSM_DOUBLE_GG_GV:
				obsModel = new ObsModelBleedThroughGG<Image> (optDoBlur);
				observation.convertRGB2GrayScale();	
#endif
				break;
			case OBSM_SINGLE_GAUSS_GV:	
				obsModel = new ObsModelGauss1VCorr<Image> (optDoBlur);
				observation.convertRGB2GrayScale();	
				break;
			case OBSM_DOUBLE_NONPAR_GV:
				obsModel = new ObsModelBleedThroughNonParametric<Image> (optDoBlur);
				observation.convertRGB2GrayScale();	
				break;
			case OBSM_SINGLE_GAUSS_CORR_SAUVOLA:
				obsModel = new ObsModelGaussCorrectedSauvola<Image>;
				observation.convertRGB2GrayScale();	
				break;
			default:
				usage (*argv);
				cerr << "Invalid observation model!\n";
				exit (1);
		}
			
		// -----------------------------------------------------------------------
		// Background substraction and creation of the corresponding pixel mask
		// -----------------------------------------------------------------------
		
		{	/*
			Histogram<int> *hist;
			unsigned char v;
			float mean1, mean2;
			int kopt;
			
			// FIRST SEPARATION BY OTSU
			hist = buildHistogram<int>(observation);
			kopt = hist->getThresholdValueFisher (0, 256, &mean1, &mean2);
			
			// Change the optimal threshold -> we want to be sur that no text
			// pixel is classified as background
			kopt = (int) (mean1 + (optBGSubAlpha*((double)kopt - mean2)));
			mask = new Image (observation);
    		thresholdFixed (*mask, kopt);
    		DBSAVE(*mask,ofname,"_mask_firstpass.pgm","");
    		delete hist;
    		
    		reverseVideo (*mask);
    		dilate (*mask, 1, false, 255); 
    		DBSAVE(*mask,ofname,"_mask.pgm","");
    		*/		
		}
					
		// -----------------------------------------------------------------------
		// MRF optimization
		// -----------------------------------------------------------------------
		
		// Choose the prior model	
		switch (optPRIORM)
		{
			case PRIORM_SINGLE_LOGISTIC:
				if (optNOSingleClasses==2)
					ms = new MRFSegmenterAL8N<Image> (&observation, &labRecto, obsModel,
						 optForceMRFParams);
				else
					ms = new MRFSegmenterMLL8N<Image> (&observation, &labRecto, obsModel, 
						optForceMRFParams, optNOSingleClasses, optDoGaussFilterBeforeKmeans);
				break;
			
			case PRIORM_SINGLE_POTTS:
				MRFSegmenterPotts<Image> *mspotts;
				ms = new MRFSegmenterPotts<Image> (&observation, &labRecto, obsModel, 
						eventualColorObs, optForceMRFParams, optNOSingleClasses, 
						optDoGaussFilterBeforeKmeans,
						optBetaAdjust, ColorSpaceCode (optColSpc));
				mspotts = (MRFSegmenterPotts<Image> *) ms;
				mspotts->doSauvolaInsteadOfKmeans = optDoSauvolaInsteadOfKmeans;
				break;
			
			case PRIORM_DOUBLE:
				ms = new MRFSegmenterBleedThrough<Image>  (
					&observation, &labRecto, &labVerso, obsModel, 
					eventualColorObs,
					optForceMRFParams, optDoUseMask,	
					optUseGV4RectoRecognition, optDoGaussFilterBeforeKmeans, 
					optBetaAdjust, optConstraintParamFactor,
					ColorSpaceCode (optColSpc));
				break;
				
			case PRIORM_DOUBLE_LP:
				ERR_THROW("line processes are currently not supported "
					"(must be adapted to the potts model).");
				
// 				ms = new MRFSegmenterBleedThroughLP<Image>  (
// 					&observation, &labRecto, &labVerso, obsModel, optDoUseMask,	
// 					optUseGV4RectoRecognition, optConstraintParamFactor);
				break;

			default:
				usage (*argv);
				cerr << "Invalid prior model!\n";
				exit (1);
		}
				
		ms->init();
		
		// Cut the suffix for the intermediate files
		cp = ofname_seg+(strlen(ofname_seg)-4);
		if ((strcmp(cp,".pgm")==0) || (strcmp(cp,".ppm")==0))
		{
			strcpy(buf,ofname_seg);
			strcat(buf+strlen(buf)-4,"_it");
		}
		else
			sprintf(buf,"%s_it",ofname_seg);
				
		noIterations=optISA;
		if (strcmp(optOptAlg,OPTALG_GRCT)==0)
			noIterations=optIExpMove;
				
		ms->iteratedConditionalEstimation (optIICE, optT1, optC, noIterations, optIICM, 
			optDoMedianFilter, !optDoDerinForAll, 
			optOptAlg, optDebugImageOptions,
			(optDebugImageOptions==DBG_WIMAGES_NO ? NULL : buf));		
			
		cerr << "Image segmented." << endl;		
	
		/***************************************************************************
		 * OUTPUT - LABELS
		 ***************************************************************************/
				
		if (optDoWriteLabelField)
		{
			ms->writeLabelField((char *) LABEL_FILE_NAME, DBG_WIMAGES_LABELS);
			cerr << "Label field written to " << LABEL_FILE_NAME << "." << endl;
		}
		
		/***************************************************************************
		 * OUTPUT - VISUALIZATION OF THE SEGMENTATION:
		 * The normal case: just call the visualization method of the 
		 * segmenter object
		 ***************************************************************************/
		 
		if ((optOBSM!=OBSM_SINGLE_GAUSS_GV) || (!optForceBleedThroughVisualization))		
		{	
			cerr << "Calling segmenter's visualization method.\n";	
			ms->writeLabelField(ofname_seg, DBG_WIMAGES_SEG);	
			ms->writeLabelField(ofname_lab, DBG_WIMAGES_LABELS);
		}
		
		/***************************************************************************
		 * OUTPUT - VISUALIZATION OF THE SEGMENTATION:
		 * A special case: a single MRF has been used to perform bleed through rem.		 
		 * Manually separate the segmented image into a recto and a verso part
		 ***************************************************************************/
		
		else
		{				
			Image *lab	  = ms->getLab(),
			      *labv   = new Image (  lab->xsize, lab->ysize, 1),
			      *labtot = new Image (2*lab->xsize, lab->ysize, 1),
			      *labout = new Image (2*lab->xsize, lab->ysize, 1),
			      *vi     = obsModel->visualizeSegmentation();
			int grayvalues[3];
			bool complete=true;																		
			
			cerr << "Manual visualization.\n";
			
			// Recognize the labels and reorder them,			
			reOrderLabels(lab, getBGLabel(lab));
			
			// Write the labels into the label output file
			ms->writeLabelField(ofname_lab, DBG_WIMAGES_LABELS);
			
			// separate them into recto and verso	
			recto2RectoVerso(lab, labv, &observation, optUseGV4RectoRecognition, -1);
			labtot->paste (*lab,0,0);
			labtot->paste (*labv,lab->xsize,0);					
						
			/* Determine the colors of the classes from the 
			 * segmented and from the visualized image
			 */
			for (int i=0; i<3; ++i)
				grayvalues[i]=-1;			 
			for (int y=0; y<labtot->ysize; ++y)
			for (int x=0; x<labtot->xsize; ++x)
			{
				unsigned char l=labtot->get(PLANE_RED,x,y);
				if (l>2)
					ERR_THROW ("Internal error in main() - (1).");
					
				// The label depends on where we are in the image
				// (the image has already been separated into
				//  the recto and verso part)
				if (x<lab->xsize)
					grayvalues[(l?CL_IND_RECTO:CL_IND_BG)]=vi->get(PLANE_RED,x,y);
				else
					grayvalues[(l?CL_IND_VERSO:CL_IND_BG)]=vi->get(PLANE_RED,x-lab->xsize,y);
				
				complete=true;
				for (int i=0; i<3; ++i)
					if (grayvalues[i]<0)
						complete=false;
				if (complete)
					break;
			} 
			if (!complete)
				ERR_THROW ("Cannot separate single MRF-segmented image: "
					"not all labels are present.");		
							
			for (int i=0; i<3; ++i)
				cerr << "label colors: " << i << " = " << grayvalues[i] << endl;				
						
			// Replace the labels by the class colors
			// flip the right (verso) part of the image
			for (int y=0; y<lab->ysize; ++y)
			for (int x=0; x<lab->xsize; ++x)				
			{
 				unsigned int label;
 				label = (labtot->get(PLANE_RED,x,y) ? CL_IND_RECTO : CL_IND_BG);
				labout->set (PLANE_RED,x,y,grayvalues[label]);
				label = (labtot->get(PLANE_RED,x+lab->xsize,y) ? CL_IND_RECTO : CL_IND_BG);
				labout->set (PLANE_RED,labtot->xsize-x-1,y,grayvalues[label]);				
				
			}
			labout->write (ofname_seg);
						
			// Clean up
			delete vi;
			delete labv;
	 		delete labtot;
	 		delete labout;
		}		
		
		/***************************************************************************
		 * OUTPUT - RESTORATION
		 ***************************************************************************/
		
		if ((optPRIORM==PRIORM_SINGLE_LOGISTIC) || 
			(optPRIORM==PRIORM_SINGLE_POTTS))
			cerr << "No image restoration with single MRF.\n";
		else
		{
			MRFSegmenterBleedThrough<Image> *ms_cast=(MRFSegmenterBleedThrough<Image> *)ms;
			HRULE;
			ms_cast->restoreFromSegmentation();
			observation.write(ofname_res);
			cerr << "Image restored to '" << ofname_res << "'." << endl;
		}					
	}
	
	// Error handling
	catch (EError &e)
	{
		cerr << e.message.c_str() << endl;
		exit (1);
	}
	
	// Cleanup 
	delete obsModel;
	delete ms;
	if (eventualColorObs!=NULL)
		delete eventualColorObs;
	
    return 0;
}
