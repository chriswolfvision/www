/***************************************************************************
 * Visualize.cpp  
 *
 * Author: Christian Wolf
 * Begin: 17.2.2003
 ***************************************************************************/

// C
#include <math.h>

// C++
#include <iostream>
#include <set>
#include <map>

// From the main module
#include <CIL.h>

// From the MATHEMATICS module
#include <Histogram.h>

// From the IMAGE PROCESSING module
#include <ImageProc.h>

// From the NUMERICAL RECIPES IN C library
#ifdef HAVE_NUM_RECIPES
#include <NumericalRecipes.h>
#endif

// From this module
#include "Visualize.h"

// The directions of the eigen vectors
enum DIRECTIONS {
	HORIZ	=0,
	VERTICAL=1,
	R_DIAG	=2,
	L_DIAG	=3,
};

// The direction orthogonal to the index direction
static int orthogonal_direction [] = {
1,0,3,2
};

// The patterns used to draw directions into an image,
// indexes are enum DIRECTIONS
static unsigned char dirpatterns[] = {
0,0,0,1,1,1,0,0,0,	// HORIZ
0,1,0,0,1,0,0,1,0,	// VERTICAL
1,0,0,0,1,0,0,0,1, 	// R_DIAG
0,0,1,0,1,0,1,0,0,	// L_DIAG
};

// *********************************************************************
// Draw the directions
// *********************************************************************

static void drawDirections (Image *im, int x, int y, int plane, int val, int directions) {
	int intensity;
	for (int j=0; j<3; ++j)
		for (int i=0; i<3; ++i) {
  			intensity = abs(val*dirpatterns[directions*9+j*3+i]);     		
     		if (intensity>255)
       			intensity=255;
			im->set (plane, x*3+i, y*3+j, intensity);
   		}
}

// *********************************************************************
// Approximate an angle by the four principal image directions.
// *********************************************************************

inline int approximate_angle (float phi) {
	if (phi<-3.0*M_PI/8.0)
	    return VERTICAL;
	else
		if (phi<-1.0*M_PI/8.0)
			return R_DIAG;
		else
			if (phi<M_PI/8.0)
				return HORIZ;
			else
				if (phi<3.0*M_PI/8.0)
					return L_DIAG;
				else
					return VERTICAL;								
}

// *********************************************************************
// Draw directions into a file
// If the intensity image is NULL, then all intensity have
// the same value (=255).
// *********************************************************************

void visualizeDirections (Image &dircodes, FloatMatrix *intensity, char *filename) {
	Image *dirimage;
 	int val;
	
	dirimage = new Image (3*dircodes.xsize, 3*dircodes.ysize, 1);
	dirimage->setPlaneToValue(PLANE_RED, 255);
	
	for (int y=0; y<dircodes.ysize; ++y) {	
		for (int x=0; x<dircodes.xsize; ++x) {

      		if (intensity==NULL)
        		val = 255;
          	else
           		val = (int) rint(intensity->get(x,y));
        		
	   		drawDirections(dirimage, x, y, PLANE_RED, val,
    	  		orthogonal_direction[dircodes.get(PLANE_RED,x,y)]);	
    	}
    }
	
	dirimage->write (filename);
	delete dirimage;
}

// *********************************************************************
// Visualizes the histogram for each class of a segmented image
// I ............ the original grayvalue image
// segI ......... the segmented image. If NULL, then the whole image
//                is considered to be a single class
// st ........... the stream where the distributions will be written
// noErosions ... perform this number of erosions on each class before
//                collecting the pixel values
// ignoreLabel .. Ignore labels with this value
// factor ....... multiplied with each bin of the histogram before print.
// *********************************************************************

void plotClassDistributions (Image *I, Image *segI, ostream &st, 
	unsigned int noErosions, int ignoreLabel, float factor) 
{
	unsigned int xs=I->xsize;
	unsigned int ys=I->ysize;
	unsigned int noClasses;
	set<unsigned char> labels;
	map<unsigned char, unsigned int> labmap;
	typedef Histogram<float> *HistogramPointer;
	HistogramPointer *h;
	Image *segISelected;
	
	if (I->nbColorPlanes()!=1)
		ERR_THROW ("visualizeClassDistributions(): color images are not yet supported!");
		
	/* -----------------------------------------
	 * The whole image is a single class 
     * -----------------------------------------
     */	
	if (segI==NULL)
	{
		Histogram<float> h (256, 0, 255, true);		
		for (unsigned int y=0; y<ys; ++y)
		for (unsigned int x=0; x<xs; ++x)
			h.add(I->get(x,y));				
		h.normalize();
		for (unsigned int g=0; g<=255; ++g)
			st << (float) h[g]*factor << endl;
		return;
	}
	
	/* -----------------------------------------
	 * multiple classes
     * -----------------------------------------
     */	
			
	// Determine the class labels
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x) 
	{
		unsigned char l=segI->get(x,y);
		if (l!=ignoreLabel)
			labels.insert(l);
	}

	//cerr << "No of labels: " << labels.size();
		
	noClasses=labels.size();
	h = new HistogramPointer [noClasses];
	
	// Assign a number to each label. The numbers begin with 1,
	// since 0 is reserved for the "background" (the image background
	// itself being actually a component).
	set<unsigned char>::iterator iter=labels.begin();
	unsigned int i=1;
	while (iter!=labels.end())
	{
		// cerr << "Label: " << (int) *iter << "," << i << endl;
		labmap[*iter] = i;
		++iter, ++i;		
	}
	
	// Travers all classes and build the histogram for each class
	for (unsigned int c=1; c<=noClasses; ++c)
	{
		unsigned char cL=0;
		
		// Get the label for this class index
		for (map<unsigned char, unsigned int>::iterator i=labmap.begin(); i!=labmap.end(); ++i)
			if (i->second==c)
				cL = i->first;
		
		// Perform the erosions if requested
		if (noErosions>0)
		{
			segISelected = new Image (xs, ys, 1);
			for (unsigned int y=0; y<ys; ++y)
			for (unsigned int x=0; x<xs; ++x)
				segISelected->set(x,y, (segI->get(x,y)==cL ? c : 0));
			erode(*segISelected, noErosions, false, c); 
		}
		else
			segISelected = segI;
	
		// Build the histogram for this class
		h[c-1] = new Histogram<float> (256, 0, 255, true);		
		for (unsigned int y=0; y<ys; ++y)
		for (unsigned int x=0; x<xs; ++x)
			if (segISelected->get(x,y)==c)
				h[c-1]->add(I->get(x,y));				
		h[c-1]->normalize();
		
		if (noErosions>0)
			delete segISelected;
	}	
	
	// Travers the gray values and print the information
	for (unsigned int g=0; g<=255; ++g)
	{
		for (unsigned int c=0; c<noClasses; ++c)
			st << (float) ((*h[c])[g])*factor << " ";
		st << endl;
	}
	
	// Clean up
	for (unsigned int c=0; c<noClasses; ++c)
		delete h[c];
	delete [] h;
}

// *********************************************************************
// Create a function plot of a gaussian distribution with 
// given parameters in the given interval
// m .............. mean
// ss ............. variance (sigma squared)
// min,max,step ... The interval
// *********************************************************************

void plotGaussian (ostream &st, float m, float ss, 
	float min, float max, float step, float factor, bool plotXAxis)
{
	float s=sqrt(ss);
	for (float g=min; g<=max; g+=step)
	{
		if (plotXAxis)
			st << g << " ";
		st << factor*1./(sqrt(2.*M_PI)*s)*exp(-0.5*(g-m)*(g-m)/ss) << endl;
	}
}

// *********************************************************************
// Create a function plot of the generalized gaussian distribution with 
// given parameters in the given interval
// m .............. mean
// ss ............. variance (sigma squared)
// p .............. shape
// min,max,step ... The interval
// *********************************************************************

#ifdef HAVE_NUM_RECIPES
void plotGeneralizedGaussian (ostream &st, float m, float ss, float p, 
	float min, float max, float step, bool plotXAxis)
{
	for (float g=min; g<=max; g+=step)
	{
		float n=sqrt(exp(numrec::gammln(3./p))/(ss*exp(numrec::gammln(1./p))));
		if (plotXAxis)
			st << g << " ";
		st << 1./(2.*exp(numrec::gammln(1./p)))*n*p*exp(-1.*pow(n*fabs((float)g-m),p)) << "\n";
	}
}
#endif
	
// *********************************************************************
// Create a function plot of the log normal distribution with 
// given parameters in the given interval
// m .............. mean
// ss ............. variance (sigma squared)
// min,max,step ... The interval
// *********************************************************************

void plotLogNormal (ostream &st, float m, float ss, 
	float min, float max, float step)
{
	float factor=1./(sqrt(2.*M_PI)*sqrt(ss));
	
	for (float g=min; g<=max; g+=step)
		st << factor/g * exp (-0.5*pow(log(g) - m, 2)/ss) << "\n";
}

