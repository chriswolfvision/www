/*********************************************************************
 * Astronomical image processing
 *
 * Christian Wolf, chriswolf@gmx.at
 *********************************************************************/

// C
#include <errno.h>
#include <math.h>

// From the MATH module
#include <MathFuncs.h>

// From the IMAGE module
#include <FloatMatrix.h>

// From our own module
#include "ImageProc.h"

/*********************************************************************
 * calculate the first order and second order moments
 *********************************************************************/
 
static Moments calcMoments (FloatMatrix &i)
{
	int xs=i.xsize;
	int ys=i.ysize;
	double mx,my,mxx,myy,mxy;
	double max, min;
	double sumI;
	int peakX,peakY;
	Moments rv;
	
	max=min=i.get(0,0);
	mx=my=mxx=myy=mxy=0;
	peakX=0,peakY=0;
	sumI=0;
	
	// First pass
	for (int x=0; x<xs; ++x)
	for (int y=0; y<ys; ++y)
	{
		double I=i.get(x,y);
		
		sumI += I;
		if (I>max) 
		{
			max=I;
			peakX=x;
			peakY=y;			
		}
		if (I<min) 
			min=I;
		
		mx += x*I;
		my += y*I;			
	}
	
	// Second pass
	for (int x=0; x<xs; ++x)
	for (int y=0; y<ys; ++y)
	{
		double I=i.get(x,y);
		int xn=x-peakX;
		int yn=y-peakY;
		mxx += xn*xn*I;
		mxy += xn*yn*I;
		myy += yn*yn*I;		 		
	}
	
	rv.x = mx / sumI;
	rv.y = my / sumI;	
	rv.xx = mxx / sumI;
	rv.xy = mxy / sumI;
	rv.yy = myy / sumI;
	rv.peakX = peakX;
	rv.peakY = peakY;
	rv.max = max;
	rv.min = min;
	
	return rv;
} 

/*********************************************************************
 * Calculate a couple of shape features from a single band 
 * galaxy image
 *********************************************************************/
 
void galaxyShapeSingleBand (FloatMatrix &i, FloatMatrix &o) 
{
	// int xs=i.xsize;
	// int ys=i.ysize;
	Moments m;
	double r;
		
	m = calcMoments(i);
	
	cerr << "x=" << m.x << " "
	     << "y=" << m.y << " "
		 << "xx=" << m.xx << " "
		 << "xy=" << m.xy << " "
		 << "yy=" << m.yy << " " << endl
		 << "min=" << m.min << " "
		 << "max=" << m.max << " "
		 << "peakX=" << m.peakX << " "
		 << "peakY=" << m.peakY << endl
		 ;
		 
	// -------------------------------------------------------------
	// Draw the ellipse
	
	o.copy(i);	
	
	r=(m.xx+m.xy+m.yy);
	r=r*sqrt(r);
	r=m.xx*m.xx;
	cerr << "radius=" << sqrt(r) << endl;
	
	drawEllipseFromMoments (o, m, r, m.max);
	
	// -------------------------------------------------------------
	// Calculate the radial profile
	
	/*
	Matrix<float> covinv(2,2);
	Vector<float> center(2), p(2);
	
	// invert the covariance matrix
	covinv.set (0,0, m.xx);
	covinv.set (1,1, m.yy);
	covinv.set (0,1, m.xy);
	covinv.set (1,0, m.xy);
	// covinv.invert();
	
	// the mode
	center[0] = m.peakX;
	center[1] = m.peakY;
	
	for (int x=0; x<xs; ++x)
	for (int y=0; y<ys; ++y)
	{	
		double dist;
		
		// Get the radius of this point, which can be 
		// calculated using the Mahalanobis distance:
		p[0] = x;
		p[1] = y;
		// dist = Mahalanobis (center, p, covinv);
		
		cout << dist << " " << i.get(x,y) << endl;	
	}	
	*/	 	
		 
	/*
	// Filter the image
	FloatMatrix ffm (xs, ys);
	FilterMask mask (FLTTPE_GAUSS_3X3);
    filterFloatMatrix (i, ffm, mask);	
	filterFloatMatrix (ffm, i, mask);	
	
	m = calcMoments(i);
	
	cerr << "x=" << m.x << " "
	     << "y=" << m.y << " "
		 << "xx=" << m.xx << " "
		 << "xy=" << m.xy << " "
		 << "yy=" << m.yy << " "
		 << "max=" << m.max << " "
		 << "peakX=" << m.peakX << " "
		 << "peakY=" << m.peakY << endl
		 ;
	*/
}
