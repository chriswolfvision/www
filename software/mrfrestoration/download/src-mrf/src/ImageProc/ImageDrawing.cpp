/*********************************************************************
 * Drawing routines
 *
 * Christian Wolf, chriswolf@gmx.at
 *********************************************************************/

// C
#include <errno.h>
#include <math.h>

// From the IMAGE module
#include <FloatMatrix.h>

// From our own module
#include "ImageProc.h"

/*********************************************************************
 * Draw an ellipse into an image
 * i .... Image to draw into
 * m ... the moments
 * r .... radius
 * v .... the value to draw
 *********************************************************************/
 
#define MARK(x,y,i)		{ double I=i.get(x,y); i.set(x,y,((I-m.min) > (m.max-I))?m.min:m.max); }

void drawEllipseFromMoments (FloatMatrix &i, Moments &m, double r_square, double v)
{
	int xs=i.xsize;
	int ys=i.ysize;
	int x,ry,y;
	double t1,t2;
		
	for (int rx=0; rx<xs; ++rx)
	{
		x=rx-m.peakX;
	
		t1 = 0.5/m.xx;
		t2 = 4.0*m.xy*m.xy*x*x-4.0*m.xx*(m.yy*x*x-r_square);
		if (t2<0)
			continue;			
						
		t2=sqrt(t2);
		
		y = (int) rint(t1 * (2.0*m.xy*x + t2));
		ry = (int) rint(y + m.peakY);
		if (ry>0 && ry<ys)
			MARK (rx, ry, i);
		
		y = (int) rint(t1 * (2.0*m.xy*x - t2));
		ry = (int) rint(y + m.peakY);
		if (ry>0 && ry<ys)
			MARK (rx, ry, i);		
	}	
}

