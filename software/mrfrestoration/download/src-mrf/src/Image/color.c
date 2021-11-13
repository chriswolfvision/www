/**********************************************************
 * color.c
 * Functions for color calculations
 *
 * Christian Wolf
 * Beginn: 9.7.1997
 **********************************************************/

// C
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include <visualbug.h>

#include "color.h"

/************************************************************************
 * Static help variables
 ************************************************************************/

static double MATRIX_XYZ_TO_RGB[] =	
{
	 3.240479, -1.537150, -0.498535,
	-0.969256,	1.875992,  0.041556,
	 0.055648, -0.204043,  1.057311
};

static double MATRIX_RGB_TO_XYZ[] =	
{
	0.412453, 0.357580,	0.180423,
	0.212671, 0.715160,	0.072169,
	0.019334, 0.119193,	0.950227
};

static double REF_WHITE[] =	{0.950456, 1, 1.088755};
static double un_ref, vn_ref;
double maxDE;

/************************************************************************
 * Help functions
 ************************************************************************/

static double u_from_xyz (double x, double y, double z) 
{
	double tmp;
	tmp	= x	+ 15*y + 3*z;

	if (tmp==0)
		return 0;
	else
		return 4*x/tmp;
}

static double v_from_xyz (double x, double y, double z) 
{
	double tmp;
	tmp	= x	+ 15*y + 3*z;

	if (tmp==0)
		return 0;
	else
		return 9*y/tmp;
}

static double f_ab(double I) 
{
	if (I>0.008856)
		return pow(I,1.0/3.0);
	else
		return 7.787*I + 16.0/116.0;
}



/************************************************************************
 * The euclidean color difference
 ************************************************************************/

double diffColorEuclidean (double l1, double u1, double v1, double l2, double u2, double v2) 
{
	return sqrt	((l1-l2)*(l1-l2)+(u1-u2)*(u1-u2)+(v1-v2)*(v1-v2));
}

/************************************************************************
 * Compare colors
 ************************************************************************/
double cmpColorS (double l1, double	u1,	double v1, double l2, double u2, double v2) 
	{
	return diffColorEuclidean	(l1,u1,v1,l2,u2,v2);
}

/************************************************************************
 * Initialize the color system
 ************************************************************************/

void init_color	() 
{
   un_ref =	u_from_xyz (REF_WHITE[0],REF_WHITE[1],REF_WHITE[2]);
   vn_ref =	v_from_xyz (REF_WHITE[0],REF_WHITE[1],REF_WHITE[2]);

   maxDE = diffColorEuclidean	(MAXVALUE_L,MAXVALUE_U,MAXVALUE_V,
		MINVALUE_L,MINVALUE_U,MINVALUE_V);
}

/************************************************************************
 * Convert RGB to grayscale
 ************************************************************************/

byte rgb2grayscale (byte r, byte g, byte b) 
{
	int gl  = (int)(rint (0.299*(double) r +
	                      0.587*(double) g +
	                      0.114*(double) b)/1);
	return gl>255 ? 255 : gl;
}

/************************************************************************
 * Converts	a color	from the RGB model to the LUV Model
 ************************************************************************/

void rgb2luv (double r,	double g, double b,
					double *l, double *u, double *v) {


#ifdef DUMMY

// The Algorithm taken from	the	Java library,
// to cross	check my own algo.

float M11=(float) 0.431	;
float M12=(float) 0.342	;
float M13=(float) 0.178	;
float M21=(float) 0.222	;
float M22=(float) 0.707	;
float M23=(float) 0.071	;
float M31=(float) 0.020	;
float M32=(float) 0.130	;
float M33=(float) 0.939	;
float Un=(float) 0.19793943	;  // 4*0.951/(0.951+15+1.089*3)
float Vn=(float) 0.46831096	;  // 9.0/(0.951+15+1.089*3)

	float x, y,	z ;
	float tu,tv	;
	float tmp ;

	r*=255;
	g*=255;
	b*=255;

	x=(M11*r+M12*g+M13*b)/(float)255.0 ;
	y=(M21*r+M22*g+M23*b)/(float)255.0 ;
	z=(M31*r+M32*g+M33*b)/(float)255.0 ;
	tmp=x+15*y+3*z ;
	if(tmp==0.0) tmp=(float)1.0	;
	tu=(float)4.0*x/tmp	;
	tv=(float)9.0*y/tmp	;
	if(y>0.008856) *l=116*(float)pow(y,1.0/3)-16;
	else *l=(float)903.3*y;
	*u=13*(*l)*(tu-Un) ;
	*v=13*(*l)*(tv-Vn) ;
#endif


	double x,y,z,uc,vc;

	// values RGB 0-1 are expected!!!!!

	#define	A(x,y)		MATRIX_RGB_TO_XYZ[y*3+x]
	x =	A(0,0)*r + A(1,0)*g	+ A(2,0)*b;
	y =	A(0,1)*r + A(1,1)*g	+ A(2,1)*b;
	z =	A(0,2)*r + A(1,2)*g	+ A(2,2)*b;
	#undef A

	if (x+y+z == 0)	{
		*l = 0;
		*u = 0;
		*v = 0;
		return;
	}

   /* From my color	*/
   uc =	u_from_xyz (x,y,z);
   vc =	v_from_xyz (x,y,z);

   /* Computing	the	lightness */
   if (y/REF_WHITE[1] <= 0.008856)
		*l = 903.3 * (y/REF_WHITE[1]);
	else
		*l = (116 *	pow	((y	/ REF_WHITE[1]),(1.0/3.0)))	- 16;

   *u =	13 * *l	* (uc -	un_ref);
   *v =	13 * *l	* (vc -	vn_ref);
}

/************************************************************************
 * Converts	a color	from the LUV Model,	where just * the U and V
 * coordinates are used, to	the	RGB	Model.
 ************************************************************************/

void luv2rgbscaled (double l, double u,	double v, byte *br,	byte *bg, byte *bb)	
{
	double r,g,b;

	luv2rgb	(l,u,v,&r,&g,&b);

	r *= MAX_RED;
	g *= MAX_GREEN;
	b *= MAX_BLUE;

	if (r<0) r=0;
	if (g<0) g=0;
	if (b<0) b=0;
	if (r>255) *br=255;	else *br = (int) r/1;
	if (g>255) *bg=255;	else *bg = (int) g/1;
	if (b>255) *bb=255;	else *bb = (int) b/1;
}

void luv2rgb (double l,	double u, double v,	double *r, double *g, double *b) 
{
	double x,y,z,uc,vc;

	// values RGB 0-1 are generated!!!!!

	if (l <= 7.9996248)
		y =	(l * REF_WHITE[1]) / 903.3;
	else
		y =	REF_WHITE[1] * pow(((l+16)/116),3);

	if (l==0) {
		x=y=z=0;
	}
	else {

		uc = u/(13*l)+un_ref;
		vc = v/(13*l)+vn_ref;

		x =	(9*y*uc)/(4*vc);

		/* z = (4*x-uc*x-15*uc*y)/(3*uc);				Inge Nr. 1 */
		/* z = (4*x-9*y*uc*uc/(4*vc)-15*y*uc)/(3*uc);	Inge Nr. 2 */

		z =	(9*y/vc	- x	- 15*y)/3;
	}

	#define	A(x,y)		MATRIX_XYZ_TO_RGB[y*3+x]
	*r = A(0,0)*x +	A(1,0)*y + A(2,0)*z;
	*g = A(0,1)*x +	A(1,1)*y + A(2,1)*z;
	*b = A(0,2)*x +	A(1,2)*y + A(2,2)*z;
	#undef A

	if (*r > 1)	*r = 1;
	if (*g > 1)	*g = 1;
	if (*b > 1)	*b = 1;

	if (*r < 0)	*r = 0;
	if (*g < 0)	*g = 0;
	if (*b < 0)	*b = 0;
}


/************************************************************************
 * Converts	a color	from the RGB model to the Lab Model
 * INPUTVALUES in [0,1]
 ************************************************************************/

#define labf(x)	((x>0.008856)?(pow(x,1./3.)):7.787*x+16./116.)

void rgb2lab (double r,	double g, double blue, double *l, double *a, double *bstar) 
{
	double x,y,z;

	#define	A(x,y)		MATRIX_RGB_TO_XYZ[y*3+x]
	x =	A(0,0)*r + A(1,0)*g	+ A(2,0)*blue;
	y =	A(0,1)*r + A(1,1)*g	+ A(2,1)*blue;
	z =	A(0,2)*r + A(1,2)*g	+ A(2,2)*blue;
	#undef A

	if (x+y+z == 0)	
	{
		*l = 0;
		*a = 0;
		*bstar = 0;
		return;
	}

// OLD VERSION (not correct)
// 	*l = 116.0*f_ab(y / REF_WHITE[2]) - 16.0;
//    	*a = 500.0 * (f_ab(x/REF_WHITE[1])-f_ab(y/REF_WHITE[2]));
//    	*b_ = 500.0 * (f_ab(y/REF_WHITE[2])-f_ab(z/REF_WHITE[3]));

	
	if (y/REF_WHITE[1] > 0.008856)
		*l = 116. * pow((y/REF_WHITE[1]),1./3.)  - 16.;
	else
		*l = 903.3 * y/REF_WHITE[1];
	*a = 500.*(labf(x/REF_WHITE[0]) - labf(y/REF_WHITE[1]));
	*bstar = 200.*(labf(y/REF_WHITE[1]) - labf(z/REF_WHITE[2]));		
}

/************************************************************************
 * OUTPUT RGB VALUES in [0,1]
 ************************************************************************/

void lab2rgb (double l,	double a, double bstar, double *r, double *g, double *blue) 
{
	double x,y,z;
	double p = (l + 16.) / 116.;
	
	x = REF_WHITE[0] * pow((p+a/500.),3.);
	y = REF_WHITE[1] * p*p*p;
	z = REF_WHITE[2] * pow((p - bstar/200.),3.);
   
	#define	A(x,y)		MATRIX_XYZ_TO_RGB[y*3+x]
	*r    = A(0,0)*x +	A(1,0)*y + A(2,0)*z;
	*g    = A(0,1)*x +	A(1,1)*y + A(2,1)*z;
	*blue = A(0,2)*x +	A(1,2)*y + A(2,2)*z;
	#undef A

	if (*r >    1)	*r = 1;
	if (*g >    1)	*g = 1;
	if (*blue > 1)	*blue = 1;

	if (*r <    0)	*r = 0;
	if (*g <    0)	*g = 0;
	if (*blue < 0)	*blue = 0;
}

/************************************************************************
 * The CIE94 Color difference formula
 ************************************************************************/

#define KL 1.0
#define KC 1.0
#define KH 1.0

double diffLabCIE94 (double l1, double a1, double b1, double l2, double a2, double b2) 
{
	double Cab1, Cab2, Cab, H1, H2, SL, SC, SH,
		   t1, t2, t3;
	
	Cab1 = sqrt(a1*a1+b1*b1);
	Cab2 = sqrt(a2*a2+b2*b2);
	Cab = sqrt(Cab1*Cab2);
	if (a1!=0)
		H1 = atan(b1/a1);
	else
		if (b1>0)
			H1 = M_PI/2.0;
		else
			H1 = M_PI/-2.0;
	if (a2!=0)
		H2 = atan(b2/a2);
	else
		if (b2>0)
			H2 = M_PI/2.0;
		else
			H2 = M_PI/-2.0;
	SL=1.0;
	SC=1.0 + 0.045*Cab;
	SH=1.0 + 0.015*Cab;

	t1=(l1-l2)/(KL*SL);
	t2=(Cab1-Cab2)/(KC*SC);
	t3=(H1-H2)/(KH*SH);
	
	return sqrt(t1*t1 + t2*t2 + t3*t3);	
}

/************************************************************************
 * Converts	a color	to the HSV Model
 * This routine has been taken (and modified) from 
 * http://www.cs.rit.edu/~ncs/color/t_convert.html
 *
 * THE RGB VALUES ARE SUPPOSED TO BE BETWEEN 0 AND 1!
 * OUTPUT VALUES ARE
 * H=[0,360]
 * S=[0,1]
 * V=[0,1]
 *
 ************************************************************************/

void rgb2hsv (double r,	double g, double b, double *h, double *s, double *v)
{
	float min, max, delta;
	char argmax, argmin;

	min=r;
	argmin='r';
	if (g<min) 
	{
		min=g;
		argmin='g';
	}
	if (b<min) 
	{
		min=b;
		argmin='b';
	}
	max=r;
	argmax='r';
	if (g>max) 
	{
		max=g;
		argmax='g';
	}
	if (b>max) 
	{
		max=b;
		argmax='b';
	}
	
	*v = max;				// v
	delta = max - min;
	
// 	printf ("RGB values: %f, %f, %f, %f, %f, %f\n",r,g,b,255*r, 255*g, 255*b);
// 	printf ("Max = %f, min = %f\n",max, min);

	if( max != 0 )
		*s = delta / max;		// s
	else 
	{
		// r = g = b = 0		// s = 0, v is undefined
		*s = 0;
		*h = -1;
		return;
	}

	if(argmax=='r')
		*h = ( g - b ) / delta;		// between yellow & magenta
	else 
		if(argmax=='g')
			*h = 2 + ( b - r ) / delta;	// between cyan & yellow
		else
			*h = 4 + ( r - g ) / delta;	// between magenta & cyan
		

	*h *= 60;				// degrees
	if( *h < 0 )
		*h += 360;
}

/************************************************************************
 * INPUT VALUES ARE
 * H=[0,360]
 * S=[0,1]
 * V=[0,1]
 * OUTPUT VALUES ARE [0,255]
 ************************************************************************/

void hsv2rgbscaled (double h, double s,	double v, byte *br,	byte *bg, byte *bb)	
{
	double r,g,b;

	hsv2rgb	(h,s,v,&r,&g,&b);

	r *= MAX_RED;
	g *= MAX_GREEN;
	b *= MAX_BLUE;

	if (r<0) r=0;
	if (g<0) g=0;
	if (b<0) b=0;
	if (r>255) *br=255;	else *br = (int) r/1;
	if (g>255) *bg=255;	else *bg = (int) g/1;
	if (b>255) *bb=255;	else *bb = (int) b/1;
}

/************************************************************************
 * INPUT VALUES ARE
 * H=[0,360]
 * S=[0,1]
 * V=[0,1]
 * OUTPUT VALUES ARE [0,1]
 ************************************************************************/

void hsv2rgb(double h, double s, double v,double *r, double *g, double *b)
{
	int i;
	double f, p, q, t;

	if( s == 0 ) {
		// achromatic (grey)
		*r = *g = *b = v;
		return;
	}

	h /= 60;			// sector 0 to 5
	i = floor( h );
	f = h - i;			// factorial part of h
	p = v * ( 1 - s );
	q = v * ( 1 - s * f );
	t = v * ( 1 - s * ( 1 - f ) );

	switch( i ) {
		case 0:
			*r = v;
			*g = t;
			*b = p;
			break;
		case 1:
			*r = q;
			*g = v;
			*b = p;
			break;
		case 2:
			*r = p;
			*g = v;
			*b = t;
			break;
		case 3:
			*r = p;
			*g = q;
			*b = v;
			break;
		case 4:
			*r = t;
			*g = p;
			*b = v;
			break;
		default:		// case 5:
			*r = v;
			*g = p;
			*b = q;
			break;
	}
}

/************************************************************************
 * Gamma correction	for	a value	in the range 0-1
 ************************************************************************/

double gammaCorrection4Float (double val, double gamma)	
{
	val	= pow(val,1.0/gamma);
	if (val<0) val=0;
	if (val>1) val=1;
	return val;
}

/************************************************************************
 * Gamma correction	for	a grayvalue
 ************************************************************************/

byte gammaCorrection4GrayValue (byte gv, double	gamma) 
{
	double val;
	int	n;

	val	=  (double)	gv / 255.0;
	val	= 255.0	*  pow(val,1.0/gamma);
	n =	(int) val /	1;
	if (n<0) n=0;
	if (n>255) n=255;

	return n;
}

/************************************************************************
 * Scales a	value from one range to	another	one
 ************************************************************************/

double scale_value (double val,	double omin, double	omax, double nmin, double nmax)	
{
	double color=val;
	
	color -= omin;
	color =	(color*(nmax-nmin))/(omax-omin);
	color += nmin;

	if (color>nmax)
		color =	nmax;
	else
		if (color <	nmin)
			color =	nmin;

	return color;
}

