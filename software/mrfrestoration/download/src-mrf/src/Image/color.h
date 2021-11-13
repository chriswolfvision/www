/***************************************************************************
                          color.h  -  description
                             -------------------
    begin                : Sometime in 1995
    email                : wolf@rfv.insa-lyon.fr
 ***************************************************************************/

#ifndef	_WOLF_COLOR_H_
#define	_WOLF_COLOR_H_

#define	MAX_RED			255
#define	MAX_GREEN		255
#define	MAX_BLUE		255

// From	this threshold of the L	part of	a LUV color, the color
// is considered as	dark. We suppose that the range	has
// been	scaled to 0-255!!!
#define	THR_L_SCALED_DARK		120

// the same, unscaled, Range MINVALUE_L	- MAXVALUE_L
#define	THR_L_UNSCALED_DARK		50

typedef	unsigned char byte;

#ifdef __cplusplus
extern "C" {
#endif

/**************************************************************************	
 *	A couple of definitions
 **************************************************************************/
 
enum ColorSpaceCode
{
	COLSPC_RGB=0,
	COLSPC_HSV,
	COLSPC_LUV,
	COLSPC_LAB,
	COLSPC_LAB_CIE94			// Different distance function
};

void init_color	();

// Conversion functions for the L*u*v* color model
void rgb2luv (double r,	double g, double b,	double *l, double *u, double *v);
void luv2rgb (double l,	double u, double v,	double *r, double *g, double *b);
void luv2rgbscaled (double l, double u,	double v, byte *br,	byte *bg, byte *bb);

// Conversion functions for the L*a*b* color model
void lab2rgb (double l,	double a, double bstar, double *r, double *g, double *blue); 
void rgb2lab (double r,	double g, double blue, double *l, double *a, double *bstar);

// Conversion functions for the HSV color model
void rgb2hsv (double r,	double g, double b,	double *h, double *s, double *v);
void hsv2rgb (double h,	double s, double v,	double *r, double *g, double *b);
void hsv2rgbscaled (double h, double s,	double v, byte *br,	byte *bg, byte *bb);

// The CIE94 color difference
double diffLabCIE94 (double l1, double a1, double b1, double l2, double a2, double b2);
double diffColorEuclidean (double l1, double u1, double v1, double l2, double u2, double v2);

byte rgb2grayscale (byte r, byte g, byte b);



double gammaCorrection4Float (double val, double gamma);
byte gammaCorrection4GrayValue (byte gv, double	gamma);
double scale_value (double val,	double omin, double	omax, double nmin, double nmax);


#ifdef __cplusplus
}
#endif

#define	MINVALUE_L	000.0
#define	MAXVALUE_L	100.0
#define	MINVALUE_U	-084.0
#define	MAXVALUE_U	176.0
#define	MINVALUE_V	-131.0
#define	MAXVALUE_V	109.0
#define	MINVALUE_A	-120
#define	MAXVALUE_A	+120
#define	MINVALUE_B	-120
#define	MAXVALUE_B	+120

#endif

