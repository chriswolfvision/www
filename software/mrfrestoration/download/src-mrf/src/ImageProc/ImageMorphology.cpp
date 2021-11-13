// C
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

// From the MATHEMATICS library
#include <Matrix.h>

// From this library
#include "ImageProc.h"

// *************************************************************
// Erosion
// gvForm is the grayvalue of the object
// The grayvalue of the background is assumed to be 0
// *************************************************************

void erode (Image &im, int Nit, bool horizOnly, unsigned char gvObject) 
{
	int	i, j, k, x,	y, n;
	unsigned char		**I, **J ;

	I =	im.R;
	J =	(unsigned char **) CREATE_IMAGE(im.ysize,im.xsize) ;

	/*********************************************
	 * Iterations
	 *********************************************/

	for	(n = 0 ; n < Nit ; n++)	
	{

		for	(i = 0 ; i < im.xsize ; i++)
		for	(j = 0 ; j < im.ysize ; j++)

			if (I[i][j]>0)	
			{
				k =	1 ;
				y =	j;
				x =	i;

				if (horizOnly) 
				{
					for	(x = i-1 ; x <=	i+1	; x++)
						if ((x >= 0) &&	(x < im.xsize))
							k *= I[x][y] ;
				}
				else  
				{
					for	(x = i-1 ; x <=	i+1	; x++)
						if ((x >= 0) &&	(x < im.xsize))
							for	(y = j-1 ; y <=	j+1	; y++)
								if ((y >= 0) &&	(y < im.ysize))
									k *= I[x][y] ;
				}

				J[i][j]	= (k !=	0 ?	gvObject : 0) ;

			}
			else
				J[i][j]	= 0	;

		for	(i = 0 ; i < im.xsize ; i++)
		for	(j = 0 ; j < im.ysize ; j++)
			I[i][j]	= J[i][j] ;
	}
	FREE_IMAGE(J);
}

// *************************************************************
// Dilation
// gvForm is the grayvalue of the object
// The grayvalue of the background is assumed to be 0
// *************************************************************

void dilate (Image &im, int Nit, bool horizOnly, unsigned char gvObject) 
{
	 int i,	j, k, x, y,	n;
	 unsigned char		**I, **J;

	 J = (unsigned char	**)	CREATE_IMAGE(im.ysize,im.xsize) ;
	 I = im.R;

	/*********************************************
	 * Iterations
	 *********************************************/

	for	(n = 0 ; n < Nit ; n++)	
	{
		/* Travers the whole image */
		for	(i = 0 ; i < im.xsize ; i++)
		for	(j = 0 ; j < im.ysize ; j++)	
		{
			if (I[i][j]	== 0) 
			{
				k =	0;
				y =	j;
				x =	i;

				if (horizOnly) 
				{
					for	(x = i-1 ; x <=	i+1	; x++)
						if ((x >= 0) &&	(x < im.xsize))
						   k +=	I[x][y]	;
				}
				else 
				{
					for	(x = i-1 ; x <=	i+1	; x++)
						if ((x >= 0) &&	(x < im.xsize))
							for	(y = j-1 ; y <=	j+1	; y++)
								if ((y >= 0) &&	(y < im.ysize))
									k += I[x][y] ;
				}

				J[i][j]	= (k ? gvObject : 0) ;
			 }
			 else
				J[i][j]	= gvObject;
		}

		for	(i = 0 ; i < im.xsize ; i++)
		for	(j = 0 ; j < im.ysize ; j++)
			I[i][j]	= J[i][j] ;
	}

	FREE_IMAGE (J);
}

void dilate	(Image &i, int iterations) 
{
	dilate (i, iterations, 0, 255);	
}

void dilateHoriz (Image &i, int iterations) 
{	
	dilate	(i, iterations, 1, 255);
}

void erode (Image &i, int	iterations)	
{
	erode (i, iterations, 0, 255);
}

void erodeHoriz	(Image &i, int iterations) 
{
	erode (i, iterations, 1, 255);
}

// *************************************************************
// Erosion - on gray level images
// *************************************************************

void erodeGray (Image &im, int Nit, bool horizOnly) {
	int	i, j, x, y, n, min;
	unsigned char **I, **J ;

	I =	im.R;
	J =	(unsigned char **) CREATE_IMAGE(im.ysize,im.xsize) ;

	/*********************************************
	 * Iterations
	 *********************************************/

	for	(n = 0 ; n < Nit ; n++)	{

		for	(i = 0 ; i < im.xsize ; i++)
		for	(j = 0 ; j < im.ysize ; j++) {

			y =	j;
			x =	i;

			if (horizOnly) {
				min = 255;
				for	(x = i-1 ; x <=	i+1	; x++)
					if ((x >= 0) &&	(x < im.xsize))
						if (I[x][y] < min)
							min = I[x][y] ;
			}
			else  {
				min = 255;
				for	(x = i-1 ; x <=	i+1	; x++)
					if ((x >= 0) &&	(x < im.xsize))
						for	(y = j-1 ; y <=	j+1	; y++)
							if ((y >= 0) &&	(y < im.ysize))
								if (I[x][y] < min)
									min = I[x][y] ;
			}

			J[i][j]	= min;
		}

		for	(i = 0 ; i < im.xsize ; i++)
		for	(j = 0 ; j < im.ysize ; j++)
			I[i][j]	= J[i][j] ;
	}
	FREE_IMAGE(J);
}

// *************************************************************
// Calculates the distance from	the	border
// *************************************************************

void levelCurve4 (Image &im) {

	int	i,j,k,min,pix[3];
	unsigned char **B;

	B =	im.R;

	for	(i=0; i<im.xsize; i++)
		 for (j=0; j<im.ysize; j++)
			 if	(B[i][j]!=0)
				B[i][j]=1;
			 else
				B[i][j]=0;


	for	(i=0; i<im.xsize; i++)
		for	(j=0; j<im.ysize; j++)
			{
			 if	(i>0) pix[0]=B[i-1][j] +1;
					  else pix[0]=1;
			 if	(j>0) pix[1]=B[i][j-1] +1;
					  else pix[1]=1;
			 if( B[i][j] !=	0)
				 {
				 min=pix[0];
				 if	(pix[1]< min) min=pix[1];
				 B[i][j]=min;
				 }
			}

	for	(i=im.xsize-1; i>=0;	i--)
		for	(j=im.ysize-1; j>=0;	j--)
			 {
			 if	(i+1<im.xsize) pix[0]=B[i+1][j] +1;
						 else  pix[0]=1;
			 if	(j+1<im.ysize) pix[1]=B[i][j+1] +1;
						else pix[1]=1;
			 pix[2]=B[i][j];
			 if(B[i][j]	!= 0)
						{
				 min=pix[0];
				 for(k=1;k<=2;k++) if (pix[k]<min) min=pix[k];
				 B[i][j]=min;
						}
			}
}

// *************************************************************
// Calculates the distance from	the	border -
// Only	horizontal
// *************************************************************

void levelCurve4Horiz (Image &im) {

	int	i,j,min,pix[3];
	unsigned char **B;

	B =	im.R;

	for	(i=0; i<im.xsize; i++)
		 for (j=0; j<im.ysize; j++)
			 if	(B[i][j]!=0)
				B[i][j]=1;
			 else
				B[i][j]=0;


	for	(i=0; i<im.xsize; i++)
		for	(j=0; j<im.ysize; j++) {
			 if	(i>0) min=B[i-1][j]	+ 1;
					  else min=1;
			 if( B[i][j] !=	0)
				 B[i][j]=min;
			}

	for	(i=im.xsize-1; i>=0;	i--)
		for	(j=im.ysize-1; j>=0;	j--) {
			 if	(i+1<im.xsize) pix[0]=B[i+1][j] +1;
						 else  pix[0]=1;
			 pix[1]=B[i][j];
			 if(B[i][j]	!= 0){
				 min=pix[0];
				 if	(pix[1]< min) min=pix[1];
				 B[i][j]=min;
			}
		 }
}

// *************************************************************
// Calculates and thresholds the ferret	diameter
// *************************************************************

void ferret (Image &im, unsigned int minX, unsigned int maxX,
	unsigned int minY, unsigned	int	maxY) {

	int	i,j;
	unsigned int **buff;

	buff=(unsigned int **)malloc(im.xsize*sizeof(void *));
	for	(j=0; j<im.xsize; j++) {
		buff[j]=(unsigned int *)malloc(im.ysize*sizeof(int));
		memset(buff[j],0,im.ysize*sizeof(int));
	}

	for	(i=1; i<im.ysize; i++)
		for	(j=1; j<im.xsize; j++)
			  buff[j][i]=im.R[j][i];

	/* diametre	de ferret vertical */
	for	(i=1; i<im.ysize; i++) {
		for	(j=1; j<im.xsize; j++) {
			 if	(im.R[j][i]>0)
				 im.R[j][i]=im.R[j][i-1]+1;
		}
	}

	for	(i=im.ysize-2; i>0; i--)	{
		for	(j=1; j<im.xsize; j++) {
			 if	(im.R[j][i]>0)
				if (im.R[j][i]<im.R[j][i+1]) im.R[j][i]=im.R[j][i+1];
		}
	}

	/* diametre	de ferret horizontal */
	for	(j=1; j<im.xsize; j++) {
		for	(i=1; i<im.ysize; i++) {
			 if	(buff[j][i]>0)
				 buff[j][i]=buff[j-1][i]+1;
		}
	 }

	for	(j=im.xsize-2; j>0; j--)	{
		for	(i=1; i<im.ysize; i++) {
			 if	(buff[j][i]>0)
				if (buff[j][i]<buff[j+1][i]) buff[j][i]=buff[j+1][i];
		}
	}

	for	(i=1; i<im.ysize; i++) {
		for	(j=1; j<im.xsize; j++)
		  if (im.R[j][i]>0) {
			if ( ((unsigned	int) (im.R[j][i]*3)>maxY)&&
				 ((unsigned	int) (buff[j][i]*8)<maxX))
				im.R[j][i]=0;
			else
				if ( ((unsigned	int)(im.R[j][i]*3)<minY)||
					 ((unsigned	int)(buff[j][i]*8)<minX))
					im.R[j][i]=0;
		}
	}

	for	(j=0; j<im.xsize; j++)
		free(buff[j]);
	free(buff);
}

