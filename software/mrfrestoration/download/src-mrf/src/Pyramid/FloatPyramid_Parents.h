/***************************************************************************
                          FloatPyramid_Parents.h  -  description
                             -------------------
    begin                : Thu May 30 2002
    copyright            : (C) 2002 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 *
 ***************************************************************************
 * These functions calculate the parent pixel from a given child pixel.
 * Depending on the filter kernel, the reduction function and the position
 * of the pixel on the grid, a different number of parent pixels have to
 * to be combined with different weigths.
 *
 * For all the following accessParents*() functions, see
 * research notebook 29.5.2002, pp 122-123
 *
 * There is not sophisticated border treatment: If one of the indices
 * is out of range, then the center pixel is taken (i,j)
 *
 * Parameters: The Image of the PARENT level, and the coordinates of the
 * 				pixel on the CHILD level.
 ****************************************************************************/

// *******************************************************************************
// The macros responsible for the border treatment
// *******************************************************************************

// for 3x3, 5x3 and 3x5 gaussians
#define CHECK_I  	if(i<=0) return X->get(i,j); else
#define CHECK_J  	if(j<=0) return X->get(i,j); else
#define CHECK_IJ  	if(i<=0||j<=0) return X->get(i,j); else
#define CHECK_IB  	if(i<=0||i>=X->xsize-1) return X->get(i,j); else
#define CHECK_JB  	if(j<=0||j>=X->ysize-1) return X->get(i,j); else
#define CHECK_IB_J 	if(i<=0||i>=X->xsize-1||j<=0) return X->get(i,j); else
#define CHECK_I_JB 	if(i<=0||j<=0||j>=X->ysize-1) return X->get(i,j); else
#define CHECK_IB_JB if(i<=0||i>=X->xsize-1||j<=0||j>=X->ysize-1) return X->get(i,j); else

#define XS			X->xsize
#define YS			X->ysize
#define BORDER_ELSE	return X->get(i,j); else

// *******************************************************************************
// A 3x3 gaussian filter kernel with 2x2 reduction function
// *******************************************************************************

inline float accessParents33red22 (FloatMatrix *X, int x, int y) {
	int i=x/2;
	int j=y/2;
	if (x%2==0) {
		if (y%2==0) {
			CHECK_IJ return (X->get(i-1,j-1) + X->get(i,j-1) + X->get(i-1,j) + X->get(i,j)) / 4.0;
		}
		else  {
			CHECK_IJ return (X->get(i-1,j-1) + X->get(i,j-1)) / 2.0;					
		}
	}
	else {
		if (y%2==0) {
			CHECK_J return (X->get(i,j-1) + X->get(i,j)) / 2.0;			
		}
		else		
			return X->get(i,j);
	}
}

// *******************************************************************************
// A 5x3 gaussian filter kernel with 2x2 reduction function
// *******************************************************************************

inline float accessParents53red22 (FloatMatrix *X, int x, int y) {
	int i=x/2;
	int j=y/2;
	if (x%2==0) {
		if (y%2==0) {
			CHECK_IJ return (X->get(i-1,j-1) + X->get(i,j-1) + X->get(i-1,j) + X->get(i,j)) / 4.0;
		}
		else {
			CHECK_I return (X->get(i-1,j) + X->get(i,j)) / 2.0;			
		}
	}
	else {
		if (y%2==0) {
			CHECK_IB_J
			return (X->get(i-1,j-1) + 6.0*X->get(i,j-1) + X->get(i+1,j-1) +
					X->get(i-1,j  ) + 6.0*X->get(i,j  ) + X->get(i+1,j)) / 16.0;
		}					
		else {
			CHECK_IB return (X->get(i-1,j) + 6.0*X->get(i,j)+ X->get(i+1,j)) / 8.0;
		}
	}
}

// *******************************************************************************
// A 3x5 gaussian filter kernel with 2x2 reduction function
// *******************************************************************************

inline float accessParents35red22 (FloatMatrix *X, int x, int y) {
	int i=x/2;
	int j=y/2;
	if (x%2==0) {
		if (y%2==0) {
			CHECK_IJ return (X->get(i-1,j-1) + X->get(i,j-1) + X->get(i-1,j) + X->get(i,j)) / 4.0;
		}
		else {
			CHECK_I_JB return (X->get(i-1,j-1) + X->get(i,j-1) + X->get(i-1,j) +
								X->get(i,j) + X->get(i-1,j+1) + X->get(i,j+1)) /16.0;
		}
	}
	else {
		if (y%2==0) {
			CHECK_J return (X->get(i,j-1) + X->get(i,j)) / 2.0;
		}					
		else {
			CHECK_JB return (X->get(i,j-1) + 6.0*X->get(i,j)+ X->get(i,j+1)) / 8.0;
		}
	}
}

// *******************************************************************************
// A 7x3 gaussian filter kernel with 2x2 reduction function
// *******************************************************************************

inline float accessParents73red22 (FloatMatrix *X, int x, int y) {
	int i=x/2;
	int j=y/2;
	if (x%2==0) {
		if (y%2==0) {
			if (i<2||i>=XS-1||j<=0) BORDER_ELSE
				return (X->get(i-2,j-1)+15.0*X->get(i-1,j-1)+15.0*X->get(i,j-1)+
					X->get(i+1,j-1)+X->get(i-2,j)+15.0*X->get(i-1,j)+
					15.0*X->get(i,j)+X->get(i+1,j))/64.0;
		}
		else {
			if (i<2||i>=XS-1) BORDER_ELSE
			return (X->get(i-2,j) + 15.0*X->get(i-1,j) + 15.0*X->get(i,j) +
					X->get(i+1,j)) / 32.0;
			
		}
	}
	else {
		if (y%2==0) {
			CHECK_IB_J return (3.0*X->get(i-1,j-1)+10.0*X->get(i,j-1)+3.0*X->get(i+1,j-1)+
					3.0*X->get(i-1,j)+10.0*X->get(i,j)+3.0*X->get(i+1,j)) / 32.0;			
		}					
		else {
			CHECK_IB return (3.0*X->get(i-1,j)+10.0*X->get(i,j)+3.0*X->get(i+1,j))/16.0;
		}
	}
}

// *******************************************************************************
// A 3x7 gaussian filter kernel with 2x2 reduction function
// *******************************************************************************

inline float accessParents37red22 (FloatMatrix *X, int x, int y) {
	int i=x/2;
	int j=y/2;
	if (x%2==0) {
		if (y%2==0) {
			if (i<=0||j<2||j>=YS-1) BORDER_ELSE
				return (X->get(i-1,j-2) + X->get(i,j-2)+
					15.0*X->get(i-1,j-1) + 15.0*X->get(i,j-1)+
					15.0*X->get(i-1,j) + 15.0*X->get(i,j)+
					X->get(i-1,j+1) + X->get(i,j+1))/64.0;
		}
		else {
			CHECK_I_JB	return (3.0*X->get(i-1,j-1) + 3.0*X->get(i,j-1) +
					10.0*X->get(i-1,j) + 10.0*X->get(i,j) +
					3.0*X->get(i-1,j+1) + 3.0*X->get(i,j+1)
					) / 32.0;			
		}
	}
	else {
		if (y%2==0) {
			if (j<2||j>=YS-1) BORDER_ELSE
				return (X->get(i,j-2) + 15.0*X->get(i,j-1)+
					15.0*X->get(i,j) + X->get(i,j+1)) / 32.0;			
		}					
		else {
			CHECK_JB return (3.0*X->get(i,j-1)+10.0*X->get(i,j)+3.0*X->get(i,j+1))/16.0;
		}
	}
}

// *******************************************************************************
// A 7x7 R-diagonal gaussian filter kernel with 2x2 reduction function
// *******************************************************************************

inline float accessParents77RDred22 (FloatMatrix *X, int x, int y) {
	int i=x/2;
	int j=y/2;
	if (x%2==0) {
		if (y%2==0) {
			if (i<2||i>=XS-1||j<1||j>=YS-1) BORDER_ELSE
				return (1.0*X->get(i-2,j+0)+ 8.0*X->get(i-1,j-1)+ 22.0*X->get(i-1,j+0)+
					1.0*X->get(i-1,j+1)+1.0*X->get(i+0,j-2)+22.0*X->get(i+0,j-1)+
					8.0*X->get(i+0,j+0)+1.0*X->get(i+1,j-1))/64.0;
		}
		else {
			if (i<2||i>=XS-1||j<1||j>=YS-1) BORDER_ELSE
				return(1.0*X->get(i-2,j+1)+1.0*X->get(i-1,j-1)+23.0*X->get(i-1,j+0)+
					7.0*X->get(i-1,j+1)+7.0*X->get(i+0,j-1)+23.0*X->get(i+0,j+0)+
					1.0*X->get(i+0,j+1)+1.0*X->get(i+1,j-1))/64.0;					
		}
	}
	else {
		if (y%2==0) {
			if (i<1||i>=XS-1||j<2||j>=YS-1) BORDER_ELSE
				return (1.0*X->get(i-1,j-1)+7.0*X->get(i-1,j+0)+1.0*X->get(i-1,j+1)+
					23.0*X->get(i+0,j-1)+23.0*X->get(i+0,j+0)+1.0*X->get(i+1,j-2)+
					7.0*X->get(i+1,j-1)+1.0*X->get(i+1,j+0))/64.0 ;
		}					
		else {
			CHECK_IB_JB return (4.0*X->get(i-1,j+0)+3.0*X->get(i-1,j+1)+4.0*X->get(i+0,j-1)+
				42.0*X->get(i+0,j+0)+4.0*X->get(i+0,j+1)+3.0*X->get(i+1,j-1)+
				4.0*X->get(i+1,j+0))/64.0;				
		}
	}
}

// *******************************************************************************
// A 7x7 L-diagonal gaussian filter kernel with 2x2 reduction function
// *******************************************************************************

inline float accessParents77LDred22 (FloatMatrix *X, int x, int y) {
	int i=x/2;
	int j=y/2;
	if (x%2==0) {
		if (y%2==0) {
			if (i<2||i>=XS-1||j<2||j>=YS-1) BORDER_ELSE
				return (1.0*X->get(i-2,j-1)+1.0*X->get(i-1,j-2)+22.0*X->get(i-1,j-1)+
					8.0*X->get(i-1,j+0)+8.0*X->get(i+0,j-1)+22.0*X->get(i+0,j+0)+
					1.0*X->get(i+0,j+1)+1.0*X->get(i+1,j+0))/64.0;
		}
		else {
			if (i<2||i>=XS-1||j<1||j>=YS-1) BORDER_ELSE
				return (1.0*X->get(i-2,j-1)+7.0*X->get(i-1,j-1)+23.0*X->get(i-1,j+0)+
					1.0*X->get(i-1,j+1)+1.0*X->get(i+0,j-1)+23.0*X->get(i+0,j+0)+
					7.0*X->get(i+0,j+1)+1.0*X->get(i+1,j+1))/64.0;
		}
	}
	else {
		if (y%2==0) {
			if (i<1||i>=XS-1||j<2||j>=YS-1) BORDER_ELSE
				return(1.0*X->get(i-1,j-2)+7.0*X->get(i-1,j-1)+1.0*X->get(i-1,j+0)+
					23.0*X->get(i+0,j-1)+23.0*X->get(i+0,j+0)+1.0*X->get(i+1,j-1)+
					7.0*X->get(i+1,j+0)+1.0*X->get(i+1,j+1))/64.0;
		}					
		else {
			CHECK_IB_JB return 	(3.0*X->get(i-1,j-1)+4.0*X->get(i-1,j+0)+4.0*X->get(i+0,j-1)+
				42.0*X->get(i+0,j+0)+4.0*X->get(i+0,j+1)+4.0*X->get(i+1,j+0)+
				3.0*X->get(i+1,j+1))/64.0;
		}
	}
}


