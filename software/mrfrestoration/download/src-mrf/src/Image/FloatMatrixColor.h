/**********************************************************
 * FloatMatrixColor.h
 *
 * Christian Wolf, e9226297@stud1.tuwien.ac.at
 * Beginn: 12.6.2000
 **********************************************************/

#ifndef	_WOLF_FloatMatrixCOLOR_H_
#define	_WOLF_FloatMatrixCOLOR_H_

#include "FloatMatrix.h"

class Image;

class FloatMatrixColor :	public FloatMatrix {

	public:

		// Constructors
		FloatMatrixColor	();
		FloatMatrixColor	(int sx, int sy) { _alloc (sx, sy);	}
		FloatMatrixColor (const FloatMatrixColor &other);
		FloatMatrixColor	(Image &other);

		// destructor
		virtual	~FloatMatrixColor ();

		// Access functions
		inline float get (int p, int x,	int	y) const;
		inline void	set	(int p,	int	x, int y, float	v);
		inline void	inc	(int p,	int	x, int y);

		// Access and change the whole data plane
		inline float ** getPlane(int plane);
		inline void changePlane (int plane, float **data);

		// IO
		virtual	void read (char	*filename, int sx, int sy);
		virtual	void write (char *filename);

		// Conversions
		void convertRGB2LUV();
		void convertRGB2LAB();
		void convertLUV2RGB();
		void convertRGB2GCM(); // Conversion function for the Gaussian color model
							  // After JM Geusebroek et al, "Color Invariants",
							  // PAMI 23(12)2001 pp1338-1350

		// Math & data operations
		virtual	void setZero ();
		void combine (int plane, double	w1,	double w2, double w3);
		void setZero (int plane);
		void add(Image &i);
		void sub(Image &i);
		void add(FloatMatrixColor &other, int xoffset, int yoffset);
		void scalePlane	(int plane,	double omin, double	omax,
			double nmin, double	nmax);

		void debugPrint();

	protected:	// DATA

		float **data2, **data3;

	private:

		void _alloc	(int sx, int sy);

		inline double sobelOnPixel (int	plane, int i, int j) {
			return
				 get(plane,i+1,j+1)	-	get(plane,i-1,j+1)
			+	get(plane,i+1,j-1)	-	get(plane,i-1,j-1)
			+2*( get(plane,i+1,j)	-	get(plane,i-1,j))
		  ;
		}

};

// **************************************************************
// Access methods
// **************************************************************

inline void	FloatMatrixColor::set (int p, int x, int y, float v) {
    switch (p) {
    	case 1:
    		data1[x][y] = v;
    		break;
    	case 2:
    		data2[x][y] = v;
    		break;
    	case 3:
    		data3[x][y] = v;
    		break;
    	default:
    		break;
    }
}

inline float FloatMatrixColor::get (int p, int x, int y) const {
	switch (p) {
		case 1:
			return data1[x][y];
		case 2:
			return data2[x][y];
		case 3:
			return data3[x][y];
		default:
			return 0;
	}
}

inline void	FloatMatrixColor::inc (int p, int x, int y) {
    switch (p) {
    	case 1:
    		++data1[x][y];
    		break;
    	case 2:
    		++data2[x][y];
    		break;
    	case 3:
    		++data3[x][y];
    		break;
    	default:
    		break;
    }
}

// **************************************************************
// Access and change the whole data plane
// **************************************************************

inline float ** FloatMatrixColor::getPlane(int plane) {
	switch (plane) {
		case 1:
			return data1;
		case 2:
			return data2;
		case 3:
			return data3;
#ifdef CHECK_CODE			
		default:
			cerr << "Internal error in FloatMatrixColor::getPlane()\n";
			CORE_DUMP;
#endif
	}		
}

inline void FloatMatrixColor::changePlane (int plane, float **d) {
	float **p;
	switch (plane) {
		case 1:
			p=data1;
			data1=d;
			break;
		case 2:
			p=data2;
			data2=d;
			break;
		case 3:
			p=data3;
			data3=d;
			break;
#ifdef CHECK_CODE
		default:
			cerr << "Internal error in FloatMatrixColor::getPlane()\n";
			CORE_DUMP;
#endif
	}
	if (p!=NULL) {
		float *p2=*p;
		delete [] p2;
		delete [] p;
	}
}
	


#endif


