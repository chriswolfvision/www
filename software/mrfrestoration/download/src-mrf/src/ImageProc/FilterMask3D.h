#ifndef	_WOLF_FILTERMASK3D_H_
#define	_WOLF_FILTERMASK3D_H_

#include "FilterMask.h"

class FilterMask3D:	public FilterMask {

	public:

		FilterMask3D(int x,	int	y, int z);

		void normalize();

		int	getIndex (int x, int y,	int	z) {
			return z*ysize*xsize+y*xsize+x;
		}

		void clear() { for (int	i=0; i<xsize*ysize*zsize; ++i) M[i]=0; }
		double get(int x, int y, int z)	{ return M[getIndex(x,y,z)]; }
		void   set(int x, int y, int z,	double v) {	M[getIndex(x,y,z)]=v; }

		int	xsize, ysize, zsize;

	private:

		double *M;

};

#endif

