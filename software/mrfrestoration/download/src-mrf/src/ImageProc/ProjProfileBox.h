#ifndef	_WOLF_PROJPROFILEBOX_H_
#define	_WOLF_PROJPROFILEBOX_H_

#include "ProjProfile.h"

class ProjProfileBox {

	public:

		// Constructors
		ProjProfileBox () {	horiz =	vert = NULL; }

		// Destructor
		~ProjProfileBox	() {
			if (horiz!=NULL) {
				delete horiz;
				delete vert;
			}
		}

		// The bounding	box
		int	getTop ()	{ return (int) vert->getMin(); }
		int	getBottom () { return (int) vert->getMax(); }
		int	getLeft	()	{ return (int) horiz->getMin(); }
		int	getRight ()	{ return (int) horiz->getMax(); }

		bool check4Text(Image *im, int vbeg, int vend);

	private:

		ProjProfile	*horiz,	*vert;

};

#endif

