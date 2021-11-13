#ifndef	_WOLF_CCOMPONENT_LIST_H_
#define	_WOLF_CCOMPONENT_LIST_H_

// C++ STL Headers
#include <algorithm>
#include <list>

// From the MATHEMATICS module
#include "Matrix.h"

// From this module
#include "CComponent.h"

class RectText;
class RectTextList;
class TextBoxList;

class CComponentList {

	public:

		// Destructor
		~CComponentList	() 
		{
			 for (iterator iter=begin(); iter!=end(); ++iter)
				delete (*iter);
		}

		// The iterator
		typedef	list<CComponent	*>::iterator iterator;
		iterator begin() { return clist.begin(); }
		iterator end() { return	clist.end(); }

		void add (CComponent * comp) { clist.push_back (comp); }
		int	size () { return clist.size();	}
		
		// Get a Matrix of the size of the original image (the size has to be
		// specified), where each element (pixel) is either NULL if it does not
		// belong to any cluster or a pointer to the respective component
		Matrix<CComponent *> * getPointerMatrix (int xsize, int ysize);

		// To debug...
		void print (FILE *fp);

		// DATA

		list<CComponent	*> clist;
};

#endif

