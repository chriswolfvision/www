#ifndef	_WOLF_CCOMPONENT_H_
#define	_WOLF_CCOMPONENT_H_

// C
#include <stdlib.h>
#include <stdio.h>

// From the IMAGE library
#include <Rect.h>

// From this library
#include "CPixel.h"

class Image;

class CComponent 
{

   public:

		// Constructors
		CComponent () 
		{	
			count=0; collected=false; pixels=NULL;
			boundingBoxDirty =	true;
		}

		// Destructor
		~CComponent	();

		// Basic operations
		void insertPixel (int x, int y);
		void merge (CComponent *other);

		// Geometry
		Rect boundingBox ();
		double fillRate() {	return ((double) count / ((double) boundingBox().area())); };
				
		// To debug	...
		void print (FILE *fp);

		// data
		int	count;
		bool collected;
		CPixel *pixels;

	private:

		bool boundingBoxDirty;
		Rect box;
};

#endif

