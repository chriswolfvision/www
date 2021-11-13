#include <stdlib.h>
#include <stdio.h>

#include "CComponent.h"

/************************************************************************
 * Destructor
 ************************************************************************/

CComponent::~CComponent	() 
{
	CPixel *pp=pixels, *save_next;
	while (pp!=NULL) {
		save_next =	pp->next;
		delete pp;
		pp = save_next;
	}
}

/************************************************************************
 * Inserts a pixel
 ************************************************************************/

void CComponent::insertPixel (int x, int y)	
{
	CPixel *myPixel	= new CPixel (x,y);

	myPixel->next =	this->pixels;
	this->pixels = myPixel;
	++count;

	boundingBoxDirty = true;
}

/************************************************************************
 * Merges two regions, i.e.	appends	the	pixel list of the second region
 * to the first	region and frees the structure of the second region.
 ************************************************************************/

void CComponent::merge (CComponent *other) 
{
	CPixel *cur, *prec;

	cur	= this->pixels;
	prec = NULL;
	while (cur!=NULL) 
	{
		prec = cur;
		cur	= cur->next;
	}

	/* first region	is empty */
	if (prec==NULL)
		this->pixels = other->pixels;
	else
		prec->next = other->pixels;

	this->count	+= other->count;

	// delete other;
	other->pixels =	NULL;
	delete other;
}

/************************************************************************
 * Calculates the bounding box of a	connected component
 ************************************************************************/

Rect CComponent::boundingBox ()	
{

	int	lt = 999999, ll	= 999999;
	int	lb = 0,	lr = 0;
	CPixel *p;

	if (!boundingBoxDirty) {
		return box;
	}

	p =	pixels;
	while (p!=NULL)	
	{

		if (p->y < lt) lt =	p->y;
		if (p->x < ll) ll =	p->x;
		if (p->y > lb) lb =	p->y;
		if (p->x > lr) lr =	p->x;

		p =	p->next;
	}

	box.top	= lt;
	box.bottom = lb;
	box.left = ll;
	box.right =	lr;

	boundingBoxDirty = false;

	return box;
}

/************************************************************************
 * Calculates the bounding box of a	connected component
 ************************************************************************/

void CComponent::print (FILE *fp) 
{
	Rect b;

	b =	boundingBox	();
	fprintf	(fp, "[Component: %d pixels	box: %d,%d-%d,%d area: %d shape: %6.2f "
		"fr: %6.2f]\n",
		count, b.left, b.top, b.right, b.bottom,
		b.area(), b.shape(), fillRate());
}

