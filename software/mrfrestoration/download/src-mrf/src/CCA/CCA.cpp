/**********************************************************
 * Image_Segments.cc
 * The connected components	algorithm
 * Originally written in C, then ported to C++
 *
 * Christian Wolf
 * chriswolf@gmx.at
 * Begin: 1997 at PRIP, TU-Wien
 * Major rewrite: 9.6.2005
 **********************************************************/

// C++
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

// From the IMAGE module
#include <Image.h>

// From this module
#include "CCA.h"

/************************************************************************
 * A locally only used class
 ************************************************************************/

class CComponentContainer
{
	public:
		CComponentContainer (CComponent *c) {  component=c; pixels=NULL; }
	
	public: // data
		CComponent *component;
		CComponentContainer *next;
		CPixel *pixels;
};

/************************************************************************
 * Travers the global container	list and replace components
 ************************************************************************/

/*
static void	replaceComponents (CComponent *replace, CComponent *to,
	CComponentContainer	*root) 
{
	CComponentContainer	*p = root;
	while (p!=NULL)	
	{
		if (p->component==replace)
			p->component=to;
		p=p->next;
	}
}
*/

/************************************************************************
 * Travers the global container	list and replace components
 ************************************************************************/

static void	merge (CComponentContainer *willBeDeleted, CComponentContainer *willStay,
	Matrix<CComponentContainer *> &m) 
{
	// Replace the component pointers of the pixels in component 1
	// with the pointer to container 2
	CPixel *p=willBeDeleted->component->pixels;
	while (p!=NULL)	
	{
		m(p->y,p->x) = willStay;
		p=p->next;
	}
	
	// Merge the components
	willStay->component->merge(willBeDeleted->component);
	
	// The container is not used anymore
	// We don't delete it because the list traversal would take too
	// much overhead. We just mark the component as unused,
	// the memory will be freed at the end when the data is collected
	willBeDeleted->component=NULL;
}

/************************************************************************
 * Splits up a cluster and performs	all	the	spatial	stuff (connected
 * components analysis)
 ************************************************************************/

#define NBSET(i,x,y)		{ nbx[i]=x; nby[i]=y; } else nbx[i]=0;

CComponentList * connectedComponents (Image &im, CCANeighborhoodType nhType, bool ignoreBG)
{
	int	i,x,y;
	CComponentContainer	*containerRoot, *cur, *last;
	CComponentList *ret = new CComponentList;
	Matrix<CComponentContainer *> cca(im.ysize, im.xsize, NULL);

	// Travers all pixels
	containerRoot = NULL;
	for	(y=0; y<im.ysize; ++y) // { cerr << "[" << y << "]"; cerr.flush();
	for	(x=0; x<im.xsize; ++x) 
	{
		int nbx[4], nby[4];
		int	maxnb;
		unsigned int label = im.get(x,y);
		CComponentContainer **e = cca.getPointerYX(y,x);

		// The current pixel belongs to	the	cluster
		if (ignoreBG && label==0) 
			continue;

		// Which neighborhood do we use?
		if (nhType==CCA_NBHD_4)
			maxnb =	2;
		else
			maxnb =	4;
			
		// Border treatment: 
		// Check for the availability of each neighbor
		if (y>0) 						NBSET(0,x,y-1);
		if (x>0) 						NBSET(1,x-1,y);
		if (nhType==CCA_NBHD_8)
		{
			if (y>0	&& x>0)				NBSET(2,x-1,y-1);
			if (y>0	&& x<(im.xsize-1)) 	NBSET(3,x+1,y-1);
		}

		// Travers all neighbors and merge them together to one
		// component, if their labels are equal.
		cur=NULL;
		for	(i=0; i<maxnb; ++i)	
		{	
			CComponentContainer *n;
							
			// The neighbor is not available
			if (nbx[i]==0)
				continue;			
						
			// The neighbor has not the same grayvalue
			if (label != im.get(nbx[i],nby[i]))
				continue;
				
			n = cca(nby[i],nbx[i]);

			// Match:
			// The first one: use it as	initial	component
			if (cur==NULL)
				cur	= n;

			// Match:
			// Not the first one: append it
			else 
			{
				if (n->component != cur->component) 
				{
					/*
					replaceComponents (n->component, cur->component, containerRoot);
					cur->component->merge(n->component);					
					*/
					merge (n, cur, cca);
				}
			}
		}

		// No component	found, i.e.	no pixel in	the	neigborhood
		// has the same graylavel. Create a	new component
		if (cur==NULL) 
		{
			*e = new CComponentContainer (new CComponent);
			

			// insert into container list
			(*e)->next = containerRoot;
			containerRoot=(*e);
		}

		// One or more components found. After the ev. merging
		// step	it's only one, so we just add the pixel	to it.			 
		else
			(*e)=cur;			

		(*e)->component->insertPixel (x,y);
	}
	// }

	/********************************************************************
	 * Build the CComponentList we will return
	 ********************************************************************/

	cur	= containerRoot;
	last = NULL;
	while (cur != NULL)	
	{
		// Check if the container has a component attached.
		if (cur->component)
		{
			// Check if we already collected this component
			if (!cur->component->collected)	
			{
				cur->component->collected=true;
				ret->add (cur->component);
			}
		}

		last = cur;
		cur	= cur->next;

		// delete the container
		delete last;
	}

	return ret;
}

// *************************************************************
// Component Info
// Stores information about the connected components in
// some optional matrices
// *************************************************************

void componentInfo (Image im, Matrix<int> *height, Matrix<int> *top, Matrix<int> *compid, 
					Matrix<int> *cnt,
					CCANeighborhoodType nhType, bool ignoreBG) 
{
	
	CComponentList *clist;
	CPixel *pp;
	int component_id;

	/*
	{	Image tmp (im);
		im.write ("x_component_all.pgm");
		for (int y=0; y<im.ysize; ++y)
		for (int x=0; x<im.xsize; ++x)
			tmp.set(x,y, 100*tmp.get(x,y));
		tmp.write ("x_component_all_brightened.pgm");
	}
	*/
	
	// Perform a connected components analysis	
	clist = connectedComponents (im, nhType, ignoreBG);
	
	// Travers all the components, calculate their bounding boxes
	// and project the infos back into the matrices
	component_id = 0;
	for (CComponentList::iterator iter=clist->begin(); iter!=clist->end(); ++iter) 
	{	
		/*
		{	Image foo (im.xsize, im.ysize, 1);
			foo.setPlaneToValue (1, 255);
			
			pp = (*iter)->pixels;
			while (pp!=NULL) 
			{
				foo.set(pp->x, pp->y, 0);
				pp = pp->next;
			}
			
			char buf[1000];
			sprintf (buf, "x_component%4.4d", component_id);	
			foo.write(buf);
		}	
		*/
		
		// Fill the height information for this component into the height matrix
		if (height != NULL) 
		{
			int heightv = (*iter)->boundingBox().height();
		
			pp = (*iter)->pixels;
			while (pp!=NULL) 
			{
				height->set (pp->x, pp->y, heightv);
				pp = pp->next;
			}
		}
		
		// Fill the top value for this component into the top matrix		
		if (top != NULL) 
		{
			int topv = (*iter)->boundingBox().top;
		
    		pp = (*iter)->pixels;
    		while (pp!=NULL) 
			{
    			top->set (pp->x, pp->y, topv);
    			pp = pp->next;
    		}
    	}
    	
    	// Fill the component id
		if (compid != NULL) 
		{
    		pp = (*iter)->pixels;
    		while (pp!=NULL) 
			{
    			compid->set (pp->x, pp->y, component_id);				
    			pp = pp->next;
    		}
    	}   
    	 	
    	// Fill the pixel cnt
		if (cnt != NULL) 
		{
			int thiscnt = (*iter)->count;
			
    		pp = (*iter)->pixels;
    		while (pp!=NULL) 
			{
    			cnt->set (pp->x, pp->y, thiscnt);
    			pp = pp->next;
    		}
    	} 
    	
    	++component_id;
	}
	
	delete clist;
}









