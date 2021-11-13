/***********************************************************************
 * A hierarchical pyramid structure used to replace the verso pixels
 * in a degraded document image
 *
 * Author: Christian Wolf, chriswolf@gmx.at
 * Begin: 28.10.2005
 ***********************************************************************/

#ifndef	_WOLF_VERSOREPLACEMENTPYRAMID_H_
#define	_WOLF_VERSOREPLACEMENTPYRAMID_H_

// C
#include <math.h>

// C++
#include <vector>
#include <queue>
#include <set>

// From the main module
#include <CIL.h>

// From the IMAGE module
#include <FloatMatrix.h>
#include <FloatMatrix.h>
#include <Image.h>

// From this module
#include "CommonBleedThrough.h"

/*********************************************************************************
 * Configuration
 *********************************************************************************/

#define FSIZE				3
#define PARENT_FIFO_SIZE	30
#define PIXEL_ARRAY_SIZE	30

/*********************************************************************************
 * Two helper classes
 *********************************************************************************/

class Node
{
	public: 
		Node () 							{}
		Node (unsigned int xl, unsigned int xx, unsigned int xy, unsigned int xc) 	
											{ l=xl; x=xx; y=xy; cnt=xc;} 	
		unsigned int l,x,y,cnt;
};

class Node4Med : public Node
{
	public:
		Node4Med(const Node &n, unsigned char xo) 
			: Node (n.l,n.x,n.y,n.cnt) { obs=xo; }
		bool operator < (const Node4Med &ot) const { return obs<ot.obs; }
		unsigned char obs;					// the observation		
};

class FIFO
{
	public:
		FIFO()								{ cnt=0; }
		void push (Node &n)					{ q.push(n); cnt+=n.cnt;
											  // cerr << "+(" << n.x << "," << n.y << "|" << n.l << "#"
											  //	     	 << n.cnt << ")"; 
											}
		Node pop()							{ Node rv; rv=q.front(); q.pop(); cnt-=rv.cnt; return rv; 
											  // cerr << "-";
											}
		unsigned int size()					{ return q.size(); }
		
		queue<Node> q;
		unsigned int cnt;
		
		void print()						{ 	while (q.size()>0) 
												{ 
													Node n=q.front(); q.pop();
													cerr << "(" << n.x << "," << n.y << "|" << n.l << "#"
											  	     	 << n.cnt << ")"; 
											  	}
											}
};

/*********************************************************************************
 * CLASS
 *********************************************************************************/
 
template <class TI>
class VersoReplacementPyramid 
{
	public:
	
		typedef typename TI::PixelType PixelType;

		// Constructor and destructor
		VersoReplacementPyramid	(TI *obs, Image *lab, 
			unsigned char labSrc, unsigned char labDst, unsigned int nLevels, 
			unsigned int xMinCount, unsigned char xBGClassMean);
		~VersoReplacementPyramid ();
				
		void restore();
		void debugOutput (char *filename);

	protected:	// METHODS
	
		unsigned int chi2par     (unsigned int i) { return i/2; }
		void addParent (int x, int y, unsigned int l);
		void getParents(unsigned int x, unsigned int y, unsigned int l);
		void combineNodes (FIFO &fifo, unsigned char *result);
		void searchBGPixelsAndReplace (FIFO &fifo, unsigned int x, unsigned int y); 

		void initMask();
		void build (int begLevel=0, int endLevel=-1);
 		void replace();
 				
	private:	// DATA
	
		vector<TI *> levelVector;	
		vector<Image *>maskVector;
		Image *lab, *labv;
		unsigned int minCount;
		unsigned char labSrc, labDst, BGClassMean;
		
		// help data for the replace() function
		vector<Node> parents;
		unsigned int parCount;
};

// ***************************************************************
// Constructor
// ***************************************************************

template <class TI> 
VersoReplacementPyramid<TI>::VersoReplacementPyramid (TI *obs, Image *xlab,
	unsigned char xLabSrc, unsigned char xLabDst, 
	unsigned int nLevels, unsigned int xMinCount, unsigned char xBGClassMean) 
{
	cerr << "VersoReplacementPyramid.\n";
	for (unsigned int i=0;i<nLevels; ++i)
	{
		levelVector.push_back(NULL);
		maskVector.push_back(NULL);
	}
	levelVector[0]=obs;
	lab=xlab;	
	labSrc=xLabSrc;
	labDst=xLabDst;
	minCount = xMinCount;
	parents.resize(minCount);
	BGClassMean = xBGClassMean;
}

// ***************************************************************
// Initialise the first level of the mask.
// ***************************************************************

template <class TI>
void VersoReplacementPyramid<TI>::initMask () 
{
	unsigned int xs=levelVector[0]->xsize;
	unsigned int ys=levelVector[0]->ysize;
	
	if (maskVector[0]==NULL)
		maskVector[0] = new Image (xs, ys);
		
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		// If we have a background pixel, enable the mask pixel
		if (lab->get(x,y)==labSrc)
			maskVector[0]->set(x,y, 1);
		else
			maskVector[0]->set(x,y, 0);
	}
}

// ***************************************************************
// Destructor
// ***************************************************************

template <class TI> 
VersoReplacementPyramid<TI>::~VersoReplacementPyramid () 
{
	// The first level is a _POINTER_ to the input image, the other ones
	// are created.
	for (unsigned int i=1; i<levelVector.size(); ++i) 
	{
		if (levelVector[i]!=NULL)
			delete levelVector[i];
	}	
	for (unsigned int i=0; i<levelVector.size(); ++i) 
	{
		if (maskVector[i]!=NULL)
			delete maskVector[i];			
	}	
}

/****************************************************************************
 * Calculate the mean or the median
 * median: sort the data with a	variation of bubble	sort (is fastest with
 * low data	amount,	and	we don't need to sort the whole array			 
****************************************************************************/

template <class T>
T combineValues (T *arr, unsigned int noChildren)
{	
	if (noChildren<1)
		return 255;
	
#ifdef RESTAURATION_WITH_MEDIAN		
	// Calculate the median
	T foo;
	for	(unsigned int curSort=0; curSort<=noChildren; ++curSort)
	for	(unsigned int i=curSort; i<noChildren; ++i) 
	{
		if (arr[curSort] > arr[i]) 
		{
			foo	= arr[curSort];
			arr[curSort] = arr[i];
			arr[i] = foo;
		}
	}	
	return arr[noChildren/2];
#else
	// Calculate the mean
	float mean=0;	
	for	(unsigned int i=0; i<noChildren; ++i) 
		mean += arr[i];
	mean=rint(mean/(float)noChildren);
	TRIMGRAY(mean);
	return (T) mean;
#endif	
}	
		
/****************************************************************************
 * Build the pyramid
 ****************************************************************************/
 
#define ADD(i,j)	{ unsigned int np=CHIM->get(i,j); \
					  if (np) { arr[noChildren]=CHI->get(i,j); \
                                noPixels += np; \
                                ++noChildren; } }
template <class TI>
void VersoReplacementPyramid<TI>::build (int begLevel, int endLevel)	
{
	// the array holding the child observations
	PixelType arr[FSIZE*FSIZE];
	unsigned int noChildren, noPixels;
	
	if (levelVector[0]->nbColorPlanes()>1)
		ERR_THROW ("Color replacement not yet implemented");

	if (begLevel<0)
		begLevel=0;
	if (endLevel==-1 || endLevel>=(int)levelVector.size())
		endLevel=levelVector.size()-1;
		
	if (levelVector[begLevel] == NULL) 
		ERR_THROW("Internal error in VersoReplacementPyramid<TI>::build()\n");
		
	// Travers the levels
	for(int l=begLevel+1;l<=endLevel;l++) 
	{
		unsigned int endx, endy;
		TI *CHI, *PAR;
		Image *CHIM, *PARM;
				
		// The border
		unsigned int bottomysize=levelVector[l-1]->ysize;
		unsigned int bottomxsize=levelVector[l-1]->xsize;
		
		cerr << "Level l=" << l << "\txs=" << bottomxsize
			 << "\tys=" << bottomysize << endl;

		// Allocation of the level
		if (levelVector[l] == NULL)
			levelVector[l] = new TI (chi2par(bottomxsize), chi2par(bottomysize));		
		if (maskVector[l] == NULL)
			maskVector[l]  = new Image (chi2par(bottomxsize), chi2par(bottomysize));		
					
		PAR = levelVector[l];
		CHI = levelVector[l-1];
		PARM= maskVector[l];
		CHIM= maskVector[l-1];		
						
		/* --------------------------------------------
		 * The main part of the image - not the borders
		 * does not check for borders, is faster
		 * -------------------------------------------- */
		 
		endx = (bottomxsize%2 ? bottomxsize-1 : bottomxsize-2);
		endy = (bottomysize%2 ? bottomysize-1 : bottomysize-2);				
		
		for (unsigned int y=1; y<=endy; y+=2) 
		for (unsigned int x=1; x<=endx; x+=2)  
		{
			/* Calculate the filter	boundaries */	
			unsigned int bx = x-(FSIZE/2);
			unsigned int ex = x+(FSIZE/2);
			unsigned int by = y-(FSIZE/2);    			
			unsigned int ey = y+(FSIZE/2);
			
			noChildren=noPixels=0;
			for	(unsigned int ly=by; ly<=ey; ++ly)
			for	(unsigned int lx=bx; lx<=ex; ++lx)
				ADD(lx,ly);
			
			// Collect the results
			PARM->set(chi2par(x), chi2par(y), (noPixels>255?255:noPixels));
			PAR->set (chi2par(x), chi2par(y), combineValues(arr,noChildren));			
		}
		
		/* --------------------------------------------
		 * The border treatment
		 * Slower than the main part because it does
		 * the necessaryu border checking.
		 * We ignore the unavailable values
		 * -------------------------------------------- */
		 
		if (endx<bottomxsize-1)
		{
			//cerr << "special: column x=" << bottomxsize-1 << " (parent: " 
			//     << (bottomxsize-1)/2 << ")" << endl;
			unsigned int x = bottomxsize-1;
			for (unsigned int y=1; y<=endy; y+=2) 
			{							
				/* Calculate the filter	boundaries */	
				unsigned int bx = x-(FSIZE/2);
				unsigned int ex = x;
				unsigned int by = y-(FSIZE/2);    			
				unsigned int ey = y+(FSIZE/2);
				
				noChildren=noPixels=0;
				for	(unsigned int ly=by; ly<=ey; ++ly)
				for	(unsigned int lx=bx; lx<=ex; ++lx)
					ADD(lx,ly);
				
				// Collect the results
				PARM->set(chi2par(x), chi2par(y), (noPixels>255?255:noPixels));
				PAR->set (chi2par(x), chi2par(y), combineValues(arr,noChildren));
			}
		}
		
		if (endy<bottomysize-1)
		{
			//cerr << "special: row y=" << bottomysize-1 << " (parent: " << (bottomysize-1)/2 << ")" << endl;
			unsigned int y = bottomysize-1;
			for (unsigned int x=1; x<=endx; x+=2) 
			{							
				/* Calculate the filter	boundaries */	
				unsigned int bx = x-(FSIZE/2);
				unsigned int ex = x+(FSIZE/2);
				unsigned int by = y-(FSIZE/2);			
				unsigned int ey = y;
				
				noChildren=noPixels=0;
				for	(unsigned int ly=by; ly<=ey; ++ly)
				for	(unsigned int lx=bx; lx<=ex; ++lx)
					ADD(lx,ly);				
				
				// Collect the results
				PARM->set(chi2par(x), chi2par(y), (noPixels>255?255:noPixels));
				PAR->set (chi2par(x), chi2par(y), combineValues(arr,noChildren));
			}			
		}
		
		if (endx<bottomxsize-1 && endy<bottomysize-1)
		{							
			//cerr << "special: pixel " << bottomxsize-1 << "," << bottomysize-1 
			//	 << " (parent: " << (bottomxsize-1)/2 << "," << (bottomysize-1)/2 << ")" << endl;
			/* Calculate the filter	boundaries */	
			unsigned int  x = bottomxsize-1;
			unsigned int  y = bottomysize-1;
			unsigned int bx = x-(FSIZE/2);
			unsigned int ex = x;
			unsigned int by = y-(FSIZE/2);			
			unsigned int ey = y;
			
			noChildren=noPixels=0;
			for	(unsigned int ly=by; ly<=ey; ++ly)
			for	(unsigned int lx=bx; lx<=ex; ++lx)
				ADD(lx,ly);
			
			// Collect the results
			PARM->set(chi2par(x), chi2par(y), (noPixels>255?255:noPixels));
			PAR->set (chi2par(x), chi2par(y), combineValues(arr,noChildren));			
		}
	}   
}
#undef ADD

/***************************************************************************
 * Get the indices of the parent for a given child index
 * take into account that for levels of odd size, the first line/column
 * is duplicated.
 ***************************************************************************/

template <class TI> 
void VersoReplacementPyramid<TI>::addParent(int x, int y, unsigned int l)
{
	if ((x>=0) && (x<levelVector[l]->xsize) && (y>=0) && (y<levelVector[l]->ysize))
	{
		parents[parCount].x=x;
		parents[parCount].y=y;
		parents[parCount].l=l;
		parents[parCount].cnt=maskVector[l]->get(x,y);
		++parCount; 
	}
}

template <class TI> 
void VersoReplacementPyramid<TI>::getParents(unsigned int x, unsigned int y, unsigned int l) 
{
	int i=(int)x/2;
	int j=(int)y/2;
	parCount=0;
	
	// x: odd
	if (x%2)
	{
		// x & y odd
		if (y%2)
		{
			addParent((int)i,(int)j,l+1);
		}
		
		// x odd, y even
		else
		{
			addParent((int)i,(int)j-1,l+1);
			addParent((int)i,(int)j,  l+1);
		}		
	}
	
	// x: even
	else
	{
		// x even, y odd
		if (y%2)
		{
			addParent((int)i-1,(int)j,l+1);
			addParent((int)i,  (int)j,l+1);
		}
		
		// x and y even
		else
		{
			addParent((int)i-1, (int)j-1,l+1);
			addParent((int)i,   (int)j-1,l+1);
			addParent((int)i-1, (int)j,  l+1);
			addParent((int)i,   (int)j,  l+1);
		}
	}	
}
#undef ADD

// ***************************************************************
// Combine the pixels which are currently in the queue
// ***************************************************************

template <class TI>
void VersoReplacementPyramid<TI>::combineNodes (FIFO &fifo, unsigned char *result) 
{ 		
	unsigned int nbPlanes = levelVector[0]->nbColorPlanes();
#ifdef RESTAURATION_WITH_MEDIAN	
	// A color image: we calculate the mean
	if (nbPlanes>1)
	{
#endif	
		unsigned int cnt=0;
		float mean[3];
		
		DO_COLORS_IF(nbPlanes)
			mean[col-1]=0;
		while (fifo.size()>0)
		{
			Node n=fifo.pop();
			DO_COLORS_IF(nbPlanes)
 				mean[col-1]+=n.cnt*(levelVector[n.l]->get(col,n.x,n.y));
 			cnt+=n.cnt;
		}	
		DO_COLORS_IF(nbPlanes)
		{
			mean[col-1] = rint(mean[col-1]/(float)cnt);
			TRIMGRAY(mean[col-1]);
			result[col-1]=(unsigned char)mean[col-1];		
		}		
			
		
#ifdef RESTAURATION_WITH_MEDIAN			
	}		
	// A gray scale image: we calculate the mean or the median
	else
	{	
		// Insert all the nodes plus their observed values into a sorted set
		set<Node4Med> nodeset;
		unsigned int midIndex=0,curCnt;
		int oldVal,curVal,oldCnt;
		set<Node4Med>::iterator iter;
		
		while (fifo.size()>0)
		{
			Node n=fifo.pop();				
 			nodeset.insert(Node4Med(n,levelVector[n.l]->get(n.x,n.y))); 		
			midIndex+=n.cnt;
		}
		
		// Get the middle of the set
		midIndex/=2;
		iter=nodeset.begin();
		oldCnt=curCnt=0;
		oldVal=curVal=-1;
		while (iter!=nodeset.end())
		{
			curVal=iter->obs;			
						
			if (curCnt>=midIndex)
				break;
										
			oldCnt=curCnt;
			curCnt+=iter->cnt;
			oldVal=curVal;
			++iter;			
		}	
					
		// Let's see whether this value or the previous value are the good one.
		if ((midIndex-oldCnt<=curCnt-midIndex) && (oldVal!=-1))
			return oldVal;
		else
			return curVal;	
	}
#endif						
}

// ***************************************************************
// Each site s has possibly more than one parent, so we can't just
// walk up the pyramid and search for the first site having enough
// pixels. We have to combine the values of the different parents.
// ***************************************************************

template <class TI>
void VersoReplacementPyramid<TI>::searchBGPixelsAndReplace (FIFO &fifo,
	unsigned int x, unsigned int y) 
{
	Node curNode;
	unsigned char result[3];
	
	while (fifo.cnt < minCount)
	{		
		if (fifo.size()<=0)
			ERR_THROW ("VersoReplacementPyramid: fifo is empty!!");
			
		// Get the first node from the FIFO	
		curNode = fifo.pop();	
		
		// Are we at the last level? Quite improbable but possible if the 
		// verso area is really large. We just take the mean value for this class		
		if (curNode.l>=levelVector.size()-1)
		{	
			if (levelVector[0]->nbColorPlanes()>1)
				ERR_THROW ("VersoReplacementPyramid<TI>::searchBGPixelsAndReplace():"
					<< "Color not yet implemented.\n");
		
			levelVector[0]->set(x,y,BGClassMean);
			return;
		}		
		
		// Search its parents and store them in the (parents[],ParCount) structure		
		getParents(curNode.x, curNode.y, curNode.l);				
			
		// Add all the parents to the fifo 
		for (unsigned int pi=0; pi<parCount; ++pi)
			fifo.push(parents[pi]);
	}
	
	// We have enough pixels		
	combineNodes(fifo, result);
	DO_COLORS_IF(levelVector[0]->nbColorPlanes())
		levelVector[0]->set(col,x,y,result[col-1]);						
}

// ***************************************************************
// Replace the verso pixels in the bottom image. 
// Search using a bottom up algorithm.
// ***************************************************************

template <class TI>
void VersoReplacementPyramid<TI>::replace () 
{
	unsigned int xs=levelVector[0]->xsize;
	unsigned int ys=levelVector[0]->ysize;
	FIFO fifo;	
	
	// Traverse the image
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		// Is this a verso pixel?
		// Start the breadth-first bottom up process	
		if (lab->get(x,y)==labDst)
		{			
			Node n(0,x,y,maskVector[0]->get(x,y));		
			while (fifo.size()>0) 
				fifo.pop();
			// cerr << "BEGIN: ";
			fifo.push(n);			
			searchBGPixelsAndReplace(fifo,x,y);
			// cerr << "END;\n";
		}
	}
}

// ***************************************************************
// Do the whole thing
// ***************************************************************

template <class TI>
void VersoReplacementPyramid<TI>::restore () 
{
	initMask();
	build();	
	replace();
}

// ***************************************************************
// Debug output of the whole pyramid
// ***************************************************************

template <class TI>
void VersoReplacementPyramid<TI>::debugOutput (char *filename) 
{	
	for (unsigned int l=0; l<levelVector.size(); ++l)
	{
		ostringstream s;
		Image tmp (levelVector[l]->xsize, levelVector[l]->ysize);
		s << filename << "_obs_lev" << l << ".pgm";
		for (int y=0; y<tmp.ysize; ++y)
		for (int x=0; x<tmp.xsize; ++x)
			if (maskVector[l]->get(x,y)<=0)
				tmp.set(x,y,255);
			else
				tmp.set(x,y,levelVector[l]->get(x,y));				
		tmp.write(s.str().c_str());
		
		ostringstream s2;
		s2 << filename << "_mask_lev" << l << ".pgm";
		maskVector[l]->write(s2.str().c_str());
	}
}

#endif
