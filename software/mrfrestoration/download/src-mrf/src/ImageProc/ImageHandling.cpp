/***********************************************************************
 * Image Handling
 *
 * Author: Christian Wolf, christian.wolf@insa-lyon.fr
 * Begin: 16.6.2006
 ***********************************************************************/
 
// C++
#include <map> 
 
// From the IMAGE PROCESSING module
#include "ImageProc.h"

/**************************************************************
 * A class used in the SwapColor() function
 **************************************************************/ 

class SwapColor
{
	public:
	
		SwapColor () { Depth=0; }
		
		void load (unsigned char v) 	
		{
			Depth=1;
			vals[0]=v;
			vals[1]=v;
			vals[2]=v;
		}
		void load (unsigned char r, unsigned char b, unsigned char g)
		{
			Depth=3;
			vals[0]=r;
			vals[1]=g;
			vals[2]=b;
		}		
		void add (unsigned char v)
		{
			vals[Depth++]=v;
			if (Depth>3)
				ERR_THROW ("Syntax error: too many color bands!");
		}
		
		bool operator == (const SwapColor &o) const
		{
			for (unsigned int i=0; i<Depth; ++i)
			{
				if (vals[i]!=o.vals[i])
					return false;
			}
			return true;			
		}
		
		bool operator < (const SwapColor &o) const
		{
			for (unsigned int i=0; i<Depth; ++i)
			{
				if (vals[i]<o.vals[i])
					return true;
				else
					if (vals[i]>o.vals[i])
						return false;					
			}
			return false;			
		}
	
		unsigned char vals[3];
		unsigned char Depth;
};

/**************************************************************
 * Parse a line in the file specifying the rules for 
 * swapping colors
 **************************************************************/ 

#define EAT			{ while (isspace(*cp) && *cp!='\0' && *cp != ':' && *cp!=';')\
					  ++cp; }

static void parseLine (char *cp, SwapColor &src, SwapColor &dst)
{		
	char *ep;
		
	src.Depth=0;	
	do 
	{
		EAT;
		if (*cp==':' || *cp==0)
			break;
		src.add(strtol (cp, &ep, 10));
		if (ep==cp)
			ERR_THROW ("syntax error: number expected!"); 
		cp=ep;	
	} while (1);
			
	if (*cp!=':')
		ERR_THROW ("syntax error: ':' expected");	
	++cp;
		
	dst.Depth=0;
	do 
	{
		EAT;
		if (*cp==';' || *cp==0)
			break;
		dst.add(strtol (cp, &ep, 10));
		if (ep==cp)
			ERR_THROW ("syntax error: number expected: " << cp); 
		cp=ep;	
	} while (1);
		
	if (*cp!=';')
		ERR_THROW ("syntax error: ';' expected");	
}

/**************************************************************
 * Change the colors or gray values in an image according to
 * the specification in a text file
 **************************************************************/ 
  
Image * ColorSwap (Image *srcim, char *fname) 
{	
	char buf[2048], *cp;
	SwapColor src,dst,*dstp;
	int srcDepth, dstDepth;
	map<SwapColor, SwapColor> rulemap;
	map<SwapColor, SwapColor>::iterator iter;
	Image *dstim;
	bool firstLine=true;
	
	ifstream st (fname,  ios::in);
   	if (!st.good())
   		ERR_THROW ("Cannot open file " << fname << "for read!!\n");
    		    		
	while (st.getline (buf,2047)) 
	{
		cp=buf;
		EAT;
	
		// empty line?
		if (*cp==0)
			continue;
	
		parseLine (cp, src, dst);
		
		if (firstLine)
		{
			firstLine=false;
			srcDepth=src.Depth;
			dstDepth=dst.Depth;
			if (srcDepth!=srcim->nbColorPlanes())
				ERR_THROW ("Image depth does not correspond to rule depth!");			
			dstim = new Image (srcim->xsize, srcim->ysize, dstDepth);
		}
		
		if (srcDepth!=src.Depth || dstDepth!=dst.Depth)
			ERR_THROW ("Rules files does not have consistent color depths!");
					
		rulemap[src]=dst;
	}
	st.close();			
	
	cerr << "ColorSwap: " << rulemap.size() << " swap definitions.\n";
	
	// Travers the image and apply the rules
	for (int y=0; y<srcim->ysize; ++y)
	for (int x=0; x<srcim->xsize; ++x)
	{
		if (srcDepth==3)
			src.load (srcim->get(PLANE_RED,x,y),
				 	  srcim->get(PLANE_GREEN,x,y),
				 	  srcim->get(PLANE_BLUE,x,y));
		else
			src.load (srcim->get(PLANE_RED,x,y));
							 
		iter = rulemap.find (src);
		if (iter==rulemap.end())
			dstp=&src;
		else
		{
			dst=iter->second;
			dstp=&dst;
		}
				
		dstim->set(PLANE_RED,x,y,dstp->vals[0]);
		if (dstDepth==3)
		{					
			dstim->set(PLANE_GREEN,x,y,dstp->vals[1]);
			dstim->set(PLANE_BLUE,x,y,dstp->vals[2]);
		}		
	}
	
	return dstim;
}

