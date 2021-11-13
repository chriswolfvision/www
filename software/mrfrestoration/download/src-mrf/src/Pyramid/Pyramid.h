#ifndef	_WOLF_PYRAMID_H_
#define	_WOLF_PYRAMID_H_

// C++
#include <vector>

// From the IMAGE library
#include <Image.h>
#include <FloatMatrix.h>

/*********************************************************************************
 * A Pyramid
 * Base level (finest resolution): 0
 * Top (coarsest resolution): levels()-1
 *********************************************************************************/
 
class RectTextList;

class Pyramid 
{

	public:

		// Constructors
		Pyramid	(Image *inputimage, int nLevels);
  		Pyramid (Image *inputimage, int maxxsize, int maxysize, int minxsize, int minysize);

		// Destructor
		~Pyramid ();
		
		void changeImage(Image *inputimage);

		// The iterators
		typedef	vector<Image *>::iterator iterator;
		iterator begin() { return levelVector.begin(); }
		iterator end()	{ return levelVector.end(); }

		// Access to the levels
		Image *& operator[]	(int level)	{ return levelVector[level]; }
		Image &	baseImage()	{ return *levelVector[0]; }
		
		// Return the image of the level and unlink it,
		// i.e. do not destroy it when freeing the pyramid
		Image * getAndUnlinkLevel (int level);

		// Build the Pyramides
		void gaussian (int begLevel, int endLevel);		

		int	levels() { return levelVector.size(); }

		void write	(char	*filename);

		RectTextList * locateText(bool displayResults);

	protected:	// METHODS

 		void _alloc (Image *inputimage, int nLevels);	
		
	protected:	// DATA

		int imagetype;
		vector<Image *> levelVector;	
		
	private:	// DATA
	
		FloatMatrix *reductionBuffer;			

};

#endif

