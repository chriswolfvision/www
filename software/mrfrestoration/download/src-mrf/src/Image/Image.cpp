// *************************************************************
// Image.cc
// The basic methods for the image class
//
// author: Christian Wolf
// *************************************************************

// C
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#ifdef UNIX_ARCHITECTURE
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

// C++
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <string.h>

// From the main module
#include <CIL.h>

#include <visualbug.h>

using namespace std;

// From the main module
#include <CIL.h>

// From this module
#include "Image.h"
#include "FloatMatrixColor.h"
#include "Rect.h"

#define	abs(x)			( x	> 0	? x	: -(x))
#define	MAX_ANZAHL		4096

// *************************************************************
// The newline for this Architecture

#ifdef UNIX_ARCHITECTURE
#define ARCH_NEWLINE	"\n"
#else
static unsigned char windows_newline [] = {10,0};
#define ARCH_NEWLINE	windows_newline
#endif

// *************************************************************
// Constructor - load from file
// *************************************************************

Image::Image (const char *filename) {
	read(filename);
}

// *************************************************************
// Constructor - load from file
// *************************************************************

Image::Image (const char *filename, bool color) {
	if (color)
		readColor (filename);
	else
		readGray (filename);
}

// *************************************************************
// Constructor - new image
// *************************************************************

Image::Image (int xs, int ys, int colorPlanes) {
	_alloc (xs,	ys,	colorPlanes);
}

// *************************************************************
// Constructor - new image
// *************************************************************

Image::Image (int xs, int ys) {
	_alloc (xs,	ys,	1);
}

// *************************************************************
// Constructor - from FloatMatrix
// *************************************************************

Image::Image (FloatMatrix &other) {
	_alloc (other.xsize, other.ysize, 1);
	_copyFloat (other);
}

// *************************************************************
// Constructor - from color	FloatMatrix
// *************************************************************

Image::Image (FloatMatrixColor &other) {
	_alloc (other.xsize, other.ysize, 3);
	_copyFloat (other);
}

// *************************************************************
// Constructor - from a	one-dimensional	array
// *************************************************************

Image::Image (int xs, int ys, byte *oned, int nrPlanes)	{
	byte *cur;

	_alloc (xs,	ys,	nrPlanes);

	cur	= oned;

	// Grayvalue images
	if (nrPlanes ==	1) {
		for	(int y = 0 ; y < ysize ; y++) {
			for	(int x = 0 ; x < xsize ; x++) {
				set(1,x,y,*cur);
				++cur;
			}
		}
	}

	// Color images
	else {
		for	(int y = 0 ; y < ysize ; y++) {
			for	(int x = 0 ; x < xsize ; x++) {
				set(1,x,y,*cur); ++cur;
				set(2,x,y,*cur); ++cur;
				set(3,x,y,*cur); ++cur;
			}
		}
	}
}

// *************************************************************
// Constructor - from three	one-dimensional	arrays
// *************************************************************

Image::Image (byte *oned[],	long int xsize[], long int ysize[])	{
	byte *cur;
	int	xs,	ys;

	xs = xsize[0];
	ys = ysize[0];


	_alloc (xs,	ys,	3);

	for	(int plane=0; plane<3; ++plane)	{

		cur	= oned[plane];

		for	(int y = 0 ; y < ys; y++) {
			for	(int x = 0 ; x < xs; x++) {
				set(plane,x,y,*cur); ++cur;
			}
		}
	}
}

// *************************************************************
// A list of files - paste them all together
// either (1) vertically
//        (2) horizontally
//		  (3) vertically, more than one tile per line
//        (4) Not used anymore in this function.
//            Use Image::CreateOCRPages() instead
// *************************************************************

Image::Image (const char *filename, StackMethod method, int border, byte bordervalue, bool do_order) {
	istream *st=NULL;
	ifstream *fst=NULL;
	int full_type = 2;
	int full_x=0, full_y=0;
	int run_x, run_y;
	int size=0;
	bool first=true;
	bool isVertical, fixedWidth;

	// Different data structures containing the tiles
	multiset<Image *, image_xless> images_xo;
	multiset<Image *, image_ymore> images_yo;
	multiset<Image *, image_xless>::iterator xiter;
	multiset<Image *, image_ymore>::iterator yiter;
	vector<Image *> images;

	// Set the properties of the different methods
	if ((method==STACKMETHOD_VERTICAL) || (method==STACKMETHOD_VERT_TILED))
		isVertical = true;
	else
		isVertical = false;

	if ((method==STACKMETHOD_VERT_TILED))
		fixedWidth = true;
	else
		fixedWidth = false;

	// Stdin
	if (strcmp(filename,"-")==0)
		st = &cin;

	// A regular file
	else {
    	fst = new ifstream (filename,  ios::in);
    	if (!st->good()) {
			ERR_THROW("can't open file " << filename << "for reading!\n");
			/*
			ostringstream s;
    		s << "can't open file " << filename << "for reading!\n";
    		ERR_THROW(s.str());
			*/
    	}
    	st = fst;
	}

	if ((method==STACKMETHOD_VERT_TILED))
		do_order = true;

	// read the input stream and read all the images.
	// put them into an array.
	while (!st->eof()) {
		string s;

		*st >> s;

		if (s.length()>0) {
			Image *im = new Image (s.c_str());
			if (do_order) {
			 	if (isVertical)
					images_xo.insert(im);
				else
					images_yo.insert(im);
			}
			else
				images.push_back(im);
			++size;
			if (im->type==3)
				full_type = 3;

			// Calculate the full image's size for vertical stacking
			if (isVertical) {
				if (im->xsize > full_x)
					full_x = im->xsize;
				full_y += im->ysize;
				if (!first)
					full_y += border;
			}

			// Calculate the full image's size for horizontal stacking
			else {
				if (im->ysize > full_y)
					full_y = im->ysize;
				full_x += im->xsize;
				if (!first)
					full_x += border;
			}
			first = false;
		}

		if (st->eof())
			break;
		if (st->fail()) {
			ERR_THROW ("*** ERROR: Syntax error in the data!\n");
		}
	}

	// For the tiled methods we always use a fixed page width. If the text
	// boxes are smaller then we repeat them
	if (fixedWidth) {
		if (full_x < OCR_PAGE_WIDTH)
			full_x = OCR_PAGE_WIDTH;
	}

	// Allocate the new stacked image
	_alloc (full_x, full_y, full_type);
	setPlaneToValue (PLANE_RED, bordervalue);
	if (type==3) {
		setPlaneToValue (PLANE_GREEN, bordervalue);
		setPlaneToValue (PLANE_BLUE, bordervalue);
	}

	int i;
	switch (method) {

    	// ----------------------------------------------
    	// Paste all the images into the full image
    	// THE VERTICAL TILED METHOD
    	// ----------------------------------------------

		case STACKMETHOD_VERT_TILED:

        	run_y = 0;
        	yiter = images_yo.begin();
        	for (i=0; i<size; ++i) {

    			Image *im = *yiter;
    			++yiter;

    			run_x = 0;
    			do {
    				paste (*im, run_x, run_y);
    				run_x += im->xsize + border;
    			} while (run_x + im->xsize <= OCR_PAGE_WIDTH);

        		run_y += im->ysize + border;
        		delete im;
        	}
        	break;

    	// ----------------------------------------------
    	// Paste all the images into the full image
    	// HORIZONTAL - VERTICAL
    	// ----------------------------------------------

    	default:

        	run_x = run_y = 0;
        	if (isVertical)
		      	xiter = images_xo.begin();
		    else
			    yiter = images_yo.begin();

        	for (i=0; i<size; ++i) {
        		Image *im;

        		if (do_order)  {
        			if (isVertical) {
	        			im = *xiter;
    	    			++xiter;
    	    		}
    	    		else {
    	    			im = *yiter;
    	    			++yiter;
    	    		}
        		}
        		else
        			im = images[i];

        		paste (*im, run_x, run_y);


        		if (method==STACKMETHOD_VERTICAL)
        			run_y += (im->ysize + border);
        		else
        			run_x += (im->xsize + border);

        		delete im;
        	}
        	break;
    }

	// Clean up
	if (fst!=NULL) {
		fst->close();
		delete fst;
	}
}


// *************************************************************
// Copy	constructor
// *************************************************************

Image::Image (const	Image &other) {
	_alloc (other.xsize, other.ysize, other.type);
	_copy (other);
}

// *************************************************************
// Assignment operator
// *************************************************************

Image &	Image::operator= (const	Image &other){

	if (this==&other)
		return *this;

	// Delete the old object
	_free();
	_alloc (other.xsize, other.ysize, other.type);
	_copy (other);

	return *this;
}

// *************************************************************
// Assignment operator for FloatMatrixColor
// *************************************************************

Image &	Image::operator= (const	FloatMatrixColor	&other){

	// Delete the old object
	_free();
	_alloc (other.xsize, other.ysize, 3);
	_copyFloat (other);

	return *this;
}


// *************************************************************
// Destructor
// *************************************************************

Image::~Image () {
	_free();
}

// *************************************************************
// Help	function for the constructors
// *************************************************************

void Image::_alloc (int	x, int y, int planes) {
	xsize =	x;
	ysize =	y;
	type = planes;

#ifdef CHECK_CODE
	if (x<=0 || y<=0) {
		ERR_THROW ("Invalid image sizes in Image::_alloc():\n"
			 	   << "xsize: " << x << ", ysize: " << y << endl);
	}
#endif

	if ((xsize==0)||(ysize==0))
		return;

	switch (planes)	{

		case 3:
			G =	CREATE_IMAGE (y,x);
			B =	CREATE_IMAGE (y,x);
			/* on purpose no break */

		case 2:
		case 1:
			R =	CREATE_IMAGE (y,x);
			break;

		default:
			ERR_THROW ("Error	in Image::_alloc()\n");
	}
}

// *************************************************************
// Help	function for the desctructor
// *************************************************************

void Image::_free () {

	if ((xsize==0) || (ysize==0))
		return;

	FREE_IMAGE(	R);

	if (type ==	3) {
		FREE_IMAGE(	G);
		FREE_IMAGE(	B);
	}
}

// *************************************************************
// Help	function for the constructors
// *************************************************************

void Image::_copy (const Image &other) {
	xsize =	other.xsize;
	ysize =	other.ysize;
	type = other.type;

	if ((xsize==0)||(ysize==0))
		return;

	memcpy (*R, other.R[0], xsize*ysize*sizeof(byte));

	if (type==3){
		memcpy (*G, other.G[0], xsize*ysize*sizeof(byte));
		memcpy (*B, other.B[0], xsize*ysize*sizeof(byte));
	}
}

// *************************************************************
// Help	function for the constructors
// *************************************************************

void Image::_copyFloat (const FloatMatrix &other) 
{
	float v;
	for	(int y=0; y<ysize; ++y)
	for	(int x=0; x<xsize; ++x)
	{
		v =	other.get(x,y);
		if (v>255.0) v = 255;
		if (v<0.0) v = 0;
		set(1,x,y, (int) (v/1)	);
	}
}

// *************************************************************
// Help	function for the constructors
// *************************************************************

void Image::_copyFloat (const FloatMatrixColor &other) {
	float v;

	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{

			v =	other.get(1,x,y);
			if (v>255.0) v = 255;
			if (v<0.0) v = 0;

			set(1,x,y, (int) (v/1)	);

			v =	other.get(2,x,y);
			if (v>255.0) v = 255;
			if (v<0.0) v = 0;

			set(2,x,y, (int) (v/1)	);

			v =	other.get(3,x,y);
			if (v>255.0) v = 255;
			if (v<0.0) v = 0;


			set(3,x,y, (int) (v/1)	);
		}
	}
}

// *************************************************************
// Sets	a colorplane to	zero
// *************************************************************

void Image::setZero	(int plane)	{
	byte **P = getPlane	(plane);	
	memset (P[0], 0, sizeof(byte)*ysize*xsize);
}

// *************************************************************
// resize the image, don't keep	the	contents
// *************************************************************

void Image::resizeWithDestroy (int xs, int ys, int planes) {
	if (xs==xsize && ys==ysize && planes==type)
		return;
	_free();
	_alloc (xs,	ys,	planes);
}

// *************************************************************
// Reads a grayscale image
// *************************************************************

void Image::readGray (const char *filename) {

	char *buf;
	char shortbuf[256];
	short int x, y;
	int	color, foo;
	char c;
	FILE * inpic;
	int	entete,	z;

	if ( (inpic	= fopen(filename,"r+b")) == NULL)	{
		ERR_THROW ("can't open file '" << filename << "': " << strerror(errno) << endl);
	}

	type = 2 ;
	if (fscanf(inpic,"%c%c\n",&c,&c) !=	2) {
		ERR_THROW ("Image::readGray():\n Wrong Image Format: no .ppm!!\n"
  			 << "filename: " << filename << endl);
	}

	if (c == '6')  {
		z =	3 ;
		ERR_THROW ("Image::readGray():: disabled due to bug.\n"
			"Use Image::readColor() + Image::convertRGB2GrayScale() instead\n");
	}
	else
	{
		if (c != '5')  
			ERR_THROW ("Image::readGray():: wrong image format: "
				"for .ppm only versions P5 and P6 are supported!");
		z =	1 ;
	}

	fscanf(inpic,"%c",&c) ;
	entete = 3 ;
	while (c ==	'#') {
		entete++ ;
		while (c !=	'\n') {
			entete++ ;
			fscanf(inpic,"%c",&c) ;
		}
		fscanf(inpic,"%c",&c) ;
	}

	if ( (inpic	= freopen(filename,"r+b",inpic)) == NULL)	{
		fprintf(stderr,"can't open file '%s': %s\n",
			filename, strerror(errno));
		exit(1);
	}
	fread(shortbuf,1,entete,inpic);

	if (fscanf(inpic,"%d%d\n%d",&xsize,&ysize,&color) != 3)	{
		ERR_THROW ("Image::readGray(): Internal error (2):" << filename << endl);
	}

	fread(shortbuf,1,1,inpic) ;

	buf	= new char [z*xsize+10];

	R =	CREATE_IMAGE(ysize,xsize) ;
	for	( y	= 0	; y	< ysize	; y++) {

		if ((foo=fread(buf,1,z*xsize,inpic)) != z*xsize) {
			ostringstream s;
			s << "file " << filename << ":\nrow " << y << " input failure: "
				<< "got " << foo << " instead of " << z*xsize << " bytes!\n";
			
			if (!feof(inpic))
				s << "No ";			
			s << "EOF occured.\n";
			if (!ferror(inpic))
				s << "No ";			
			s << "error in the sense of ferror() occured.\n";
			ERR_THROW (s.str());
		}
		else {
			if (z == 1)	{
				for	( x	= 0	; x	< xsize	; x++)
					R[x][y]	= buf[x] ;
			}
			else {
				for	( x	= 0	; x	< z*xsize ;	x += z )
					R[x/z][y] =	(int)(.299*(float)buf[x] + 0.587*(float)buf[x+1]
						+ 0.114*(float)buf[x+2]);
			}
		}
	}
	fclose (inpic);
	delete [] buf;
}

// *************************************************************
// Writes a	grayscale image
// *************************************************************

void Image::writeGray(const char *filename) {

	// We write to stdout
	if (strcmp(filename,"-")==0) {
		_writeGray (stdout);
	}
	else {
    	FILE *fp;
    	if ((fp=fopen(filename,"w+b"))==NULL) 
    	{
    		ERR_THROW ("Cannot create output file '" << filename << "': " << strerror(errno) << "!\n");
    	}
    	_writeGray (fp);
    	fclose(fp);
 	}
}

// *************************************************************
// Writes a	grayscale image.
// Two versions, one working with the C library,
// one using UNIX file descriptors (necessary for the
// reliable temp file code
// *************************************************************

#ifdef UNIX_ARCHITECTURE
void Image::_writeGray(int fd) {
	char *buf;
	short int y, x;

	buf = new char [xsize+10];

	sprintf(buf,"P5%s%d	%d%s255%s",ARCH_NEWLINE,xsize,ysize,ARCH_NEWLINE,ARCH_NEWLINE)	;
	x =	strlen(buf);
	if (::write(fd,buf,x)!=x) {
		ERR_THROW ("Could not write image to file (Image::writeGray())!\n");
	}

	for	( y	= 0	; y	< ysize	; y++)	{
		for	( x	= 0	; x	< xsize	; x++ )	{
			buf[x] = R[x][y];
		}

		if (::write(fd,buf,xsize) != xsize ) {
			ERR_THROW ("Could not write image to file (Image::writeGray())!\n");
		} 
	}
	delete [] buf;
}
#endif
void Image::_writeGray(FILE *fp) {
	char *buf;
	short int y, x;

	buf = new char [xsize+10];

	sprintf(buf,"P5%s%d	%d%s255%s",ARCH_NEWLINE,xsize,ysize,ARCH_NEWLINE,ARCH_NEWLINE)	;
	x =	strlen(buf);
	clearerr(fp);
	fwrite(buf,1,x,fp);
	if (ferror(fp)) 
		ERR_THROW ("Could not write image to file (Image::writeGray())!\n");

	for	( y	= 0	; y	< ysize	; y++)	{
		for	( x	= 0	; x	< xsize	; x++ )	{
			buf[x] = R[x][y];
		}

		clearerr(fp);
		fwrite(buf,1,xsize,fp);
		if (ferror(fp))
			ERR_THROW ("Could not write image to file (Image::writeGray())!\n");
	}
	delete [] buf;
}

// *************************************************************
// Reads an image, either color or gray, depending what is
// stored
// *************************************************************
			
void Image::read(const char *filename) {
	char shortbuf[256];
	FILE * inpic;

	if ((inpic	= fopen(filename,"r")) == NULL)
		ERR_THROW ("Cannot open file " << filename << ": " << strerror(errno));
	
	// Read the type flag
	fgets(shortbuf,250,inpic);
	fclose (inpic);

	if (shortbuf[0] != 'P') {
		bool couldload=false;
#ifdef HAVE_LIBJPEG
		couldload=readJPEG (filename);
#endif
		if (couldload)
			return;
		else {
			ERR_THROW ("File '" << filename << "': unsupported file format!\n");
		}
	}

	// The image is in PPM or PGM format
	if (shortbuf[1]=='6')
		readColor (filename);
	else
		readGray (filename);
}


// *************************************************************
// Writes an image, either color or gray, depending what is
// stored
// *************************************************************

void Image::write (const char *filename) {
	if (type==3)
		writeColor (filename);
	else
		writeGray (filename);
}

// *************************************************************
// Writes an image, either color or gray, depending what is
// stored,
// INTO a temporary file whose name pattern is stored in filename.
// the substring XXXXX is replaced by something uniq
// *************************************************************
#ifdef UNIX_ARCHITECTURE
void Image::writeToTempFile (char *filename) {
	int fd;
	if ((fd=mkstemp(filename))==-1)	{
		ERR_THROW ("can't create temp file with pattern " << filename 
			<< ": " << strerror(errno) << endl);
	}

	if (type==3)
		_writeColor (fd);
	else
		_writeGray (fd);

	close(fd);
}
#endif



// *************************************************************
// Reads a color image
// *************************************************************

void Image::readColor(const char *filename) {

	char shortbuf[256];
	char *buf;
	short int	 x, y ;
	FILE * inpic;

	if ((inpic	= fopen(filename,"r+b")) == NULL)	{
		ERR_THROW ("can't open file '" << filename << "': "
			<< strerror(errno)  << endl);
	}

	type = 3;

	fgets(shortbuf,250,inpic);

	if (shortbuf[0] != 'P') {
		bool couldload=false;
		fclose (inpic);

#ifdef HAVE_LIBJPEG
		couldload=readJPEG (filename);		
#endif
		if (couldload)
			return;
		else {
			ERR_THROW ("File '" << filename << "': unsupported file format!\n");
		}
	}
	
	// The image is a grayscale image. Let's read it using readGray()
	// and copy the colour planes afterwards.
	if (shortbuf[1]	== '5')	
	{			
		fclose (inpic);
		readGray (filename);
		G =	CREATE_IMAGE(ysize,xsize) ;
		B =	CREATE_IMAGE(ysize,xsize) ;
		copyPlane (2,1);
		copyPlane (3,1);
		type = 3;		
		return;		
	}
	if (shortbuf[1] != '6')
		ERR_THROW ("Image::readColor():: wrong image format: "
				"for .ppm only versions P5 and P6 are supported!");

	do {
		if (fgets(shortbuf,MAX_ANZAHL,inpic) !=	NULL);
	} while	((*shortbuf=='#')||(*shortbuf=='\n'));

	sscanf(shortbuf,"%d	%d",&xsize,&ysize);
	fgets(shortbuf,MAX_ANZAHL,inpic);

	R =	CREATE_IMAGE(ysize,xsize) ;
	G =	CREATE_IMAGE(ysize,xsize) ;
	B =	CREATE_IMAGE(ysize,xsize) ;


	buf	= new char [3*xsize+10];

	for	( y	= 0	; y	< ysize	; y++) {

		if (fread(buf,1,3*xsize,inpic) != (size_t) 3*xsize)	{
			ERR_THROW ("row	" << y << " input failure!\n");
		}
		else {
			for	( x	= 0	; x	< xsize	; x++ )	{
				R[x][y]	= buf[3*x] ;
				G[x][y]	= buf[3*x+1] ;
				B[x][y]	= buf[3*x+2] ;
			}
		}
	}
	fclose (inpic);
	delete [] buf;
}

// *************************************************************
// Writes a	color image:
// WRAPPER FUNCTION
// *************************************************************

void Image::writeColor(const char *filename) {
	
	// We write to stdout
	if (strcmp(filename,"-")==0) {
		_writeColor (stdout);
	}

	// We write to disk file
	else {
		FILE *fp;
    	if ((fp=fopen(filename,"w+b"))==NULL)	{
    		ERR_THROW ("Can't create file " << filename <<
    			": " << strerror(errno) << endl);       		
    	}
    	_writeColor (fp);
    	fclose(fp);
 	}
}

// *************************************************************
// Writes a	color image
// Two versions, one working with the C library,
// one using UNIX file descriptors (necessary for the
// reliable temp file code
// *************************************************************

#ifdef UNIX_ARCHITECTURE
void Image::_writeColor(int fd) {
    char shortbuf[1024];
	char *buf;
	int	len;

	// Write the PPM Header
	sprintf(shortbuf,"P6%s%d %d%s255%s",ARCH_NEWLINE,xsize,ysize,ARCH_NEWLINE,ARCH_NEWLINE);
	len	= strlen(shortbuf);
	if (::write(fd,shortbuf,len)!=len) {
		ERR_THROW ("Could not write image to file (Image::writeColor())!\n");
	}

	buf = new char [3*xsize+10];

	// Write the data
	for	(int y = 0 ; y < ysize ; y++){
		for	(int x = 0 ; x < xsize ; x++ ) {


			buf[3*x+0] = R[x][y];

			// Real	color image
			if (type ==	3) {
				buf[3*x+1] = G[x][y];
				buf[3*x+2] = B[x][y];
			}

			// Grayscale image:	copy grayplane into	red	and	blue
			else {
				buf[3*x+1] = R[x][y];
				buf[3*x+2] = R[x][y];
			}
		}

		if (::write(fd,buf,3*xsize) != 3*xsize )
			ERR_THROW ("Could not write image to file (Image::writeColor())!\n");
			
	}
	delete [] buf;
}
#endif
void Image::_writeColor(FILE *fp) {
    char shortbuf[1024];
	char *buf;
	int	len;

	// Write the PPM Header
	sprintf(shortbuf,"P6%s%d %d%s255%s",ARCH_NEWLINE,xsize,ysize,ARCH_NEWLINE,ARCH_NEWLINE);
	len	= strlen(shortbuf);
	clearerr(fp);
	fwrite(shortbuf,1,len,fp);
	if (ferror(fp))
		ERR_THROW ("Could not write image to file (Image::writeColor())!\n");
		
	buf = new char [3*xsize+10];

	// Write the data
	for	(int y = 0 ; y < ysize ; y++){
		for	(int x = 0 ; x < xsize ; x++ ) {



			buf[3*x+0] = R[x][y];

			// Real	color image
			if (type ==	3) {
				buf[3*x+1] = G[x][y];
				buf[3*x+2] = B[x][y];
			}

			// Grayscale image:	copy grayplane into	red	and	blue
			else {
				buf[3*x+1] = R[x][y];
				buf[3*x+2] = R[x][y];
			}
		}

		clearerr(fp);
		fwrite(buf,1,3*xsize,fp);
		if (ferror(fp))
			ERR_THROW ("Could not write image to file (Image::writeColor())!\n");
	}
	delete [] buf;
}

// *************************************************************
// INTERNAL: Creates an	internal Image structure, used in the
// Image class to store	a color	plane.
// *************************************************************

unsigned char **CREATE_IMAGE (int ysize, int xsize)	{

	unsigned char ** im;
	unsigned char *big;
	typedef byte * pbyte;

	im = new pbyte [xsize];
	big	= new byte [xsize*ysize];

	/*
	memset (im,	0, xsize*sizeof(byte));
	memset (big, 0,	xsize*ysize*sizeof(byte));
	*/


	for	(int i = 0 ; i < xsize ; i++)
		im[i] =	big	+ i*ysize;	

	return (im);
}

// *************************************************************
// INTERNAL: FREEs an internal Image structure
// *************************************************************

void FREE_IMAGE	(byte **im)	
{
	delete [] im[0];
	delete [] im;
}

// *************************************************************
// INTERNAL: Clears	an internal	Image structure
// *************************************************************


void CLEAR_IMAGE (int xsize, int ysize,	byte **im) 
{
	memset (im[0], 0, xsize*ysize*sizeof(byte));
}

// *********************************************************
// Mark	a pixel	with high contrast
// *********************************************************


void Image::mark (int x, int y)	
{

	R[x][y]	= R[x][y]<128 ?	255	: 0;

	if (type==3) 
	{
		G[x][y]	= G[x][y]<128 ?	255	: 0;
		B[x][y]	= B[x][y]<128 ?	255	: 0;
	}
}

// *********************************************************
// Mark	a pixel	with high contrast in a given plane, and
// mark a cross if wished.
// *********************************************************

void Image::mark (int plane, int x, int y, bool docross)	
{

	set (plane, x, y, get(plane,x,y) < 128 ? 255 : 0);


	if (docross) 
	{
    	if (x>0) set (plane, x-1, y, get(plane,x-1,y) < 128 ? 255 : 0);
    	if (y>0) set (plane, x, y-1, get(plane,x,y-1) < 128 ? 255 : 0);
    	if (x<xsize-1) set (plane, x+1, y, get(plane,x+1,y) < 128 ? 255 : 0);
    	if (y<ysize-1) set (plane, x, y+1, get(plane,x,y+1) < 128 ? 255 : 0);
    }
}

// *********************************************************
// Draw a horizontal dashed line with a fixed count of
// dashes. Used to debug image algorithms (If a number needs
// to written in the pixel, write a count of dashes instead)
// *********************************************************

void Image::markDashes (int xb, int yb, int count, int dashlength) 
{
	int u;

	cerr << "[Dashes: " << count << "]";

	for (int i=0; i<count; ++i) 
	{	
		for (u=0; u<dashlength; ++u) 
		{
			mark (xb, yb);
			++xb;
		}
		for (u=0; u<dashlength; ++u)

			++xb;
	}			
}

// *********************************************************
// draws text (digits only) into the image
// Returns the length in pixels (width) of the space which
// has been used.
// *********************************************************

int Image::drawText (int xb, int yb, const char *s, bool isBig) {
	int curpos, slen=strlen(s);
	Image *digits;
	int *widths;
	int distance;
	
	// Load the digits file
	if (isBig) {
	    digits = GetDigitsBig();
	    widths = Image::digitwidthsBig;
	    distance = Image::digitdistanceBig;

	}	
	else {
		digits = GetDigits();
		widths = Image::digitwidths;
		distance = Image::digitdistance;
	}

    curpos = 0;
	for (int i=0; i<slen; ++i) {	
	
	    char c=s[i];
	    int index = (c-'0');

		// Special characters:
	    if (c<'0'||c>'9') {
	    	switch(c) {
			  	case '-':
	     			index=11;
	       			break;
	       		default:
	    			continue;
	     	}
	    }
	    	
	    // We touch the right border
	    if (xb+curpos+widths[index] >= xsize)
	    	return curpos;
	    	  	
		paste (*digits, xb+curpos, yb,
			index*distance, 0, widths[index], digits->ysize);
			
		curpos += widths[index];
		if (isBig)
			curpos += (distance / 10);
	}
	
	if (isBig)
		curpos -=  (distance / 10);
	
	return curpos;
}

// *********************************************************
// draws text (digits only) into the image
// takes an integer value as input
// *********************************************************


int Image::drawText (int xb, int yb, int value, bool isBig) {
	char buf[100];
	sprintf (buf,"%d",value);
	return drawText (xb, yb, buf, isBig);
}

// *************************************************************
// A small routine which just draws	the	borders	of the
// components
// *************************************************************

void Image::componentBorders () {

	int	i, j, k, l,	m ;
	unsigned char **I, **J;

	J =	(unsigned char **) CREATE_IMAGE(ysize,xsize) ;
	I =	R;

	for	(i = 0 ; i < xsize ; i++)
	for	(j = 0 ; j < ysize ; j++)
		J[i][j]	= 0	;

	for	(i = 1 ; i < xsize-1 ; i++)
	for	(j = 1 ; j < ysize-1 ; j++)
	if (I[i][j]	== 255)	{
		m =	0 ;
		for	(k = -1	; k	<= 1 ; k++)
		for	(l = -1	; l	<= 1 ; l++)
			if (I[i+k][j+l]	!= 255)	m++	;

		if ((m < 8)	&& (m >	1))	J[i][j]	= 255 ;
	}

	R =	J;
}

// *************************************************************
// DEBUG PRINT
// *************************************************************

void Image::print (FILE	*fp) {
	int	i,j;


	for	(i = 0 ; i < xsize ; i++)
	for	(j = 0 ; j < ysize ; j++)
		fprintf	(fp, "(%d)",R[i][j]);
}

// *************************************************************

// Cuts a sub image out of an image
// *************************************************************

Image *Image::cutSubImage (Rect r, int growX,	int	growY) 
{
	Image *im;
	Rect b;
	int width, height;

	b =	r;
	b.growAndClip (growX, growY, xsize, ysize);

	// Check the size
	width = b.width();

	height = b.height();
	if ((width==0) || (height==0))

		return NULL;
		
#ifdef CHECK_CODE
	if (r.right>=xsize || r.bottom>=ysize) {
		cerr << "Warning in cutSubImage(), has been corrected\n"
			 << "rect: " << r
			 << ", image size: " << xsize << "x" << ysize << endl;

		r.right=xsize-1;
		r.bottom=ysize-1;
	}
#endif		
	
	// Create the image		
	im = new Image (width, height, type);
	im->setZero();
	
	// Copy	the	part of	the	big	original image
	for	(int y=0; y<b.height();	++y) {
		for	(int x=0; x<b.width(); ++x)	{

			im->set	(1,x,y,	get(1,b.left+x,b.top+y));

			// Color images
			if (im->type ==	3) {
				im->set	(2,x,y,	get(2,b.left+x,b.top+y));
				im->set	(3,x,y,	get(3,b.left+x,b.top+y));
			}
		}
	}
	
	return im;
}			
