// C
#include "errno.h"
#include "string.h"

// C++
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <functional>

#ifndef TELESUN_INSA
#include <sstream>
#endif

using namespace std;

// From the main module
#include <CIL.h>

// From this module
#include "Image.h"

#define MAX_ANZAHL	200

// *********************************************************************
// STATIC CLASS VARIABLES
// Mainly necessary for the text drawing functions
// *********************************************************************

int Image::digitwidths[] = {6, 5, 6, 5, 6, 7, 6, 7, 6, 7, 0, 6};	
int Image::digitdistance = 10;

int Image::digitwidthsBig[] = {65, 43, 66, 59, 68, 61, 64, 61, 58, 63, 0, 39};	
int Image::digitdistanceBig = 100;

Image *Image::_Digits 	 = NULL;
Image *Image::_DigitsBig = NULL;

// *********************************************************************
// Get the directory name of the application
// *********************************************************************

char *Image::GetAppDir() {
#ifdef ENVVAR_ECAVDIR
	return ENVVAR_ECAVDIR;
#else
	char *dir=getenv("ECAVDIR");
		
    if (dir==NULL) {
     	cerr << "Please specify the ECAV main directory in the environment \n"
     		 << "variable ECAVDIR!\n"
     		 << "E.g. export ECAVDIR=/home/toto/ECAV\n";
     	exit (1);
    }
    return dir;
#endif
}

// *********************************************************************
// The library contains a digits image which is used to draw text
// into an image. It is loaded on demand.
// *********************************************************************

Image * Image::GetDigits () {
	
	// Not yet loaded
    if (_Digits==NULL) {
		char fname[1024];
		char *dir = GetAppDir();
		
        sprintf (fname, "%s/data/digits.ppm", dir);
		_Digits = new Image (fname);		
	}
	
	return _Digits;
}

Image * Image::GetDigitsBig () {
	
	// Not yet loaded
    if (_DigitsBig==NULL) {
		char fname[1024];
		char *dir = GetAppDir();

        sprintf (fname, "%s/data/digits_big.ppm", dir);
		_DigitsBig = new Image (fname);				
	}	
	return _DigitsBig;
}

// *************************************************************
// Determines if a stored image is a color Image
// *************************************************************

bool Image::isColorImage (const char *filename) 
{

	char shortbuf[256];
	FILE * inpic;
	
	char c;

	if ((inpic	= fopen(filename,"r")) == NULL)	
		ERR_THROW ("Cannot open file " << filename << ": " << strerror(errno));
		
	// Read the type flag
	fgets(shortbuf,MAX_ANZAHL,inpic);	
	c=shortbuf[1];	
	fclose (inpic);
	
	return (c == '6');
}	

#ifndef TELESUN_INSA

// *************************************************************
// Creates a vector of pages out of a list of tiles.
// Stacks the tiles vertically, but begins each line with
// marker images which are recognized by a OCR Program
// *************************************************************
void Image::CreateOCRPages (const char *filename,
						 	char *name_template,
						 	bool order,
						 	int border,
						 	byte bordervalue,
						 	bool doTile) {
	istream *st=NULL;
	ifstream *fst=NULL;
	int full_type = 2;
	int run_x, run_y, count_x;
	int no_tiles=0;
	int tile_nr;
	int no_page;
	bool lineIsMarker;
	Image *marker_template;
	Image *cur_page=NULL;
	int marker_x, marker_y;
	int digits[3];
	int xoff;

	// Different data structures containing the tiles
	multiset<Image *, image_ymore> images_yo;
	multiset<Image *, image_ymore>::iterator yiter;
	vector<Image *> images;

	// Load the marker and seperator tokens
	char *ecavdir = GetAppDir();
	ostringstream str1;
	str1 << ecavdir << "/data/token.ppm";
	marker_template = new Image (str1.str().c_str());

	// Stdin
	if (strcmp(filename,"-")==0)
		st = &cin;

	// A regular file
	else {
    	fst = new ifstream (filename,  ios::in);
    	if (!st->good()) {
    		cerr << "Cannot open file " << filename << "for read!!\n";
    		exit (1);
    	}
    	st = fst;
	}

	// read the input stream and read all the images.
	// put them into an array.
	while (!st->eof()) {
		string s;

		*st >> s;

		if (s.length()>0) {
			Image *im = new Image (s.c_str());
			if (order)
				images_yo.insert(im);
			else
				images.push_back(im);
			++no_tiles;
			if (im->type==3)
				full_type = 3;
		}

		if (st->eof())
			break;
		if (st->fail()) {
			cerr << "*** ERROR: Syntax error in the data!\n";
			exit (1);
		}
	}

	// -----------------------------------------------------------
	// Travers the tiles
	// -----------------------------------------------------------

	run_y = 0;
	no_page = 0;
	tile_nr = 0;
	lineIsMarker = false;
	yiter = images_yo.begin();
	while (tile_nr < no_tiles) {
		Image *im=NULL;		

        // ----------------------------------------------------------
        // Access the tile/marker
        
        // This line is a marker line, create the marker
		if (lineIsMarker) {
			digits[0] = tile_nr%10;
    		digits[1] = (tile_nr/10)%10;
    		digits[2] = (tile_nr/100)%10;
    		marker_y = marker_template->ysize;
    		marker_x = marker_template->xsize + 60 +
    				   Image::digitwidthsBig[digits[2]] +
    				   Image::digitwidthsBig[digits[1]] +
    				   Image::digitwidthsBig[digits[0]];
    		im = new Image (marker_x, marker_y, 3);
    		im->setPlaneToValue(PLANE_RED, bordervalue);
    		im->setPlaneToValue(PLANE_GREEN, bordervalue);
    		im->setPlaneToValue(PLANE_BLUE, bordervalue);
    		im->paste (*marker_template,0,0);
    		xoff = marker_template->xsize + 30;
    		for (int d=2; d>=0; --d) {
    			im->drawText (xoff, 0, digits[d], true);
    			xoff += Image::digitwidthsBig[digits[d]] + 10;
    		}
		}
		else {
			// Access the new tile
			if (order) {
        		im = *yiter;
        		++yiter;
    		}
    		else
    			im = images[tile_nr];
    		++tile_nr;
		}

		// ----------------------------------------------------------
		// The page is full, create a new one

		
		if ((no_page==0) || (run_y + im->ysize > OCR_PAGE_HEIGHT)) {

			// Save finished page in the return set of pages
			if (no_page!=0) {
				ostringstream str;
				str << name_template << setfill ('0') << setw(2) << no_page << ".ppm";
				cur_page->write (str.str().c_str());
				delete cur_page;
			}

			// Allocate and init the new page
			cur_page = new Image (OCR_PAGE_WIDTH, OCR_PAGE_HEIGHT, full_type);
        	cur_page->setPlaneToValue (PLANE_RED, bordervalue);
        	if (full_type==3) {
        		cur_page->setPlaneToValue (PLANE_GREEN, bordervalue);
        		cur_page->setPlaneToValue (PLANE_BLUE, bordervalue);
        	}

			// Set counters to zero
			run_y = 0;
        	++no_page;

        	cout << "# Page " << no_page << ":\n" << "P " << no_page << endl;
        }

		// Paste the tile/marker into the line
		run_x = count_x = 0;
		while (run_x + im->xsize <= OCR_PAGE_WIDTH) {
			cur_page->paste (*im, run_x, run_y);
			run_x += im->xsize + border;
			++count_x;

			// If we do not tile, i.e. copy only one tile on the page, then exit
			if (!doTile)
				break;
		}

		run_y += im->ysize + border;

		cout << (lineIsMarker ? "M " : "T ") << tile_nr << " " << count_x << endl;
		lineIsMarker = ( lineIsMarker ? false : true );

		delete im;
	}

	// Save the last page
	ostringstream str;
	str << name_template << setfill ('0') << setw(2) << no_page << ".ppm";
	cur_page->write (str.str().c_str());
	// delete cur_page;

	// Clean up
	if (fst!=NULL) {
		fst->close();
		delete fst;
	}
	
	// delete marker_template;
}

#endif
