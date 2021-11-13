/**************************************************************************
                          RandomNumberGenerator.h  -  description
                             -------------------
    begin                : Thu Aug 23 2001
    copyright            : (C) 2001 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 ***************************************************************************/

#ifndef _WOLF_RANDOMNUMERGENERATOR_H_
#define _WOLF_RANDOMNUMERGENERATOR_H_

// C
#include <math.h>

// C++
#include <iostream>

// From the main library
#include <CIL.h>

#define MAX_RI				256*256*256*256
#define STATEARRAYLENGTH	256

class RandomNumberGenerator {

	public:
	
		RandomNumberGenerator (bool useRandomDevice, double min, double max, int buffer_size);
		~RandomNumberGenerator ();
		
		double next();
		double next16();
		unsigned char nextByte();
		void getGaussian16(double &z1, double &z2);

    private:

    	void fillBuf();

    private: 	// DATA

    	unsigned int buffersize;
    	unsigned int bufpointer;

    	// Used for the kernel device generator
    	unsigned char *kbuf;
     	
    	// Used for the C library generator
    	double *cbuf;
    	
    	double max, min, delta;
    	
    	// File descriptor for the random device
    	bool usedevice;
    	int fd;	
    	
    	char statearray[STATEARRAYLENGTH];  	    	
};


// **********************************************************************
// Get the next random number
// A double value calculated from 32 bit integer precision
// **********************************************************************

inline double RandomNumberGenerator::next() 
{	
	if (usedevice) 
	{		
		// throw EError ("RandomNumberGenerator::Next(): Not implemented!\n");
		
		if (bufpointer+4>=buffersize)
			fillBuf();

    	unsigned char *cp = kbuf+bufpointer;
    	double ri = cp[0] +
    		 cp[1] * 256. +
    		 cp[2] * 256.*256. +
    		 cp[3] * 256.*256.*256.;
    	bufpointer+=4;

    	return  (ri / (256.*256.*256.*256.-1.)) * delta + min;
	}

	else 
	{
		if (bufpointer+sizeof(double)>=buffersize)
			fillBuf();
	
		double rv = cbuf[bufpointer] / (double) RAND_MAX;
		rv = rv*delta+min;
		++bufpointer;
		return rv;
	}
}

// **********************************************************************
// Get the next random number
// A double value calculated from 16 bit integer precision
// **********************************************************************

inline double RandomNumberGenerator::next16() 
{	
	double rv;
	if (usedevice) 
	{		
		if (bufpointer+2>=buffersize)
			fillBuf();

    	unsigned char *cp = kbuf+bufpointer;
    	double ri = cp[0] + cp[1] * 256.;
    	bufpointer+=2;
    	rv = (ri / (256.*256.-1.)) * delta + min;
#ifdef PEDANTIC_CHECK_CODE
		if (!finite(rv))
			ERR_THROW ("Numerical error in RandomNumberGenerator::next16: "
				<< " ri=" << ri
				<< " delta=" << delta
				<< " min=" << min);
#endif    	
    	return  rv;
	}

	else 
	{
		ERR_THROW ("next16() not yet implemented for this device!");
	}
}

// **********************************************************************
// Get 2 samples from a zero mean unit veriance Gaussian distribution
// A double value calculated from 16 bit integer precision
// **********************************************************************

inline void RandomNumberGenerator::getGaussian16(
	double &z1, double &z2)
{ 
	double r1,r2,u1,u2;
	do { r1 = next16(); } while (r1<=0);
	do { r2 = next16(); } while (r2<=0);
	u1 = -2.*log(r1);
	u2 = -2.*log(r2);
	z1 = sqrt(u1)*cos(2.*M_PI*u2);
	z2 = sqrt(u2)*sin(2.*M_PI*u2);
#ifdef PEDANTIC_CHECK_CODE
		if (!(finite(z1) && finite(z2)))
			ERR_THROW ("Numerical error in RandomNumberGenerator::getGaussian: "
				<< " z1=" << z1
				<< " z2=" << z2
				<< " u1=" << u1
				<< " u2=" << u2
				<< " r1=" << r1
				<< " r2=" << r2
				);
#endif  
 }

// **********************************************************************
// Get the next random byte
// **********************************************************************

inline unsigned char RandomNumberGenerator::nextByte() 
{
#ifdef UNIX_ARCHITECTURE
	unsigned char rv;

	if (bufpointer+sizeof(float)>=buffersize)
		fillBuf();

	if (usedevice) 
	{

		rv = kbuf[bufpointer];
		++bufpointer;
		return rv;
	}
	else 
		ERR_THROW ("RandomNumberGenerator::NextByte(): Not implemented!\n");
#else
	ERR_THROW ("Internal error: RandumNumberGenerator may only be used on UNIX architectures!\n");
#endif
}

#endif


