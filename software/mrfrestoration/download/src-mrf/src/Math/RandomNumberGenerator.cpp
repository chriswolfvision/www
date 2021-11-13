/***************************************************************************
                          RandomNumberGenerator.cc  -  description
                             -------------------
    begin                : Thu Aug 23 2001
    copyright            : (C) 2001 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 ***************************************************************************/

// C
#include <stdlib.h>
#include <stdio.h>
#ifdef UNIX_ARCHITECTURE
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>	
#include <unistd.h>
#include <time.h>
#endif

// From the main module
#include <CIL.h>

// Own classes
#include "RandomNumberGenerator.h"

using namespace std;

// **********************************************************************
// Constructor
// **********************************************************************

RandomNumberGenerator::RandomNumberGenerator (bool useDev, double mi, double ma, int bufsize) 
{
	min = mi;
	max = ma;
	delta = max-min;
	buffersize = bufsize;
	usedevice = useDev;	
	bufpointer=buffersize+1;

#ifdef UNIX_ARCHITECTURE	

	if (useDev) 
	{
		kbuf = new unsigned char [buffersize];
			
		// Open	the	random number device
        if ((fd = open ("/dev/urandom", O_RDONLY))<0) 
        	ERR_THROW ("Internal error in RandomNumberGenerator::RandomNumberGenerator()!\n"
        		"Could not open	random number device!\n");
	}
	
	// Initialize the C Library random number system	
	else 
	{
		ERR_THROW ("RandumNumberGenerator: C library not supported for the moment (unfound bug)");
	
		unsigned int seed = time(NULL);
		cbuf = new double [buffersize];
		// srand(seed);
		initstate(seed, statearray, STATEARRAYLENGTH);
		srandom(seed);
	}	

#else
	ERR_THROW ("Internal error: RandumNumberGenerator may only be used on UNIX architectures!\n");
#endif

}

// **********************************************************************
// Destructor
// **********************************************************************

RandomNumberGenerator::~RandomNumberGenerator () 
{
#ifdef UNIX_ARCHITECTURE
	if (usedevice) 
	{
		close (fd);
		delete [] kbuf;
	}
	else
		delete [] cbuf;
#else
	ERR_THROW ("Internal error: RandumNumberGenerator may only be used on UNIX architectures!\n");
#endif

}

// **********************************************************************
// Fill the buffer with random bytes from the
// random number kernel device
// **********************************************************************

void RandomNumberGenerator::fillBuf () 
{

#ifdef UNIX_ARCHITECTURE

	// Fill the buffer using the kernel device
	if (usedevice) 
	{

    	if (read (fd, kbuf, buffersize)<0)	
    		ERR_THROW ("Internal error in RandomNumberGenerator::RandomNumberGenerator()!\n"
    			 "Could not read from random number device!\n");
	}
	
	// Fill the buffer using the C Library functions
	else 
	{
		for (unsigned int i=0; i<buffersize; ++i) 
		{
			cbuf[i] = random();
		}
				
	}	
	bufpointer=0;
#else
	ERR_THROW ("Internal error: RandumNumberGenerator may only be used on UNIX architectures!\n");
#endif

}	

