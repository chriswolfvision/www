/***************************************************************************
                          DataDiffAccumulator.h  -  description
                             -------------------
    begin                : Sat Apr 12 2003
    copyright            : (C) 2003 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 ***************************************************************************/

#ifndef _WOLF_DataDiffAccumulator_H_
#define _WOLF_DataDiffAccumulator_H_

// C
#include <math.h>

using namespace std;

class DataDiffAccumulator {

	public:

		// Constructors
		inline DataDiffAccumulator ()		{ clear(); }

		// Add and/or remove value with immediate recalculation of
		// the statistics
		inline void add (float v);
		inline void clear();

		// Access functions		
		inline float get();
		
	private:

		// DATA
		float sum, last_value;
		unsigned int count;
};


// *******************************************************************************
// Clear - remove all values
// *******************************************************************************

inline void DataDiffAccumulator::clear ()	{
	sum = last_value = 0;
	count = 0;
}

// *******************************************************************************
// Add a value with recalculation of the statistics
// *******************************************************************************

inline void DataDiffAccumulator::add (float v) {
	if (count>0)
		sum += fabs(last_value-v);
	++count;
	last_value = v;	
}

// *******************************************************************************
// Add a value with recalculation of the statistics
// *******************************************************************************

inline float DataDiffAccumulator::get () {
	if (count<2)
		return 0;
	else	
		return sum / (float) (count-1);
}

#endif
