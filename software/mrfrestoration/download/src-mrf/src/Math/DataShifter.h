/******************************************************************************
 * DataShifter
 *
 * As a data serie, but uses incremental calculation
 * Not all statistics can be calculated, but the available
 * ones are optimized for fast re-calculation after a single
 * value has been added and another one at the beginning has been removed.
 *
 * Begin: Christian Wolf chriswolf@gmx.at,
 * March. 2002
 ******************************************************************************/

#ifndef _WOLF_DataShifterINCREMENTAL_H_
#define _WOLF_DataShifterINCREMENTAL_H_

// C
#include <math.h>

// From the main module
#include <CIL.h>

using namespace std;

class DataShifter {

	friend inline ostream & operator << (ostream & os, DataShifter &d);

	public:

		// Constructors
		inline DataShifter ()		{ clear(); }

		// Add and/or remove value with immediate recalculation of
		// the statistics
		inline void add (float v);
		inline void add_sub (float a);

		// Access functions
		inline float operator[] (int i) const;
		inline int size() const { return values.size(); }
		inline void clear();

		// get the stats
		inline float getMin()		{ return min; }
		inline float getMax()		{ return max; }
		inline float getSum()		{ return sum; }
		inline float getMean();
		inline float getVariance();

	private:

		// DATA
		vector<float> values;
		unsigned int index_afterend;
		float sum, sum_sq, min, max;
};

// *******************************************************************************
// The access method
// *******************************************************************************

inline float DataShifter::operator[] (int i) const {
#ifdef PEDANTIC_CHECK_CODE
	if (i<0 || i>=(int)values.size()) {
		ERR_THROW ("Out of bounds in DataShifter::operator[]() ["  << values.size()
			 << "]: i=" << i << endl);
	}
#endif
	return values[i];
}

// *******************************************************************************
// Clear - remove all values
// *******************************************************************************

inline void DataShifter::clear ()	{
	values.clear();
	index_afterend=0;
	sum=sum_sq=min=max=0;
}

// *******************************************************************************
// Add a value with recalculation of the statistics
// *******************************************************************************

inline void DataShifter::add (float v) {
	if (v<min || values.size()==0) min=v;
	if (v>max || values.size()==0) max=v;
	values.push_back (v);
	++index_afterend;
	sum += v;
	sum_sq += v*v;
}

// *******************************************************************************
// Add a value and remove another value with recalculation of the statistics
// *******************************************************************************

inline void DataShifter::add_sub (float a) {
	float r;

	if (index_afterend < values.size()) {
#ifdef PEDANTIC_CHECK_CODE
    	if (index_afterend<0 || index_afterend>=values.size()) {
    		ERR_THROW ("Out of bounds in DataShifter::add_sub() ["  << values.size()
    			 << "]: index_afterend=" << index_afterend << endl);
    	}
#endif
		r=values[index_afterend];
		values[index_afterend]=a;
		++index_afterend;
	}
	else {
		r=values[0];
		values[0]=a;
		index_afterend=0;
	}

	//
	sum -= r;
	sum += a;
	sum_sq -= r*r;
	sum_sq += a*a;

	if (a<min)
		min=a;
	else {
		// we removed the minimum -> recalculate the new minimum
    	if (r<=min) {
     		min=values[0];
       		for (unsigned int i=1; i<values.size(); ++i) {
				float v=values[i];
				if (v<min) min=v;
			}
     	}
	}

	if (a>max)
		max=a;
	else {
		// we removed the minimum -> recalculate the new minimum
    	if (r>=max) {
     		max=values[0];
       		for (unsigned int i=1; i<values.size(); ++i) {
				float v=values[i];
				if (v>max) max=v;
			}
     	}
	}
}

// *******************************************************************************
// get the mean value
// *******************************************************************************

inline float DataShifter::getMean() {
	return sum / (float) values.size();
}

// *******************************************************************************
// get the variance
// *******************************************************************************

inline float DataShifter::getVariance() {
	float s=(float) values.size();
	return (sum_sq - (sum*sum)/s)/s;
}

// *******************************************************************************
// Output to a stream
// *******************************************************************************

inline ostream & operator << (ostream & os, DataShifter &d) {
	for (int i=0; i<d.size(); ++i)
		os << d.values[i] << " ";
	os << endl;
	return os;
}

#endif
