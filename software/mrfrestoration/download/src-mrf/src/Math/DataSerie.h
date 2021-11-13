/******************************************************************************
 * DataSerie
 *
 * Stores a data serie, i.e. a serie of double values.
 * It can be used to calculate statistics on these values (mean, median,
 * variance etc.), or to do some correlation measures.
 *
 * Begin: Christian Wolf chriswolf@gmx.at, Dec. 2000
 ******************************************************************************/

#ifndef _WOLF_DATASERIE_H_
#define _WOLF_DATASERIE_H_

// C
#include <math.h>
#include <string.h>

// C++
#include <iostream>
#include <fstream>
#include <vector>

// From the main module
#include <CIL.h>

using namespace std;

#include <visualbug.h>

class FilterMask;

class DataSerie {

	public:
	
		// Constructors
		inline DataSerie ();	
		inline DataSerie (istream &s);
		inline DataSerie (char *filename);
		inline DataSerie (int lborder, int rborder);

		// The iterator		
		typedef	vector<double>::iterator iterator;
		typedef	vector<double>::const_iterator const_iterator;
		iterator begin() 						{ return values.begin(); }
		iterator end()							{ return values.end(); }
		const_iterator begin() const 			{ return values.begin(); }
		const_iterator end()   const 			{ return values.end(); }
			
		// Output of the data
		inline void write (const char *filename);
	
		// Add a value
		inline void add (double v);
		
		// Access functions
		inline double operator[] (int i) const 	{ return values[i];}
		inline void set (int i, double val) 	{ values[i] = val; }
		inline int size() const 				{ return values.size(); }
		inline void setZero() 					{ for(int i=0; i<size(); ++i) values[i]=0;}
		inline void clear();

		inline void resize(int n) { values.resize(n); haveStats=haveMedian=false; }
		
		// get the stats
		inline double getSum();
		inline double getMean();
		inline double getVariance();
		inline double getMin();
		inline double getMax();
		inline double getMedian();
		inline double getMeanOfNonZero ();
		
		DataSerie *FFTMagnitude (DataSerie &im_serie);
		DataSerie *FFTPhase     (DataSerie &im_serie);
		void reOrderAfterFFT ();
		
		// Filter the data with a mask (convolution)
		void normalize();
		void filter (DataSerie &fm);
		
		// Trim non-zero valued borders
		void trim ();

		// The earth mover's distance
		double EMD (DataSerie &other);
		double L2  (DataSerie &other);
		
		// Get the number of peaks
		int getPeaksAndCuts 		(double k,          vector<int> &out_couts);
		int getPeaksAndCutsWindowed (double k, int win, vector<int> &out_couts);
		
		DataSerie * correlation (DataSerie &other, int maxsigma) const;
		
		void debugOut (ostream & os);
		
		friend ostream & operator << (ostream & os, DataSerie &d);
		
		// HELP FUNCTIONS	
	
		void calcStats();
		void calcMedian();
		void _alloc (istream &s);
	
		// DATA
	
		vector<double> values;
			
		double sum, mean, variance, min, max, median;
		
		int left_border, right_border;
		bool haveStats, haveMedian;

};

// *******************************************************************************
// Constructor
// *******************************************************************************    	

inline DataSerie::DataSerie () {
	haveStats = haveMedian = false;	
	left_border = right_border = 0;
}

// *******************************************************************************
// Constructor
// *******************************************************************************    	

inline DataSerie::DataSerie (int lborder, int rborder) {
	haveStats = haveMedian = false;	
	left_border = lborder;
	right_border = rborder;
}

// *******************************************************************************
// Constructor
// *******************************************************************************    	

inline DataSerie::DataSerie (char *filename) {

	// Stdin
	if (strcmp(filename,"-")==0)
		_alloc (cin);
		
	// A regular file
	else {
    	ifstream st (filename,  ios::in);
    	if (!st.good()) {
    		ERR_THROW ("Cannot open file " << filename << "for read!!\n");
    	}
    	_alloc (st);	
    	st.close();	
	}
}

// *******************************************************************************
// Output of the data
// *******************************************************************************    	

inline void DataSerie::write (const char *filename) {

	// Stdout
	if (strcmp(filename,"-")==0)
		cout << *this;
		
	// A regular file
	else {
    	ofstream st (filename,  ios::out);
    	if (!st.good()) {
    		ERR_THROW ("Cannot open file " << filename << "for write!!\n");
    	}
    	st << *this;
    	st.close();	
	}
}

// *******************************************************************************
// Constructor
// *******************************************************************************    	

inline DataSerie::DataSerie (istream &st) {
	_alloc (st);	
}

// *******************************************************************************
// Clear - remove all values
// *******************************************************************************    	

inline void DataSerie::clear ()	{
	values.clear();
	haveStats = false;
	haveMedian = false;
}

// *******************************************************************************
// Add a value
// *******************************************************************************    	

inline void DataSerie::add (double v) {
	values.push_back (v);
	haveStats = haveMedian = false;
}

// *******************************************************************************
// Get the mean value
// *******************************************************************************    	

inline double DataSerie::getMean() {
	if (!haveStats)
		calcStats();
	return mean;	
}

// *******************************************************************************
// Get the variance
// *******************************************************************************    	

inline double DataSerie::getVariance() {
	if (!haveStats)
		calcStats();
	return variance;	
}

// *******************************************************************************
// Get the median
// *******************************************************************************    	

inline double DataSerie::getMedian() {
	if (!haveMedian)
		calcMedian();
	return median;	
}

// *******************************************************************************
// Get the minimum
// *******************************************************************************    	

inline double DataSerie::getMin() {
	if (!haveStats)
		calcStats();
	return min;	
}

// *******************************************************************************
// Get the minimum
// *******************************************************************************    	

inline double DataSerie::getMax() {
	if (!haveStats)
		calcStats();
	return max;	
}

// *******************************************************************************
// Get the sum
// *******************************************************************************    	

inline double DataSerie::getSum() {
	if (!haveStats)
		calcStats();
	return sum;	
}   	  	

// ***************************************************************
// Compute the mean	value but just of
// the values which	are	non	zero.
// ***************************************************************

double DataSerie::getMeanOfNonZero ()	{
	int	vcount = 0;
	double mean=0;
	int s=size();

	for	(int i=0; i<s; ++i) {
		if (values[i]>0) {
			mean +=	values[i];
			++vcount;
		}
	}
	return mean	/ (double) vcount;
}

#endif

