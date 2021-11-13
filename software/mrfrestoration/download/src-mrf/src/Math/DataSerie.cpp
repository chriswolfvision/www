// C
#include <math.h>

// C++
#include <algorithm>
#include <iostream>

// From the main module
#include <CIL.h>

// Own Classes
#include "DataSerie.h"

#include "visualbug.h"
using namespace std;

// *******************************************************************************
// Help function for the constructor - read values from a stream
// *******************************************************************************    	

void DataSerie::_alloc (istream &st) {
	double val;

	while (!st.eof()) {
		
		st >> val;
		
		if (st.eof())
			break;
		if (st.fail())
			ERR_THROW ("*** ERROR: Syntax error in the data!\n");
			
		add(val);
	}
}

// *******************************************************************************
// Calculate the magnitude of two series containing the real
// and imaginary FF transformed data
// *******************************************************************************    	

DataSerie *DataSerie::FFTMagnitude (DataSerie &im_serie) {
	DataSerie *rv = new DataSerie;
	int s=size();
	for	(int x=0; x<s; ++x)	{
		float re = (*this)[x];
		float im = im_serie[x];
		rv->add(sqrt(re*re+im*im));
	}
	return rv;
}

// *******************************************************************************
// Calculate the phase of two series containing the real
// and imaginary FF transformed data
// *******************************************************************************    	

DataSerie *DataSerie::FFTPhase (DataSerie &im_serie) {
	DataSerie *rv = new DataSerie;
	int s=size();
	for	(int x=0; x<s; ++x)	{
		float re = (*this)[x];
		float im = im_serie[x];
		if (re==0)
			rv->add(0);
		else
			rv->add(atan(im/re));
	}
	return rv;
}

/****************************************************************
 * Order the coefficents after a FFT transform
 ****************************************************************/

void DataSerie::reOrderAfterFFT ()	{
	int s = size();
	int	offset = s / 2;
	DataSerie copy	(*this);

	for	(int x=0; x<s/2; ++x) {
		values[x] =			copy[x+offset];
		values[x+offset] =	copy[x];
	}
}

// *******************************************************************************
// Calculate the statistics
// *******************************************************************************    	

void DataSerie::calcStats () 
{
	int size, sizerange;
	int lbcount, rbcount;
	double sum_sq, v;
	bool first;
	vector<double>::iterator itleft, itright;
		
	// The borders of the data
	size = values.size();
	lbcount = (size * left_border) / 100;
	rbcount = (size * right_border) / 100;	
	sizerange = size - lbcount - rbcount;
			
	// We have elements
	if (sizerange>0) 
	{
		
		itleft = values.begin();
		itleft += lbcount;
		itright = values.end();
		itright -= rbcount;							
									
		// Travers the data
		sum = sum_sq = 0;
		first = true;

		for (int i=lbcount; i<size-rbcount; ++i) 
		{
			v = values[i];
			sum    += v;
			sum_sq += v*v;
			if (first) 
			{				
				min = v;
				max = v;				
				first = false;
			}
			else 
			{
				if (v < min) min=v;
				if (v > max) max=v;				
			}			
		}		
		
		mean = (sum / (double) sizerange);
		variance = (sum_sq - (sum*sum)/sizerange) / (double) sizerange;
    }

    // No elements
    else 
    {
#ifdef CHECK_CODE  	
		if (size>0)
			ERR_THROW ("Internal error in DataSerie::calcStats()\n");
#endif
    	min = max = mean = variance = 0;

    }
	
	haveStats = true;
}

// *******************************************************************************
// Calculate the Median
// *******************************************************************************    	

void DataSerie::calcMedian () {
	int size, sizerange;
	int lbcount, rbcount;
	vector<double>::iterator itleft, itright;
	
	// Create a copy of the data
	vector<double> copy (values);
		
	// The borders of the data
	size = copy.size();
	lbcount = (size * left_border) / 100;
	rbcount = (size * right_border) / 100;	
	sizerange = size - lbcount - rbcount;
			
	// We have elements
	if (sizerange>0) {
		
		// Sort the vector in the valid field - needed for the median
		itleft = copy.begin();
		itleft += lbcount;
		itright = copy.end();
		itright -= rbcount;							
		sort (itleft, itright);
									
		median = copy[(lbcount + size - rbcount)/2];	
    }

    // No elements
    else
    	median = 0;
	
	haveMedian = true;
}

// *******************************************************************************
// Trim non-zero valued borders
// *******************************************************************************    	

void DataSerie::trim () {
	int s=size();
	int left, right;
	DataSerie new_serie;
	
	// Search the borders
	for (left=0; left<s; ++left)
		if (values[left]!=0)
			break;
	for (right=s-1; right>=0; --right)
		if (values[right]!=0)
			break;
	
	for (int i=left; i<=right; ++i)
		new_serie.add (values[i]);
		
	*this = new_serie;
}

// *******************************************************************************
// Calculate the correlation between two sets of data
// *******************************************************************************    	

DataSerie * DataSerie::correlation (DataSerie &other, int maxsigma) const {		
	DataSerie *corrserie;
	int st, so, corrlength;
	double corr, sumt, sumo, norm;
	
	corrserie = new DataSerie;
	
	// Calculate borders and overlap area
	st = this->size();
	so = other.size();	
	corrlength = st-maxsigma < so ? st-maxsigma : so;
	
	if (corrlength<1) {
		ERR_THROW ("**** ERROR in DataSerie::correlation()\n"
			 << "     maxsigma too big for the length of this data!\n"
			 <<	"	  corrlength: " << corrlength << endl);
	}
	
	// Travers the sigmas = shift between the two series
	for (int sig=0; sig<maxsigma; ++sig) {

		// Travers the overlap area of the two series and compare	
		corr=sumt=sumo=0;						
		for (int i=0; i<corrlength; ++i) {
			double t = (*this)[i];
			double o = other[i];
			
			corr += (*this)[i+sig] * o;			
			sumt += t*t;
			sumo += o*o;
		}
		corr /= (double) corrlength;
		norm = sqrt (sumt * sumo) / (double) corrlength;
		corrserie->add (corr / norm);
	}	
		
	return corrserie;			
}

/************************************************************************
 * Normalize the filter
 ************************************************************************/

void DataSerie::normalize(){
	double sum = 0;

	for	(unsigned int x=0; x<values.size(); ++x)	{
		sum	+= values[x];
	}

	for	(unsigned int x=0; x<values.size(); ++x)	{
		values[x] = values[x]/sum;
	}
}

// *************************************************************
// Filter the data with a given filter mask.
// The mask must of course be 1D
// *************************************************************

void DataSerie::filter (DataSerie &fm) 
{
	vector<double> *FIL;
	int	ix, bx, fx;
	double newVal;
	int s=size();
	int fmx = fm.size();
	
	// Make a copy of the data
	FIL = new vector<double> (values);

	// Treat the borders
	for (int i=0; i<fmx/2; ++i) 
	{
		(*FIL)[i] = values[i];
		(*FIL)[s-1-i] = values[s-1-i];
	}

	// Travers the data
	for	(int x=0+fmx/2; x<(s-fmx/2); ++x ) 
	{
		// the left border
		bx = x-(fmx/2);

    	// convolve
    	newVal = 0;
        for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx)
			newVal += (((double)values[ix])*fm[fx]);

		(*FIL)[x]=	(int) rint (newVal);
	}

	// Copy result and clean up
	values = *FIL;
	haveStats = haveMedian = false;
	delete FIL;
}

// *******************************************************************************
// Search the peaks in the data - Global version
// - using a method which ressembles Niblack's adaptive thresholding method
//
// Returns the number of peaks
// Stores in the vector "out_cuts", the points where the data has to be cut.
// The first cut to the left of the first peak will be added as well as the last
// cut on the right of the last peak.
// I.e. If one peak has been found, the cut vector has 2 entries (to trim the
// boxes ...)
// *******************************************************************************

int DataSerie::getPeaksAndCuts (double k, vector<int> &out_cuts) {
	int xsize=size();
	double sum, m, s, th;
	double R;
	int curval, oldval, peaks, posofmin, posoflastpeak=0;
	int begin_x, end_x;
	double min;	
	
	// No data
	if (xsize<1)
		return 0;
	
	// Compute the dynamics of the data
	R = getMax() - getMin();
	
	// Search beginning and end of the data
	begin_x=0;
	for (int i=0; i<xsize; ++i)
		if (values[i]>0) {
			begin_x=i;
			break;
		}
	end_x=xsize-1;
	for (int i=xsize-1; i>=0; --i)
		if (values[i]>0) {
			end_x=i;
			break;
		}
			

   	// Calculate the statistics and the threshold
   	sum = 0;
   	for	(int w=begin_x; w<=end_x; w++)
   		sum += values[w];
   	m  = sum / (double) (end_x-begin_x+1);
   	s = 0;
   	for	(int w=begin_x; w<=end_x; w++) {
		double d  = m - values[w];
		s += d*d;
	}
	s = sqrt (s/(double)(end_x-begin_x+1));
	th = m * (1.0 + k*(s/R - 1.0));
   	
   	// Travers the data, threshold it and count the number
   	// of ranges of 1s.
   	oldval = 0;
   	peaks = 0;
   	min = values[0];
   	posofmin = 0;
	for (int i=0; i<xsize; ++i) {
		curval = (values[i] >= th ? 1 : 0);		
		if (curval && !oldval) {
			++peaks;					
			posoflastpeak = i;
			out_cuts.push_back(posofmin);						
			min = values[i];
		}			
		
		if (values[i] <= min) {
			min = values[i];	
			posofmin = i;
		}	
		oldval = curval;
	}   		
	
	// Search the last minimum on the right border and add the position
	// to the cuts
	min = values[xsize-1];
	posofmin = xsize-1;
	for (int i=xsize-1; i>posoflastpeak; --i) {
		if (values[i] <= min) {
			min = values[i];	
			posofmin = i;
		}			
	}
	out_cuts.push_back(posofmin);
	
	cerr << "th: " << th << endl;	   			
	
   	return peaks;
}

// *******************************************************************************
// Search the peaks in the data - Windowed version
// - using a method which ressembles Niblack's adaptive thresholding method
//
// Returns the number of peaks
// Stores in the vector "out_cuts", the points where the data has to be cut.
// The first cut to the left of the first peak will be added as well as the last
// cut on the right of the last peak.
// I.e. If one peak has been found, the cut vector has 2 entries (to trim the
// boxes ...)
// *******************************************************************************

int DataSerie::getPeaksAndCutsWindowed (double k, int win, vector<int> &out_cuts) {
	int xsize=size();
	int wh=win/2;
	DataSerie *th_serie;
	double sum, m, s, th;
	double R;
	int curval, oldval, peaks, posofmin, posoflastpeak=0;
	double min;	
	
	// No data
	if (xsize<1)
		return 0;
	
	// Allocate the threshold serie
	th_serie = new DataSerie;
	
	// Compute the dynamics of the data
	R = getMax() - getMin();
	
	// Functor, gets the standard deviation of the window
    // --------------------------------------------------
    	
    class _getsd { public:
    	double operator()(double mean, int left) {
    		double d, ls2 = 0;
    		for	(int x=0 ; x<lwin; x++) {
    			d  = mean - (*lvalues)[left+x];
    			ls2 += d*d;
    		}
    		return sqrt (ls2/(double)lwin);
    	}
    		_getsd(vector<double> *xv, int xwin)
    		: lvalues (xv), lwin(xwin) {};
    	vector<double> *lvalues;
    	int lwin;		
    } getsd (&values, win);
    	
    	
    // Functor, calculates the treshold
    // --------------------------------------------------
    	
    class _thresh { public:
    	double operator()(double mean, double sd) {
    		return mean * (1.0 + lk*(sd/lR - 1.0));
    	}
    	_thresh(double xk, double xR) : lk(xk), lR(xR){};
    	double lk, lR;
    } thresh (k, R);
			
	// Travers all lines
	// -----------------------------

   	// Calculate the initial window at the beginning of the line
   	sum = 0;
   	for	(int w=0 ; w<win; w++)
   		sum += values[w];
   	m  = sum / (double) (win);
   	s  = getsd (m, 0);
   	th = thresh(m, s);
   	
   	// Left border treatment
   	for (int i=0; i<wh+1; ++i)
   		th_serie->add(th);
   	   	
   	// Shift the window, add and remove	new/old values to the histogram
   	for	(int i=1 ; i <= xsize-win; i++) {
		
   		// Remove the left old column and add the right new column
  		sum -= values[i-1];
   		sum += values[i+win-1];
   		m  = sum / (double) win;
   		s  = getsd (m, i);
   		th = thresh(m, s);
   		th_serie->add (th);						
   	}
   	
   	// Right border treatment
   	for (int i=0; i<wh; ++i)
   		th_serie->add(th);   		
   	
   	// Travers the data, threshold it and count the number
   	// of ranges of 1s.
   	oldval = 0;
   	peaks = 0;
   	min = values[0];
   	posofmin = 0;
	for (int i=0; i<xsize; ++i) {
		curval = (values[i] >= (*th_serie)[i] ? 1 : 0);		
		if (curval && !oldval) {
			++peaks;					
			posoflastpeak = i;
			out_cuts.push_back(posofmin);						
			min = values[i];
		}			
		
		if (values[i] <= min) {
			min = values[i];	
			posofmin = i;
		}	
		oldval = curval;
	}   		
	
	// Search the last minimum on the right border and add the position
	// to the cuts
	min = values[xsize-1];
	posofmin = xsize-1;
	for (int i=xsize-1; i>posoflastpeak; --i) {
		if (values[i] <= min) {
			min = values[i];	
			posofmin = i;
		}			
	}
	out_cuts.push_back(posofmin);
	
	// Clean up
	delete th_serie;
	
   	return peaks;
}

// *********************************************************************
// The L2 norm
// *********************************************************************

double DataSerie::L2 (DataSerie &other) {
 	double rv=0;
  	int s=size();
	for (int i=0; i<s; ++i) {
		double t=values[i]-other.values[i];
		rv += (t*t);
	}
	return rv;
}

// *********************************************************************
// The Earth Mover's Distance
// *********************************************************************

double DataSerie::EMD (DataSerie &other) {
	double rv;

	double Hi[256], Hj[256];
	double masseToMove, masse,si, sj;
	double x[256][256];		// the amount of earth to be moved

	int	w[256],	y[256],
		cij[256][256],			// the cost	matrix
		c, k, l, i,	j, s ,	i_min, i_max,
		j_min, j_max,
		N;

	if (size()!=other.size()) 
		ERR_THROW ("Internal error in Histogram<T>::EMD():\n"
			 "Different sizes!\n");
		
	if (size()>256)
		ERR_THROW ("Internal error in DataSerie::EMD():\n"
			 "DataSeries with sizes > 256 are not yet supported!\n");
		
	N  = size();
	si = getSum();
	sj = other.getSum();
	masseToMove	= masse	= (si <	sj ? si	: sj);

	if (masseToMove	== 0)
		return 0;

	 //	Initialise tables
	 //	CW could be	optimised: doesn't need to be done every time	...
	 for (i	= 0	; i	< N	; i++) {
		Hi[i] =	values[i];
		Hj[i] =	other.values[i];

		for	(j = 0 ; j < N ; j++) {
			x[i][j]	= 0	;
			cij[i][j] =	(i > j ? i-j : j-i)	;
		}
	}

	// Search the boundries	where the histograms begin to be non-zero
	i_min =	0 ;
	while (Hi[i_min] ==	0) i_min++ ;
	i_max =	N-1	;
	while (Hi[i_max] ==	0) i_max-- ;
	j_min =	0 ;
	while (Hj[j_min] ==	0) j_min++ ;
	j_max =	N-1	;
	while (Hj[j_max] ==	0) j_max-- ;

	while (masse > 0.0001) {

		// Step	1
		for	(i = i_min ; i <= i_max	; i++)
			if (Hi[i] != 0)
				if (i <	(j_max+j_min)/2)
					w[i] = j_max - i ;
				else
					w[i] = i - j_min ;

		for	(j = j_min ; j <= j_max	; j++)
			if (Hj[j] != 0)
				if (j <	(i_max+i_min)/2)
					y[j] = i_max - j ;
				else
					y[j] = j - i_min ;

		// Step	2 -	Determine source and destination bin
		k =	l =	s =	-1 ;
		for	(i = i_min ; i <= i_max	; i++) {
			if (Hi[i] != 0)	{
				for	(j = j_min ; j <= j_max	; j++) {
					if (Hj[j] != 0)	{
						c =	w[i] + y[j]	- cij[i][j]	;
						if (c >	s) {
							s =	c ;
							k =	i ;
							l =	j ;
						}
					}
				}
			}
		}

		// Step	3 -	How	much earth do we move?
		x[k][l]	= (Hi[k] < Hj[l] ? Hi[k] : Hj[l]) ;

		// Step	4 -	move the earth
		Hi[k] -= x[k][l] ;
		Hj[l] -= x[k][l] ;
		masse -= x[k][l] ;

		if (masse != 0)	{

			// Adjust the boundries, if	the	non-zero part of the
			// histogram got smaller.
			if (Hi[k] == 0)	{
				if (k == i_min)
					while (Hi[i_min] ==	0) i_min++ ;

				if (k == i_max)
					while (Hi[i_max] ==	0) i_max-- ;
			}

			if (Hj[l] == 0)	{
				if (l == j_min)
					while (Hj[j_min] ==	0) j_min++ ;
				if (l == j_max)
					while (Hj[j_max] ==	0) j_max-- ;
			}
		}

	}

	// Take	all	the	amounts	we moved, multiply them	by the costs,
	// and normalize it	by the total mass.
	rv = 0.	;
	for	(i = 0 ; i < N ; i++)
		for	(j = 0 ; j < N ; j++) {
			rv += x[i][j] *	cij[i][j] ;
		}

	rv = rv/masseToMove;

	return rv;
}

// *******************************************************************************
// Output to a stream
// *******************************************************************************

ostream & operator << (ostream & os, DataSerie &d) {
	for (int i=0; i<d.size(); ++i)
		os << d[i] << endl;
	return os;
}

// ***************************************************************
// Write all statistics to a stream
// ***************************************************************

void DataSerie::debugOut (ostream & os) {
	os << "Mean:        " << getMean() << endl
  	 << "Variance:    " << getVariance() << endl
     << "Std.Dev.:    " << sqrt(getVariance()) << endl
     << "Median:      " << getMedian() << endl
     << "Minimum:     " << getMin() << endl
     << "Maximum:     " << getMax() << endl						
     << "Max-Min:     " << getMax() - getMin() << endl;		
}

