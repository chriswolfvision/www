/***************************************************************************
 * VectorDistances.h
 ***************************************************************************/

#include "MathFuncs.h"

inline double Bhattacharyya (double x1, double x2, double v1, double v2) 
{
	double var_comb = (v1+v2)/2.0;
	
	// Mean shift term
    return 	(x1-x2)*(x1-x2)/(var_comb*8.0) +	
   	// Covariance shift term
    		0.5*log(fabs(var_comb)/( (sqrt(fabs(v1))*sqrt(fabs(v2)))));    		    		
}

template <class T>
double MahalanobisDiff (Vector<T> &diff_vec, Matrix<T> &invcov) 
{
	unsigned int s = diff_vec.size;
	Vector<T> t(s);
	T sum;
	
	// The first multiplikation: The diff vector (left) with the
	// inverse of the covariance matrix
	for (unsigned int i=0; i<s; ++i) 
	{
		sum = 0;
		for (unsigned int j=0; j<s; ++j)
			sum += (diff_vec[j]*invcov.get(j,i));
		t.set(i, sum);
	}
	
	// The second multiplication: The scalar product of the first
	// result with the diff vector (right)
	sum = 0;
	for (unsigned int i=0; i<s; ++i)
		sum += (t[i]*diff_vec[i]);
	
	return sum;	
}

template <class T>
double Mahalanobis (Vector<T> &x, Vector<T> &y, Matrix<T> &invcov) 
{
	Vector<T> diff_vec = x;
	diff_vec -= y;	
	return MahalanobisDiff (diff_vec, invcov);
}