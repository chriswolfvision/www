/**********************************************************
 * EMClustererGauss1V.h
 * A class for the Expectaction-Maximization clustering 
 * algorithm
 * Model: univariate Gaussian distribution
 * See research notebook 25.10.2005
 *
 * Christian Wolf
 * Begin: 25.10.2005
 **********************************************************/

#ifndef	_WOLF_EMCLUSTERER1V_H_
#define	_WOLF_EMCLUSTERER1V_H_

// C
#include <math.h>

// C++
#include <vector>

// From the main module
#include <CIL.h>

// From this module
#include <Matrix.h>
#include <Histogram.h>

// T is the scalar type
template <class T>
class EMClustererGauss1V
{
	public:
			
		// Constructor & Destructor
		EMClustererGauss1V (unsigned int nC);
		EMClustererGauss1V (Image &im, unsigned int nC);
		~EMClustererGauss1V();
		
		// Acessor
		unsigned int size()									{ return data.size(); }
		unsigned int operator[] (unsigned int i) const 		{ return labels[i]; }
		
		// Add data
		void add (T value, float weight);	
		
		// Initialize clusters
		void initComponent(unsigned int c, T prior, T mean, T sigma_square);
				
		// Perform the clustering
		void doCluster(int maxIterations);
		
	private:	// METHODS
	
		void _constructor(unsigned int nC);
		void prepare();
		void calcDataLabels();
		void visualize();
		void debugOutMatrices(unsigned int iteration);
	
	private: 	// DATA
		
		unsigned int noClasses;
		vector<T> data;
		vector<float> weights;
		vector<unsigned int> labels;
		
		vector<T> priors, means, sigma_squares;
		bool didInit;
		
		Matrix<T> *p_xn_j; 
		Matrix<T> *p_j_xn; 
		
		// For visualization only
		Histogram<float> *histo;	
};

// ******************************************************
// Constructor: empty data
// ******************************************************

template <class T>
EMClustererGauss1V<T>::EMClustererGauss1V (unsigned int nC) 	
{ 	
	_constructor(nC);
}

// ******************************************************
// Constructor: get data from an integer image
// ******************************************************

template <class T>
EMClustererGauss1V<T>::EMClustererGauss1V (Image &im, unsigned int nC)
{
	_constructor(nC);
	histo = buildHistogram<float> (im);		
	for (unsigned int g=0; g<256; ++g)
		add(g,(*histo)[g]);
	histo->normalize();
}

// ******************************************************
// Help function for the constructor
// ******************************************************

template <class T>
void EMClustererGauss1V<T>::_constructor (unsigned int nC)
{
	noClasses=nC;
	didInit=false;
	p_j_xn = p_xn_j = NULL;
	histo=NULL;
}

// ******************************************************
// Destructor
// ******************************************************

template <class T>
EMClustererGauss1V<T>::~EMClustererGauss1V()
{
	if (p_j_xn!=NULL)
	{
		delete p_j_xn;
		delete p_xn_j;
	}
	if (histo!=NULL)
		delete histo;
}

// ******************************************************
// Initialize the model
// ******************************************************

template <class T>
void EMClustererGauss1V<T>::initComponent(unsigned int c, T prior, T mean, T sigma_square)
{
	if (priors.size() <= c)
	{
		priors.resize(c+1);
		means.resize(c+1);
		sigma_squares.resize(c+1);
	}
	priors[c]=prior;
	means[c]=mean;
	sigma_squares[c]=sigma_square;
	didInit=true;
}

// ******************************************************
// Add data
// ******************************************************

template <class T>
inline void EMClustererGauss1V<T>::add (T value, float weight=1.)
{
	data.push_back(value);	
	weights.push_back(weight);
	labels.push_back(0);
}

// ******************************************************
// Prepare the two matrices storing the 
// likelihoods and the posterior probabilities
// ******************************************************

template <class T>
void EMClustererGauss1V<T>::prepare()
{
	unsigned int N=data.size();
	
	// Create the matrices if they don't exist or 
	// if their size changed.
	if (p_j_xn==NULL)
	{
		p_j_xn = new Matrix<T> (noClasses, N);
		p_xn_j = new Matrix<T> (N, noClasses);
	} 
	else
	{
		if ((p_j_xn->rows()!=noClasses) || (p_j_xn->columns()!=data.size()))
		{
			delete p_j_xn;
			delete p_xn_j;
			p_j_xn = new Matrix<T> (noClasses, N);
			p_xn_j = new Matrix<T> (N, noClasses);
		}
	}

	// Prepare the matrix of likelihoods	
	for (unsigned int j=0; j<noClasses; ++j)
	for (unsigned int n=0; n<N; ++n)
	{
		float d=means[j]-data[n];
		(*p_xn_j)(n,j) = 1./sqrt(2*M_PI*sigma_squares[j])*exp(-0.5*d*d/sigma_squares[j]);
	}
	
	// Prepare the matrix of posterior probabilities
	for (unsigned int n=0; n<N; ++n)
	{
		float p_xn=0;
		for (unsigned int j=0; j<noClasses; ++j)	
			p_xn += (*p_xn_j)(n,j)*priors[j];
		for (unsigned int j=0; j<noClasses; ++j)	
			(*p_j_xn)(j,n) = (*p_xn_j)(n,j)*priors[j]/p_xn;
	}
}

// ******************************************************
// Label each data vector according to the 
// highest posterior probablity
// ******************************************************

template <class T>
void EMClustererGauss1V<T>::calcDataLabels()
{
	T max=0;
	unsigned int argmax=0;
	
	prepare();
	
	for (unsigned int n=0; n<data.size(); ++n)
	{
		for (unsigned int j=0; j<noClasses; ++j)
		{
			if (j==0 || (*p_j_xn)(j,n)>max)
			{
				argmax=j;
				max=(*p_j_xn)(j,n);				
			}
		}
		labels[n]=argmax;
	}
}

// ******************************************************
// Perform the clustering
// ******************************************************

template <class T>
void EMClustererGauss1V<T>::visualize ()
{
#ifdef HAVE_VISUALIZE
	static unsigned int funcCallCount=0;
	
	// Create a new file
	ostringstream s; 
	s << "likelihoodplot_em_intern_" << setfill ('0') << setw(3) 
	  << funcCallCount++ << ".txt"; 
	ofstream st (s.str().c_str(), ios::out);
	if (!st.good())
		ERR_THROW ("Cannot open file " << s.str() << "for writing!\n");

	// Plot the histogram		
	st << *histo << endl;

	// Debug: plot the class histograms as well as the fitted distributions
	for (unsigned int c=0; c<noClasses; ++c)
	{							 												 
		plotGaussian (st, means[c], sigma_squares[c], 0, 255, 1, priors[c], true);
		st << endl;		
	}
	st.close();	
#else
	ERR_THROW ("Visualization has not been compiled into the code!\n";
#endif
}

// ******************************************************
// Perform the clustering
// ******************************************************

template <class T>
void EMClustererGauss1V<T>::doCluster (int maxIter)
{	
	int it=0;	
	unsigned int N=data.size();
	
	cerr << "EM-Algo Clustering.\n"; cerr.flush();
	
	if (!didInit)
		throw EError ("EMClustererGauss1V::doCluser(): you need to initialize first!");		
		
	cerr << "EM: initialization::\n";
	for (unsigned int j=0; j<noClasses; ++j)
		cerr << "class " << j << ": m=" << means[j] 
				<< "\tsigma=" << sqrt(sigma_squares[j]) 
				<< "\tP(j)=" << priors[j] << endl;

#ifdef WRITE_DEBUG_IMAGES
	visualize();
#endif			
	
	// Iterate
	do 
	{							
		prepare();
		
#ifdef WRITE_DEBUG_IMAGES
		debugOutMatrices(it);
#endif
		
		// Update the parameters of the mixture	
		for (unsigned int j=0; j<noClasses; ++j)
		{
			float z,res,NWeighted;			
			
			z=0;
			for (unsigned int n=0; n<N; ++n)
				z += weights[n]*(*p_j_xn)(j,n);
				
			// Update the mean
			res=0;
			for (unsigned int n=0; n<N; ++n)
				res += weights[n]*(*p_j_xn)(j,n)*data[n];				
			means[j]=res/z;
		
			// Update the variance
			res=0;
			for (unsigned int n=0; n<N; ++n)
			{
				float d=data[n]-means[j];
				res += weights[n]*(*p_j_xn)(j,n)*d*d;
			}
			sigma_squares[j]=res/z;	
			
			// Update the prior
			res=NWeighted=0;
			for (unsigned int n=0; n<N; ++n)
			{
				res += weights[n]*(*p_j_xn)(j,n);
				NWeighted+=weights[n];
			}
			priors[j]=res/NWeighted;					
		}	
			
		cerr << "EM: result after " << it+1 << " iterations:\n";
		
		// Update the parameters of the mixture	
		for (unsigned int j=0; j<noClasses; ++j)
			cerr << "class " << j << ": m=" << means[j] 
				 << "\tsigma=" << sqrt(sigma_squares[j]) 
				 << "\tP(j)=" << priors[j] << endl;
					 
#ifdef WRITE_DEBUG_IMAGES
		visualize();
#endif			
		++it;
	} while (it<maxIter);
			
	cerr << endl;	
}

#endif


// ******************************************************
// debugOutput
// ******************************************************

template <class T>
void EMClustererGauss1V<T>::debugOutMatrices(unsigned int iteration)
{
	unsigned int N=data.size();
	
	{ 	
		ostringstream s; 
		s << "p_xn_j_" << iteration << ".txt"; 
		ofstream st (s.str().c_str(), ios::out);
		if (!st.good())
			ERR_THROW ("Cannot open file '" << s.str().c_str() << "' for writing!\n");
		
		for (unsigned int j=0; j<noClasses; ++j)
		{
			for (unsigned int n=0; n<N; ++n)
				st << n << " " << (*p_xn_j)(n,j) << endl;
			st << endl;
		}
		
		st.close();
	}
	
	{ 	
		ostringstream s; 
		s << "p_j_xn_" << iteration << ".txt"; 
		ofstream st (s.str().c_str(), ios::out);
		if (!st.good())
			ERR_THROW ("Cannot open file '" << s.str().c_str() << "' for writing!\n");
		
		for (unsigned int j=0; j<noClasses; ++j)
		{
			for (unsigned int n=0; n<N; ++n)
				st << n << " " << (*p_j_xn)(j,n) << endl;
			st << endl;
		}
		
		st.close();
	}
}
