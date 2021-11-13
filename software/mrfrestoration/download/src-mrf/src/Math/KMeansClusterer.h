/**********************************************************
 * KMeansClusterer.h
 * A class for the K-Means clustering algorithm
 *
 * Christian Wolf
 * Beginn: 30.5.2005
 **********************************************************/

#ifndef	_WOLF_KMEANS_H_
#define	_WOLF_KMEANS_H_

// C++
#include <vector>

// From the main module
#include <CIL.h>

// From this module
#include "Vector.h"
#include "RandomNumberGenerator.h"

// ******************************************************
// Class definition
// T... the datatype storing a scalar, not a vector!!!!
// ******************************************************

// This class can be given as an argument to provide the
// clusterer with a custom distance function
class ColorDistanceFunctional
{
	public:
		virtual float distance (Vector<float> &, Vector<float> &)=0;
};

// The clusterer itself
template <class T>
class KMeansClusterer
{
	public:
	
		// typedef Vector<T> VectorType;		
		typedef Vector<float> VectorTypeFloat;
		
		// Constructor & Destructor
		KMeansClusterer (unsigned int d, unsigned int nC);
		~KMeansClusterer();
		
		// Acessor
		unsigned int operator[] (unsigned int i) const 	{ return labels[i]; }
		unsigned int size()								{ return labels.size(); }
		Vector<float> getCenter(unsigned int i) const 	{ return centers[i]; }
		
		// Add data
		void add (VectorTypeFloat v);	
		
		// Initialize clusters
		void initRandom();
		void initZero();
		void init(unsigned int i, T v);		
		void init(unsigned int i, Vector<float> &v);
		
		// Perform the clustering
		void doCluster(int maxIterations=20);
		
		// define a custom distance function between data points
		void setDistance(ColorDistanceFunctional *d) 		{ distance=d; }
		
		void printMeans();
		void printData ();
		
	private:	// METHODS
	
		void  calcDataLabels();
		float calcClusterCenters();
		void _initRandom();
		void getRange();
	
	protected: 	// DATA
		
		unsigned int noClasses;
		unsigned int dim;
		vector<VectorTypeFloat> data;
		vector<unsigned int> labels;
		
		// The centers are always floats, even
		// if we are working with integer data
		vector<VectorTypeFloat> centers;		
		
		bool didInit;
		Vector<float> *min, *max;
		RandomNumberGenerator *rngen;
		
		// The distance function used 
		ColorDistanceFunctional *distance;
};

// ******************************************************
// Constructor
// ******************************************************

template <class T>
KMeansClusterer<T>::KMeansClusterer (unsigned int d, unsigned int nC) 	
{ 	
	dim=d; 
	noClasses=nC; 
	didInit=false;
	distance=NULL;
	min = new Vector<float>(dim);
	max = new Vector<float>(dim);
	rngen = new RandomNumberGenerator (true, 0, 1, 32);
}

// ******************************************************
// Destructor
// ******************************************************

template <class T>
KMeansClusterer<T>::~KMeansClusterer () 	
{ 	
	delete min;
	delete max;
	delete rngen;
}

// ******************************************************
// Add data
// ******************************************************

template <class T>
inline void KMeansClusterer<T>::add (VectorTypeFloat v)
{
	data.push_back(v);	
	labels.push_back(0);
}

// ******************************************************
// Get the range of the data
// ******************************************************

template <class T>
inline void KMeansClusterer<T>::getRange() 
{
	unsigned int no=data.size();
	
	for (unsigned int d=0; d<dim; ++d)
		(*min)[d]=data[0][d];
	*max=*min;
	for (unsigned int i=1; i<no; ++i)
	{
		for (unsigned int d=0; d<dim; ++d)
		{
			if (data[i][d]<(*min)[d])
				(*min)[d]=data[i][d];
			if (data[i][d]>(*max)[d])
				(*max)[d]=data[i][d];
		}
	}	
}

// ******************************************************
// Initialize to zero clusters
// ******************************************************

template <class T>
void KMeansClusterer<T>::initZero()
{
	centers.resize(noClasses);
	didInit=true;
}

// ******************************************************
// Initialize a cluster center (manually)
// ******************************************************

template <class T>
void KMeansClusterer<T>::init(unsigned int i, Vector<float> &v)
{
	centers[i]=v;
}

// template <class T>
// void KMeansClusterer<T>::init(unsigned int i, VectorType &v)
// {
// 	Vector<float> vf;
// 	for (unsigned int j=0; j<v.size(); ++j)
// 		vf[j] = (float) v[j];
// 	centers[i]=vf;
// }

// This one only works for one dimensional data....
template <class T>
void KMeansClusterer<T>::init(unsigned int i, T v)
{
	Vector<float> vf(1);
	vf[0] = (float) v;
	centers[i]=vf;
}

// ******************************************************
// Initialize: use random clusters
// ******************************************************

template <class T>
inline void KMeansClusterer<T>::initRandom ()
{
	getRange();
	_initRandom();
}

template <class T>
inline void KMeansClusterer<T>::_initRandom ()
{
	
	unsigned int no=data.size();	
	
	if (no<1)
		throw ("KMeansClusterer::initRandom(): no data!");	
	
	// Create the cluster variables
	centers.resize(noClasses, Vector<float>::createZero(dim));
	
	// Create the clusters in the range
	for (unsigned int c=0; c<noClasses; ++c)
		for (unsigned int d=0; d<dim; ++d)
			centers[c][d] = rngen->next()*((*max)[d]-(*min)[d])+(*min)[d];

	didInit=true;
}

// ******************************************************
// Label each data vector according to the nearest
// cluster center
// ******************************************************

template <class T>
void KMeansClusterer<T>::calcDataLabels()
{
	unsigned int no=data.size();
	
	// ---------------------------------------------------
	// The standard euclidean distance between the 
	// data points
	
	if (distance==NULL)
	{	
		// Traverse the data
		for (unsigned int i=0; i<no; ++i)
		{	
			float dist, minDist;
			unsigned int minIndex;
				
			// For each point determine the nearest cluster center
			minDist = centers[0].euclideanDistanceSquared(data[i]);
			minIndex = 0;
			for (unsigned int c=1; c<noClasses; ++c)
			{
				dist = centers[c].euclideanDistanceSquared(data[i]);
				if (dist < minDist)
				{
					minDist=dist;
					minIndex=c;
				}			
			}
			labels[i]=minIndex;
		}
	}
		
	// ---------------------------------------------------
	// A user defined distance between the data points
		
	else
	{	
		// Traverse the data
		for (unsigned int i=0; i<no; ++i)
		{	
			float dist, minDist;
			unsigned int minIndex;
				
			// For each point determine the nearest cluster center
			minDist = distance->distance(centers[0],data[i]);
			minIndex = 0;
			for (unsigned int c=1; c<noClasses; ++c)
			{
				dist = distance->distance(centers[c],data[i]);
				if (dist < minDist)
				{
					minDist=dist;
					minIndex=c;
				}			
			}
			labels[i]=minIndex;
		}
	}
}

// ******************************************************
// Calculate the cluster centers from the
// labeling
// Return the sum of the distances between the old and the
// new centers.
// or -1 if one of the centers vanished.
// ******************************************************

template <class T>
float KMeansClusterer<T>::calcClusterCenters()
{
	vector<VectorTypeFloat> newcenters;		
	vector<unsigned int> noDataPerCenter;
	unsigned int no=data.size();
	float dist, distTot; 
	
	newcenters.resize(noClasses, VectorTypeFloat::createZero(dim));
	noDataPerCenter.resize(noClasses);
	for (unsigned int c=0; c<noClasses; ++c)
	{
		newcenters[c] = VectorTypeFloat::createZero(dim);
		noDataPerCenter[c] = 0;
	}
		
	// Traverse the data
	for (unsigned int i=0; i<no; ++i)
	{
		newcenters[labels[i]] += data[i];		
		++noDataPerCenter[labels[i]];
	}
	
	for (unsigned int c=0; c<noClasses; ++c)
	{
		// A center vanished: we stop the process
		// and inform the calling function which 
		// will reinitialize the process
		if (noDataPerCenter[c]<=0)
			return -1;
		newcenters[c] /= noDataPerCenter[c];
	}
		
	// Calculate the amount the centers moved
	// and replace the centers
	distTot=0;
	for (unsigned int c=0; c<noClasses; ++c)
	{		
		dist = newcenters[c].euclideanDistance(centers[c]);
		distTot += dist;
#ifdef CHECK_CODE
		if (isnan(dist))
			ERR_THROW ("Internal error in KMeansClusterer<T>::calcClusterCenters\n"
				<< "eucl. distance is nan for center c" << c 
				<< ": old=" << centers[c][0] 
				<< ", new=" << newcenters[c][0]);
#endif		
		centers[c] = newcenters[c];		
	}	
	return distTot;
}

// ******************************************************
// Perform the clustering
// ******************************************************

template <class T>
void KMeansClusterer<T>::doCluster (int maxIter)
{	
	int it=0, 
		noCenterDisappearances=0;
	float dist;
	
	if (!didInit)
		ERR_THROW ("KMeansClusterer::doCluser(): you need to initialize first!");
	
	// Iterate
	do 
	{	
		do
		{
			// printMeans();
			
			calcDataLabels();				
			dist = calcClusterCenters();
				
			if (dist<0)
			{			
				if (it==0)
				{	
					cerr << "A cluster center disappeared. Re-initializing\n";
					initRandom();
					++noCenterDisappearances;
				}
				else
					ERR_THROW ("K-means clustering: a cluster center disappeared, "
						"something is fishy here. Aborting.");
			}
			
		} while (dist<0 && noCenterDisappearances<5);
		
		if (dist<0)
			ERR_THROW ("K-means clustering: Cannot initialize the cluster labels w/o "
				" causing a cluster center to disappear after the first iteration. "
				" Something is wrong with the data, aborting.");
				
		cerr << "[" << dist << "]";		
		++it;
			
	} while (it<maxIter && dist>0);
	
	cerr << endl;	
}

// ******************************************************
// Print the means to stdout
// ******************************************************

template <class T>
void KMeansClusterer<T>::printMeans ()
{
	for (unsigned int c=0; c<noClasses; ++c)
		cout << "Mean 1: " << centers[c];
}

// ******************************************************
// Print the data to stdout
// ******************************************************

template <class T>
void KMeansClusterer<T>::printData ()
{
	for (unsigned int i=0; i<data.size(); ++i)
	{
		VectorTypeFloat v=data[i];
		cout << "Clustering-data:";
		for (unsigned int d=0; d<dim; ++d)
			cout << " " << (double) v[d];
		cout << endl;
	}
}

#endif

