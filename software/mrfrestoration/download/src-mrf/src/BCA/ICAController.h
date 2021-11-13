/***********************************************************************
 * ICAController.h
 *
 * Bleed through separation using "quasi-linear" ICA 
 * (Independent Component Analysis)
 *
 * Author: Christian Wolf
 *         christian.wolf@insa-lyon.fr
 * 
 * Changelog:
 * 25.04.2006 cw: First version
 *
 * 1 tab = 4 spaces
 ***********************************************************************/
 
#ifndef _WOLF_ICACONTROLLER_H_
#define _WOLF_ICACONTROLLER_H_

// From the main module
#include <CIL.h>

// From the NUMERICAL RECIPES IN C library
#include <NumericalRecipes.h>

/***********************************************************************
 * PARAMETERS
 ***********************************************************************/
 
// To keep the values in a certain range
#define ADJUST_FACTOR			1.0

// Stopping criterion for the minimization algorithm
#define FRACT_TOLERANCE 		0.1 
//#define GRAD_TOLERANCE 			0.1 

/***********************************************************************/
 
// C++
#include <vector>
#include <iostream>

// From the main module
#include <CIL.h>

// From the MATH module
#include <Vector.h>
#include <Matrix.h>
#include <MathFuncs.h>

// From the IMAGE module
#include <Image.h>

// From this module
#include "BCAObjectiveFuncData.h"

enum LabelType
{
	LAB_RECTO=0,
	LAB_VERSO,
	LAB_BG
};

typedef float ScalarT;
typedef Vector<ScalarT> VectorT;
typedef Matrix<ScalarT> MatrixT;

/***********************************************************************
 * Template parameters:
 * Tim ..... The image type
 ***********************************************************************/

template <class Tim>
class ICAController
{
	public:
	
		// Constructor & Destructor
		ICAController (Tim *input);
		ICAController(){}
		
		void prepareData (int stepx, int stepy);
		void freeData();
		void solve();		
		void simulate (Tim *output);
		void runInverseSystem (Tim *out_recto, Tim *out_verso);
		
	protected:
	
		void initialEstimate ();
		void printParameters();
		
	protected: // DATA
	
		unsigned int dim;
		MatrixT A, W;
		VectorT p, pInv;
		
		Tim *inputImage;
		
	public:		// CLASS MEMBERS
		
		
 		// The data used by the objective function must be globally accessible
		static BCAObjectiveFuncData *ObjFuncData;
};


/***********************************************************************
 * Initialisation of the static class members
 ***********************************************************************/

template <class T>
BCAObjectiveFuncData *ICAController<T>::ObjFuncData = NULL;

/***********************************************************************
 * The constructor
 ***********************************************************************/
 
template <class Tim> 
ICAController<Tim>::ICAController (Tim *input)
{
	inputImage=input;
	
	// The dimension of the problem equals 
	// the number of colors/bands in the image
	dim = input->nbColorPlanes();	
	A.resize(dim,dim);
	p.resize(dim);
	
	// Manually set up the parameters:
	// The weight vector
	p[LAB_RECTO]=100;
	p[LAB_VERSO]=2;
	p[LAB_BG]=1;
	
	// The mixing matrix for image2_200x200.ppm
	A(0,0)=46; 	// The color of the recto side
	A(1,0)=25;
	A(2,0)=12;
	A(0,1)=150;	// The color of the verso side
	A(1,1)=113;
	A(2,1)=74;
	A(0,2)=216;	// The color of the back ground
	A(1,2)=183;
	A(2,2)=138;
		
	// Invert the mixing matrix 
 	// and the weight vector 
 	// (In the paper we invert a triangular weight matrix,
 	// but it's more efficient to store that data in a vector
	W = A;
	W.invert();	
	pInv = p;
	pInv.elInverse();
	printParameters();
}

/***********************************************************************
 * Print the model parameters
 ***********************************************************************/
 
template <class Tim> 
void ICAController<Tim>::printParameters ()
{
	HRULE;
	cerr << "Mixing matrix (A):" << endl << A;
	cerr << "Inverse mixing matrix (W):" << endl << W;
	cerr << "Weights (p):" << endl << p;
	cerr << "Inverse weights matrix (pInv):" << endl << pInv;	
	HRULE;
}

/***********************************************************************
 * Take the parameters (vectors and matrices) and put them into a 
 * single 1D vector which can be used in optimization algorithms
 ***********************************************************************/
 
inline void serializeParameters (MatrixT &W, VectorT &pInv, float *theta, 
	unsigned int dim, unsigned int dim_theta)
{
#warning NO ESTIMATION OF THE WEIGHTS!!!
	for (unsigned int i=0; i<dim; ++i)
	for (unsigned int j=0; j<dim; ++j)
		theta[i*dim+j]=W(i,j);
	/*
	for (unsigned int i=0; i<dim-1; ++i)
		theta[i]=pInv[i];
	for (unsigned int i=0; i<dim; ++i)
	for (unsigned int j=0; j<dim; ++j)
		theta[dim-1+i*dim+j]=W(i,j);
	*/
	
	// Adjust the range of the values
	for (unsigned int i=0; i<dim_theta; ++i)
		theta[i] *= ADJUST_FACTOR;
}

/***********************************************************************
 * The inverse operation: recreate the parameters (vectors and matrices) 
 * from the single 1D vector which can be used in optimization algorithms
 ***********************************************************************/
 
inline void deserializeParameters (float *theta, MatrixT &W, VectorT &pInv, 
	unsigned int dim, unsigned int dim_theta)
{
#warning NO ESTIMATION OF THE WEIGHTS!!!
	for (unsigned int i=0; i<dim; ++i)
	for (unsigned int j=0; j<dim; ++j)
		W(i,j)=theta[i*dim+j];
	pInv[0]=1./10.;
	pInv[1]=1./2.;
	pInv[2]=1.;
	/*
	for (unsigned int i=0; i<dim-1; ++i)
		pInv[i]=theta[i];
	for (unsigned int i=0; i<dim; ++i)
	for (unsigned int j=0; j<dim; ++j)
		W(i,j)=theta[dim-1+i*dim+j];
	*/
	
	for (unsigned int i=0; i<dim_theta; ++i)
		theta[i] /= ADJUST_FACTOR;
}

/***********************************************************************
 * run the inverse system on a single vector
 ***********************************************************************/

inline void runInverseOnVector (MatrixT &W, VectorT &pInv, VectorT &p, VectorT &v)
{							 				 
	v *= v.dotProduct(p);
	v = W*v;
	v = v.elProduct(pInv);		
}

inline void denormalizeVector (VectorT &v, unsigned int dim, char *debugStr)
{	
	float max_abs,signed_of_max_abs;
	
	if (debugStr!=NULL)
	{
		cerr << debugStr << "s= " << v;
	}

	signed_of_max_abs=v[0];
	max_abs=fabs(signed_of_max_abs);		
	if (debugStr!=NULL)
		cerr << debugStr << "[max=" << max_abs << "]\n";
	for (unsigned int i=1; i<dim; ++i)
	{
		if (debugStr!=NULL)
			cerr << debugStr << "[max=" << max_abs << " test=" << fabs(v[i]) << "]\n";
		if (fabs(v[i])>max_abs)
		{
			max_abs=fabs(v[i]);
			signed_of_max_abs=v[i];
		}
	}
	if (debugStr!=NULL)
		cerr << debugStr << "[max=" << max_abs << "]\n";
	
	v *= (1/signed_of_max_abs);
	
	if (debugStr!=NULL)
	{
		cerr << debugStr << "sdenorm= " << v;
		cerr << debugStr << "div=" << signed_of_max_abs << endl;			
	}
}

/***********************************************************************
 * Calculate an initial estimate of p and A
 ***********************************************************************/

template <class Tim> 
void ICAController<Tim>::initialEstimate ()
{		
	p[0]=10;	// the weight for recto
	p[1]=2;		// the weight for verso
	p[2]=1;		// the weight for the background
#warning THE INITIAL ESTIMATE IS DONE MANUALLY WITH SIMULATED CORRECT VALUES!!!	
	A(0,0)=46; 	// The color of the recto side
	A(1,0)=25;
	A(2,0)=12;
	A(0,1)=150;	// The color of the verso side
	A(1,1)=113;
	A(2,1)=74;
	A(0,2)=216;	// The color of the back ground
	A(1,2)=183;
	A(2,2)=138;
	W = A;
	W.invert();	
	pInv = p;
	pInv.elInverse();
}

/***********************************************************************
 * Prepares the data which is used by the objective function, which 
 * is used by the solver
 ***********************************************************************/
 
template <>
void ICAController<Image>::prepareData (int stepx, int stepy)
{		
	ObjFuncData = new BCAObjectiveFuncData;
	VectorT *v;
	
	for (int y=0; y<inputImage->ysize; y+=stepy)	
	for (int x=0; x<inputImage->xsize; x+=stepx)
	{
		v = new VectorT(3);
		(*v)[0]=inputImage->get(PLANE_RED,x,y);
		(*v)[1]=inputImage->get(PLANE_GREEN,x,y);
		(*v)[2]=inputImage->get(PLANE_BLUE,x,y);
		
		ICAController<Image>::ObjFuncData->data.push_back(v);
	}	
	
	cerr << "Loaded " << ICAController<Image>::ObjFuncData->size() 
		 << " input vectors." << endl;
}

/***********************************************************************
 * Free the data used by the objective function
 ***********************************************************************/
 
template <class Tim> 
void ICAController<Tim>::freeData ()
{			
	delete ICAController<Tim>::ObjFuncFata;	
	ICAController<Tim>::ObjFuncFata=NULL;
}

/***********************************************************************
 * The function which is to be minimized  
 * It must be a global function, not a method.
 * See Research Notebook 30.4.2006, p. 56
 ***********************************************************************/
 
// Access to the different parts inside the parameter vector
// (a vector and a matrix) 
#define param_p(p,i)		((p)[i])
#define param_W(p,dim,i,j)	((p)[dim-1+(i)*(dim)+(j)])
 
template <class Tim>
float objectiveFunction (float *theta)
{	
	unsigned int dim=ICAController<Tim>::ObjFuncData->dim;
	unsigned int cntUsedData;
	VectorT s(dim),sn(dim), *px;
	ScalarT angle,min,mean_angle,s_norm;
	ScalarT regValue;
	MatrixT _W(dim,dim);
	VectorT _pInv(dim), _p(dim);
	vector<unsigned int> cnts;
		
	cerr << "OBJFUNC";
	for (unsigned int i=0; i<ICAController<Tim>::ObjFuncData->dim_theta; ++i)		
		cerr << " " << theta[i];
	cerr << endl;
		
	// Restore the parameters into "W" and "pInv"
	deserializeParameters(theta,_W,_pInv,dim,ICAController<Tim>::ObjFuncData->dim_theta);
	_p = _pInv;
	_p.elInverse();
	
	/* DEBUG: a trivial objective function to test
	 * the optimization routine
	mean_angle=0;
	MatrixT m(dim,dim);
	VectorT v(dim);
	for (unsigned int i=0; i<dim; ++i)
	{
		mean_angle+=pow(i-_pInv[i],2);
		for (unsigned int j=0; j<dim; ++j)
			mean_angle+=pow(i*dim+j-_W(i,j),2);
	}		
	*/
	
	mean_angle=0;
	cntUsedData=0;
	for (unsigned int i=0; i<ICAController<Tim>::ObjFuncData->noTargets(); ++i)
		cnts.push_back(0);
	
	// Travers the input data	
	for (unsigned int j=0; j<ICAController<Tim>::ObjFuncData->size(); ++j)
	{
		// Calculate the inverse s for a given x
		px=(*ICAController<Tim>::ObjFuncData)[j];
		
		s=*px;
		runInverseOnVector(_W,_pInv,_p,s);
		s_norm=s.norm();
		
		// If the norm is zero, we do not use this test vector
#warning ZERO TEST SHOULD PERHAPS BE DONE WITH A THRESHOLD		
		if (s_norm==0)
			continue;
		++cntUsedData;
		
		// ------------------------------------------------------------------------------
		// The 4 terms corresponding to the 4 binary components
		for (unsigned int comp=0; comp<ICAController<Tim>::ObjFuncData->noTargets(); ++comp)
		{
			VectorT *pb=ICAController<Tim>::ObjFuncData->targets[comp];
			
#ifdef CHECK_CODE			
			if (s_norm==0)
			{
				cerr << "Norm of s is zero! s=" << s << ", x=" << *px << endl;
				exit(1);
			}
#endif			

			// #warning DIFFERENCE INSTEAD OF ANGLE!!!!
			// angle=(s-(*pb)).norm();
			
			angle=acos(s.dotProduct(*pb)/
				(s_norm*ICAController<Tim>::ObjFuncData->target_norms[comp]));
			
			if (comp==0)
				min=angle;
			else
				min=NeuroFuzzyMinimum(min,angle,1.);
		}
		mean_angle += min;
		
		// ------------------------------------------------------------------------------
		// The regularisation term
		
		sn = s;
		denormalizeVector(sn,dim,NULL);
		for (unsigned int i=0; i<dim; ++i)
			if (sn[i]>0.75) 
				sn[i]=1;
			else
				sn[i]=0;
						
		// Count the number of pixels for each target vector
		for (unsigned int i=0; i<ICAController<Tim>::ObjFuncData->noTargets(); ++i)
		{
			if (sn==*(ICAController<Tim>::ObjFuncData->targets[i]))
			{
				++cnts[i];
				break;
			}			
		}				
	}
	mean_angle/= (ScalarT) cntUsedData;
	
	cerr << "CLASS COUNTS:\n";
	for (unsigned int i=0; i<ICAController<Tim>::ObjFuncData->noTargets(); ++i)
		cerr << cnts[i] << " px for class: " << *(ICAController<Tim>::ObjFuncData->targets[i]);
	
	// Add the regularisation term
	regValue=1;
	DataSerie ds;
	for (unsigned int comp=0; comp<ICAController<Tim>::ObjFuncData->noTargets(); ++comp)
	{
		// We do not use the RECTO&VERSO pixels for regularisation. There aren't 
		// necessarily one of those, and 1 zero in the counts penalizes to the maximum!
		if (comp!=TARGET_IND_RECTO_VERSO)
		{	
			ds.add((ScalarT) cnts[comp] / (ScalarT) cntUsedData);
			cerr << "reg: " << ((ScalarT) cnts[comp] / (ScalarT) cntUsedData) << endl;
		}
	}
	regValue = 10*(ds.getVariance()/ds.getMean());
		
	cerr << "VALUE=" << mean_angle << " + " << regValue << " = " << mean_angle+regValue << endl;
	
	return mean_angle+regValue;
} 

/***********************************************************************
 * Solve for the parameters p and A
 ***********************************************************************/
 
#define UNITOFFSET(x)		((x)-1)
 
template <class Tim> 
void ICAController<Tim>::solve ()
{		
	int noIterations;
	float min_func_value;
	float *theta; 	
	float **xi;
	
	if (ObjFuncData==NULL)
		ERR_THROW ("prepareData() must be called before solve()!\n");
	
	/* The parameter vector theta:
	 * 0..dim-2:              the first dim-1 values of inv(p) (the last value is supposed to be 1!!)
	 * dim-1..dim*dim+dim-2:  W, row by row
	 */
	theta = new float [ObjFuncData->dim_theta];
	
	// Estimate an initial parameters and 
	// copy them into the parameter vector
	initialEstimate();		
	serializeParameters(W, pInv, theta, dim, ObjFuncData->dim_theta);
	
	cerr << "The initial estimate:\n";
	printParameters();
	
	cerr << "Obj. function on the initial estimate: " << objectiveFunction<Image>(theta) << endl;
	theta[3] += 0.005;
	cerr << "After 1. change: " << objectiveFunction<Image>(theta) << endl;
				
	// ---------------------------------------------------------------
	// Minimize the objective function with Powell's method from the 
	// book "Numerical Recipes in C"
	// ---------------------------------------------------------------
	
	// Create the matrix with the initial directions (in the columns)
	// Attention! the numrec library works with matrices whose 
	// indices begin with 1!!!!
	xi = numrec::matrix(1,ObjFuncData->dim_theta,1,ObjFuncData->dim_theta);
	for (unsigned int r=1; r<=ObjFuncData->dim_theta; ++r)
	for (unsigned int c=1; c<=ObjFuncData->dim_theta; ++c)
	{
		if (c==r)
			xi[r][c]=1.;
		else
			xi[r][c]=0;
	}	
	
	// ATTENTION!!!
	// p = [1..dim_theta]
	// xi = [1..dim_theta,1..dim_theta]
	numrec::powell(UNITOFFSET(theta), xi, ObjFuncData->dim_theta, FRACT_TOLERANCE, &noIterations, 
		&min_func_value, objectiveFunction<Image>);
	
	cerr << "Optimal value: " << min_func_value << endl
		 << "Reached in " << noIterations << " iterations.\n";
	
	// Get the parameters from the vector
	deserializeParameters(theta,W,pInv,ObjFuncData->dim, 
		ICAController<Tim>::ObjFuncData->dim_theta);
		
	A = W;
	A.invert();	
	p = pInv;
	p.elInverse();
	cerr << "Parameters of the solved system: " << endl;
	printParameters();
		
	// clean up
	numrec::free_matrix(xi,1,ObjFuncData->dim_theta,1,ObjFuncData->dim_theta);	
	delete theta;
}

/***********************************************************************
 * Run the inverse system
 ***********************************************************************/
 
inline unsigned char trim (ScalarT v)
{
	if (v<0)
		v*=-1;
	if (v>255)
		return 255;
	else
		return (unsigned char) rint(v);
}
 
template <class Tim> 
void ICAController<Tim>::runInverseSystem (Tim *out_recto, Tim *out_verso)
{
	unsigned int xs = inputImage->xsize;
	unsigned int ys = inputImage->ysize;	
	
	if (out_recto==NULL || out_verso==NULL)
		ERR_THROW ("internal error in ICAController<TIM>::runInverseSystem(): "
			"no output image!\n");
		
		
	cerr << "Running the restoration process...\n"; cerr.flush();
		
	for (unsigned y=0; y<ys; ++y)
	for (unsigned x=0; x<xs; ++x)
	{
		VectorT v ((ScalarT) inputImage->get(PLANE_RED, x, y),
				   (ScalarT) inputImage->get(PLANE_GREEN, x, y),
				   (ScalarT) inputImage->get(PLANE_BLUE, x, y));						 				 
						 						 
		if ((x==22 && y==126) || (x==212 && y==140) || (x==232 && y==148) || (x==222 && y==215))
		{
			char buf[100];
			sprintf (buf,"(%d,%d):",x,y);			
			cerr << buf << "x=" << v;						
			runInverseOnVector (W,pInv,p,v);
			denormalizeVector (v,dim,buf);
		}
		else
		{		
			runInverseOnVector (W,pInv,p,v);
			denormalizeVector (v,dim,NULL);
		}
						 				 
		v *= 150.;					
				
		TRIMGRAY(v[0]);
		TRIMGRAY(v[1]);
		out_recto->set(PLANE_RED, x, y, (unsigned char) rint(v[0]));
		out_verso->set(PLANE_RED, x, y, (unsigned char) rint(v[1]));		
	}
	
	cerr << "Image restored.\n"; cerr.flush();
}

/***********************************************************************
 * Simulate the mixing process
 ***********************************************************************/
 
template <class Tim> 
void ICAController<Tim>::simulate (Tim *output)
{
	unsigned int xs = inputImage->xsize;
	unsigned int ys = inputImage->ysize;
	VectorT w;
	
	if (output==NULL)
		ERR_THROW ("internal error in ICAController<TIM>::simulate(): "
			"no output image!\n");
		
	for (unsigned y=0; y<ys; ++y)
	for (unsigned x=0; x<xs; ++x)
	{
		VectorT v ((ScalarT) inputImage->get(PLANE_RED, x, y),
				   (ScalarT) inputImage->get(PLANE_GREEN, x, y),
				   (ScalarT) inputImage->get(PLANE_BLUE, x, y));
						 				 
		w = v.elProduct(p);						 				 
		w = A*w;
		w /= v.dotProduct(p);
			
		output->set(PLANE_RED,   x, y, (unsigned char) w[0]);
		output->set(PLANE_GREEN, x, y, (unsigned char) w[1]);
		output->set(PLANE_BLUE,  x, y, (unsigned char) w[2]);
	}	
}

#endif


