/***********************************************************************
 * BCAObjectiveFuncData.h
 *
 * Data needed by the objective function with is minimized in the 
 * framework of the BCA "Binary Components Analysis" algorithm
 *
 * Author: Christian Wolf
 *         christian.wolf@insa-lyon.fr
 * 
 * Changelog:
 * 30.04.2006 cw: First version
 *
 * 1 tab = 4 spaces
 ***********************************************************************/
 
#ifndef _WOLF_BCAOBJECTIVEFUNCDATA_H_
#define _WOLF_BCAOBJECTIVEFUNCDATA_H_

// C++
#include <vector>

typedef float ScalarT;
typedef Vector<ScalarT> VectorT;

#define TARGET_COUNT				4
#define TARGET_IND_RECTO			0
#define TARGET_IND_VERSO			1
#define TARGET_IND_BG				2
#define TARGET_IND_RECTO_VERSO		3

class BCAObjectiveFuncData
{
	public:
		BCAObjectiveFuncData();
		~BCAObjectiveFuncData();
		
		// Accessors
		VectorT * operator [] (int x) 			{ return data[x]; }
		unsigned int size () 					{ return data.size(); }
		unsigned int noTargets()				{ return targets.size(); }
		
	public: // DATA
		
		unsigned int dim, dim_theta;
		vector<VectorT *> targets;
		vector<VectorT *> data;
		vector<ScalarT> target_norms;
};

/***********************************************************************
 * The constructor
 ***********************************************************************/

BCAObjectiveFuncData::BCAObjectiveFuncData()
{
	dim=3;
	dim_theta=dim*dim+dim-1;	
	
#warning NO ESTIMATION OF WEIGHTS!!!	
	dim_theta=dim*dim;
	
	targets.resize(TARGET_COUNT);
	targets[TARGET_IND_RECTO] = 		new VectorT (1, 0, 0);	
	targets[TARGET_IND_VERSO] = 		new VectorT (0, 1, 0);	
	targets[TARGET_IND_BG] = 			new VectorT (0, 0, 1);	
	targets[TARGET_IND_RECTO_VERSO] = 	new VectorT (1, 1, 0);			
	
	// Precalculated the norms of the targets
	for (unsigned int i=0; i<noTargets(); ++i)
		target_norms.push_back(targets[i]->norm());
}

/***********************************************************************
 * The destructor
 ***********************************************************************/

BCAObjectiveFuncData::~BCAObjectiveFuncData()
{
	for (vector<VectorT *>::iterator iter=data.begin(); iter!=data.end(); ++iter)
		delete *iter;
	for (unsigned int i=0; i<targets.size(); ++i)
		delete targets[i];		
}

#endif
