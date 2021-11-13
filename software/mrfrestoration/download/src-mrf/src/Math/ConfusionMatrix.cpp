/***********************************************************************
 * Confusion matrix
 * or: how to evaluate
 *
 * Author: Christian Wolf, christian.wolf@insa-lyon.fr
 * Begin: 14.6.2006
 ***********************************************************************/
 
// C++
#include <map> 
#include <set>
#include <vector>
  
// From this module
#include "ConfusionMatrix.h"

/**************************************************************
 * Add another matrix
 * See research notebook 20.9.2006, page 87
 **************************************************************/

void ConfusionMatrix::operator += (ConfusionMatrix &o)
{
	unsigned char gs,gg;
	int gs_index, gg_index=0;
	Matrix<float> tmp;
	
	// Travers all entries of the first matrix
	for (unsigned int y=0; y<rows(); ++y)	
	for (unsigned int x=0; x<columns(); ++x)	
	{
		// Search the gray values in the legend vectors
		// of the second matrix
		gs=gv_seg[y];
		gg=gv_gt[x];
		gs_index=-1;
		for (unsigned int i=0; i<o.rows(); ++i)
		{
			if (gs==o.get_gv_seg(i))
			{
				for (unsigned int j=0; j<o.columns(); ++j)
				{
					// Found both gray values
					if (gg==o.get_gv_gt(j))
					{
						gs_index=i;
						gg_index=j;
					}
				}
			}
		}		
		
		if (gs_index!=-1)
		{	
			mat(y,x) += o(gs_index,gg_index);
			o(gs_index,gg_index)=0;
		}
		else
			cerr << "Warning: adding two conf-matrices: gray value not found in (2).\n";
	}
		
	// Travers all entries of the second matrix 
	// and search for non zero entries
	for (unsigned int i=0; i<o.rows(); ++i)	
	for (unsigned int j=0; j<o.columns(); ++j)	
	{
		if (o(i,j)>0)
		{
			// Search whether this gray value exists
			gs=o.get_gv_seg(i);
			gs_index=-1;			
			for (unsigned int y=0; y<rows(); ++y)
			{
				if (gs==gv_seg[y])
					gs_index=y;
			}
			
			// Gray value not found
			// increase vector and matrix
			if (gs_index==-1)
			{	
				cerr << "Warning: adding two conf-matrices: gray value not found in (1a).\n"
					 << "         adding label " << (int) gs << " to result image.\n";	
				gv_seg.push_back(gs);
				gs_index=gv_seg.size()-1,
				mat.resizeNonDestructiveZero(gv_seg.size(), gv_gt.size());
			}
			
			// Search whether this gray value exists
			gg=o.get_gv_gt(j);
			gg_index=-1;			
			for (unsigned int x=0; x<columns(); ++x)
			{
				if (gg==gv_gt[x])
					gg_index=x;
			}
			
			// Gray value not found
			// increase vector and matrix
			if (gg_index==-1)
			{	
				cerr << "Warning: adding two conf-matrices: gray value not found in (1b).\n"
				     << "         adding label " <<  (int) gg << " to ground truth image.\n";	
				gv_gt.push_back(gg);
				gg_index=gv_gt.size()-1,
				mat.resizeNonDestructiveZero(gv_seg.size(), gv_gt.size());								
			}
					
			// Add the entry to first matrix
			mat(gs_index,gg_index)+=o(i,j);
			o(i,j)=0;
		}		
	}
		
	/* This code does not allow unequal matrices
	if ((gv_seg.size()!=o.gv_seg.size()) ||
	    (gv_gt.size() !=o.gv_gt.size()))
		ERR_THROW ("Confusion matrices are not compatible: "
			<< "#seg-labels="
			<< gv_seg.size() << "," << o.gv_seg.size() << "; "
			<< "#gt-labels="
	    	<< gv_gt.size() << "," << o.gv_gt.size());
	
	mat += o.mat;	
	*/	
}

/***********************************************************************
 * Search among all permutations of result labels the one
 * with highest accuracy
 * + Changes the confusion matrix
 * + Returns the permutation vector which has been applied: 
 *   each element transforms a column (a ground truth label) into 
 *   another column
 ***********************************************************************/ 

vector<int> * ConfusionMatrix::searchMostProbableLabelPermutation()
{
 	ConfusionMatrix *pCM, CMmin;
 	float acc, curAccuracy, maxAccuracy;
	int pos,cnt=mat.columns();
	int next_val;
	set<int> used_labels;
	vector<int> permutation(cnt), 
				last_tested(cnt),
				*minPerm;
	            	
	minPerm    = new vector<int> (cnt);
	curAccuracy=maxAccuracy=getAccuracy();	
	pos=0;
	for (int i=0; i<cnt; ++i)
	{
		(*minPerm)[i]=i;
		last_tested[i]=-1;
	}

	while (pos>=0)
	{			
		// get the next valid value for this position
		next_val=-1;
		for (int x=last_tested[pos]+1; x<cnt; ++x)
		{
			if (used_labels.find(x)==used_labels.end())
			{
				next_val=x;
				break;
			}
		}
		
		/*
		cerr << "pos=" << pos << ": ";
		for (int i=0; i<pos; ++i)
			cerr << " " << permutation[i];
		cerr << "; ";
		for (int i=0; i<cnt; ++i)
			cerr << " " << last_tested[i];
		cerr << " next=" << next_val 
			 << " #used=" << used_labels.size() << endl;
		*/
				
		if (next_val==-1)
		{
			used_labels.erase(permutation[pos-1]);
			--pos;
		}
		else
		{
			permutation[pos]=next_val;
			last_tested[pos]=next_val;
			for (int i=pos+1; i<cnt; ++i)
				last_tested[i]=-1;
			if (pos>=cnt-1)
			{
				// We found a valid permutation, test it				
				pCM = getPermutatedColumns(permutation);
				acc=pCM->getAccuracy();				
				if (acc>maxAccuracy)
				{
					maxAccuracy=acc;
					CMmin=*pCM;
					for (int i=0; i<cnt; ++i)
						(*minPerm)[i]=permutation[i];					
				}				
				delete pCM;
				
				/*
				cerr << "NEW PERMUTATION:";
				for (int i=0; i<=pos; ++i)
					cerr << " " << permutation[i];
 				cerr << "; ACC=" << acc << endl;
 				*/
				
				used_labels.erase(permutation[pos-1]);
				--pos;
			}
			else
			{
				used_labels.insert(next_val);
				++pos;
			}			
		}
	}
	
	if (maxAccuracy>curAccuracy)
		*this = CMmin;
			
	return minPerm;
}

/***********************************************************************
 * Apply a permutation to the columns of the matrix
 ***********************************************************************/ 

ConfusionMatrix * ConfusionMatrix::getPermutatedColumns (vector<int> permutation)
{
	ConfusionMatrix *rv=new ConfusionMatrix (*this);
	for (unsigned int y=0; y<rows(); ++y)
	for (unsigned int x=0; x<columns(); ++x)
		rv->mat(y,permutation[x]) = mat(y,x);		
	return rv;
}

/***********************************************************************
 * Get accuracy
 ***********************************************************************/ 

float ConfusionMatrix::getAccuracy()
{
	float corr=0, all=0;
	for (unsigned int y=0; y<mat.rows(); ++y)
	for (unsigned int x=0; x<mat.columns(); ++x)
	{
		if (x==y)
			corr += mat(y,x);
		all+=mat(y,x);
	}
	return corr / all;
}

/***********************************************************************
 * If one of the labels is designated as the object class label,
 * then precision and recall values can be calculated
 ***********************************************************************/

float ConfusionMatrix::getRecall(int objectClassNumber)
{
	float allgt=0;
	float corrdet=mat(objectClassNumber, objectClassNumber);
	for (unsigned int y=0; y<mat.rows(); ++y)
		allgt+=mat(y,objectClassNumber);
	if (allgt==NULL)
		return 1.;
	else
		return corrdet/allgt;
}

float ConfusionMatrix::getPrecision(int objectClassNumber)
{
	float alldet=0;
	float corrdet=mat(objectClassNumber, objectClassNumber);
	for (unsigned int x=0; x<mat.columns(); ++x)
		alldet+=mat(objectClassNumber,x);
	if (alldet==NULL)
		return 1.;
	else
		return corrdet/alldet;
}

/**************************************************************
 * Divide the matrix by a scalar
 **************************************************************/

void ConfusionMatrix::normalize () 
{			
	mat /= mat.sum();
	mat *= 100;
}

// *********************************************************************
// Output
// width and precision controle the formatting.
// if width<0 then there is no formatting.
// *********************************************************************

void ConfusionMatrix::print (ostream &os, bool latex, int width, int precision)
{
	if (!latex)
	{
		os << "CONFUSION MATRIX"<<endl;
		os << "labels in result: ";
		for (unsigned int i=0; i<gv_seg.size(); ++i) 
			os << (int) gv_seg[i] << " ";
		os << "\nlabels in ground truth: ";
		for (unsigned int i=0; i<gv_gt.size(); ++i) 
			os << (int) gv_gt[i] << " ";
		os << endl;		
	}
	
	if (latex)
		mat.LaTeXOut(os, width, precision);
	else
		os << mat;
}

ostream & operator << (ostream &os, ConfusionMatrix &m) 
{
	m.print (os, false);
	return os;
}
