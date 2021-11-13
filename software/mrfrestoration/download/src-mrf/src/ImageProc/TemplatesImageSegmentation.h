/**********************************************************
 * ImageSegmentation.h
 * Routines for Image Segmentation
 * Templates
 *
 * Christian Wolf
 * Beginn: 30.5.2005
 **********************************************************/
 
// From the MATH module 
#include <Vector.h> 
#include <KMeansClusterer.h> 

// From this module
#include "ImageProc.h"
#include "color.h"

/*******************************************************
 * Use the K-Means Algorithm
 * If the last argument is NULL, then it must be cast
 * since the argument is a template:
 * (vector<PixelType> *) NULL
 *
 * if out_segmented_visualize is not NULL, then
 * a pointer to an image is returned which contains
 * the segmented image with pixel values being set to the 
 * cluster centers.
 *
 * Works only if the parameter T is an image (as opposed
 * to a float image), else it must be set to NULL!
 *******************************************************/
 
// A possible custom color distance replacing the euclidean
// distance between two colors: The CIE94 color distance
class LabCie94Distance : public ColorDistanceFunctional
{
	public:
		float distance (Vector<float> &x, Vector<float> &y)
		{
			return diffLabCIE94 (x[0], x[1], x[2], y[0], y[1], y[2]);
		}
};
 
template <class T>
void segmentKMeans(T &im, Image &o, unsigned int noClasses, ColorSpaceCode colspc, 
	vector<float> *initValues, 
	Image **out_segmented_visualize)
{		
	unsigned int dim=im.nbColorPlanes();
	unsigned int xs=im.xsize;
	unsigned int ys=im.ysize;	
	unsigned int pixnum;	
	double pl, pu, pv;
	Vector<float> v(dim);		
	KMeansClusterer<float> kmc(dim, noClasses);
	LabCie94Distance *customDistance;
	
	if ((dim!=3) && (colspc!=COLSPC_RGB))
		ERR_THROW ("Incompatible colorspace and dimension in segmentKMeans()\n");
		
	// Collect the pixels and put them into the
	// K-Means clusterer
	switch (colspc)
	{
		case COLSPC_RGB:
			for (unsigned int y=0; y<ys; ++y)
			for (unsigned int x=0; x<xs; ++x)
			{
				for (unsigned int d=0; d<dim; ++d)
					v.set(d, im.get(d+1,x,y));
				kmc.add (v);
			}
			break;
			
		case COLSPC_LUV:
			for (unsigned int y=0; y<ys; ++y)
			for (unsigned int x=0; x<xs; ++x)
			{
				rgb2luv (im.get(PLANE_RED,x,y)/255., 
						 im.get(PLANE_GREEN,x,y)/255., 
						 im.get(PLANE_BLUE,x,y)/255.,
						 &pl, &pu, &pv);
				v.set(0, pl);
				v.set(1, pu);
				v.set(2, pv);
				kmc.add (v);
				
			}
			break;	
			
		case COLSPC_HSV:
			for (unsigned int y=0; y<ys; ++y)
			for (unsigned int x=0; x<xs; ++x)
			{
				rgb2hsv (im.get(PLANE_RED,x,y)/255., 
						 im.get(PLANE_GREEN,x,y)/255., 
						 im.get(PLANE_BLUE,x,y)/255.,
						 &pl, &pu, &pv);
				pl = pl*255./360.;
				pu = pu*255.;
				pv = pv*255.;
				v.set(0, pl);
				v.set(1, pu);
				v.set(2, pv);		
				kmc.add (v);
			}
			break;
			
		case COLSPC_LAB:
		case COLSPC_LAB_CIE94:
			for (unsigned int y=0; y<ys; ++y)
			for (unsigned int x=0; x<xs; ++x)
			{
				rgb2lab (im.get(PLANE_RED,x,y)/255., 
						 im.get(PLANE_GREEN,x,y)/255., 
						 im.get(PLANE_BLUE,x,y)/255.,
						 &pl, &pu, &pv);
				v.set(0, pl);
				v.set(1, pu);
				v.set(2, pv);
				kmc.add (v);				
			}
			break;
			
		default:
			ERR_THROW ("Unknown color space in segmentKMeans()\n");
	}
	
	// Do we use a color space with non euclidean distance?
	if (colspc==COLSPC_LAB_CIE94)
	{
		customDistance = new LabCie94Distance;
		kmc.setDistance(customDistance);	
	}
	
	// We initialize manually
	if (initValues!=NULL)
	{
		kmc.initZero();
		for (unsigned int i=0; i<noClasses; ++i)
			kmc.init(i,(*initValues)[i]);
	}
			
	// We initialize randomly
	else
		kmc.initRandom();
	
	kmc.doCluster();	
// 	kmc.printData();
// 	kmc.printMeans(); 	
	
	// Collect the results
	pixnum=0;
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		o.set(x,y, kmc[pixnum]);
		++pixnum;
	}
	
	// Output the cluster centers
	for (int i=0; i<dim; ++i)
	{
		double dr,dg,db;
		unsigned char r,g,b;
		Vector<float> center=kmc.getCenter(i);
		
		cerr << "CLUSTER-CENTER ";
		switch (colspc)
		{
			case COLSPC_RGB:
				r=(byte) TRIMMEDGRAY(rint(center[0]));
				g=(byte) TRIMMEDGRAY(rint(center[1]));
				b=(byte) TRIMMEDGRAY(rint(center[2]));
				break;
		
			case COLSPC_LUV:
				luv2rgbscaled(center[0], center[1], center[2], &r,&g,&b);
				cerr << "LUV: " << center[0] << " " << center[1] << " " << center[2] << " ";
				break;
				
			case COLSPC_HSV:			
				center[0] = center[0]*360./255.;	
				center[1] = center[1]/255.;
				center[2] = center[2]/255.;
				hsv2rgbscaled(center[0], center[1], center[2], &r,&g,&b);						
				cerr << "HSV: " << center[0] << " " << center[1] << " " << center[2] << " ";
				break;
				
			case COLSPC_LAB:
			case COLSPC_LAB_CIE94:
				lab2rgb(center[0], center[1], center[2], &dr,&dg,&db);
				r = (byte) TRIMMEDGRAY(rint(255.*dr));
				g = (byte) TRIMMEDGRAY(rint(255.*dg));
				b = (byte) TRIMMEDGRAY(rint(255.*db));
				cerr << "LUV: " << center[0] << " " << center[1] << " " << center[2] << " ";
				break;
		}	
		
		cerr << "RGB: " << (int) r << " " << (int) g << " " << (int) b << endl;
	}
	
	if (out_segmented_visualize!=NULL)
	{
		Image *rv=new Image (im.xsize, im.ysize, im.nbColorPlanes());
		Vector<float> center;
		
		for (unsigned int y=0; y<ys; ++y)
		for (unsigned int x=0; x<xs; ++x)
		{
			double dr,dg,db;
			unsigned char r,g,b;
			center=kmc.getCenter(o.get(x,y));
			
			switch (colspc)
			{
				case COLSPC_RGB:
					r=(byte) TRIMMEDGRAY(rint(center[0]));
					g=(byte) TRIMMEDGRAY(rint(center[1]));
					b=(byte) TRIMMEDGRAY(rint(center[2]));
					break;
			
				case COLSPC_LUV:
					luv2rgbscaled(center[0], center[1], center[2], &r,&g,&b);
					break;
					
				case COLSPC_HSV:
					center[0] = center[0]*360./255.;	
					center[1] = center[1]/255.;
					center[2] = center[2]/255.;
					hsv2rgbscaled(center[0], center[1], center[2], &r,&g,&b);				
					break;
					
				case COLSPC_LAB:
				case COLSPC_LAB_CIE94:
					lab2rgb(center[0], center[1], center[2], &dr,&dg,&db);
					r = (byte) TRIMMEDGRAY(rint(255.*dr));
					g = (byte) TRIMMEDGRAY(rint(255.*dg));
					b = (byte) TRIMMEDGRAY(rint(255.*db));
					break;
			}
			
			rv->set(PLANE_RED, x, y, r);
			if (im.nbColorPlanes()>1)
			{
				rv->set(PLANE_GREEN, x, y, g);
				rv->set(PLANE_BLUE, x, y, b);
			}
		}
		*out_segmented_visualize=rv;
	}
	
	// Clean up
	if (colspc==COLSPC_LAB_CIE94)
		delete customDistance;
}

/*******************************************************
 * Use the K-Means Algorithm, the _PAIRWISE_ version:
 *
 * One data vector of the K-Means Algorithm is a pair
 * of neighboring pixels.
 *
 * If the last argument is NULL, then it must be cast
 * since the argument is a template:
 * (vector<PixelType> *) NULL
 *
 * if out_segmented_visualize is not NULL, then
 * a pointer to an image is returned which contains
 * the segmented image with pixel values being set to the 
 * cluster centers.
 *
 * Works only if the parameter T is an image (as opposed
 * to a float image), else it must be set to NULL!
 *******************************************************/

class PairPositions
{
	public:
		vector<unsigned int> pp;
};

class LabelCounter
{
	public:
		LabelCounter(unsigned int noClasses)	{ cnt.resize(noClasses); }
		void incLabel(unsigned int l)			{ ++cnt[l]; }
		void clear()	
		{ 
			for (unsigned int i=0; i<cnt.size(); ++i) cnt[i]=0; 
		}
		unsigned int getMaxLabel() 
		{
			unsigned int max=cnt[0];
			unsigned int argmax=0;
			for (unsigned int i=1; i<cnt.size(); ++i)
				if (cnt[i]>max)
				{
					max=cnt[i];
					argmax=i;
				}
			return argmax;
		}	
	
	private:
		vector<unsigned int> cnt;
};


template <class T>
void segmentKMeansPairWise(T &i, Image &o, unsigned int noClasses, 
	vector<float> *initValues, 
	Image **out_segmented_visualize)
{		
	typedef typename T::PixelType PixelType;
	unsigned int dims=i.nbColorPlanes();
	unsigned int dimp=2*dims;
	unsigned int xs=i.xsize;
	unsigned int ys=i.ysize;	
	unsigned int pixnum;	
	Vector<float> v(dimp);			
	KMeansClusterer<PixelType> kmc(dimp, noClasses);
	Matrix<PairPositions> M(ys,xs);
		
	// Collect the pixel pairs and put them into the
	// K-Means clusterer
	pixnum=0;
	for (unsigned int y=0; y<ys; ++y)
	{
		cerr << "[" << y << "]"; cerr.flush();
	for (unsigned int x=0; x<xs; ++x)
	{
		// cerr << "\nPIXEL " << x << "," << y << endl;	
		for (unsigned int wy=y-1; wy<=y+1; ++wy)
		for (unsigned int wx=x-1; wx<=x+1; ++wx)
		{					
			if ((wy>=0) && (wx>=0) && (wy<ys) && (wx<xs) && !((wx==x) && (wy==y)))
			{
				// cerr << "(" << wx << "," << wy << ")";
				
				for (unsigned int d=0; d<dims; ++d)
					v.set(     d, i.get(d+1,x,y));
				for (unsigned int d=0; d<dims; ++d)
					v.set(dims+d, i.get(d+1,wx,wy));
				kmc.add (v);
				M(y,x).pp.push_back(pixnum);
				++pixnum;
			}
		}
	}
	}
	
	// We initialize manually
	if (initValues!=NULL)
	{
		kmc.initZero();
		for (unsigned int i=0; i<noClasses; ++i)
			kmc.init(i,(*initValues)[i]);
	}
			
	// We initialize randomly
	else
		kmc.initRandom();
	
	kmc.doCluster();
	//kmc.printMeans(); 	
	
	// Collect the results	
	LabelCounter lc(noClasses);
	for (unsigned int y=0; y<ys; ++y)
	for (unsigned int x=0; x<xs; ++x)
	{
		vector<unsigned int> &pair_positions = M(y,x).pp;
		
		lc.clear();
		for (unsigned int i=0; i<pair_positions.size(); ++i)
			lc.incLabel(kmc[pair_positions[i]]);
			
		o.set(x,y, lc.getMaxLabel());
	}
		
	if (out_segmented_visualize!=NULL)
	{
		Image *rv=new Image (i.xsize, i.ysize, i.nbColorPlanes());
		Vector<float> center;
		
		for (unsigned int y=0; y<ys; ++y)
		for (unsigned int x=0; x<xs; ++x)
		{
			center=kmc.getCenter(o.get(x,y));
			rv->set(PLANE_RED, x, y, (byte) TRIMMEDGRAY(rint(center[0])));
			if (i.nbColorPlanes()>1)
			{
				rv->set(PLANE_GREEN, x, y, (byte) TRIMMEDGRAY(rint(center[1])));
				rv->set(PLANE_BLUE, x, y, (byte) TRIMMEDGRAY(rint(center[2])));
			}
		}
		*out_segmented_visualize=rv;
	}
}
