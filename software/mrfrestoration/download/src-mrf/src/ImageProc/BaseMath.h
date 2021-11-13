#ifndef	_WOLF_BASEMATH_H_
#define	_WOLF_BASEMATH_H_

int	** MAT_INT(int t_l,	int	t_c);
float *	VECT_FLOAT(int l);
float ** MAT_FLOAT(int t_l,	int	t_c);
float *** PRISM_FLOAT(int n,int	t_l,int	t_c);
unsigned char ** CREATE_IMAGE(int t_l, int t_c);

void FREE_VECT_FLOAT (float	*f);
void FREE_MAT_FLOAT	(int ysize,	float **M);

#endif

