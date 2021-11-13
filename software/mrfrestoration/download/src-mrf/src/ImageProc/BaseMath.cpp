// ***********************************************************
// ***********************************************************

#include		<stdio.h>
#include		<stdlib.h>
#include		<errno.h>

// **************************************************************

int	** MAT_INT(int t_l,	int	t_c) {
  short	int		i;
  int ** im;

  /*printf("Allocate %d	x %d\n",taille,taille) ;*/
  if ((im =	(int **) calloc( t_c , sizeof(int *))) == NULL)
   printf("Probleme	allocation\n") ;;
  for (	i =	0 ;	i <	t_c	; i	++)
   if ((im[i] =	(int *)	calloc(	t_l	, sizeof(int)))	== NULL)
   printf("Probleme	allocation\n") ;;

  fflush(stdout) ;
return (im);
}

// **************************************************************

float *	VECT_FLOAT(int l) {
  float	* im;

  /*printf("Allocate %d	\n",l) ;*/
   fflush(stdout) ;
  if ((im =	(float *) calloc( l	, sizeof(float *)))	== NULL)
   printf("Probleme	allocation\n") ;;

return (im);
}

void FREE_VECT_FLOAT (float	*f)	{
	free (f);
}

// **************************************************************

float ** MAT_FLOAT(int t_l,	int	t_c) {
  short	int		i;
  float	** im;

  /*printf("Allocate %d	x %d\n",t_c,t_l) ;
   fflush(stdout) ;*/
  if ((im =	(float **) calloc( t_c , sizeof(float *))) == NULL)
   printf("Probleme	allocation\n") ;;
   fflush(stdout) ;
  for (	i =	0 ;	i <	t_c	; i	++)
   if ((im[i] =	(float *) calloc( t_l ,	sizeof(float)))	== NULL)
   printf("Probleme	allocation\n") ;;

  fflush(stdout) ;
return (im);
}

void FREE_MAT_FLOAT	(int xsize,	float **M) {
   for (int	i =	0 ;	i <	xsize ;	i ++)
		free (M[i]);
	free (M);
}


// **************************************************************

float *** PRISM_FLOAT(int n, int t_l, int t_c) {

  short	int		i, j;
  float			***im;

  if ((im =	(float ***)	calloc(	n, sizeof(float	**))) == NULL)
   printf("Probleme	allocation\n") ;;
  for (j = 0 ; j < n ; j++)
  {	if ((im[j] = (float	**)	calloc(	t_c, sizeof(float *))) == NULL)
	 printf("Probleme allocation\n") ;;
	for	( i	= 0	; i	< t_c ;	i ++)
	 if	((im[j][i] = (float	*) calloc( t_l , sizeof(float))) ==	NULL)
	 printf("Probleme allocation\n") ;;
  }

return (im);
}

