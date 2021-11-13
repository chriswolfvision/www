/***********************************************************************
 * MRFLineProcesses.h
 * Different macros relevant to MRF line processes
 *
 * Author: Christian Wolf
 * Begin: 7.2.2007
 ***********************************************************************/
 
#ifndef _WOLF_MRFLINEPROCESSES_H_ 
#define _WOLF_MRFLINEPROCESSES_H_

// Not all elements of the line process matrix are actually 
// line process variables. Some of them are part of the 
// intensity process and others are just virtual labels.
// These macros tells us who is who
#define IS_ELEMENT_OF_MRFLP(x,y)	(((x)%2)!=((y)%2))
#define IS_HORIZ_EDGE_OF_MRFLP(x,y)	((((x)%2)==0)&&(((y)%2)==1))
#define IS_VERT_EDGE_OF_MRFLP(x,y)	((((x)%2)==1)&&(((y)%2)==0))
#define IS_VIRTUAL_LABEL_MRFLP(x,y)	(((x)%2==1)&&((y)%2==1))

#define CL_MASK_NS			(1+2+4+8 +32+64+128+256)

#define CL_MASK_LP_4NO		(2^12)
#define CL_MASK_LP_4SO		(2^17)
#define CL_MASK_LP_4WE		(2^14)
#define CL_MASK_LP_4EA		(2^15)
#define CL_MASK_LP_4NW_N	(2^9)
#define CL_MASK_LP_4NW_W	(2^11)
#define CL_MASK_LP_4NW_E	(2^12)
#define CL_MASK_LP_4NW_S	(2^14)
#define CL_MASK_LP_4NE_N	(2^10)
#define CL_MASK_LP_4NE_W	(2^12)
#define CL_MASK_LP_4NE_E	(2^13)
#define CL_MASK_LP_4NE_S	(2^15)
#define CL_MASK_LP_4SW_N	(2^14)
#define CL_MASK_LP_4SW_W	(2^16)
#define CL_MASK_LP_4SW_E	(2^17)
#define CL_MASK_LP_4SW_S	(2^19)
#define CL_MASK_LP_4SE_N	(2^15)
#define CL_MASK_LP_4SE_W	(2^17)
#define CL_MASK_LP_4SE_E	(2^18)
#define CL_MASK_LP_4SE_S	(2^20)


#endif
