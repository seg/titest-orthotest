/*
 * Copyright (c) 2005 by the Society of Exploration Geophysicists.
 * For more information, go to http://software.seg.org/2005/0001 .
 * You must read and accept usage terms at:
 * http://software.seg.org/disclaimer.txt before use.
 *
 * Revision history:
 * Original SEG version by Joe Dellinger, BP EPTG, July 2005.
 */

#ifndef INCLUDE_CMAT_H
#define INCLUDE_CMAT_H
/*
 * Master include file for titest and orthotest routines.
 */

/*
 * Comment this out to do calculations in
 * single precision instead.
 */
#define DOUBLE_PRECISION

#ifdef DOUBLE_PRECISION
#define FLT_DBL	double
#else
#define FLT_DBL	float
#endif

/* Pi / 180. */
#define DEGTORAD (3.14159265358979323846264338327950288419716939937511/180.)

/*
 * Table for transforming between 6x6 Voigt compressed notation
 * and 3x3x3x3 Tensor notation. See rotate_tensor.c.
 */
extern int      extern_voigt[3][3];

/*
 * Various stiffness arrays... we have room for 3 at a time.
 *
 * Each array (cc, cc1, or cc2) can be indexed in one of 3 different ways.
 *
 * CC(A,B)  uses "C" Voigt notation, with indices running from 0 to 5.
 * CCF(A,B) uses "Fortran" Voigt notation, with indices running from 1 to 6.
 * CCT(A,B,C,D) uses 3x3x3x3 tensor notation, with indices running from 0 to 2.
 *
 * All 3 notation methods index into the same 1D C array.
 */

#define CC(A,B)			cc[(B)+6*(A)]
#define CCF(A,B)		cc[((B)-1)+6*((A)-1)]
#define CCT(A,B,C,D)		CC( (extern_voigt[(A)][(B)]), (extern_voigt[(C)][(D)]) )

#define CC1(A,B)		cc1[(B)+6*(A)]
#define CCF1(A,B)		cc1[((B)-1)+6*((A)-1)]
#define CCT1(A,B,C,D)		CC1( (extern_voigt[(A)][(B)]) , (extern_voigt[(C)][(D)]) )

#define CC2(A,B)		cc2[(B)+6*(A)]
#define CCF2(A,B)		cc2[((B)-1)+6*((A)-1)]
#define CCT2(A,B,C,D)		CC2( (extern_voigt[(A)][(B)]) , (extern_voigt[(C)][(D)]) )


/*
 * 3x3 rotation matrices
 */
#define RMAT(A,B)	rmat[(B)+3*(A)]
#define RMAT1(A,B)	rmat1[(B)+3*(A)]
#define RMAT2(A,B)	rmat2[(B)+3*(A)]
#define RMAT3(A,B)	rmat3[(B)+3*(A)]

/*
 * Subroutines
 */
void            rotate_tensor (FLT_DBL *, FLT_DBL *, FLT_DBL *);
void            make_rotation_matrix (FLT_DBL, FLT_DBL, FLT_DBL, FLT_DBL *);
void            transpose_matrix (FLT_DBL *, FLT_DBL *);
void            quaternion_to_matrix (FLT_DBL *, FLT_DBL *);
void            matrix_times_vector (FLT_DBL *, FLT_DBL *, FLT_DBL *);
void            matrix_times_matrix (FLT_DBL *, FLT_DBL *, FLT_DBL *);
void            print_matrix_6x6 (FLT_DBL *);
void            format_print_matrix_6x6 (char *format, FLT_DBL *);
void            print_matrix_3x3 (FLT_DBL *);
void            read_matrix_6x6 (FLT_DBL *);
FLT_DBL         ti_distance (FLT_DBL *, FLT_DBL *);
FLT_DBL         ortho_distance (FLT_DBL *, FLT_DBL *);
FLT_DBL         norm_matrix_6x6 (FLT_DBL *);
void            vector_to_angles (FLT_DBL v[3], FLT_DBL *, FLT_DBL *);
FLT_DBL         find_ti (FLT_DBL * cc, FLT_DBL * theta_best, FLT_DBL * phi_best);
FLT_DBL         find_ortho (FLT_DBL * cc, FLT_DBL * rmat);

/*
 * Author Joe Dellinger, February 1997
 */
#endif			/* INCLUDE_CMAT_H */
