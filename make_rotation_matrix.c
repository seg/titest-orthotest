/*
 * Copyright (c) 2005 by the Society of Exploration Geophysicists.
 * For more information, go to http://software.seg.org/2005/0001 .
 * You must read and accept usage terms at:
 * http://software.seg.org/disclaimer.txt before use.
 *
 * Revision history:
 * Original SEG version by Joe Dellinger, BP EPTG, July 2005.
 */

#include <math.h>
#include "cmat.h"

/*
 * Create a rotation matrix specified by Euler angles.
 *
 * Input: degz1, degx, degz2 are the input Euler angles
 *     Rotate first by degz1 DEGREES counterclockwise about the +Z axis,
 *     then by degx DEGREES counterclockwise about the +X axis,
 *     then by degz2 DEGREES counterclockwise about the +Z axis.
 *
 * On Output:
 *     rmat is the 3x3 rotation matrix that applies this rotation.
 *
 * Author Joe Dellinger, Amoco TTC, 19 Feb 1997.
 */

void
make_rotation_matrix (FLT_DBL degz1, FLT_DBL degx, FLT_DBL degz2, FLT_DBL * rmat)
{
int             ii, jj;
int             iii, jjj, kkk;
FLT_DBL         rmat1[9], rmat2[9], rmat3[9];

    degz1 *= DEGTORAD;
    degx *= DEGTORAD;
    degz2 *= DEGTORAD;


    for (ii = 0; ii < 3; ii++)
	for (jj = 0; jj < 3; jj++)
	{
	    RMAT (ii, jj) = 0.;
	}

/*
 * Start with the identity
 */
    for (ii = 0; ii < 3; ii++)
	RMAT (ii, ii) = 1.;

/*
 * Accumulate three rotations
 */
/*
 * Rotate first by degz1 DEGREES counterclockwise about the +Z axis.
 */
    iii = 0;
    jjj = 1;
    kkk = 2;
    for (ii = 0; ii < 3; ii++)
	for (jj = 0; jj < 3; jj++)
	{
	    RMAT1 (ii, jj) = 0.;
	}
    RMAT1 (iii, iii) = cos (degz1);
    RMAT1 (jjj, jjj) = cos (degz1);
    RMAT1 (iii, jjj) = sin (degz1);
    RMAT1 (jjj, iii) = -sin (degz1);
    RMAT1 (kkk, kkk) = 1.;

    matrix_times_matrix (rmat2, rmat1, rmat);


/*
 * Rotate by degx DEGREES counterclockwise about the +X axis.
 */
    iii = 1;
    jjj = 2;
    kkk = 0;
    for (ii = 0; ii < 3; ii++)
	for (jj = 0; jj < 3; jj++)
	{
	    RMAT1 (ii, jj) = 0.;
	}
    RMAT1 (iii, iii) = cos (degx);
    RMAT1 (jjj, jjj) = cos (degx);
    RMAT1 (iii, jjj) = sin (degx);
    RMAT1 (jjj, iii) = -sin (degx);
    RMAT1 (kkk, kkk) = 1.;

    matrix_times_matrix (rmat3, rmat1, rmat2);


/*
 * Rotate by degz2 DEGREES counterclockwise about the +Z axis.
 */
    iii = 0;
    jjj = 1;
    kkk = 2;
    for (ii = 0; ii < 3; ii++)
	for (jj = 0; jj < 3; jj++)
	{
	    RMAT1 (ii, jj) = 0.;
	}
    RMAT1 (iii, iii) = cos (degz2);
    RMAT1 (jjj, jjj) = cos (degz2);
    RMAT1 (iii, jjj) = sin (degz2);
    RMAT1 (jjj, iii) = -sin (degz2);
    RMAT1 (kkk, kkk) = 1.;

    matrix_times_matrix (rmat, rmat1, rmat3);

    return;
}

/*
 * Multiply one matrix by another:
 * rmat3 = rmat1 * rmat2
 *
 * Input: 3x3 matrices rmat1, rmat2
 *
 * Output: 3x3 matrix rmat3
 */
void
matrix_times_matrix (FLT_DBL * rmat3, FLT_DBL * rmat1, FLT_DBL * rmat2)
{
int             ii, jj, ll;
double          temp;

    for (ii = 0; ii < 3; ii++)
	for (jj = 0; jj < 3; jj++)
	{
	    temp = 0.;

	    for (ll = 0; ll < 3; ll++)
		temp += RMAT1 (ll, jj) * RMAT2 (ii, ll);

	    RMAT3 (ii, jj) = temp;
	}

    return;
}

/*
 * Transpose a 3x3 matrix:
 * rmat2 = transp(rmat1)
 *
 * Input: 3x3 matrix rmat1
 *
 * Output: 3x3 matrix rmat2
 */

void
transpose_matrix (FLT_DBL * rmat2, FLT_DBL * rmat1)
{
int             ii, jj;

    for (ii = 0; ii < 3; ii++)
	for (jj = 0; jj < 3; jj++)
	{
	    RMAT2 (ii, jj) = RMAT1 (jj, ii);
	}

    return;
}
