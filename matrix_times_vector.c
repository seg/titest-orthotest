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
 * Multiply a matrix times a vector:
 * outvec = rmat * invec
 *
 * Input: 3-vector invec
 *        3x3 matrix rmat
 * Output: 3-vector outvec
 *
 * Author Joe Dellinger, Amoco TTC, 19 Feb 1997.
 */

void
matrix_times_vector (FLT_DBL * outvec, FLT_DBL * rmat, FLT_DBL * invec)
{
int             ii, jj;
double          temp;

    for (ii = 0; ii < 3; ii++)
    {
	temp = 0.;

	for (jj = 0; jj < 3; jj++)
	{
	    temp += RMAT (jj, ii) * invec[jj];
	}

	outvec[ii] = temp;
    }

    return;
}
