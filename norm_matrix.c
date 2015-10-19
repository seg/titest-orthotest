/*
 * Copyright (c) 2005 by the Society of Exploration Geophysicists.
 * For more information, go to http://software.seg.org/2005/0001 .
 * You must read and accept usage terms at:
 * http://software.seg.org/disclaimer.txt before use.
 *
 * Revision history:
 * Original SEG version by Joe Dellinger, BP EPTG, July 2005.
 */

#include "cmat.h"
#include <math.h>

/*
 * Calculate the Federov norm of a stiffness matrix.
 * (Federov, Theory of Elastic Waves in Crystals, 1968, Plenum Press, New York.)
 *
 * Input:
 * cc1 is an arbitrary 6x6 elastic stiffness matrix in Voigt notation.
 *
 * Return value:
 * Federov's 3x3x3x3 tensor norm is returned.
 *
 * Author Joe Dellinger, Amoco TTC, 19 Feb 1997.
 */

FLT_DBL
norm_matrix_6x6 (FLT_DBL * cc1)
{
int             ii, jj, kk, ll;
double          temp1, temp2;

    temp2 = 0.;
    for (ii = 0; ii < 3; ii++)
	for (jj = 0; jj < 3; jj++)
	    for (kk = 0; kk < 3; kk++)
		for (ll = 0; ll < 3; ll++)
		{
		    temp1 = CCT1 (ii, jj, kk, ll);
		    temp2 += temp1 * temp1;
		}

    return (FLT_DBL) sqrt (temp2);
}
