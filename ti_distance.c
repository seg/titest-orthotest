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
 * Find the nearest VTI medium to an elastic stiffness matrix.
 *
 * Input:
 * 	cc1 is an input 6x6 elastic stiffness matrix in Voigt notation.
 *
 * Output:
 * 	cc2 is then the nearest VTI matrix (transversely isotropic with a
 * 	Z axis of symmetry) to cc1.
 *
 * Return value: the distance between the VTI matrix cc2 and
 *      the input matrix cc1.
 *
 * Author Joe Dellinger, Amoco TTC, 19 Feb 1997.
 */

FLT_DBL
ti_distance (FLT_DBL * cc2, FLT_DBL * cc1)
{
int             ii, jj, kk, ll;
FLT_DBL         c33, c13, c55, c11, c66;
double          temp1, temp2;


/*
 * Find the nearest VTI medium.
 * Watch out! For some reason, this equation has been particularly
 * susceptible to typos in the literature.
 * This is a minimal set of 5 independent stiffness constants required
 * to define a VTI medium.
 */
    c33 = CCF1 (3, 3);

    c13 = (CCF1 (1, 3) + CCF1 (2, 3)) / 2.;

    c55 = (CCF1 (4, 4) + CCF1 (5, 5)) / 2.;

    c11 =
     (3. * CCF1 (1, 1) + 3. * CCF1 (2, 2) + 4. * CCF1 (6, 6) +
      2. * CCF1 (1, 2)) / 8.;

    c66 =
     (CCF1 (1, 1) + CCF1 (2, 2) + 4. * CCF1 (6, 6) - 2. * CCF1 (1, 2)) / 8.;


/*
 * Fill out the complete Voigt-form VTI matrix from those 5.
 * Most elements are zero.
 */
    for (ii = 0; ii < 6; ii++)
	for (jj = 0; jj < 6; jj++)
	    CC2 (ii, jj) = 0.;

/* The non-zero elements */

    CCF2 (1, 1) = c11;
    CCF2 (2, 2) = c11;

    CCF2 (3, 3) = c33;

    CCF2 (4, 4) = c55;
    CCF2 (5, 5) = c55;

    CCF2 (6, 6) = c66;

    CCF2 (1, 3) = CCF2 (3, 1) = c13;
    CCF2 (2, 3) = CCF2 (3, 2) = c13;

    CCF2 (1, 2) = CCF2 (2, 1) = c11 - 2. * c66;


/*
 * Now that we have the best-approximating VTI matrix,
 * find the distance between it and the input matrix.
 * This is more easily done using tensor (hence "T") notation.
 * The 6x6 array CCF1 and the 3x3x3x3 array CCT1
 * index into the same memory. Ditto for CCF2 and CCT2.
 * This is a straightforward implementation of Federov's distance formula.
 */

    temp2 = 0.;
    for (ii = 0; ii < 3; ii++)
	for (jj = 0; jj < 3; jj++)
	    for (kk = 0; kk < 3; kk++)
		for (ll = 0; ll < 3; ll++)
		{
		    temp1 = (CCT2 (ii, jj, kk, ll) - CCT1 (ii, jj, kk, ll));
		    temp2 += temp1 * temp1;
		}

    return (FLT_DBL) sqrt (temp2);
}
