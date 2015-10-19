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
 * Find the nearest canonically oriented Orthorhombic medium
 * to an elastic stiffness matrix.
 *
 * Input:
 * cc1 is an arbitrary elastic matrix.
 *
 * Output:
 * cc2 is then the nearest Orthorhombic matrix with X,Y,Z principal axes.
 *
 * Return value:
 * The distance between cc1 and cc2 is returned.
 *
 * Author Joe Dellinger, Amoco TTC, 19 Feb 1997.
 */

FLT_DBL
ortho_distance (FLT_DBL * cc2, FLT_DBL * cc1)
{
int             ii, jj, kk, ll;
FLT_DBL         c11, c12, c13, c22, c23, c33, c44, c55, c66;
double          temp1, temp2;


/*
 * Find the nearest Orthorhombic medium. For the Orthorhombic case,
 * it's trivial: you simply throw away the elastic constants that
 * "should be zero". Here we pick out the elastic constants we
 * aren't going to throw away from the input stiffness matrix.
 */
    c11 = CCF1 (1, 1);
    c12 = CCF1 (1, 2);
    c13 = CCF1 (1, 3);
    c22 = CCF1 (2, 2);
    c23 = CCF1 (2, 3);
    c33 = CCF1 (3, 3);
    c44 = CCF1 (4, 4);
    c55 = CCF1 (5, 5);
    c66 = CCF1 (6, 6);


/*
 * Fill out the complete Voigt stiffness matrix
 * of the nearest Orthorhombic medium. First,
 * we zero out everything.
 */
    for (ii = 0; ii < 6; ii++)
	for (jj = 0; jj < 6; jj++)
	    CC2 (ii, jj) = 0.;

/*
 * Then, we set all the ones that are not always zero.
 */
    CCF2 (1, 1) = c11;
    CCF2 (1, 2) = CCF2 (2, 1) = c12;
    CCF2 (1, 3) = CCF2 (3, 1) = c13;
    CCF2 (2, 2) = c22;
    CCF2 (2, 3) = CCF2 (3, 2) = c23;
    CCF2 (3, 3) = c33;
    CCF2 (4, 4) = c44;
    CCF2 (5, 5) = c55;
    CCF2 (6, 6) = c66;

/*
 * Find the distance between the input and that matrix.
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
