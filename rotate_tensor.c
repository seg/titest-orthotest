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

/*
 * This routine calculates the tensor product of the rotation matrix
 * rmat and the stiffness tensor cc1 and returns the result in the
 * stiffness tensor cc2.
 *
 * cc2 = rmat * cc1
 *
 * Input:
 *	cc1 is a 3x3x3x3 tensor in compressed 6x6 Voigt notation.
 *	rmat is a 3x3 rotation matrix.
 *
 * Output:
 *	cc2 is a 3x3x3x3 tensor in compressed 6x6 Voigt notation.
 *
 * Author Joe Dellinger, Amoco TTC, 19 Feb 1997.
 */


/*
 * The "*_voigt" arrays are used so that the same stiffness matrix can be
 * referred to using either 6x6 Voigt notation or 3x3x3x3 tensor notation,
 * whichever is more convenient.
 *
 * Macros are provided to take care of the details:
 *
 * Tensor notation has a "T" and 4 subscripts:
 *                   CCT1 (ii, jj, kk, ll)
 * Voigt notation has no "T" and 2 subscripts:
 *                   CC1  (ii, jj)
 *
 * Both are macros indexing into the same 1-D C array "cc1".
 */


/* This array converts from tensor notation to compressed Voigt notation */
int             extern_voigt[3][3] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};

/*
 * These two arrays convert from compressed Voigt notation to tensor notation.
 * static_voigtl is for the left subscript, static_voigtr the right.
 */
static int      static_voigtl[6] = {0, 1, 2, 1, 0, 0};
static int      static_voigtr[6] = {0, 1, 2, 2, 2, 1};


void
rotate_tensor (FLT_DBL * cc2, FLT_DBL * cc1, FLT_DBL * rmat)
{
int             ij, kl;
int             pp, qq, rr, ss;
int             ii, jj, kk, ll;
double          temp;

/*
 * How to rotate a tensor:
 *
 * Cijkl = Rip Rjq Rkr Rls Cpqrs            (Einstein summation convention)
 *
 * Where R is a rotation matrix.
 *
 */


    /*
     * ij and kl are in Voigt notation
     */
    for (ij = 0; ij < 6; ij++)
	for (kl = 0; kl <= ij; kl++)
	{
	    temp = 0.;

	    /* Uncompress ij to ii and jj, kl to kk and ll */
	    ii = static_voigtl[ij];
	    jj = static_voigtr[ij];
	    kk = static_voigtl[kl];
	    ll = static_voigtr[kl];

	    /*
	     * Einstein summation over tensor indices p, q, r, s.
	     */
	    for (pp = 0; pp < 3; pp++)
		for (qq = 0; qq < 3; qq++)
		    for (rr = 0; rr < 3; rr++)
			for (ss = 0; ss < 3; ss++)
			{
			    temp +=
			     RMAT (pp, ii) *
			     RMAT (qq, jj) *
			     RMAT (rr, kk) *
			     RMAT (ss, ll) * CCT1 (pp, qq, rr, ss);
			}

	    CC2 (kl, ij) = CC2 (ij, kl) = temp;
	}

    return;
}
