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
#include <stdio.h>

/*
 * Read a stiffness matrix from standard input.
 *
 * Output:
 * 	cc is a 6x6 elastic stiffness matrix read from standard input.
 *
 * Author Joe Dellinger, Amoco TTC, 19 Feb 1997.
 */

void
read_matrix_6x6 (FLT_DBL * cc)
{
int             ii, jj;

    for (ii = 0; ii < 6; ii++)
	for (jj = 0; jj < 6; jj++)
	{
#ifdef DOUBLE_PRECISION
	    scanf ("%lf", &CC (ii, jj));
#else
	    scanf ("%f", &CC (ii, jj));
#endif
	}

    return;
}
