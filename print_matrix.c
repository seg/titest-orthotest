/*
 * Copyright (c) 2005 by the Society of Exploration Geophysicists.
 * For more information, go to http://software.seg.org/2005/0001 .
 * You must read and accept usage terms at:
 * http://software.seg.org/disclaimer.txt before use.
 *
 * Revision history:
 * Original SEG version by Joe Dellinger, BP EPTG, July 2005.
 */

#include <stdio.h>
#include "cmat.h"

/*
 * Print out a 6x6 elastic matrix
 *
 * Input: 6x6 elastic matrix cc
 *
 * Author Joe Dellinger, Amoco TTC, 19 Feb 1997.
 */

void
print_matrix_6x6 (FLT_DBL * cc)
{
int             ii, jj;

    for (ii = 0; ii < 6; ii++)
    {
	for (jj = 0; jj < 6; jj++)
	{
	    printf ("%11.4g ", (float) CC (ii, jj));
	}
	printf ("\n");
    }

    return;
}

/*
 * Same as print_matrix_6x6, but allow the format to be specified.
 *
 * Input: 6x6 elastic matrix cc
 *        printf-style format string format
 */

void
format_print_matrix_6x6 (char *format, FLT_DBL * cc)
{
int             ii, jj;

    for (ii = 0; ii < 6; ii++)
    {
	for (jj = 0; jj < 6; jj++)
	{
	    printf (format, (float) CC (ii, jj));
	}
	printf ("\n");
    }

    return;
}

/*
 * Print out a 3x3 matrix
 *
 * Input: 3x3 elastic matrix rmat
 */

void
print_matrix_3x3 (FLT_DBL * rmat)
{
int             ii, jj;

    for (ii = 0; ii < 3; ii++)
    {
	for (jj = 0; jj < 3; jj++)
	{
	    printf ("%12.6g ", (float) RMAT (ii, jj));
	}
	printf ("\n");
    }

    return;
}
