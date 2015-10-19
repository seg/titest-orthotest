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
 * Express an orientation axis, as given by a unit vector,
 * in spherical coordinates.
 *
 * Input:
 * 	vec is an input (X,Y,Z) vector.
 *
 * Output:
 * 	phi and theta give the spherical coordinates of the direction
 * 	it points in, in _DEGREES_.
 *
 * phi=0 is the Z axis
 * phi=90 theta=0 is the X axis
 * phi=90 theta=90 is the Y axis
 *
 * Author Joe Dellinger, Amoco TTC, 19 Feb 1997.
 */

void
vector_to_angles (FLT_DBL vec[3], FLT_DBL * phi, FLT_DBL * theta)
{

    *theta = atan2 (vec[0], vec[1]) / DEGTORAD;
    *phi =
     atan2 (sqrt (vec[0] * vec[0] + vec[1] * vec[1]), vec[2]) / DEGTORAD;

    return;
}
