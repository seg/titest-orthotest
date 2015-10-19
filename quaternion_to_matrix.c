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
 * Express a quaternion as a rotation matrix.
 *
 * qq is a quaternion (see EDN March 2, 1995 for more details).
 * "Quaternions quickly transform coordinates without error buildup",
 * author "Do-While Jones" (a pen name of R. David Pogge). The
 * article is available online at:
 * http://www.reed-electronics.com/ednmag/archives/1995/030295/05df3.htm
 * or
 * http://www.edn.com/archives/1995/030295/05df3.htm
 * This subroutine implements "Algorithm 1" from that paper.
 *
 * Input:
 * 	qq is a quaternion.
 *
 * Output:
 * 	rmat is the equivalent 3x3 rotation matrix.
 *
 * Author Joe Dellinger, Amoco TTC, 19 Feb 1997.
 */

void
quaternion_to_matrix (FLT_DBL * qq, FLT_DBL * rmat)
{
double          tx, ty, tz, tq, tk;

    tx = qq[1] * qq[1];
    ty = qq[2] * qq[2];
    tz = qq[3] * qq[3];
    tq = ty + tz;

    if (tq + tx + qq[0] * qq[0] != 0.)
	tk = 2. / (tq + tx + qq[0] * qq[0]);
    else
	tk = 0.;

    RMAT (0, 0) = 1. - tk * tq;
    RMAT (1, 1) = 1. - tk * (tx + tz);
    RMAT (2, 2) = 1. - tk * (tx + ty);
    tx = tk * qq[1];
    ty = tk * qq[2];
    tq = tk * qq[3] * qq[0];
    tk = tx * qq[2];
    RMAT (0, 1) = tk - tq;
    RMAT (1, 0) = tk + tq;
    tq = ty * qq[0];
    tk = tx * qq[3];
    RMAT (0, 2) = tk + tq;
    RMAT (2, 0) = tk - tq;
    tq = tx * qq[0];
    tk = ty * qq[3];
    RMAT (1, 2) = tk - tq;
    RMAT (2, 1) = tk + tq;

    return;
}
