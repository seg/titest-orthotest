/*
 * Copyright (c) 2005 by the Society of Exploration Geophysicists.
 * For more information, go to http://software.seg.org/2005/0001 .
 * You must read and accept usage terms at:
 * http://software.seg.org/disclaimer.txt before use.
 *
 * Revision history:
 * Original SEG version by Joe Dellinger, BP EPTG, July 2005.
 */

/*
 * Given a set of elastic constants, the routine "ortho_distance" finds the
 * distance between that set of elastic constants and the nearest
 * Orthorhombic medium with the XY, XZ, and YZ planes as symmetry planes.
 *
 * We want to find the nearest Orthorhombic medium regardless of the
 * orientation of the symmetry planes. This subroutine scans over all possible
 * unique orientations, and finds a trial orthorhombic orientation that
 * produces a minimum distance. It then refines the search in that vicinity
 * until it converges to an optimal answer.
 *
 * On input:
 *	cc is a 6x6 array of Voigt-notation elastic stiffness constants.
 *
 * On output:
 *	rmat is the 3x3 array that will rotate the input stiffness tensor
 *	into the canonical coordinate system of the best-approximating
 *      orthorhombic medium, ie, with the XYZ axes in the orthorhombic symmetry
 *	planes. The ambiguity of which axes should map to "X", "Y", and "Z"
 *      is resolved by choosing "Z" to be the one of the three that best
 *      serves as a transversely isotropic symmetry axis. The "Y" axis is
 *      then the second-best choice of the three.
 *
 * Return value:
 *      The distance between the nearest transversely orthorhombic medium
 *      and the input medium, in absolute units (not normalized).
 */

#include <math.h>
#include "cmat.h"

/* How much to refine after each successive search */
#define SUBDIVIDE	5
/* How much to subdivide the rotation around the axis */
#define SUB_ROT		29
/* How much to subdivide the choice of axis */
#define SUB_POS		5

#ifdef DOUBLE_PRECISION
#define END_RES		(1.e-9)
#else
#define END_RES		(1.e-6)
#endif

#define NOT_SET_YET	1000.
#define NO_NORM		-1.

FLT_DBL
find_ortho (FLT_DBL * cc, FLT_DBL * rmat)
{
int             kk;
FLT_DBL         ccrot[6 * 6];
FLT_DBL         ccortho[6 * 6];
FLT_DBL         ccrot2[6 * 6];
FLT_DBL         ccti[6 * 6];
FLT_DBL         rmat_transp[9];
FLT_DBL         rmat_temp[9];
FLT_DBL         rmat_temp2[9];
FLT_DBL         vec[3];
FLT_DBL         vec2[3];
FLT_DBL         dist;
FLT_DBL         dista[3];
FLT_DBL         qq[4], qq_best[4];
double          center[4];
double          range[4];
int             count[4];
int             qindex[4];
double          inc[4];
FLT_DBL         phi, theta;
FLT_DBL         temp;
FLT_DBL         dist_best;


/*
 * Keep the compiler from complaining that this may be uninitialized.
 */
    dist_best = NO_NORM;

/*
 * Search over all possible orientations.
 *
 * Any orientation in 3-space can be specified by a unit vector,
 * giving an axis to rotate around, and an angle to rotate about
 * the given axis. The orientation is given with respect to some
 * fixed reference orientation.
 *
 * Since rotating by theta degrees about (A,B,C) produces the same
 * result as rotating -theta degrees about (-A,-B,-C), we only
 * need to consider 180 degrees worth of angles, not 360.
 *
 * In this application, we are finding the orientation of an orthorhombic
 * medium. Orthorhombic symmetry has three orthogonal symmetry planes,
 * so any one octant defines the whole. We thus only need to search
 * over rotation axes within one octant.
 *
 * Following the article in EDN, March 2, 1995, on page 95, author
 * "Do-While Jones" (a pen name of R. David Pogge),
 * "Quaternions quickly transform coordinates without error buildup",
 * we use quaternions to express the rotation. The article can be read
 * online here:
 * http://www.reed-electronics.com/ednmag/archives/1995/030295/05df3.htm
 *
 * If (A,B,C) is a unit vector to rotate theta degrees about, then:
 *
 * q0 = Cos (theta/2)
 * q1 = A * Sin(theta/2)
 * q2 = B * Sin(theta/2)
 * q3 = C * Sin(theta/2)
 *
 * so that q0^2 + q1^2 + q2^2 + q3^2 = 1. (A unit magnitude quaternion
 * represents a pure rotation, with no change in scale).
 *
 * For our case, taking advantage of the orthorhombic symmetry to
 * restrict the search space, we have:
 * 0 <= A <= 1
 * 0 <= B <= 1
 * 0 <= C <= 1
 * 0 <= theta <= 180 degrees.
 * The rotation axis direction is limited to within one octant,
 * and the rotation about that axis is limited to half of the full circle.
 *
 * In terms of quaternions, this bounds all four elements between 0 and 1,
 * inclusive.
 */

    /*
     * How much to subdivide each quaternion axis in the original scan. These
     * were somewhat arbitrarily chosen. These choices appear to be overkill,
     * but that ensures we won't accidentally miss the correct result by
     * insufficient sampling of the search space.
     */
    /*
     * We sample the rotation angle more finely than the rotation axis.
     */
    count[0] = SUB_ROT;
    count[1] = SUB_POS;
    count[2] = SUB_POS;
    count[3] = SUB_POS;

    /*
     * Between 0. and 1. for all 4 Q's   (That is .5 +- .5.)
     */
    for (kk = 0; kk < 4; kk++)
    {
	range[kk] = .5;
	center[kk] = .5;
	/*
	 * A number meaning "not set yet", to get us through the loop the
	 * first time. Needs to be much bigger than END_RES.
	 */
	inc[kk] = NOT_SET_YET;
    }


    while (inc[0] > END_RES && inc[1] > END_RES &&
	   inc[2] > END_RES && inc[3] > END_RES)
    {
	/*
	 * Update inc to reflect the increment for the current search
	 */
	for (kk = 0; kk < 4; kk++)
	{
	    inc[kk] = (2. * range[kk]) / (FLT_DBL) (count[kk] - 1);
	}

	/*
	 * Start the 4-dimensional search. Keep track of the best result
	 * found so far. The distance must be non-negative; we use -1 to mean
	 * "not set yet".
	 */
	dist_best = NO_NORM;

	for (qindex[3] = 0; qindex[3] < count[3]; qindex[3]++)
	    for (qindex[2] = 0; qindex[2] < count[2]; qindex[2]++)
		for (qindex[1] = 0; qindex[1] < count[1]; qindex[1]++)
		    for (qindex[0] = 0; qindex[0] < count[0]; qindex[0]++)
		    {
			/*
			 * Calculate the quaternion for this search point.
			 */
			for (kk = 0; kk < 4; kk++)
			{
			    /*
			     * The term in parenthesis ranges from -1 to +1,
			     * inclusive, so qq ranges from (-range+center)
			     * to (+range + center).
			     */
			    qq[kk] =
			     range[kk] *
			     (((FLT_DBL)
			       (2 * qindex[kk] -
				(count[kk] -
				 1))) / ((FLT_DBL) (count[kk] - 1))) +
			     center[kk];
			}

			/*
			 * Convert from a quaternion to a rotation matrix.
			 * The subroutine also takes care of normalizing the
			 * quaternion.
			 */
			quaternion_to_matrix (qq, rmat);
			/*
			 * Apply the rotation matrix to the elastic stiffness
			 * matrix.
			 */
			rotate_tensor (ccrot, cc, rmat);

			/*
			 * Find the distance of the rotated medium from
			 * orthorhombic aligned with the coordinate axes.
			 */
			dist = ortho_distance (ccortho, ccrot);

			/*
			 * If it's the best found so far, or the first time
			 * through, remember it.
			 */
			if (dist < dist_best || dist_best < 0.)
			{
			    dist_best = dist;
			    for (kk = 0; kk < 4; kk++)
				qq_best[kk] = qq[kk];
			}
		    }

	/*
	 * Refine for the next, finer, search. To avoid any possible problem
	 * caused by the optimal solution landing at an edge, we search over
	 * twice the distance between the two search points from the previous
	 * iteration.
	 */
	for (kk = 0; kk < 4; kk++)
	{
	    center[kk] = qq_best[kk];
	    count[kk] = SUBDIVIDE;
	    range[kk] = inc[kk];
	}

	/*
	 * We keep refining and searching the ever finer grid until we
	 * achieve the required accuracy, at which point we fall out the
	 * bottom of the loop here.
	 */
    }

    /*
     * We've got the answer to sufficient resolution... clean it up a bit,
     * then output it.
     */

    /*
     * Convert the best answer from a Quaternion back to a rotation matrix
     */
    quaternion_to_matrix (qq_best, rmat);

    /*
     * To make the order of the axes unique, we sort the principal axes
     * according to how well they work as a TI symmetry axis.
     * 
     * Specifically, since after rotation the medium is canonically oriented,
     * with the X, Y, and Z axes the principal axes, the INVERSE rotation
     * must take the X, Y, and Z axes to the original arbitrarily oriented
     * principal axes. So we first inverse-rotate a coordinate axis back to a
     * principal axis. We then use vector_to_angles to give us the Euler
     * angles theta and phi for the principal axis. make_rotation_matrix then
     * constructs a rotation matrix that rotates that principal axis to +Z.
     * We then use that matrix to rotate the tensor. We then measure its
     * distance from VTI, and remember that distance.
     */

    /*
     * First we need to find the inverse (the same as the transpose, because
     * it's _unitary_) of the rotation matrix rmat.
     */
    transpose_matrix (rmat_transp, rmat);

    /* Test the X axis */
    vec[0] = 1.;
    vec[1] = 0.;
    vec[2] = 0.;
    matrix_times_vector (vec2, rmat_transp, vec);
    vector_to_angles (vec2, &phi, &theta);
    make_rotation_matrix (theta, phi, 0., rmat_temp);
    rotate_tensor (ccrot2, cc, rmat_temp);
    dista[0] = ti_distance (ccti, ccrot2);

    /* Test the Y axis */
    vec[0] = 0.;
    vec[1] = 1.;
    vec[2] = 0.;
    matrix_times_vector (vec2, rmat_transp, vec);
    vector_to_angles (vec2, &phi, &theta);
    make_rotation_matrix (theta, phi, 0., rmat_temp);
    rotate_tensor (ccrot2, cc, rmat_temp);
    dista[1] = ti_distance (ccti, ccrot2);

    /* Test the Z axis */
    vec[0] = 0.;
    vec[1] = 0.;
    vec[2] = 1.;
    matrix_times_vector (vec2, rmat_transp, vec);
    vector_to_angles (vec2, &phi, &theta);
    make_rotation_matrix (theta, phi, 0., rmat_temp);
    rotate_tensor (ccrot2, cc, rmat_temp);
    dista[2] = ti_distance (ccti, ccrot2);


    /*
     * See which axis best functions as a TI symmetry axis, and make that one
     * the Z axis.
     */
    if (dista[2] <= dista[1] && dista[2] <= dista[0])
    {
	/* The Z axis is already the best. No rotation needed. */
	make_rotation_matrix (0., 0., 0., rmat_temp);
    }
    else if (dista[1] <= dista[2] && dista[1] <= dista[0])
    {
	/* Rotate Y to Z */
	make_rotation_matrix (0., 90., 0., rmat_temp);
	temp = dista[2];
	dista[2] = dista[1];
	dista[1] = temp;
    }
    else
    {
	/* Rotate X to Z */
	make_rotation_matrix (90., 90., -90., rmat_temp);
	temp = dista[2];
	dista[2] = dista[0];
	dista[0] = temp;
    }

    /*
     * Accumulate this axis-relabeling rotation (rmat_temp) onto the original
     * rotation (rmat).
     */
    matrix_times_matrix (rmat_temp2, rmat_temp, rmat);

    /*
     * Now find the next-best TI symmetry axis and make that one the Y axis.
     */
    if (dista[1] <= dista[0])
    {
	/* Already there; do nothing. */
	make_rotation_matrix (0., 0., 0., rmat_temp);
    }
    else
    {
	/* Rotate X to Y */
	make_rotation_matrix (90., 0., 0., rmat_temp);
	temp = dista[1];
	dista[1] = dista[0];
	dista[0] = temp;
    }

    /*
     * Accumulate the new axis relabeling rotation (rmat_temp) onto the
     * combined previous rotation matrix (rmat_temp2) to produce the final
     * desired result, rmat. The axes should now be in sorted order.
     */
    matrix_times_matrix (rmat, rmat_temp, rmat_temp2);

    return dist_best;
}
