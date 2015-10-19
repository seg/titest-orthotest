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
 * Given a set of elastic constants, the routine "ti_distance" finds the
 * distance between that set of elastic constants and the nearest
 * TI medium that has a Z axis of symmetry.
 *
 * We want to find the nearest TI medium regardless of the direction of the
 * symmetry axis. This subroutine scans over a hemisphere of possible symmetry
 * axis directions, finds the trial symmetry axis direction that produced
 * the minimum distance, and then refines the search in that vicinity until
 * it converges to the optimal answer.
 *
 * On input:
 * 	FLT_DBL * cc
 *		A 6x6 array of Voigt-notation elastic stiffness constants.
 *
 * On output:
 *      FLT_DBL * theta_best, FLT_DBL * phi_best
 *
 *		The axi-symmetry direction of the nearest TI approximation.
 *		   phi and theta are defined as follows:
 * 		   phi=0 is the +Z axis.
 * 		   phi=90 theta=0 is the +X axis.
 * 		   phi=90 theta=90 is the +Y axis.
 *
 * Return value:
 *	The distance between the nearest transversely isotropic medium
 *	and the input medium, in absolute units (not normalized).
 */

#include <math.h>
#include "cmat.h"

/*
 * How much the search grid is reduced in scale after each successive level
 * of refinement is complete.
 */
#define SUBDIVIDE	4

/*
 * Initial search every 5 degrees
 */
#define DEG_INC		5.

/*
 * END_RES sets at what grid-interval scale we stop refining
 * and declare victory.
 */
#ifdef DOUBLE_PRECISION
/* FLT_DBL is "double" */
#define END_RES		(1.e-9)
#else
/* FLT_DBL is "float" */
#define END_RES		(1.e-6)
#endif


FLT_DBL
find_ti (FLT_DBL * cc, FLT_DBL * theta_best, FLT_DBL * phi_best)
{
int             ii, jj, kk;
FLT_DBL         ccrot[6 * 6];
FLT_DBL         ccti[6 * 6];
FLT_DBL         rmat[9];
FLT_DBL         vec[3];
FLT_DBL         dist;
FLT_DBL         theta, phi, dist_best;
FLT_DBL         phi_min, phi_max, phi_inc;
FLT_DBL         theta_min, theta_max, theta_inc;
FLT_DBL         v0[3], v1[3], v2[3], vv[3];


/*
 * Begin the first symmetry-axis scan, spanning a hemisphere.
 * (By symmetry, the other hemisphere is equivalent, so a search over
 * a hemisphere is sufficient.)
 */

/*
 * Search in latitude from pole to equator. Proceed in increments of
 * 5 (DEG_INC) degrees.
 */
    phi_min = 0.;
    phi_max = 90.01;
    phi_inc = DEG_INC;

/* At each latitude, search all 360 degrees of longitude. */
    theta_min = 0.;
    theta_max = 360.;
/*
 * Keep track of the best so far. The norm must be non-negative, so
 * a norm of -1 indicates that we haven't got any value yet.
 */
    dist_best = -1.;

/* Loop to scan over latitude */
    for (phi = phi_min; phi < phi_max; phi += phi_inc)
    {
/*
 * Calculate a longitude increment that has the same size as the
 * latitude increment. Add a small amount of fuzz to prevent the
 * calculation from going singular at the pole.
 */
	theta_inc = phi_inc / sin ((fabs (phi) + .01) * DEGTORAD);

/* Loop to scan over longitude */
	for (theta = theta_min; theta < theta_max; theta += theta_inc)
	{
	    /*
	     * rmat is the rotation matrix that rotates the current trial
	     * symmetry axis, as defined by theta and phi, to the +Z axis.
	     */
	    make_rotation_matrix (theta, phi, 0., rmat);
	    /* Rotate the elastic constants using rmat */
	    rotate_tensor (ccrot, cc, rmat);
	    /*
	     * Find the distance of the rotated constants from VTI:
	     * transversely isotropic with a vertical (+Z) symmetry axis.
	     */
	    dist = ti_distance (ccti, ccrot);

	    /*
	     * Is it better than the best we have found so far, or is it the
	     * first time through the loop?
	     */
	    if (dist < dist_best || dist_best < 0.)
	    {
		dist_best = dist;
		*phi_best = phi;
		*theta_best = theta;
	    }
	}
    }

/*
 * We now have an approximate global answer.
 * Now progressively refine the search grid around the best point
 * we found in the previous global search. Keep looping, subdividing
 * the search grid by a factor of SUBDIVIDE each time, until the
 * resolution is smaller than END_RES, which defines the minimal
 * acceptable resolution.
 */
    while (phi_inc > END_RES)
    {
	/*
	 * Calculate the "current best" symmetry axis vector.
	 */
	/*
	 * rmat is the rotation matrix that takes the +Z axis back to this
	 * symmetry axis candidate. This is the inverse of what we needed
	 * before; hence the minus signs on phi and theta.
	 */
	make_rotation_matrix (0., -(*phi_best), -(*theta_best), rmat);

	/* Here is the +Z axis */
	vec[0] = 0.;
	vec[1] = 0.;
	vec[2] = 1.;
	/* Apply the rotation. */
	matrix_times_vector (v0, rmat, vec);
	/* v0 is now the "current best" symmetry-axis vector */

	/*
	 * Now construct two vectors perpendicular to the current best
	 * symmetry axis. We will use these to construct a small 2D grid on
	 * the surface of the sphere, with the grid centered on the best
	 * symmetry-axis candidate found so far.
	 */
	/* Rotate the +X vector */
	vec[0] = 1.;
	vec[1] = 0.;
	vec[2] = 0.;
	matrix_times_vector (v1, rmat, vec);

	/* Rotate the +Y vector */
	vec[0] = 0.;
	vec[1] = 1.;
	vec[2] = 0.;
	matrix_times_vector (v2, rmat, vec);

	/*
	 * Do a search over this small 2D grid. Keep track of the best so
	 * far. A negative distance means we don't have an answer yet.
	 */
	dist_best = -1.;

	/*
	 * Loop over a (4*SUBDIVIDE + 1)^2 grid centered on the current
	 * optimal point. The grid spacing for this search is phi_inc /
	 * SUBDIVIDE, where phi_inc was the grid spacing of the previous
	 * search. We make the search grid twice as big in each direction as
	 * we would have to to search the entire grid cell area from the
	 * previous search, so that we avoid any problems that might be
	 * caused by the optimal value lying near the edge of our search
	 * grid.
	 */

	/* The loop over basis vector v1 */
	for (ii = -2 * SUBDIVIDE; ii <= 2 * SUBDIVIDE; ii++)
	    /* The loop over basis vector v2 */
	    for (jj = -2 * SUBDIVIDE; jj <= 2 * SUBDIVIDE; jj++)
	    {
		/*
		 * Calculate the search vector's X, Y, and Z components. v0
		 * is the center of the grid; v1 and v2 are the two
		 * orthogonal basis vectors used to perturb v0.
		 */
		for (kk = 0; kk < 3; kk++)
		{
		    vv[kk] = v0[kk] +
		     tan (phi_inc * DEGTORAD) *
		     ((FLT_DBL) ii / (FLT_DBL) SUBDIVIDE) * v1[kk] +
		     tan (phi_inc * DEGTORAD) *
		     ((FLT_DBL) jj / (FLT_DBL) SUBDIVIDE) * v2[kk];
		}

		/*
		 * Convert the direction vector vv to spherical coordinates.
		 * This also normalizes it back to being on the unit sphere.
		 */
		vector_to_angles (vv, &phi, &theta);

		/*
		 * We now have a current trial symmetry direction given by
		 * phi and theta. Find the corresponding rotation matrix, and
		 * apply it to the input elastic stiffness constants to
		 * rotate that trial symmetry axis to the +Z direction.
		 */
		make_rotation_matrix (theta, phi, 0., rmat);
		rotate_tensor (ccrot, cc, rmat);

		/* Find the distance from VTI */
		dist = ti_distance (ccti, ccrot);

		/* Keep track of the best candidate found so far */
		if (dist < dist_best || dist_best < 0.)
		{
		    dist_best = dist;
		    *phi_best = phi;
		    *theta_best = theta;
		}
	    }

	/*
	 * We now have a new best candidate. Refine the grid and keep going
	 * until it's fine enough.
	 */
	phi_inc /= (FLT_DBL) SUBDIVIDE;
    }

    /* theta_best and phi_best are returned set. */
    return dist_best;
}
