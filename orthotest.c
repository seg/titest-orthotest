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
 * Usage:
 *
 * orthotest < elastic_constants
 *
 * elastic constants is a file with 36 numbers in it,
 * usually 6 numbers on each of 6 lines.
 * For example:
 *
 *   331.325      128.029      112.309     -1.30380     -23.3328     -1.92204
 *   128.029      339.374      108.716     -9.83459     -4.08399     -1.99410
 *   112.309      108.716      226.191     0.447454      1.10140      1.74841
 *  -1.30380     -9.83459     0.447454      56.8929      1.27023     -9.88887
 *  -23.3328     -4.08399      1.10140      1.27023      59.5035     -3.66209
 *  -1.92204     -1.99410      1.74841     -9.88887     -3.66209      103.658
 *
 * orthotest finds the best-fitting Orthorhombic medium to this,
 * and outputs:
 * 0) the input elastic constants,
 * 1) the elastic matrix rotated so that the best-fitting Orthorhombic
 * matrix has the X, Y, and Z axes as principal axes. The axes are ordered
 * so that the Z axis is the closest to being a TI axis of symmetry,
 * then Y, then X. Thus, if the input medium is actually TI then the Z axis
 * will be the axis of symmetry. (If the input medium is arbitrary, there is
 * no reason to expect that the Z axis found by this program will coincide
 * with the best-fitting TI axis found by titest!)
 * 2) the orthorhombic approximation to the rotated matrix,
 * 3) the orthorhombic approximation in the original coordinates,
 * 4) the percent difference between the input stiffness matrix and the
 * best-fitting orthorhombic approximation in the original coordinates
 * (normalized by dividing by the norm of the input stiffness matrix),
 * 5) the percent error from Orthorhombic, and
 * 6) the coordinates of the 3 principal axes, in both cartesian and
 * theta-phi notation.
 *
 * Phi and Theta are defined as follows:
 *  phi=0 is the +Z axis
 *  phi=90 theta=0 is the +X axis
 *  phi=90 theta=90 is the +Y axis
 *
 * For more about what "best-fitting" means, see
 * SEG extended abstracts 1991 page 1534, the article by
 * Arts, Helbig, and Rasolofosaon,
 * "General Anisotropic Elastic Tensor in Rocks: Approximation,
 * Invariants, and Particular Directions".
 * Note, however, that equation 3 in their paper is incorrect.
 * Their explanation of how to define the distance between
 * elastic stiffness matrices is correct, however.
 *
 * Author Joe Dellinger, Amoco TTC, 19 Feb 1997.
 */

#include <stdio.h>
#include <math.h>
#include "cmat.h"

int
main ()
{
int             ii, jj;
FLT_DBL         cc[6 * 6];
FLT_DBL         cc1[6 * 6];
FLT_DBL         cc2[6 * 6];
FLT_DBL         ccrot[6 * 6];
FLT_DBL         ccortho[6 * 6];
FLT_DBL         ccrot2[6 * 6];
FLT_DBL         ccti[6 * 6];
FLT_DBL         rmat[9];
FLT_DBL         rmat_transp[9];
FLT_DBL         rmat_temp[9];
FLT_DBL         vec[3];
FLT_DBL         vec2[3];
FLT_DBL         dist;
FLT_DBL         norm;
FLT_DBL         dist_best;
FLT_DBL         phi, theta;

/*
 * Read in the elastic constants
 */
    read_matrix_6x6 (cc);

    printf ("Input C matrix:\n");
    print_matrix_6x6 (cc);
    printf ("\n\n");

/*
 * Calculate the norm of this matrix, so that we can calculate
 * normalized deviations from symmetry for output.
 */
    norm = norm_matrix_6x6 (cc);

/*
 * Find the best-approximating orthorhombic medium.
 */
    dist_best = find_ortho (cc, rmat);

    transpose_matrix (rmat_transp, rmat);

    /*
     * Write out the result!
     */
    rotate_tensor (ccrot, cc, rmat);
    printf ("Rotated C matrix:\n");
    print_matrix_6x6 (ccrot);
    printf ("\n");

    dist = ortho_distance (ccortho, ccrot);
    printf ("Orthorhombic approximation:\n");
    print_matrix_6x6 (ccortho);
    printf ("\n");

    rotate_tensor (cc2, ccortho, rmat_transp);
    printf ("Orthorhombic approximation in original coordinates:\n");
    print_matrix_6x6 (cc2);
    printf ("\n");

    for (ii = 0; ii < 6; ii++)
	for (jj = 0; jj <= ii; jj++)
	{
	    CC1 (ii, jj) = CC1 (jj, ii) =
	     (CC (ii, jj) - CC2 (ii, jj)) * 100. / norm;
	}
    printf
     ("Normalized deviation from Orthorhombic in original coordinates, in percent:\n");
    format_print_matrix_6x6 ("%11.4f ", cc1);
    printf ("\n");

    printf ("Distance from Orthorhombic = %.3f percent\n",
	    100. * dist / norm);
    printf ("\n");

    /*
     * And write out the canonically ordered principal axes, and how well
     * each functions as a TI symmetry axis.
     */
    vec[0] = 1.;
    vec[1] = 0.;
    vec[2] = 0.;
    matrix_times_vector (vec2, rmat_transp, vec);
    printf ("X axis: (%.4f, %.4f, %.4f)  ", vec2[0], vec2[1], vec2[2]);
    vector_to_angles (vec2, &phi, &theta);
    make_rotation_matrix (theta, phi, 0., rmat_temp);
    rotate_tensor (ccrot2, cc, rmat_temp);
    dist = ti_distance (ccti, ccrot2);
    printf ("theta=%.3f, phi=%.3f, TI dist=%.3f%%\n",
	    theta, phi, 100. * dist / norm);


    vec[0] = 0.;
    vec[1] = 1.;
    vec[2] = 0.;
    matrix_times_vector (vec2, rmat_transp, vec);
    printf ("Y axis: (%.4f, %.4f, %.4f)  ", vec2[0], vec2[1], vec2[2]);
    vector_to_angles (vec2, &phi, &theta);
    make_rotation_matrix (theta, phi, 0., rmat_temp);
    rotate_tensor (ccrot2, cc, rmat_temp);
    dist = ti_distance (ccti, ccrot2);
    printf ("theta=%.3f, phi=%.3f, TI dist=%.3f%%\n",
	    theta, phi, 100. * dist / norm);

    vec[0] = 0.;
    vec[1] = 0.;
    vec[2] = 1.;
    matrix_times_vector (vec2, rmat_transp, vec);
    printf ("Z axis: (%.4f, %.4f, %.4f)  ", vec2[0], vec2[1], vec2[2]);
    vector_to_angles (vec2, &phi, &theta);
    make_rotation_matrix (theta, phi, 0., rmat_temp);
    rotate_tensor (ccrot2, cc, rmat_temp);
    dist = ti_distance (ccti, ccrot2);
    printf ("theta=%.3f, phi=%.3f, TI dist=%.3f%%\n",
	    theta, phi, 100. * dist / norm);

    return 0;
}
