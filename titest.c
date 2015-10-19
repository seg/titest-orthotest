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
 * titest < elastic_constants
 *
 * titest reads from standard input a fully general anisotropic
 * stiffness matrix in the form of 6 numbers on each of 6 lines of input.
 * It finds the best-fitting transversely isotropic (TI) medium,
 * and outputs:
 * 0) the input matrix,
 * 1) the input elastic constants rotated so that the
 * best-fitting TI axis is the Z axis,
 * 2) the TI approximation to the rotated matrix,
 * 3) the TI approximation in the original unrotated coordinate system,
 * 4) the percent difference between the input stiffness matrix and the
 * best-fitting TI approximation, in the original coordinates (normalized
 * by dividing each difference by the norm of the input stiffness matrix),
 * 5) the percent error from TI (0 error means the medium is exactly TI; 100%
 * is the maximum possible error, which is only possible in extreme cases,
 * for example if c16=1 and all the other elastic constants are 0.), and
 * 6) the coordinates of the axis vector, in both cartesian and theta-phi
 * notation.
 *
 * phi and theta are defined as follows:
 * phi=0 is the +Z axis
 * phi=90 theta=0 is the +X axis
 * phi=90 theta=90 is the +Y axis
 *
 * For more about what "best-fitting" means for elastic stiffness matrices, see
 * the article by Arts, Helbig, and Rasolofosaon in the SEG extended abstracts
 * for 1991, page 1534:
 * "General Anisotropic Elastic Tensor in Rocks: Approximation,
 * Invariants, and Particular Directions".
 *
 * The following stiffness matrix is TI (transversely isotropic),
 * but this fact is not obvious because it has been
 * rotated to have a symmetry axis pointing in the direction
 * phi=12.345 and theta=67.890 degrees.
 * Try it as input to titest:
 *
 *  331.325      128.029      112.309     -1.30380     -23.3328     -1.92204
 *  128.029      339.374      108.716     -9.83459     -4.08399     -1.99410
 *  112.309      108.716      226.191     0.447454      1.10140      1.74841
 *  -1.30380     -9.83459     0.447454      56.8929      1.27023     -9.88887
 *  -23.3328     -4.08399      1.10140      1.27023      59.5035     -3.66209
 *  -1.92204     -1.99410      1.74841     -9.88887     -3.66209      103.658
 *
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
FLT_DBL         ccti[6 * 6];
FLT_DBL         rmat[9];
FLT_DBL         rmat_transp[9];
FLT_DBL         vec[3];
FLT_DBL         vec_sym[3];
FLT_DBL         dist;
FLT_DBL         norm;
FLT_DBL         theta_best, phi_best, dist_best;


/*
 * Read in the input elastic constants.
 */
    read_matrix_6x6 (cc);

    printf ("Input C matrix:\n");
    print_matrix_6x6 (cc);
    printf ("\n\n");

/*
 * Find the norm of this set of elastic constants.
 * We will use this for calculating the percent error later.
 */
    norm = norm_matrix_6x6 (cc);

/*
 * Find the best-approximating TI medium.
 */
    dist_best = find_ti (cc, &theta_best, &phi_best);

    /*
     * Output the results
     */

    /*
     * Rotate the input elastic constants so that the symmetry axis is +Z, so
     * that we can recognize the symmetry, and output it.
     */
    make_rotation_matrix (theta_best, phi_best, 0., rmat);
    rotate_tensor (ccrot, cc, rmat);
    printf ("Rotated C matrix:\n");
    print_matrix_6x6 (ccrot);
    printf ("\n");

    /*
     * Output the nearest VTI elastic constants.
     */
    dist = ti_distance (ccti, ccrot);
    printf ("TI approximation:\n");
    print_matrix_6x6 (ccti);
    printf ("\n");

    /* Transpose the rotation matrix so that it rotates the other way. */
    transpose_matrix (rmat_transp, rmat);
    /*
     * Rotate the best VTI approximation back to the original symmetry axis
     * direction.
     */
    rotate_tensor (cc2, ccti, rmat_transp);
    printf ("TI approximation in original coordinate system:\n");
    print_matrix_6x6 (cc2);
    printf ("\n");

    /*
     * Calculate the deviation from TI element by element.
     */
    for (ii = 0; ii < 6; ii++)
	for (jj = 0; jj <= ii; jj++)
	{
	    CC1 (ii, jj) = CC1 (jj, ii) =
	     (CC (ii, jj) - CC2 (ii, jj)) * 100. / norm;
	}
    printf
     ("Normalized deviation from TI in original coordinate system, in percent:\n");
    format_print_matrix_6x6 ("%11.4f ", cc1);
    printf ("\n");

    /* Output the distance from TI, normalized to a percentage. */
    printf ("distance from TI = %.3f percent\n", 100. * dist / norm);

    /* Output the symmetry axis direction. */
    vec[0] = 0.;
    vec[1] = 0.;
    vec[2] = 1.;
    matrix_times_vector (vec_sym, rmat_transp, vec);
    printf ("Symmetry axis: (%.4f, %.4f, %.4f)\n",
	    vec_sym[0], vec_sym[1], vec_sym[2]);

    printf ("theta = %.3f,   phi = %.3f\n", theta_best, phi_best);

    return 0;
}
