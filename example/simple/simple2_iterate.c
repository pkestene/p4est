/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*
 * Usage: p4est_simple <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit square.
 *        o three     Refinement on a forest with three trees.
 *        o evil      Check second round of refinement with np=5 level=7
 *        o evil3     Check second round of refinement on three trees
 *        o pillow    Refinement on a 2-tree pillow-shaped domain.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *        o cubed     Refinement on a 6-tree cubed sphere surface.
 *        o disk      Refinement on a 5-tree spherical disk.
 *        o periodic  Refinement on the unit square with all-periodic b.c.
 *        o rotwrap   Refinement on the unit square with weird periodic b.c.
 *        o brick     Refinement on a brick
 *        o ring      Refinement on a ring
 *        o shell2d   Refinement on a 2d shell
 */

#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#include <p4est_iterate.h>

#include "ring_connectivity.h"
#include "shell2d_connectivity.h"

typedef enum
{
  P4EST_CONFIG_NULL,
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_EVIL,
  P4EST_CONFIG_EVIL3,
  P4EST_CONFIG_PILLOW,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR,
  P4EST_CONFIG_CUBED,
  P4EST_CONFIG_DISK,
  P4EST_CONFIG_PERIODIC,
  P4EST_CONFIG_ROTWRAP,
  P4EST_CONFIG_BRICK,
  P4EST_CONFIG_RING,
  P4EST_CONFIG_SHELL2D
}
simple_config_t;

typedef struct
{
  simple_config_t     config;
  int                 mpisize;
  int                 level;
  unsigned            checksum;
}
simple_regression_t;

typedef struct
{
  p4est_topidx_t      a;
  double              is_outside; /* 0 means no */
}
user_data_t;

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

static int          refine_level = 0;

/* *INDENT-OFF* */
static const simple_regression_t regression[] =
{{ P4EST_CONFIG_THREE, 1, 7, 0xa8d85863U },
 { P4EST_CONFIG_THREE, 2, 7, 0xa8d85863U },
 { P4EST_CONFIG_THREE, 3, 7, 0xa8d85863U },
 { P4EST_CONFIG_THREE, 4, 7, 0x20fb58edU },
 { P4EST_CONFIG_MOEBIUS, 1, 6, 0x98ab6cb2U },
 { P4EST_CONFIG_MOEBIUS, 3, 6, 0x98ab6cb2U },
 { P4EST_CONFIG_MOEBIUS, 5, 6, 0x98ab6cb2U },
 { P4EST_CONFIG_MOEBIUS, 6, 6, 0x6d2d6d6cU },
 { P4EST_CONFIG_STAR, 5, 6, 0x38d3736fU },
 { P4EST_CONFIG_STAR, 5, 7, 0xfb97aadfU },
 { P4EST_CONFIG_CUBED, 4, 3, 0x85581649U },
 { P4EST_CONFIG_CUBED, 5, 5, 0x64a1d105U },
 { P4EST_CONFIG_DISK, 5, 4, 0x4995411dU },
 { P4EST_CONFIG_DISK, 2, 6, 0x3f758706U },
 { P4EST_CONFIG_ROTWRAP, 1, 6, 0x9dd600c5U },
 { P4EST_CONFIG_ROTWRAP, 3, 6, 0x9dd600c5U },
 { P4EST_CONFIG_NULL, 0, 0, 0 }};
/* *INDENT-ON* */

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  user_data_t        *data = (user_data_t *) quadrant->p.user_data;

  data->a = which_tree;
  data->is_outside = 0.0; /* false */
}

static int
refine_normal_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                  p4est_quadrant_t * quadrant)
{
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id (quadrant) == 3) {
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static int
refine_evil_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quadrant)
{
  if ((int) quadrant->level >= refine_level) {
    return 0;
  }
  if (p4est->mpirank <= 1) {
    return 1;
  }

  return 0;
}

static int
refine_evil3_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * quadrant)
{
  p4est_qcoord_t      u2;
  p4est_quadrant_t    ref;

  P4EST_QUADRANT_INIT (&ref);

  u2 = P4EST_QUADRANT_LEN (2);

  if (which_tree == 0) {
    ref.x = 3 * u2;
    ref.y = 2 * u2;
  }
  else if (which_tree == 1) {
    ref.x = 2 * u2;
    ref.y = 3 * u2;
  }
  ref.level = 2;

  if ((int) quadrant->level >= refine_level) {
    return 0;
  }
  if ((which_tree == 0 || which_tree == 1) &&
      (p4est_quadrant_is_equal (&ref, quadrant) ||
       p4est_quadrant_is_ancestor (&ref, quadrant))) {
    return 1;
  }

  return 0;
}

static int
coarsen_evil_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * q[])
{
  if (p4est->mpirank >= 2) {
    return 1;
  }

  return 0;
}

static void
is_outside_boundary_compute_callback(p4est_iter_face_info_t * info, void *user_data) {

  int                     i, j;
  p4est_t                *p4est = info->p4est;

  p4est_quadrant_t       *quad;
  p4est_iter_face_side_t *side[2];
  sc_array_t             *sides = &(info->sides);

  /* on the outside boundary, a face has only one side */
  if (sides->elem_count == 1) {

    /* get face side data */
    side[0] = p4est_iter_fside_array_index_int (sides, 0);

    if (!side[0]->is_hanging) { /* outside boundary cell shouldn't be hanging */

      quad = side[0]->is.full.quad;

      user_data_t        *data = (user_data_t *) quad->p.user_data;
      data->is_outside = 1.0; /*true */

    }

  }

} /* end is_outside_boundary_compute_callback */

static void
is_outside_boundary_getdata_callback(p4est_iter_volume_info_t * info, void *user_data)
{
  /* we passed the array of values to fill as the user_data in the call
     to p4est_iterate */
  double             *is_outside_data = (double *) user_data;
  p4est_t            *p4est = info->p4est;
  p4est_quadrant_t   *q = info->quad;
  p4est_topidx_t      which_tree = info->treeid;
  p4est_locidx_t      local_id = info->quadid;  /* this is the index of q *within its tree's numbering*.  We want to convert it its index for all the quadrants on this process, which we do below */
  p4est_tree_t       *tree;
  user_data_t        *data = (user_data_t *) q->p.user_data;
  p4est_locidx_t      arrayoffset;

  tree = p4est_tree_array_index (p4est->trees, which_tree);
  local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */
  arrayoffset = local_id;

  is_outside_data[arrayoffset] = data->is_outside;

} /* is_outside_boundary_getdata_callback */

/*
 * perform a face iteration, and flags cells which have an outside boundary.
 */
static void
is_outside_boundary_compute(p4est_t *p4est) {

  /* iterate to find outside quadrant */
  /* *INDENT-OFF* */
  p4est_iterate (p4est,                 /* the forest */
		 NULL,                  /* the ghost layer */
		 NULL,                  /* the synchronized ghost data */
		 NULL,                  /* volume callback */
		 is_outside_boundary_compute_callback,     /* face callback */
#ifdef P4_TO_P8
		 NULL,                  /* there is no callback for the
					   edges between quadrants */
#endif
		 NULL);                 /* there is no callback for the
					   corners between quadrants */
    /* *INDENT-ON* */


} /* is_outside_boudary_compute */

void
write_user_data(p4est_t *p4est, p4est_geometry_t * geom, const char *filename)
{

  double *is_outside_data;    /* array of cell data to write */
  p4est_locidx_t      numquads;

  numquads = p4est->local_num_quadrants;

  /* create a vector with one value per local quadrant */
  is_outside_data = P4EST_ALLOC (double, numquads);

  /* Use the iterator to visit every cell and fill in the data array */
  p4est_iterate (p4est, 
		 NULL,   /* we don't need any ghost quadrants for this loop */
                 (void *) is_outside_data,     /* pass in array pointer so that we can fill it */
                 is_outside_boundary_getdata_callback,   /* callback function to get cell data */
                 NULL,
#ifdef P4_TO_P8
                 NULL,
#endif
                 NULL);

  p4est_vtk_write_header (p4est, geom, 0.95, filename);  
  p4est_vtk_write_cell_data (p4est, geom, 1, 1,
			     1, 0, 1, 0, filename, "is_outside", is_outside_data);
  p4est_vtk_write_footer (p4est, filename);


  P4EST_FREE(is_outside_data);

} /* write_user_data */

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 wrongusage;
  unsigned            crc;
  const char         *usage;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_geometry_t   *geom;
  p4est_refine_t      refine_fn;
  p4est_coarsen_t     coarsen_fn;
  simple_config_t     config;
  const simple_regression_t *r;

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  usage =
    "Arguments: <configuration> <level>\n"
    "   Configuration can be any of\n"
    "      unit|three|evil|evil3|pillow|moebius|\n"
    "         star|cubed|disk|periodic|rotwrap|brick|ring|shell2d\n"
    "   Level controls the maximum depth of refinement\n";
  wrongusage = 0;
  config = P4EST_CONFIG_NULL;
  if (!wrongusage && argc < 3) {
    wrongusage = 1;
  }
  if (!wrongusage) {
    if (!strcmp (argv[1], "unit")) {
      config = P4EST_CONFIG_UNIT;
    }
    else if (!strcmp (argv[1], "three")) {
      config = P4EST_CONFIG_THREE;
    }
    else if (!strcmp (argv[1], "evil")) {
      config = P4EST_CONFIG_EVIL;
    }
    else if (!strcmp (argv[1], "evil3")) {
      config = P4EST_CONFIG_EVIL3;
    }
    else if (!strcmp (argv[1], "pillow")) {
      config = P4EST_CONFIG_PILLOW;
    }
    else if (!strcmp (argv[1], "moebius")) {
      config = P4EST_CONFIG_MOEBIUS;
    }
    else if (!strcmp (argv[1], "star")) {
      config = P4EST_CONFIG_STAR;
    }
    else if (!strcmp (argv[1], "cubed")) {
      config = P4EST_CONFIG_CUBED;
    }
    else if (!strcmp (argv[1], "disk")) {
      config = P4EST_CONFIG_DISK;
    }
    else if (!strcmp (argv[1], "periodic")) {
      config = P4EST_CONFIG_PERIODIC;
    }
    else if (!strcmp (argv[1], "rotwrap")) {
      config = P4EST_CONFIG_ROTWRAP;
    }
    else if (!strcmp (argv[1], "brick")) {
      config = P4EST_CONFIG_BRICK;
    }
    else if (!strcmp (argv[1], "ring")) {
      config = P4EST_CONFIG_RING;
    }
    else if (!strcmp (argv[1], "shell2d")) {
      config = P4EST_CONFIG_SHELL2D;
    }
    else {
      wrongusage = 1;
    }
  }
  if (wrongusage) {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("Usage error");
  }

  /* assign variables based on configuration */
  refine_level = atoi (argv[2]);
  if (config == P4EST_CONFIG_EVIL) {
    refine_fn = refine_evil_fn;
    coarsen_fn = coarsen_evil_fn;
  }
  else if (config == P4EST_CONFIG_EVIL3) {
    refine_fn = refine_evil3_fn;
    coarsen_fn = NULL;
  }
  else {
    refine_fn = refine_normal_fn;
    coarsen_fn = NULL;
  }

  /* create connectivity and forest structures */
  geom = NULL;
  if (config == P4EST_CONFIG_THREE || config == P4EST_CONFIG_EVIL3) {
    connectivity = p4est_connectivity_new_corner ();
  }
  else if (config == P4EST_CONFIG_PILLOW) {
    connectivity = p4est_connectivity_new_pillow ();
  }
  else if (config == P4EST_CONFIG_MOEBIUS) {
    connectivity = p4est_connectivity_new_moebius ();
  }
  else if (config == P4EST_CONFIG_STAR) {
    connectivity = p4est_connectivity_new_star ();
  }
  else if (config == P4EST_CONFIG_CUBED) {
    connectivity = p4est_connectivity_new_cubed ();
  }
  else if (config == P4EST_CONFIG_DISK) {
    connectivity = p4est_connectivity_new_disk ();
  }
  else if (config == P4EST_CONFIG_PERIODIC) {
    connectivity = p4est_connectivity_new_periodic ();
  }
  else if (config == P4EST_CONFIG_ROTWRAP) {
    connectivity = p4est_connectivity_new_rotwrap ();
  }
  else if (config == P4EST_CONFIG_BRICK) {
    int numtrees_x = 1;
    int numtrees_y = 1;

    int periodic_x = 0;
    int periodic_y = 0;

    if (argc >= 4)
      numtrees_x = atoi (argv[3]);
    if (argc >= 5)
      numtrees_y = atoi (argv[4]);
    if (argc >= 6)
      periodic_x = atoi (argv[5]);
    if (argc >= 7)
      periodic_y = atoi (argv[6]);

    connectivity = p4est_connectivity_new_brick (numtrees_x, numtrees_y,
						 periodic_x, periodic_y);
  }
  else if (config == P4EST_CONFIG_RING) {

    int num_trees_radial = 3;
    int num_trees_orthoradial = 16;
    double rMin = 1.0;
    double rMax = 2.0;

    if (argc >= 4)
      num_trees_radial = atoi (argv[3]);
    if (argc >= 5)
      num_trees_orthoradial = atoi (argv[4]);
    if (argc >= 6)
      rMin = atof (argv[5]);
    if (argc >= 7)
      rMax = atof (argv[6]);

    connectivity = p4est_connectivity_new_ring (num_trees_radial, 
						num_trees_orthoradial,
						rMin,
						rMax);
  }
  else if (config == P4EST_CONFIG_SHELL2D) {
    double rMin = 0.55;
    double rMax = 1.0;

    if (argc >= 4)
      rMin = atof (argv[3]);
    if (argc >= 5)
      rMax = atof (argv[4]);

    connectivity = p4est_connectivity_new_shell2d ();
    geom = p4est_geometry_new_shell2d (connectivity, rMax, rMin);
  }
  else {
    connectivity = p4est_connectivity_new_unitsquare ();
  }
  p4est = p4est_new_ext (mpi->mpicomm, connectivity, 15, 0, 0,
                         sizeof (user_data_t), init_fn, NULL);
  p4est_vtk_write_file (p4est, geom, "simple2_new");

  /* refinement and coarsening */
  p4est_refine (p4est, 1, refine_fn, init_fn);
  if (coarsen_fn != NULL) {
    p4est_coarsen (p4est, 1, coarsen_fn, init_fn);
  }
  p4est_vtk_write_file (p4est, geom, "simple2_refined");

  /* balance */
  p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
  p4est_vtk_write_file (p4est, geom, "simple2_balanced");
  crc = p4est_checksum (p4est);

  /* partition */
  p4est_partition (p4est, 0, NULL);
  p4est_vtk_write_file (p4est, geom, "simple2_partition");

#ifdef P4EST_ENABLE_DEBUG
  /* rebalance should not change checksum */
  p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
  P4EST_ASSERT (p4est_checksum (p4est) == crc);
#endif

  /* flag all quadrants that touch the outside boundary */
  is_outside_boundary_compute(p4est);

  /* save outside boundary data */
  write_user_data(p4est, geom, "simple2_outside");
  
  /* print and verify forest checksum */
  P4EST_GLOBAL_STATISTICSF ("Tree checksum 0x%08x\n", crc);
  if (mpi->mpirank == 0) {
    for (r = regression; r->config != P4EST_CONFIG_NULL; ++r) {
      if (r->config != config || r->mpisize != mpi->mpisize
          || r->level != refine_level)
        continue;
      SC_CHECK_ABORT (crc == r->checksum, "Checksum mismatch");
      P4EST_GLOBAL_INFO ("Checksum regression OK\n");
      break;
    }
  }

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  if (geom != NULL) {
    p4est_geometry_destroy (geom);
  }
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
