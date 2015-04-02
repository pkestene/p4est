#ifdef P4_TO_P8
#include <p8est_connectivity.h>
#include <p8est_geometry.h>
#else
#include <p4est_connectivity.h>
#include <p4est_geometry.h>
#endif

#include "shell2d_connectivity.h"

p4est_connectivity_t *
p4est_connectivity_new_shell2d (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 6;
  const p4est_topidx_t num_trees    = 8;
  const p4est_topidx_t num_ctt      = 0;
  const double         vertices[6 * 3] = {
    -1,  0,  1,
     0,  0,  1,
     1,  0,  1,
    -1,  0,  2,
     0,  0,  2,
     1,  0,  2,
  };
  const p4est_topidx_t tree_to_vertex[8 * 4] = {
    0, 1, 3, 4,
    1, 2, 4, 5,
    0, 1, 3, 4,
    1, 2, 4, 5,
    0, 1, 3, 4,
    1, 2, 4, 5,
    0, 1, 3, 4,
    1, 2, 4, 5,
  };
  const p4est_topidx_t tree_to_tree[8 * 4] = {
    1, 7, 0, 0,
    2, 0, 1, 1,
    3, 1, 2, 2,
    4, 2, 3, 3,
    5, 3, 4, 4,
    6, 4, 5, 5,
    7, 5, 6, 6,
    0, 6, 7, 7,
  };
  const int8_t        tree_to_face[8 * 4] = {
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
  };
/* *INDENT-ON* */

  p4est_connectivity_t * conn =
    p4est_connectivity_new_copy (num_vertices, num_trees, 0,
				 vertices, tree_to_vertex,
				 tree_to_tree, tree_to_face,
				 NULL, &num_ctt, NULL, NULL);
  
  printf("\nconnectivity good : %d\n",p4est_connectivity_is_valid (conn));
  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
  

  return conn;
} /* p4est_connectivity_new_shell2d */


static void
p4est_geometry_connectivity_X (p4est_geometry_t * geom,
                               p4est_topidx_t which_tree,
                               const double abc[3], double xyz[3])
{
  p4est_connectivity_t *connectivity = (p4est_connectivity_t *) geom->user;
  const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
  const double       *v = connectivity->vertices;
  double              eta_x, eta_y, eta_z = 0.;
  int                 j, k;
  p4est_topidx_t      vt[P4EST_CHILDREN];

  /* retrieve corners of the tree */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    vt[k] = tree_to_vertex[which_tree * P4EST_CHILDREN + k];
  }

  /* these are reference coordinates in [0, 1]**d */
  eta_x = abc[0];
  eta_y = abc[1];
  eta_z = abc[2];

  /* bi/trilinear transformation */
  for (j = 0; j < 3; ++j) {
    /* *INDENT-OFF* */
    xyz[j] =
           ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                  eta_x  * v[3 * vt[1] + j]) +
                                  eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                  eta_x  * v[3 * vt[3] + j]))
#ifdef P4_TO_P8
            +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                  eta_x  * v[3 * vt[5] + j]) +
                                  eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                  eta_x  * v[3 * vt[7] + j]))
#endif
           );
    /* *INDENT-ON* */
  }
}



typedef enum
{
  P4EST_GEOMETRY_BUILTIN_MAGIC = 0x65F2F8DE, /* should be different from P8EST_GEOMETRY_BUILTIN_MAGIC ? */
  P4EST_GEOMETRY_BUILTIN_SHELL2D
}
p4est_geometry_builtin_type_t;

typedef struct p4est_geometry_builtin_shell2d
{
  p4est_geometry_builtin_type_t type;
  double              R2, R1;
  double              R2byR1, R1sqrbyR2, Rlog;
}
p4est_geometry_builtin_shell2d_t;

typedef struct p4est_geometry_builtin
{
  /** The geom member needs to come first; we cast to p4est_geometry_t * */
  p4est_geometry_t    geom;
  union
  {
    p4est_geometry_builtin_type_t type;
    p4est_geometry_builtin_shell2d_t shell2d;
  }
  p;
}
p4est_geometry_builtin_t;

/* geometric coordinate transformation */
static void
p4est_geometry_shell2d_X (p4est_geometry_t * geom,
			  p4est_topidx_t which_tree,
			  const double rst[3], 
			  double xyz[3])
{
  const struct p4est_geometry_builtin_shell2d *shell2d
    = &((p4est_geometry_builtin_t *) geom)->p.shell2d;
  double              x, y, R, q;
  double              abc[3];

  xyz[2] = 0.0;

  /* transform from the reference cube into vertex space */
  p4est_geometry_connectivity_X (geom, which_tree, rst, abc);

  //printf("rst %f %f %f abc %f %f %f\n",rst[0],rst[1],rst[2],abc[0],abc[1],abc[2]);


  /* assert that input points are in the expected range */
  P4EST_ASSERT (shell2d->type == P4EST_GEOMETRY_BUILTIN_SHELL2D);
  P4EST_ASSERT (0 <= which_tree && which_tree < 8);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[2] < 2.0 + SC_1000_EPS && abc[2] >  1.0 - SC_1000_EPS);

  /* abc[1] is always 0 here ... */

  /* transform abc[0] in-place for nicer grading */
  x = tan (abc[0] * M_PI_4);

  /* compute transformation ingredients */
  R = shell2d->R1sqrbyR2 * pow (shell2d->R2byR1, abc[2]);
  q = R / sqrt (x * x + 1.);

  /* assign correct coordinates based on patch id */
  switch (which_tree / 2) {
  case 3:                      /* top */
    xyz[0] = +q;
    xyz[1] = +q * x;
    break;
  case 2:                      /* left */
    xyz[0] = -q;
    xyz[1] = -q * x;
    break;
  case 1:                      /* bottom */
    xyz[0] = -q * x;
    xyz[1] =  q;
    break;
  case 0:                      /* right */
    xyz[0] = +q * x;
    xyz[1] = -q;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
} /* p4est_geometry_shell2d_X */

p4est_geometry_t   *
p4est_geometry_new_shell2d (p4est_connectivity_t * conn, double R2, double R1)
{
  p4est_geometry_builtin_t *builtin;
  struct p4est_geometry_builtin_shell2d *shell2d;

  builtin = P4EST_ALLOC_ZERO (p4est_geometry_builtin_t, 1);

  shell2d = &builtin->p.shell2d;
  shell2d->type = P4EST_GEOMETRY_BUILTIN_SHELL2D;
  shell2d->R2 = R2;
  shell2d->R1 = R1;
  shell2d->R2byR1 = R2 / R1;
  shell2d->R1sqrbyR2 = R1 * R1 / R2;
  shell2d->Rlog = log (R2 / R1);

  builtin->geom.name = "p4est_shell2d";
  builtin->geom.user = conn;
  builtin->geom.X = p4est_geometry_shell2d_X;

  return (p4est_geometry_t *) builtin;
}

