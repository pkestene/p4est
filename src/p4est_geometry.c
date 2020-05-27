/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
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

/**
 * \file p4est_geometry.c
 * We provide the identity transformation for reference.
 * Please implement p4est_geometry_t as you see fit.
 */

#ifndef P4_TO_P8
#include <p4est_geometry.h>
#else
#include <p8est_geometry.h>
#endif

typedef enum
{
  P4EST_GEOMETRY_BUILTIN_MAGIC = 0x65F2F8DE, /* should it be different from P8EST_GEOMETRY_BUILTIN_MAGIC ? */
  P4EST_GEOMETRY_BUILTIN_ICOSAHEDRON
}
p4est_geometry_builtin_type_t;

typedef struct p4est_geometry_builtin_icosahedron
{
  p4est_geometry_builtin_type_t type;
  double                           R; /* sphere radius */
}
p4est_geometry_builtin_icosahedron_t;

typedef struct p4est_geometry_builtin
{
  /** The geom member needs to come first; we cast to p4est_geometry_t * */
  p4est_geometry_t    geom;
  union
  {
    p4est_geometry_builtin_type_t type;
    p4est_geometry_builtin_icosahedron_t icosahedron;
  }
  p;
}
p4est_geometry_builtin_t;

/* geometric coordinate transformation */
static void
p4est_geometry_icosahedron_X (p4est_geometry_t * geom,
			      p4est_topidx_t which_tree,
			      const double rst[3], 
			      double xyz[3])
{
  const struct p4est_geometry_builtin_icosahedron *icosahedron
    = &((p4est_geometry_builtin_t *) geom)->p.icosahedron;
  double              x, y, z;
  double              a = 1.0;
  double              g = (1.0+sqrt(5.0))*0.5; /* golden ratio */
  double              ga = a/g;
  double              r = sqrt( (5.0-sqrt(5.0))*0.5 ); /* sqrt(a*a+ga*ga) -> current radius */

  double              R = icosahedron->R; /* target sphere radius */
  double              radius_ratio = R/r;

  /* these are reference coordinates in [0, 1]**d */
  double              eta_x, eta_y, eta_z = 0.;
  eta_x = rst[0];
  eta_y = rst[1];
  eta_z = rst[2];

  /*
   * icosahedron node cartesian coordinates
   * used for mapping connectivity vertices to 3D nodes.
   */
  const double N[12 * 3] = {
    0  , -ga,   a, /*  N0 */
    ga , - a,   0, /*  N1 */
    a  ,   0,  ga, /*  N2 */
    0  ,  ga,   a, /*  N3 */
    -a ,   0,  ga, /*  N4 */
    -ga, - a,   0, /*  N5 */
    a  ,   0, -ga, /*  N6 */
    ga ,   a,   0, /*  N7 */
    -ga,   a,   0, /*  N8 */
    -a ,   0, -ga, /*  N9 */
    0  , -ga,  -a, /* N10 */
    0  ,  ga,  -a, /* N11 */
  };

  /*
   * tree to nodes:
   *
   * tree 0:  1  6  0  2 
   * tree 1:  2  7  0  3
   * tree 2:  3  8  0  4
   * tree 3:  4  9  0  5
   * tree 4:  5 10  0  1
   * tree 5:  6 11  2  7
   * tree 6:  7 11  3  8
   * tree 7:  8 11  4  9
   * tree 8:  9 11  5 10
   * tree 9: 10 11  1  6
   *
   */
  const int tree_to_nodes[10*4] = {
    1 ,  6,  0,  2, 
    6 , 11,  2,  7,
    2 ,  7,  0,  3,
    7 , 11,  3,  8,
    3 ,  8,  0,  4,
    8 , 11,  4,  9,
    4 ,  9,  0,  5,
    9 , 11,  5, 10,
    5 , 10,  0,  1,
    10, 11,  1,  6,
  };

  /* assert that input points are in the expected range */
  P4EST_ASSERT (icosahedron->type == P4EST_GEOMETRY_BUILTIN_ICOSAHEDRON);
  P4EST_ASSERT (0 <= which_tree && which_tree < 10);

  /* transform from the reference cube into vertex space */

  /* assign correct coordinates based on patch id */
  /* use bilinear SLERP :  spherical bilinear interpolation */
  {
    int    j;
    
    /* use tree to nodes mapping to get nodes index of current tree */
    const int i0 = tree_to_nodes[which_tree*4+0];
    const int i1 = tree_to_nodes[which_tree*4+1];
    const int i2 = tree_to_nodes[which_tree*4+2];
    const int i3 = tree_to_nodes[which_tree*4+3];

    /* get 3D cartesian coordinates of our face */
    const double n0[3] = { N[i0*3 + 0], N[i0*3 + 1], N[i0*3 + 2] };
    const double n1[3] = { N[i1*3 + 0], N[i1*3 + 1], N[i1*3 + 2] };
    const double n2[3] = { N[i2*3 + 0], N[i2*3 + 1], N[i2*3 + 2] };
    const double n3[3] = { N[i3*3 + 0], N[i3*3 + 1], N[i3*3 + 2] };
    double norme2 = n0[0]*n0[0] + n0[1]*n0[1] + n0[2]*n0[2];

    /* 1. apply slerp 
     * - between n0 and n1 
     * - between n2 and n3
     */
    double xyz01[3]; // slerp along n0 -> n1
    double xyz23[3]; // slerp along n2 -> n3
    double dot1   = n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2];
    double theta1 = acos(dot1/norme2);

    for (j=0; j<3; ++j) {
      xyz01[j] = 
	sin((1.0-eta_x)*theta1)/sin(theta1)*n0[j]+
	sin((    eta_x)*theta1)/sin(theta1)*n1[j];
      xyz23[j] = 
	sin((1.0-eta_x)*theta1)/sin(theta1)*n2[j]+
	sin((    eta_x)*theta1)/sin(theta1)*n3[j];
    }

    /* apply slerp between xyz01 and xyz23 */
    double dot2 = xyz01[0]*xyz23[0] + xyz01[1]*xyz23[1] + xyz01[2]*xyz23[2];
    double theta2 = acos(dot2/norme2);
    for (j=0; j<3; ++j) {
      xyz[j] = 
	sin((1.0-eta_y)*theta2)/sin(theta2)*xyz01[j]+
	sin((    eta_y)*theta2)/sin(theta2)*xyz23[j];
      xyz[j] *= radius_ratio; /* rescale coordinates to target radius */
    }

    /* printf("DEBUG : %g | %g %g %g | \n", */
    /* 	   xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2], */
    /* 	   xyz[0],xyz[1],xyz[2] */
    /* 	   ); */

  } /* end of bilinear slerp */

} /* p4est_geometry_icosahedron_X */

void
p4est_geometry_destroy (p4est_geometry_t * geom)
{
  if (geom->destroy != NULL) {
    geom->destroy (geom);
  }
  else {
    P4EST_FREE (geom);
  }
}

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

p4est_geometry_t   *
p4est_geometry_new_connectivity (p4est_connectivity_t * conn)
{
  p4est_geometry_t   *geom;

  P4EST_ASSERT (conn->vertices != NULL);

  geom = P4EST_ALLOC_ZERO (p4est_geometry_t, 1);

  geom->name = P4EST_STRING "_connectivity";
  geom->user = conn;
  geom->X = p4est_geometry_connectivity_X;

  return geom;
}

#ifndef P4_TO_P8

p4est_geometry_t   *
p4est_geometry_new_icosahedron (p4est_connectivity_t * conn, double R)
{
  p4est_geometry_builtin_t *builtin;
  struct p4est_geometry_builtin_icosahedron *icosahedron;

  builtin = P4EST_ALLOC_ZERO (p4est_geometry_builtin_t, 1);

  icosahedron = &builtin->p.icosahedron;
  icosahedron->type = P4EST_GEOMETRY_BUILTIN_ICOSAHEDRON;
  icosahedron->R = R;

  builtin->geom.name = "p4est_icosahedron";
  builtin->geom.user = conn;
  builtin->geom.X = p4est_geometry_icosahedron_X;

  return (p4est_geometry_t *) builtin;

} /* p4est_geometry_new_icosahedron */

#endif /* !P4_TO_P8 */
