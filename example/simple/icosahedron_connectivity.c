#ifdef P4_TO_P8
#include <p8est_connectivity.h>
#include <p8est_geometry.h>
#else
#include <p4est_connectivity.h>
#include <p4est_geometry.h>
#endif

#include "icosahedron_connectivity.h"

/*
 * Nodes numbering
 *
 *    A00   A01   A02   A03   A04 
 *   /   \ /   \ /   \ /   \ /   \
 * A05---A06---A07---A08---A09---A10
 *   \   / \   / \   / \   / \   / \
 *    A11---A12---A13---A14---A15---A16
 *      \  /  \  /  \  /  \  /  \  /
 *      A17   A18   A19   A20   A21
 *
 * Origin in A05.
 *
 * Tree numbering:
 *
 * 0  2  4  6  8
 *  1  3  5  7  9
 */
p4est_connectivity_t *
p4est_connectivity_new_icosahedron (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 22;
  const p4est_topidx_t num_trees    = 10;
  const p4est_topidx_t num_corners  =  2;
  //const p4est_topidx_t num_ctt      = 10;
  const double         vertices[22 * 3] = {
    0.0 +   cos(M_PI/3),    sin(M_PI/3),  0.0, /* vertex 00 */
    1.0 +   cos(M_PI/3),    sin(M_PI/3),  0.0, /* vertex 01 */
    2.0 +   cos(M_PI/3),    sin(M_PI/3),  0.0, /* vertex 02 */
    3.0 +   cos(M_PI/3),    sin(M_PI/3),  0.0, /* vertex 03 */
    4.0 +   cos(M_PI/3),    sin(M_PI/3),  0.0, /* vertex 04 */
    0.0,  0.0,  0.0,                           /* vertex 05 */
    1.0,  0.0,  0.0,                           /* vertex 06 */
    2.0,  0.0,  0.0,                           /* vertex 07 */
    3.0,  0.0,  0.0,                           /* vertex 08 */
    4.0,  0.0,  0.0,                           /* vertex 09 */
    5.0,  0.0,  0.0,                           /* vertex 10 */
    0.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 11 */
    1.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 12 */
    2.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 13 */
    3.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 14 */
    4.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 15 */ 
    5.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 16 */ 
    0.0 + 2*cos(M_PI/3), -2*sin(M_PI/3),  0.0, /* vertex 17 */
    1.0 + 2*cos(M_PI/3), -2*sin(M_PI/3),  0.0, /* vertex 18 */
    2.0 + 2*cos(M_PI/3), -2*sin(M_PI/3),  0.0, /* vertex 19 */
    3.0 + 2*cos(M_PI/3), -2*sin(M_PI/3),  0.0, /* vertex 20 */
    4.0 + 2*cos(M_PI/3), -2*sin(M_PI/3),  0.0, /* vertex 21 */
  };
  const p4est_topidx_t tree_to_vertex[10 * 4] = {
    5,  11,  0,  6, /* tree 0 */
    11, 17,  6, 12, /* tree 1 */
    6,  12,  1,  7, /* tree 2 */
    12, 18,  7, 13, /* tree 3 */
    7,  13,  2,  8, /* tree 4 */
    13, 19,  8, 14, /* tree 5 */
    8,  14,  3,  9, /* tree 6 */
    14, 20,  9, 15, /* tree 7 */
    9,  15,  4, 10, /* tree 8 */
    15, 21, 10, 16, /* tree 9 */
  };
  const p4est_topidx_t tree_to_tree[10 * 4] = {
    8,1,9,2, /* tree 0 */
    0,3,9,2, /* tree 1 */
    0,3,1,4, /* tree 2 */
    2,5,1,4, /* tree 3 */
    2,5,3,6, /* tree 4 */
    4,7,3,6, /* tree 5 */
    4,7,5,8, /* tree 6 */ 
    6,9,5,8, /* tree 7 */
    6,9,7,0, /* tree 8 */
    8,1,7,0, /* tree 9 */
  };
  const int8_t        tree_to_face[10 * 4] = {
    3,0,3,0, /* tree 0 */
    1,2,1,2, /* tree 1 */
    3,0,3,0, /* tree 2 */
    1,2,1,2, /* tree 3 */
    3,0,3,0, /* tree 4 */
    1,2,1,2, /* tree 5 */
    3,0,3,0, /* tree 6 */
    1,2,1,2, /* tree 7 */
    3,0,3,0, /* tree 8 */
    1,2,1,2, /* tree 9 */
  };
  const p4est_topidx_t tree_to_corner[10 * 4] = {
    -1,  -1,  0,  -1, /* tree 0 */
    -1,   1, -1,  -1, /* tree 1 */
    -1,  -1,  0,  -1, /* tree 2 */
    -1,   1, -1,  -1, /* tree 3 */
    -1,  -1,  0,  -1, /* tree 4 */
    -1,   1, -1,  -1, /* tree 5 */
    -1,  -1,  0,  -1, /* tree 6 */
    -1,   1, -1,  -1, /* tree 7 */
    -1,  -1,  0,  -1, /* tree 8 */
    -1,   1, -1,  -1, /* tree 9 */
  };
  const p4est_topidx_t ctt_offset[2+1] = {
    0,5,10,
  };
  /* for each corner, report the tree numbers it is attached to */
  const p4est_topidx_t corner_to_tree[10] = {
    0,2,4,6,8, /* corner 0 */
    1,3,5,7,9, /* corner 1 */
  };

  /* a given corner belong to multiple trees;
   for each tree, we report the index identifying the vertex location
   in the tree_to_vertex.
  e.g. here :
  - corner 0 is vertex 0
  - corner 1 is vertex 17
  For each entry in corner_to_tree, we report the location of the vertex in
  tree_to_vertex
  */
  const int8_t corner_to_corner[10] = {
    2, 2, 2, 2, 2,/* corner 0 (i.e vertex  0) */
    1, 1, 1, 1, 1,/* corner 1 (i.e vertex 17) */
  };


/* *INDENT-ON* */

  p4est_connectivity_t * conn =
    p4est_connectivity_new_copy (num_vertices, num_trees, num_corners,
				 vertices, tree_to_vertex,
				 tree_to_tree, tree_to_face,
				 tree_to_corner, ctt_offset,
				 corner_to_tree, corner_to_corner);
  
  printf("\nconnectivity good : %d\n",p4est_connectivity_is_valid (conn));
  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
  

  return conn;
} /* p4est_connectivity_new_icosahedron */


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

} /* p4est_geometry_connectivity_X */



typedef enum
{
  P4EST_GEOMETRY_BUILTIN_MAGIC = 0x65F2F8DE, /* should be different from P8EST_GEOMETRY_BUILTIN_MAGIC ? */
  P4EST_GEOMETRY_BUILTIN_ICOSAHEDRON
}
p4est_geometry_builtin_type_t;

typedef struct p4est_geometry_builtin_icosahedron
{
  p4est_geometry_builtin_type_t type;
  double                           a; /* size of an edge */
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
  double              a = 0.5*icosahedron->a; /* icosahedron half edge length */

  double              g = (1.0+sqrt(5.0))*0.5; /* *golden ratio */
  double              ga = a/g;

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
    2 ,  7,  0,  3,
    3 ,  8,  0,  4,
    4 ,  9,  0,  5,
    5 , 10,  0,  1,
    6 , 11,  2,  7,
    7 , 11,  3,  8,
    8 , 11,  4,  9,
    9 , 11,  5, 10,
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
    double theta1, theta2; /* angle for SLERP (interpolation) */
    
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
    double dot1   = n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2];
    double dot2   = n0[0]*n2[0] + n0[1]*n2[1] + n0[2]*n2[2];
    theta1 = acos(dot1/norme2);
    theta2 = acos(dot2/norme2);

    /* actual computation of bilinear slerp */
    for (j=0; j<3; ++j) {
      xyz[j] = 
	sin((1.0-eta_y)*theta2)/sin(theta2) *
	( sin((1.0-eta_x)*theta1)/sin(theta1)*n0[j]+
	  sin((    eta_x)*theta1)/sin(theta1)*n1[j]) +
	sin((    eta_y)*theta2)/sin(theta2) *
	( sin((1.0-eta_x)*theta1)/sin(theta1)*n2[j]+
	  sin((    eta_x)*theta1)/sin(theta1)*n3[j]);
    }

    if (dot1<0) {
      printf("####### %d %d %d %d | %g %g || %g %g %g || %g %g %g|\n",i0,i1,i2,i3, 
	     dot1, dot2,
	     n0[0],n0[1],n0[2],
	     n1[0],n1[1],n1[2]
	     );
    }

    printf("DEBUG : %g | %g %g %g | %g %g | %g %g | %g %g\n",norme2, 
	   xyz[0],xyz[1],xyz[2], 
	   theta1, theta2,
	   dot1, norme2,
	   eta_x, eta_y);
    //printf("rst %f %f %f\n",rst[0],rst[1],rst[2]);

  } /* end of bilinear slerp */

} /* p4est_geometry_icosahedron_X */

p4est_geometry_t   *
p4est_geometry_new_icosahedron (p4est_connectivity_t * conn, double a)
{
  p4est_geometry_builtin_t *builtin;
  struct p4est_geometry_builtin_icosahedron *icosahedron;

  builtin = P4EST_ALLOC_ZERO (p4est_geometry_builtin_t, 1);

  icosahedron = &builtin->p.icosahedron;
  icosahedron->type = P4EST_GEOMETRY_BUILTIN_ICOSAHEDRON;
  icosahedron->a = a;

  builtin->geom.name = "p4est_icosahedron";
  builtin->geom.user = conn;
  builtin->geom.X = p4est_geometry_icosahedron_X;
  //builtin->geom.X = p4est_geometry_connectivity_X;

  return (p4est_geometry_t *) builtin;

} /* p4est_geometry_new_icosahedron */

