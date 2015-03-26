#ifdef P4_TO_P8
#include <p8est_connectivity.h>
#else
#include <p4est_connectivity.h>
#endif

p4est_connectivity_t *
p4est_connectivity_new_ring (int num_trees_radial,
			     int num_trees_orthoradial,
			     double rMin,
			     double rMax)
{
  p4est_topidx_t num_vertices = (num_trees_radial+1) * (num_trees_orthoradial);
  p4est_topidx_t num_trees = num_trees_radial * num_trees_orthoradial;
  p4est_topidx_t num_ctt = 0;
  double        *vertices = (double *) malloc(num_vertices * 3 * sizeof(double));

  p4est_topidx_t *tree_to_vertex = (p4est_topidx_t *) malloc(num_trees * 4 * sizeof(p4est_topidx_t));
  p4est_topidx_t *tree_to_tree   = (p4est_topidx_t *) malloc(num_trees * 4 * sizeof(p4est_topidx_t));
  int8_t         *tree_to_face   = (int8_t *)         malloc(num_trees * 4 * sizeof(int8_t));

  int iVertex = 0;
  int iRadial;
  int iOrthoradial;
  int iTree;

  // fill vertices array
  for (iRadial=0; iRadial<num_trees_radial+1; iRadial++) {

    double radius = rMin + iRadial*(rMax-rMin)/num_trees_radial;

    for (iOrthoradial=0; iOrthoradial<num_trees_orthoradial; iOrthoradial++) {
      vertices[3*iVertex+0] = radius*cos(2*M_PI/num_trees_orthoradial*iOrthoradial);
      vertices[3*iVertex+1] = radius*sin(2*M_PI/num_trees_orthoradial*iOrthoradial);
      vertices[3*iVertex+2] = 0;

      iVertex++;

    }
  }

  // fill tree to vertex
  iTree = 0;
  for (iRadial=0; iRadial<num_trees_radial; iRadial++) {
    for (iOrthoradial=0; iOrthoradial<num_trees_orthoradial; iOrthoradial++) {
      
      p4est_topidx_t left_corner_vertex = iOrthoradial + num_trees_orthoradial*iRadial;

      tree_to_vertex[4*iTree+0] = left_corner_vertex;
      tree_to_vertex[4*iTree+1] = left_corner_vertex+  num_trees_orthoradial;
      tree_to_vertex[4*iTree+2] = left_corner_vertex+1;
      tree_to_vertex[4*iTree+3] = left_corner_vertex+1+num_trees_orthoradial;

      // need to modify this when we cross iTree % num_trees_orthoradial == 0
      if (iTree % num_trees_orthoradial == num_trees_orthoradial-1) {
	tree_to_vertex[4*iTree+2] -= num_trees_orthoradial;
	tree_to_vertex[4*iTree+3] -= num_trees_orthoradial;
      }

      /* printf("Vertex of tree %d : %d %d %d %d\n",iTree, */
      /* 	     tree_to_vertex[4*iTree+0], */
      /* 	     tree_to_vertex[4*iTree+1], */
      /* 	     tree_to_vertex[4*iTree+2], */
      /* 	     tree_to_vertex[4*iTree+3] */
      /* 	     ); */

      iTree++;
    }
  }  

  // fill tree to tree
  iTree=0;
  for (iRadial=0; iRadial<num_trees_radial; iRadial++) {
    for (iOrthoradial=0; iOrthoradial<num_trees_orthoradial; iOrthoradial++) {

      tree_to_tree[4*iTree+0] = iTree-num_trees_orthoradial < 0 ? iTree : iTree-num_trees_orthoradial;
      tree_to_tree[4*iTree+1] = iTree+num_trees_orthoradial >= num_trees ? iTree : iTree+num_trees_orthoradial;
      tree_to_tree[4*iTree+2] = iTree-1;
      tree_to_tree[4*iTree+3] = iTree+1;

      // need to modify this when we cross iTree % num_trees_orthoradial == 0
      if (iTree % num_trees_orthoradial == 0) {
	tree_to_tree[4*iTree+2] += num_trees_orthoradial;
      } 
      if (iTree % num_trees_orthoradial == num_trees_orthoradial-1) {
	tree_to_tree[4*iTree+3] -= num_trees_orthoradial;
      }

      /* printf("tree to tree %d : %d %d %d %d\n",iTree, */
      /* 	     tree_to_tree[4*iTree+0], */
      /* 	     tree_to_tree[4*iTree+1], */
      /* 	     tree_to_tree[4*iTree+2], */
      /* 	     tree_to_tree[4*iTree+3] */
      /* 	     ); */

      iTree++;
    }
  }

  // tree_to_face
  iTree=0;
  for (iRadial=0; iRadial<num_trees_radial; iRadial++) {
    for (iOrthoradial=0; iOrthoradial<num_trees_orthoradial; iOrthoradial++) {
      
      tree_to_face[4*iTree+0] = 1;
      tree_to_face[4*iTree+1] = 0;
      tree_to_face[4*iTree+2] = 3;
      tree_to_face[4*iTree+3] = 2;

      // need to modifiy this when tree is on edge (inner or outer)
      if (iTree < num_trees_orthoradial) { // inner border
	tree_to_face[4*iTree+0] = 0;
      } 
      if (iTree+num_trees_orthoradial >= num_trees) { // outer border
	tree_to_face[4*iTree+1] = 1;
      } 

      iTree++;
    }
  }

  p4est_connectivity_t * conn =
    p4est_connectivity_new_copy (num_vertices, num_trees, 0,
				 vertices, tree_to_vertex,
				 tree_to_tree, tree_to_face,
				 NULL, &num_ctt, NULL, NULL);
  
  //printf("\nconnectivity good : %d\n",p4est_connectivity_is_valid (conn));
  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

  free(vertices);
  free(tree_to_vertex);
  free(tree_to_tree);
  free(tree_to_face);

  return conn;

} /* p4est_connectivity_new_ring */

