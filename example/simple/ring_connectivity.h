#ifndef RING_CONNECTIVITY_H
#define RING_CONNECTIVITY_H

p4est_connectivity_t *
p4est_connectivity_new_ring (int num_trees_radial,
			     int num_trees_orthoradial,
			     double rMin,
			     double rMax);


#endif /* !RING_CONNECTIVITY_H */
