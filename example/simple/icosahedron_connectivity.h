#ifndef ICOSAHEDRON_CONNECTIVITY_H
#define ICOSAHEDRON_CONNECTIVITY_H

p4est_connectivity_t *
p4est_connectivity_new_icosahedron ();

p4est_geometry_t   *
p4est_geometry_new_icosahedron (p4est_connectivity_t * conn, double R);

#endif /* !ICOSAHEDRON_CONNECTIVITY_H */
