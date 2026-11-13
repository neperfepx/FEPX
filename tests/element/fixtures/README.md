# Element fixtures

`tet10_reference.txt`: FEPX's 15-point tetrahedron quadrature rule and its 10-node tetrahedron
shape functions and derivatives (`src/element/tet10_mod.f90`), and the 7-point rule and 6-node
triangle of its faces. Format: `qpts 15`, then one line per quadrature point (xi1 xi2 xi3
weight); `points 30`, then one line per point (xi1 xi2 xi3, the 10 shape functions, then the 10
derivatives with respect to xi1, xi2 and xi3, Gmsh node order); points: the 15 quadrature points,
the 10 element nodes and 5 other points. Then `tri_qpts 7` (xi1 xi2 weight) and `tri_points 16`
(xi1 xi2, 6 values, 6 d/dxi1, 6 d/dxi2, Gmsh face node order) at the 7 face quadrature points,
the 6 face nodes and 3 other points. `tests/element/test_tet10.f90` checks that the rule and the
shape functions do not change (exactly, and to 1e-15).

`quad8_reference.txt`: FEPX's 8-node quadrilateral face shape functions and derivatives of the
20-node hexahedron (`src/element/hex20_mod.f90::hex20_face_shape`, `hex20_face_shape_deriv`) on
the reference square [-1, 1]², at the 8 face nodes (Gmsh order) and 4 other points:
`quad_points 12`, then xi1 xi2, 8 values, 8 d/dxi1, 8 d/dxi2. Used by
`tests/element/test_hex20.f90` (1e-15).

`geometry_reference.txt`: FEPX's element geometry (`src/element/element_mod.f90::
element_geometry_qpt`: det J and the derivatives of the shape functions with respect to x, y, z)
at the 15 quadrature points of three 10-node tetrahedra (scaled reference, affine, curved;
coordinates written in the file, Gmsh order), for `tests/element/test_element.f90` (relative
1e-13).

The three files were generated on 2026-10-25 (ROADMAP step 9.13, D91) with `dump_tet10.f90`,
`dump_quad8.f90` and `dump_geometry.f90` (this directory), built against the objects of a Release
build of `ac211915` (tag `v2-final-comparison`; `src/` unchanged since), gfortran 15.2.0, Open MPI
5.0.10:

```
B=<build dir>
mpif90 -O2 -I$B/mod dump_tet10.f90 $B/src/libfepx_core.a -o dump_tet10
./dump_tet10 > tet10_reference.txt            (the same for dump_quad8, dump_geometry)
```

They replace the fixtures `v2_tet10.txt`, `v2_quad8.txt` and `v2_geometry.txt`, the same values
written by v2-fixed's own routines (`quadrature_init`, `t10_shape_hpar`, `t10_deriv_hpar`,
`sf2d06`, `sf2d08`, `sfder_hpar`; in v2's node order and, for the quadrilateral, on v2's
reference square [0, 1]²), against which FEPX's ported rules and shape functions were checked
from step 3.1 to step 9.13 (to the same tolerances; the geometry differs from v2's by at most
6e-15). The v2 fixtures and their programs are in the tag `v2-final-comparison`. Regenerate only
if FEPX's quadrature rules, shape functions or geometry kernels change on purpose (CLAUDE.md
§2.5).

`hex_cube.msh`: a Neper mesh of the unit cube with 125 second-order hexahedra (Gmsh type 17, msh
version 2.4), used by `tests/element/test_hex20.f90` to check the Gmsh node order of the 20-node
hexahedron (ROADMAP step 3.2, D35). Generated on 2026-10-09 with Neper 5.0.1-520-devel and Gmsh
4.15.2:

```
neper -T -n 8 -morpho cube -id 1 -o cubes
neper -M cubes.tess -elttype hex -order 2 -rcl 2 -o hex_cube
```
