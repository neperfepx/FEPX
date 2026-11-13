# Mesh reader fixtures (tests/io)

Hand-written for the parser tests (ROADMAP step 2.2); not a meaningful FE mesh (no geometry checks).

- `tiny.msh`: Neper msh file version 2.4 with two 10-node tets sharing a face (nodes 1–5 vertices,
  6–14 edge nodes, Gmsh order), preceded by a 0D and a 2D element (skipped by the reader). Elset
  ids 3 and 7 (elset indices 1 and 2), partitions 1 and 2, two nsets, one faset (face 1 2 3 of
  element id 3, as the outward-oriented Gmsh triangle `1 3 2 8 7 6`), one periodicity relation, out-of-order elset orientations (euler-bunge:active),
  groups (phases 1 and 2), physical names, a blank line between sections and an unknown section
  `$FutureSection` (skipped with a warning).
- `tiny.opt`: auxiliary value sections (`$ElsetCrss` with 1 and 2 values per line, out of order;
  `$ElementCrssSat`, out of order).
- `tiny.ori`: auxiliary `$ElementOrientations` (quaternion:passive, out of order); replaces the
  elset orientations. Preceded by `$ElsetCrySym` (triclinic, the mesh file's), as Neper 5 writes
  its orientation files (step 8.5).

The tests also read the regression meshes `tests/regression/cases/meshes/poly_part4.msh` (Neper,
`-part 4`) and `sx2.msh` (unpartitioned, 2 cells with one orientation); see that directory's `generate_meshes.sh`.
The error cases (`test_mesh_file_errors`) write modified copies of `tiny.msh`.
