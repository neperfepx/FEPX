#!/bin/bash
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Generate the meshes of the FEPX regression cases (ROADMAP step 0.3).
#
# The meshes are committed; this script documents how they were made. Regenerated on 2026-11-05
# with Neper 5.0.1-731-devel and Gmsh 4.15.2 built with Netgen (Ubuntu's gmsh package lacks Netgen,
# which Neper's 3D meshing requires), all at msh version 2.4 (step 9.20, ROADMAP D113; before, step
# 9.18 with 5.0.1-713-devel, D110), except the files named below as not generated here. 731 gives
# the same tessellations as 713 but other meshes, coarser at the same -rcl (sx2 263 tets with 713,
# poly 615, poly_periodic 363, tutorial_overwriting 5904). Run from this directory.
#
# Output:
#   sx2.tess, sx2.msh                single crystal in the unit cube meshed as a 2-cell
#                                    tessellation (-id 1), 145 second-order tets; both cells have
#                                    the orientation of Neper's random draw for -n 1 -id 1, so the
#                                    material is one crystal
#   sx2_part{2,4}.msh                the same mesh with Neper's partition table
#   poly.tess                        8-grain Voronoi polycrystal in the unit cube (-id 1)
#   poly_part{1,2,4}.msh             the same polycrystal mesh (454 second-order tets), with
#                                    Neper partition tables for 1, 2 and 4 partitions
#   poly_periodic.tess, poly_periodic_part{2,4}.msh  an
#                                    8-grain periodic Voronoi polycrystal of the unit cube
#                                    (-periodicity all: the grains cross the domain boundary,
#                                    whose faces are not planar), regularized (-reg 1), meshed
#                                    with second-order tets (239) and 2 and 4 partitions;
#                                    $Periodicity relations, read by FEPX (periodic regression
#                                    cases). Regenerated on 2026-11-03 (step 9.19, same Neper and
#                                    Gmsh) with -reg 1: without it the tessellation's tiny edges
#                                    gave sliver tets (median/min volume 2.8e5, 2149 tets) that
#                                    amplified rounding ~100x (results differed between CPUs);
#                                    with 731, median/min volume 5.6, quality
#                                    6*sqrt(2)*V/max_edge^3 >= 0.093 (median 0.271)
#   triaxial_cube.msh, dual_bcc_fcc.msh  not generated here: the meshes of v2's cases
#                                    LEGACY_triaxial_csr (= LEGACY_triaxial_clr) and
#                                    05_dual_bcc_fcc (2 grains, one phase and one partition each),
#                                    copied from the tag v2-final-comparison in step 9.13;
#                                    dual_bcc_fcc_part4.msh: its copy with 4 partitions, below
# Note: the partition tables of poly_part2/4.msh differ between runs of this script (Neper's
# partitioner is not deterministic); the geometry and the 1-partition mesh are reproducible.
set -euo pipefail

# Single crystal meshed as two cells with one orientation (the -n 1 -id 1 orientation).
neper -T -n 2 -id 1 -ori "rodrigues(0.919953402851,1.412193625247,-0.358468593221)" -o sx2
neper -M sx2.tess -order 2 -rcl 1 -o sx2
for num_parts in 2 4; do neper -M -loadmesh sx2.msh -part ${num_parts} -o sx2_part${num_parts}; done

# Polycrystal, one mesh per partition count.
neper -T -n 8 -id 1 -o poly
for num_parts in 1 2 4; do
  neper -M poly.tess -order 2 -rcl 2 -part ${num_parts} -o poly_part${num_parts}
done

# Periodic polycrystal (step 8.3b; regularized, step 9.19).
neper -T -n 8 -id 1 -periodicity all -reg 1 -o poly_periodic
for num_parts in 2 4; do
  neper -M poly_periodic.tess -order 2 -rcl 2 -part ${num_parts} -o poly_periodic_part${num_parts}
done

# Version (ROADMAP D107): the files not generated here (triaxial_cube, dual_bcc_fcc and v2's
# tutorial meshes) were written at msh version 2.3 and converted to 2.4 by
# convert_msh23_to_msh24.py (the version line and the face lines; every other byte kept).
#
# Partitions (ROADMAP D108, D110): FEPX runs a mesh on as many ranks as it has partitions. A mesh
# that runs on another number of ranks has a copy <mesh>_part<n>.msh with Neper's table; the tests
# pick it by name (cmake/FEPXTesting.cmake, fepx_mesh_for_ranks). tutorial_periodic_part4.msh
# (step 9.20) needs Neper 5.0.1-731-devel or later: earlier versions crash re-partitioning a
# periodic mesh (dev/neper_findings.md, A11). Its re-save leaves out tutorial_periodic.msh's
# $NSets and $Fasets (all empty sets; A10), which its boundary conditions do not use.
#
# Tutorial meshes (step 8.5, 2026-10-22, Neper 5.0.1-526-devel), from FEPX v2's tutorials
# (legacy/fepx2-fixed/doc/tutorials, legacy/fepx2-fixed/tutorials):
#   tutorial_simple.msh       simple_simulation and simple_rve_simulation (msh 2.2.3)
#   tutorial_periodic.msh     periodic_rve_simulation (2.2.3)
#   tutorial_triax.msh        triaxCSR and restart (the same file, 2.2.1)
#   tutorial_triax_dual.msh   triaxCLR (2.2.1, two phases in $Groups)
# are v2's own files with the two changes of msh version 2.3 for them (ROADMAP D33): the version
# line, and the swapped orientation label of $ElsetOrientations (active <-> passive, the same
# numbers). Neper's re-save (neper -M -loadmesh) writes the same nodes, elements, fasets and
# partitions but renames the boundary node sets and fasets (x0, ..., x0y0z0 become f1, ...,
# f1f2f3), which the boundary conditions name, so it is not used. v2's files are in the tag
# v2-final-comparison (step 9.13): V2_DIR is a checkout of its legacy/fepx2-fixed, e.g.
#   git worktree add /tmp/v2final v2-final-comparison; V2_DIR=/tmp/v2final/legacy/fepx2-fixed
v2_dir=${V2_DIR:?set V2_DIR to legacy/fepx2-fixed of the tag v2-final-comparison}
to_msh23() {
  sed -e '5s/^2\.2\.[13]$/2.3/' \
      -e '/^\$ElsetOrientations$/{n;s/:active$/:PASSIVE_TMP/;s/:passive$/:active/;s/:PASSIVE_TMP$/:passive/}' \
      "$1" > "$2"
}
to_msh23 ${v2_dir}/doc/tutorials/simple_simulation/simulation.msh tutorial_simple.msh
to_msh23 ${v2_dir}/doc/tutorials/periodic_rve_simulation/simulation.msh tutorial_periodic.msh
to_msh23 ${v2_dir}/tutorials/triaxCSR/simulation.msh tutorial_triax.msh
to_msh23 ${v2_dir}/tutorials/triaxCLR/simulation.msh tutorial_triax_dual.msh
# tutorial_cruciform.{tess,msh} (step 9.18; with Neper 5.0.1-731-devel, whose faces of the cut
# surfaces all point outward, dev/neper_findings.md A12): the cruciform
# specimen, a plate 1 x 1 x 0.1 with its corners cut by circular arcs and a central hole, the
# dimensions fitted to v2's file (which came without commands), run by the manual only
# (tutorials.txt tier "manual"). The cuts remove cells: -n 120 leaves 88 grains, 12737 tets
# (14986 with 5.0.1-713-devel). (-n 130 does not mesh with 713: "Meshing of face 258 failed".)
neper -T -n 120 -id 1 -domain "cube(1,1,0.1)" \
  -transform "cut(cylinder(0.05,0.05,0.05,0,0,1,0.3),cylinder(0.95,0.05,0.05,0,0,1,0.3),cylinder(0.05,0.95,0.05,0,0,1,0.3),cylinder(0.95,0.95,0.05,0,0,1,0.3),cylinder(0.5,0.5,0.05,0,0,1,0.2))" \
  -o tutorial_cruciform
neper -M tutorial_cruciform.tess -order 2 -rcl 1 -part 4 -o tutorial_cruciform

# tutorial_overwriting.{tess,msh,ori} (3835 tets): the overwriting tutorial's mesh regenerated with its own
# generate_mesh.sh commands (its shipped msh 2.2.1 file has only the face node sets, and
# v2-fixed's minimal conditions also fix corners, CHANGES.md C27). Two changes for Neper 5: the
# orientation file holds its first 50 rows (the shipped orifile has 100 for 50 grains; Neper 4.2
# took the first 50, as the shipped simulation.ori shows; Neper 5 stops with a segmentation
# fault), and -orispread replaces -oridistrib (not read by Neper 5).
head -50 ${v2_dir}/tutorials/overwriting/orifile > orifile50
neper -T -n 50 -reg 1 -group "id<=25?1:2" -ori "file(orifile50,des=rodrigues)" \
  -orispread "normal(5)" -oridescriptor "rodrigues:active" -morpho gg -o tutorial_overwriting
neper -M tutorial_overwriting.tess -order 2 -format "msh,ori" -part 4 -o tutorial_overwriting
rm -f orifile50

# The partition copies of the files above (see "Partitions").
neper -M -loadmesh dual_bcc_fcc.msh -part 4 -o dual_bcc_fcc_part4
neper -M -loadmesh tutorial_simple.msh -part 4 -o tutorial_simple_part4
neper -M -loadmesh tutorial_periodic.msh -part 4 -o tutorial_periodic_part4
