#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Write structured quarter-annulus meshes (Neper msh 2.4 layout) for the Lamé verification.

Source : FEPX verification tests (ROADMAP step 4.4, dev/theory/addendum_lame.md §3, D43 (4)).
         Development tool (CLAUDE.md §6). Neper does not yet write hex face sets (ROADMAP Q12),
         so these meshes are generated here, for tets and hexes alike.

Geometry: a <= r <= b, 0 <= theta <= pi/2, 0 <= z <= h, mapped from the unit cube of parameters
(s, t, u): r = a + (b - a) s, theta = (pi/2) t, z = h u. Cells: nr x nt x nz. Every node is the
image of a point of the half-step parameter grid, so mid-nodes (and the face and cell centres
used by the tets' edge nodes) lie on the curved geometry and are shared exactly between elements.

Elements (Gmsh order, ROADMAP D35):
  hex20 (type 17): corners 1-4 at u = 0 and 5-8 at u = 1, counterclockwise in (s, t), then the
                   edge nodes of edges 1-2, 1-4, 1-5, 2-3, 2-6, 3-4, 3-7, 4-8, 5-6, 5-8, 6-7, 7-8;
  tet10 (type 11): each cell split into the 6 tetrahedra along its diagonal from parameter corner
                   (0,0,0) to (1,1,1) (one per ordering of the three axes), corners ordered so that
                   det J > 0, then the edge nodes of edges 1-2, 2-3, 3-1, 4-1, 4-3, 4-2.

Sections written: $MeshFormat 2.2 0 8, $MeshVersion 2.4, $Nodes, $Elements (3D elements only;
tags: elset 1, elementary 1, partition 0), $NSets (x0: theta = pi/2, y0: theta = 0, z0, z1,
inner: r = a, outer: r = b), $PhysicalNames, $ElsetCrySym (triclinic), $ElsetOrientations
(one elset, rodrigues:passive 0 0 0). Coordinates with 12 decimals, as Neper.

Usage: generate_annulus.py --type tet10|hex20 --nr N --nt N --nz N [--a A --b B --h H] OUT.msh
"""

import argparse
import itertools
import math

# hex20 edges (corner pairs, 1-based, Gmsh order) and the 6 Kuhn tetrahedra of a cube (corner
# offsets in (s, t, u), the path from (0,0,0) to (1,1,1) along one ordering of the axes).
HEX_EDGES = [(1, 2), (1, 4), (1, 5), (2, 3), (2, 6), (3, 4), (3, 7), (4, 8), (5, 6), (5, 8),
             (6, 7), (7, 8)]
HEX_CORNERS = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1),
               (0, 1, 1)]
TET_EDGES = [(1, 2), (2, 3), (3, 1), (4, 1), (4, 3), (4, 2)]


def kuhn_tets():
    """The 6 tetrahedra of the unit cube along its main diagonal, as corner offsets."""
    tets = []
    for perm in itertools.permutations(range(3)):
        path = [(0, 0, 0)]
        point = [0, 0, 0]
        for axis in perm:
            point[axis] = 1
            path.append(tuple(point))
        tets.append(path)
    return tets


def main(argv=None):
    """Command-line interface."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("out")
    parser.add_argument("--type", choices=("tet10", "hex20"), required=True)
    parser.add_argument("--nr", type=int, required=True)
    parser.add_argument("--nt", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--a", type=float, default=1.0)
    parser.add_argument("--b", type=float, default=2.0)
    parser.add_argument("--h", type=float, default=0.25)
    args = parser.parse_args(argv)
    nr, nt, nz = args.nr, args.nt, args.nz

    # Node ids on the half-step grid (assigned on first use) and their coordinates.
    node_id = {}
    coords = []

    def node(i, j, k):
        """Id of the node at half-step grid point (i, j, k), created if new."""
        key = (i, j, k)
        if key not in node_id:
            s, t, u = i / (2 * nr), j / (2 * nt), k / (2 * nz)
            r = args.a + (args.b - args.a) * s
            theta = 0.5 * math.pi * t
            coords.append((r * math.cos(theta), r * math.sin(theta), args.h * u))
            node_id[key] = len(coords)
        return node_id[key]

    def mid(p, q):
        """Half-step grid point halfway between two grid points (doubled coordinates)."""
        return tuple((p[d] + q[d]) // 2 for d in range(3))

    def signed_volume(ids):
        """6 times the signed volume of the tetrahedron of four node ids."""
        x = [coords[i - 1] for i in ids]
        e = [[x[m][d] - x[0][d] for d in range(3)] for m in (1, 2, 3)]
        return (e[0][0] * (e[1][1] * e[2][2] - e[1][2] * e[2][1])
                - e[0][1] * (e[1][0] * e[2][2] - e[1][2] * e[2][0])
                + e[0][2] * (e[1][0] * e[2][1] - e[1][1] * e[2][0]))

    # Elements, cell by cell (grid points in doubled coordinates: cell corner (2i, 2j, 2k)).
    elements = []
    for k in range(nz):
        for j in range(nt):
            for i in range(nr):
                base = (2 * i, 2 * j, 2 * k)
                if args.type == "hex20":
                    corners = [tuple(base[d] + 2 * c[d] for d in range(3)) for c in HEX_CORNERS]
                    ids = [node(*p) for p in corners]
                    ids += [node(*mid(corners[e1 - 1], corners[e2 - 1])) for e1, e2 in HEX_EDGES]
                    elements.append((17, ids))
                else:
                    for tet in kuhn_tets():
                        corners = [tuple(base[d] + 2 * c[d] for d in range(3)) for c in tet]
                        ids = [node(*p) for p in corners]
                        if signed_volume(ids) < 0.0:
                            corners[1], corners[2] = corners[2], corners[1]
                            ids[1], ids[2] = ids[2], ids[1]
                        ids += [node(*mid(corners[e1 - 1], corners[e2 - 1]))
                                for e1, e2 in TET_EDGES]
                        elements.append((11, ids))

    # Node sets by position on the grid (doubled coordinates).
    nsets = {"x0": [], "y0": [], "z0": [], "z1": [], "inner": [], "outer": []}
    for (i, j, k), nid in sorted(node_id.items(), key=lambda item: item[1]):
        if j == 2 * nt:
            nsets["x0"].append(nid)
        if j == 0:
            nsets["y0"].append(nid)
        if k == 0:
            nsets["z0"].append(nid)
        if k == 2 * nz:
            nsets["z1"].append(nid)
        if i == 0:
            nsets["inner"].append(nid)
        if i == 2 * nr:
            nsets["outer"].append(nid)

    # Write the file.
    with open(args.out, "w") as out:
        out.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$MeshVersion\n2.4\n$EndMeshVersion\n")
        out.write(f"$Nodes\n{len(coords)}\n")
        for nid, (x, y, z) in enumerate(coords, start=1):
            out.write(f"{nid} {x:.12f} {y:.12f} {z:.12f}\n")
        out.write("$EndNodes\n")
        out.write(f"$Elements\n{len(elements)}\n")
        for eid, (etype, ids) in enumerate(elements, start=1):
            out.write(f"{eid} {etype} 3 1 1 0 " + " ".join(str(n) for n in ids) + "\n")
        out.write("$EndElements\n")
        out.write(f"$NSets\n{len(nsets)}\n")
        for label, ids in nsets.items():
            out.write(f"{label}\n{len(ids)}\n" + "".join(f"{n}\n" for n in ids))
        out.write("$EndNSets\n")
        out.write('$PhysicalNames\n1\n3 1 "poly1"\n$EndPhysicalNames\n')
        out.write("$ElsetCrySym\ntriclinic\n$EndElsetCrySym\n")
        out.write("$ElsetOrientations\n1 rodrigues:passive\n1 0.000000000000 0.000000000000 "
                  "0.000000000000\n$EndElsetOrientations\n")


if __name__ == "__main__":
    main()
