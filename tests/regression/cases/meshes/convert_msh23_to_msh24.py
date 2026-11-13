#!/usr/bin/env python3
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

"""Convert Neper msh files from version 2.3 to version 2.4 in place (step 9.17, ROADMAP D107).

The two versions differ in the node order of the $Fasets face lines: 2.3 writes each face's
nodes in reverse, 2.4 in Gmsh's native order. This changes the $MeshVersion line and reverses
the nodes of each face line (the element id stays first). Every other byte is kept, so nodes,
elements and partitions are unchanged, which a re-save by Neper does not promise for files it
did not write (a hand-written fixture fails to load).

Usage: convert_msh23_to_msh24.py <file.msh>...
"""
import sys


def convert(path):
    """Convert one file; it must be at version 2.3."""
    with open(path, encoding="ascii", newline="") as handle:
        lines = handle.read().split("\n")

    # The version line.
    i_version = lines.index("$MeshVersion") + 1
    if lines[i_version] != "2.3":
        sys.exit(f"{path}: version {lines[i_version]}, expected 2.3")
    lines[i_version] = "2.4"

    # The face lines, when the file has a $Fasets section: count, then label, face count, lines.
    if "$Fasets" in lines:
        i_line = lines.index("$Fasets") + 1
        num_fasets = int(lines[i_line])
        i_line += 1
        for _ in range(num_fasets):
            num_faces = int(lines[i_line + 1])
            i_line += 2
            for i_face in range(i_line, i_line + num_faces):
                fields = lines[i_face].split()
                lines[i_face] = " ".join([fields[0]] + fields[:0:-1])
            i_line += num_faces
        if lines[i_line] != "$EndFasets":
            sys.exit(f"{path}: $EndFasets not found where expected")

    with open(path, "w", encoding="ascii", newline="") as handle:
        handle.write("\n".join(lines))


for name in sys.argv[1:]:
    convert(name)
