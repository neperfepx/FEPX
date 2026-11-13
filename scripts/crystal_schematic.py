#!/usr/bin/env python3
"""Schematics of crystal cells with slip systems: hcp alone, or a bcc and an hcp cell side by side.

hcp: the hexagonal prism (a = 1, c/a = 1.75) with FEPX's three slip families,
  - a basal plane {0001} (green): the bottom hexagon;
  - a prismatic plane {10-10} (red): one side face;
  - a first-order pyramidal plane {10-11} (blue): the trapezoid from the bottom edge of a side
    face to the long diagonal of the top hexagon parallel to that edge;
  - an <a> direction <11-20>: the bottom edge shared by the basal and prismatic planes;
  - two <c+a> directions <11-23>: the two slanted edges of the pyramidal trapezoid (both lie in
    that plane: FEPX's two pyramidal <c+a> systems of each {10-11} plane).
bcc: the cubic cell (a = 1, corner atoms and the body-centre atom) with FEPX's default slip family
  {110}<111>,
  - a {110} plane (yellow): the diagonal plane (110), through the corners (1,0,0), (0,1,0),
    (0,1,1), (1,0,1) and the body centre (it faces the viewer: (1-10) is nearly edge-on in this
    view);
  - two <111> directions: the body diagonals [-11-1], from (1,0,1) to (0,1,0), and [-111], from
    (1,0,0) to (0,1,1), the two that lie in (110) (FEPX's two {110}<111> systems of that plane).
The families are those of FEPX (src/material/crystal_mod.f90: hcp {0001}<11-20>, {10-10}<11-20>,
{10-11}<11-23>; bcc {110}<111>).  Orthographic projection, the same view for both cells; no
labels.  The planes are translucent and drawn in depth order: the edges and atoms behind a plane
show through it, those in front of it are drawn over it.

Usage: crystal_schematic.py hcp|bcc+hcp <output.svg>
The website's images are made, from the website's top directory, with
  python3 scripts/crystal_schematic.py hcp hcp.svg
  inkscape hcp.svg --export-type=png --export-filename=imgs/gallery/crystal-hcp-schematic.png -w 800
(not shown on the website since the gallery's crystal series was cut to four images; kept until
the release, as the other replaced images), and
  python3 scripts/crystal_schematic.py bcc+hcp bcc-hcp.svg
  inkscape bcc-hcp.svg --export-type=png --export-filename=imgs/gallery/crystal-bcc-hcp-schematic.png -w 800
(the landing page's "Crystal behavior" card).
"""

import math
import sys

# hcp lattice: a = 1 (hexagon side = corner distance from the axis), c/a = 1.75 (exaggerated from
# ~1.6 so the cell does not look stout); bcc lattice: a = 1
hcp_a = 1.0
hcp_c = 1.75 * hcp_a
bcc_a = 1.0

# View: azimuth about the vertical (c, or z) axis and elevation above the horizontal plane [degrees]
view_azimuth = -60.0
view_elevation = 18.0

# Colours (light fills) and drawing sizes [px]
colour_basal = "#cfeec0"
colour_prismatic = "#f0b0a8"
colour_pyramidal = "#b3d7f7"
colour_bcc_plane = "#f7e3a1"
colour_arrow = "#9b2335"
colour_line = "#000000"
image_size = 800
line_width = 4.0
atom_radius = 13.5
arrow_width = 5.0
# Plane opacity (the edges and atoms behind a plane show through it)
plane_opacity = 0.75
# Number of pieces each edge is cut into for the depth order (an edge may cross a plane)
edge_pieces = 64

# Drawing extents on the 800 px image: hcp alone spans 58 % of the image (a little smaller than
# the Neper cubes of the other landing-page cards, which span about 60-65 % of theirs).  The pair:
# the hcp cell at the upper left, the bcc cell at the lower right, over it (drawn after it; the
# hcp cell shows through it, unless `bcc_silhouette`, a white silhouette that hides the hcp cell
# behind the bcc cell, is set); heights as fractions of the image, the
# bcc cell's centre offset from the hcp cell's (right, down) as fractions of the image; the pair
# is then centred and scaled to span at most `pair_span` of the image
single_span = 0.58
pair_hcp_height = 0.50
pair_bcc_height = 0.38
pair_offset = (0.36, 0.24)
pair_span = 0.76
# The bcc silhouette (off: the hcp cell shows through the bcc cell) and its white halo [px] (a gap
# between the cells where they overlap)
bcc_silhouette = False
halo_width = 10.0


def project(point):
  """Orthographic projection: (screen x, screen y, depth); larger depth = nearer the viewer."""
  azimuth = math.radians(view_azimuth)
  elevation = math.radians(view_elevation)
  x, y, z = point
  # Rotate about the vertical axis so that the view direction lies in the x-z plane
  x_rotated = x * math.cos(azimuth) - y * math.sin(azimuth)
  y_rotated = x * math.sin(azimuth) + y * math.cos(azimuth)
  # Tilt: screen up mixes z and the in-plane coordinate toward the viewer
  screen_x = y_rotated
  screen_y = z * math.cos(elevation) - x_rotated * math.sin(elevation)
  depth = x_rotated * math.cos(elevation) + z * math.sin(elevation)
  return screen_x, screen_y, depth


def hcp_cell(pyramidal_face, prismatic_face, basal_on_top):
  """hcp cell: (atoms, edges, planes as (polygon, fill, opacity), arrows as (start, end)).

  pyramidal_face  side face (0-5) from whose bottom edge the pyramidal plane rises
  prismatic_face  side face (0-5) shaded as the prismatic plane
  basal_on_top    shade the top hexagon as the basal plane (else the bottom one); the <a> arrow
                  is on the edge the basal and prismatic planes share
  """

  def corner(index, height):
    """Corner `index` (0-5) of the hexagon at `height`; corners at 30 + 60*index degrees."""
    angle = math.radians(30.0 + 60.0 * index)
    return (hcp_a * math.cos(angle), hcp_a * math.sin(angle), height)

  # Corners of the bottom (z = 0) and top (z = c) hexagons
  bottom = [corner(i, 0.0) for i in range(6)]
  top = [corner(i, hcp_c) for i in range(6)]
  # Side face k lies between corners k and k+1; its outward normal is at 60*(k+1) degrees.  The
  # pyramidal plane rises from the bottom edge of face k (corners k and k+1) to the long diagonal
  # of the top hexagon parallel to that edge (corners k-1 and k+2): (10-11) has intercepts
  # a/2*sqrt(3) along its basal normal at z = 0 and 0 at z = c.  In this view face 0 faces the
  # viewer, faces 5 and 1 are the front-left and front-right faces.
  k = pyramidal_face
  basal_plane = top if basal_on_top else bottom
  prismatic_plane = [bottom[prismatic_face], bottom[(prismatic_face + 1) % 6],
                     top[(prismatic_face + 1) % 6], top[prismatic_face]]
  pyramidal_plane = [bottom[k], bottom[(k + 1) % 6], top[(k + 2) % 6], top[(k - 1) % 6]]
  # In drawing order, back to front where they overlap on screen (the depth order of the edges
  # and atoms is handled in `draw`): viewed from above, a bottom basal plane is behind the
  # pyramidal plane, a top one in front of it; the prismatic plane does not overlap the
  # pyramidal plane for the faces used here
  basal = (basal_plane, colour_basal, plane_opacity)
  prismatic = (prismatic_plane, colour_prismatic, plane_opacity)
  pyramidal = (pyramidal_plane, colour_pyramidal, plane_opacity)
  planes = [pyramidal, prismatic, basal] if basal_on_top else [basal, prismatic, pyramidal]
  # <a> along the edge of the prismatic face on the basal plane; <c+a> along the two slanted
  # edges of the pyramidal plane, from the top corner k+2 down to the bottom corner k+1 and from
  # the top corner k-1 down to the bottom corner k
  arrows = [(basal_plane[(prismatic_face + 1) % 6], basal_plane[prismatic_face]),
            (top[(k + 2) % 6], bottom[(k + 1) % 6]),
            (top[(k - 1) % 6], bottom[k])]
  # Edges of the prism: the two hexagons and the six vertical edges
  edges = [(bottom[i], bottom[(i + 1) % 6]) for i in range(6)]
  edges += [(top[i], top[(i + 1) % 6]) for i in range(6)]
  edges += [(bottom[i], top[i]) for i in range(6)]
  return bottom + top, edges, planes, arrows


def bcc_cell():
  """bcc cell: (atoms, edges, planes as (polygon, fill, opacity), arrows as (start, end))."""
  # Corners (i, j, k) in {0, a}^3 and the body centre
  corners = [(i * bcc_a, j * bcc_a, k * bcc_a) for i in (0, 1) for j in (0, 1) for k in (0, 1)]
  centre = (0.5 * bcc_a, 0.5 * bcc_a, 0.5 * bcc_a)
  # The twelve edges: pairs of corners that differ in one coordinate
  edges = [(p, q) for p in corners for q in corners
           if p < q and sum(1 for u, v in zip(p, q) if u != v) == 1]
  # (110) plane through (1,0,0), (0,1,0), (0,1,1), (1,0,1); the body diagonals from (1,0,1) to
  # (0,1,0), [-11-1], and from (1,0,0) to (0,1,1), [-111], lie in it
  a = bcc_a
  plane = [(a, 0.0, 0.0), (0.0, a, 0.0), (0.0, a, a), (a, 0.0, a)]
  arrows = [((a, 0.0, a), (0.0, a, 0.0)), ((a, 0.0, 0.0), (0.0, a, a))]
  return corners + [centre], edges, [(plane, colour_bcc_plane, plane_opacity)], arrows


def screen_box(atoms):
  """Projected extent (x min, x max, y min, y max) of a cell's atoms."""
  projected = [project(p) for p in atoms]
  xs = [p[0] for p in projected]
  ys = [p[1] for p in projected]
  return min(xs), max(xs), min(ys), max(ys)


def draw(cell, scale, x_centre, y_centre, occlude=False):
  """SVG parts of a cell drawn at `scale` [px per unit], its projected box centred at (x, y);
  with `occlude`, preceded by its white silhouette (the convex hull of its atoms, widened by the
  atom radius and `halo_width`), which hides what was drawn before behind it.

  Depth order (painter's algorithm with translucent planes): every atom and every piece of an
  edge is drawn after the planes behind it and before the planes in front of it, where a plane is
  "in front" of a point when, at the point's screen position, the plane is nearer the viewer
  (the projection is orthographic, so a plane's depth is linear in the screen coordinates).
  Only the planes that overlap the point on screen (within the atom radius or the line width)
  count; the arrows, which lie in the planes, are drawn last.
  """
  atoms, edges, planes, arrows = cell
  x_min, x_max, y_min, y_max = screen_box(atoms)

  def to_pixels(point):
    """Pixel coordinates (x right, y down) of a 3D point."""
    screen_x, screen_y, _ = project(point)
    return (x_centre + scale * (screen_x - 0.5 * (x_min + x_max)),
            y_centre - scale * (screen_y - 0.5 * (y_min + y_max)))

  def depth(point):
    """Depth of a 3D point [px] (larger = nearer the viewer)."""
    return scale * project(point)[2]

  # Each plane: its pixel polygon and the coefficients of its depth, d = c0 + c1 x + c2 y in
  # pixels (from three of its corners)
  plane_data = []
  for polygon, fill, opacity in planes:
    pixels = [to_pixels(p) for p in polygon]
    depths = [depth(p) for p in polygon]
    (xa, ya), (xb, yb), (xc, yc) = pixels[0], pixels[1], pixels[2]
    determinant = (xb - xa) * (yc - ya) - (xc - xa) * (yb - ya)
    c1 = ((depths[1] - depths[0]) * (yc - ya) - (depths[2] - depths[0]) * (yb - ya)) / determinant
    c2 = ((xb - xa) * (depths[2] - depths[0]) - (xc - xa) * (depths[1] - depths[0])) / determinant
    c0 = depths[0] - c1 * xa - c2 * ya
    plane_data.append((pixels, (c0, c1, c2), fill, opacity))

  def near_polygon(x, y, pixels, margin):
    """True if the pixel (x, y) is inside the polygon or within `margin` px of its outline."""
    inside = False
    for i in range(len(pixels)):
      (x1, y1), (x2, y2) = pixels[i], pixels[(i + 1) % len(pixels)]
      if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
        inside = not inside
      # Distance to the side
      dx, dy = x2 - x1, y2 - y1
      t = max(0.0, min(1.0, ((x - x1) * dx + (y - y1) * dy) / (dx * dx + dy * dy)))
      if math.hypot(x - x1 - t * dx, y - y1 - t * dy) <= margin:
        return True
    return inside

  def slot(point, margin):
    """Number of planes (in drawing order) to draw before the point: one past the last plane
    that overlaps it on screen and lies behind it (a point in a plane is in front of it)."""
    x, y = to_pixels(point)
    d = depth(point)
    last_behind = -1
    for index, (pixels, (c0, c1, c2), _, _) in enumerate(plane_data):
      if near_polygon(x, y, pixels, margin) and c0 + c1 * x + c2 * y < d + 1e-6 * scale:
        last_behind = index
    return last_behind + 1

  # Edge pieces and atoms, by slot: layers[s] is drawn after the first s planes
  layers = [[] for _ in range(len(planes) + 1)]
  for start, end in edges:
    # Cut the edge into pieces; consecutive pieces in the same slot form one line
    points = [tuple(s + (e - s) * i / edge_pieces for s, e in zip(start, end))
              for i in range(edge_pieces + 1)]
    piece_slots = [slot(tuple(0.5 * (p + q) for p, q in zip(points[i], points[i + 1])),
                        line_width / 2) for i in range(edge_pieces)]
    first = 0
    for i in range(1, edge_pieces + 1):
      if i == edge_pieces or piece_slots[i] != piece_slots[first]:
        (x1, y1), (x2, y2) = to_pixels(points[first]), to_pixels(points[i])
        layers[piece_slots[first]].append(
          '<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="%.1f" '
          'stroke-linecap="round"/>' % (x1, y1, x2, y2, colour_line, line_width))
        first = i
  for point in sorted(atoms, key=depth):
    x, y = to_pixels(point)
    layers[slot(point, atom_radius)].append(
      '<circle cx="%.2f" cy="%.2f" r="%.1f" fill="#000000"/>' % (x, y, atom_radius))

  # Silhouette: convex hull of the atoms' pixels (monotone chain), filled and stroked in white
  parts = []
  if occlude:
    points = sorted(set((round(x, 3), round(y, 3)) for x, y in (to_pixels(p) for p in atoms)))

    def cross(o, a, b):
      """z component of (a - o) x (b - o)."""
      return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for point in points:
      while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
        lower.pop()
      lower.append(point)
    for point in reversed(points):
      while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
        upper.pop()
      upper.append(point)
    hull = lower[:-1] + upper[:-1]
    parts.append('<polygon points="%s" fill="#ffffff" stroke="#ffffff" stroke-width="%.1f" '
                 'stroke-linejoin="round"/>'
                 % (" ".join("%.2f,%.2f" % p for p in hull), 2.0 * (atom_radius + halo_width)))

  # Assemble: layer 0, plane 1, layer 1, plane 2, ...
  parts += layers[0]
  for index, (pixels, _, fill, opacity) in enumerate(plane_data):
    coordinates = " ".join("%.2f,%.2f" % p for p in pixels)
    parts.append('<polygon points="%s" fill="%s" fill-opacity="%.2f" stroke="none"/>'
                 % (coordinates, fill, opacity))
    parts += layers[index + 1]

  # Arrows on top, double-headed, shortened at each end to clear the atoms
  inset = atom_radius + 6
  for start, end in arrows:
    (x1, y1), (x2, y2) = to_pixels(start), to_pixels(end)
    length = math.hypot(x2 - x1, y2 - y1)
    ux, uy = (x2 - x1) / length, (y2 - y1) / length
    x1, y1, x2, y2 = x1 + inset * ux, y1 + inset * uy, x2 - inset * ux, y2 - inset * uy
    parts.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="%.1f" '
                 'marker-start="url(#head)" marker-end="url(#head)"/>'
                 % (x1, y1, x2, y2, colour_arrow, arrow_width))
  return parts


# Cells and their placement
mode = sys.argv[1]
if mode == "hcp":
  # Pyramidal plane on the front-right face's bottom edge, prismatic plane the front-left face,
  # basal plane the bottom hexagon
  hcp = hcp_cell(pyramidal_face=1, prismatic_face=5, basal_on_top=False)
  # hcp alone, its larger projected side spanning `single_span` of the image, centred
  x_min, x_max, y_min, y_max = screen_box(hcp[0])
  scale = single_span * image_size / max(x_max - x_min, y_max - y_min)
  body = draw(hcp, scale, image_size / 2, image_size / 2)
elif mode == "bcc+hcp":
  # hcp at the upper left, bcc at the lower right over it.  The hcp planes are kept away from
  # the bcc cell: the pyramidal plane on the front-left face's bottom edge, the prismatic plane
  # the front-right face (its lower part behind the bcc cell), the basal plane the top hexagon.
  # First at the nominal sizes and offset, then the union of the two projected boxes centred
  # and scaled to `pair_span`
  hcp = hcp_cell(pyramidal_face=5, prismatic_face=1, basal_on_top=True)
  bcc = bcc_cell()
  cells = (hcp, bcc)
  boxes = [screen_box(cell[0]) for cell in cells]
  heights = [pair_hcp_height * image_size, pair_bcc_height * image_size]
  scales = [h / (b[3] - b[2]) for h, b in zip(heights, boxes)]
  widths = [s * (b[1] - b[0]) for s, b in zip(scales, boxes)]
  centres = [(0.0, 0.0), (pair_offset[0] * image_size, pair_offset[1] * image_size)]
  left = min(c[0] - w / 2 for c, w in zip(centres, widths))
  right = max(c[0] + w / 2 for c, w in zip(centres, widths))
  top = min(c[1] - h / 2 for c, h in zip(centres, heights))
  bottom = max(c[1] + h / 2 for c, h in zip(centres, heights))
  factor = pair_span * image_size / max(right - left, bottom - top)
  body = []
  for index, cell in enumerate(cells):
    x = image_size / 2 + factor * (centres[index][0] - 0.5 * (left + right))
    y = image_size / 2 + factor * (centres[index][1] - 0.5 * (top + bottom))
    body += draw(cell, factor * scales[index], x, y, occlude=(bcc_silhouette and cell is bcc))
else:
  sys.exit("usage: crystal_schematic.py hcp|bcc+hcp <output.svg>")

# Assemble: white background, the arrow head marker, the cells
parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
         % (image_size, image_size, image_size, image_size),
         '<defs><marker id="head" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="4" '
         'markerHeight="4" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="%s"/>'
         '</marker></defs>' % colour_arrow,
         '<rect width="100%" height="100%" fill="#ffffff"/>']
parts += body
parts.append('</svg>')

with open(sys.argv[2], "w") as output:
  output.write("\n".join(parts) + "\n")
