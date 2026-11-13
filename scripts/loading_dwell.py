#!/usr/bin/env python3
"""Cartoon of a triaxial load-dwell-unload history: the landing page's "Loading" card.

Draws the three normal stresses sigma_11, sigma_22, sigma_33 against time over two cycles, each a
proportional ramp to a peak, a dwell (the stresses held), and an unloading to zero.  The history
is one FEPX can run, triaxial load control with 33 the control component (manual, loading.rst,
"Load and stress rates, dwells" and "Triaxial load control"); with stresses in MPa and times in
s, its loading block is

  boundary_conditions triaxial
  control_component 33
  number_of_steps 6
  #                 comp  step 1  step 2  step 3  step 4  step 5  step 6
  stress_rate       33    100.0   0.0    -100.0   100.0   0.0    -100.0
  target_stress     33    500.0   -       0.0     500.0   -       0.0
  target_duration         -       5.0     -       -       5.0     -
  target_stress     22    250.0   250.0   0.0     250.0   250.0   0.0
  target_stress     11   -150.0  -150.0   0.0    -150.0  -150.0   0.0
  number_of_incrs         10      -       10      10      -       10
  dtime                   -       0.5     -       -       0.5     -

(the values are illustrative; the figure has no numbers on its axes).  Colours: three points of
viridis, the colour map of the other cards.

Usage, from the website's top directory:
  python3 scripts/loading_dwell.py imgs/gallery/load-dwell.png
"""

import io
import sys

import matplotlib
from PIL import Image, ImageOps

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (after the backend is chosen)

# Peak stresses of the three components [MPa] (the figure shows no numbers)
peak_stress = {"33": 500.0, "22": 250.0, "11": -150.0}

# One cycle, as (time [s], fraction of the peak) corners: ramp at 100 MPa/s on 33 (5 s), dwell of
# 5 s, unloading at -100 MPa/s (5 s); the second cycle follows the first directly
ramp_time = peak_stress["33"] / 100.0
dwell_time = 5.0
cycle = [(0.0, 0.0), (ramp_time, 1.0), (ramp_time + dwell_time, 1.0),
         (2.0 * ramp_time + dwell_time, 0.0)]
cycle_length = cycle[-1][0]
corners = cycle + [(cycle_length + t, f) for t, f in cycle[1:]]
times = [t for t, _ in corners]
fractions = [f for _, f in corners]

# Colours: viridis at 0.15, 0.5, 0.85 (dark purple, teal, yellow-green) for 11, 22, 33
viridis = matplotlib.colormaps["viridis"]
colours = {"11": viridis(0.15), "22": viridis(0.5), "33": viridis(0.82)}

# Figure: square, as the other cards' images (800 x 800 px), large type for the small card
figure, axes = plt.subplots(figsize=(8, 8), dpi=100)
figure.subplots_adjust(left=0.16, right=0.95, bottom=0.2, top=0.83)

# The three stress histories (33 last, so it is drawn on top)
for component in ("11", "22", "33"):
  axes.plot(times, [f * peak_stress[component] for f in fractions], color=colours[component],
            linewidth=7, solid_joinstyle="round", solid_capstyle="round",
            label=r"$\sigma_{%s}$" % component)

# Zero line, axes as arrows (no ticks, no frame)
axes.axhline(0.0, color="#000000", linewidth=2.5, zorder=0)
for side in ("top", "right", "bottom"):
  axes.spines[side].set_visible(False)
axes.spines["left"].set_linewidth(2.5)
axes.set_xticks([])
axes.set_yticks([])
axes.set_xlim(0.0, times[-1] * 1.05)
axes.set_ylim(1.15 * peak_stress["11"], 1.15 * peak_stress["33"])
axes.set_xlabel("time", fontsize=34, labelpad=16)
axes.set_ylabel("stress", fontsize=34, labelpad=16)
axes.xaxis.set_label_coords(0.5, -0.06)

# Legend above the plot, one row
axes.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False,
            fontsize=36, handlelength=1.2, columnspacing=1.2)

# Save: the drawing (lines, labels, legend) cropped to its extent at twice the resolution, scaled
# down so that its longer side is 600 px (about the size of the other cards' drawings, the Neper
# cubes 550-610 px and the hcp prism 450-490 px), and placed on a white 800 x 800 px image with
# the plot box (from the stress axis to the end of the zero line) centred horizontally and the
# whole drawing centred vertically
save_dpi = 200
buffer = io.BytesIO()
figure.savefig(buffer, format="png", dpi=save_dpi, facecolor="#ffffff")
buffer.seek(0)
drawing = Image.open(buffer).convert("RGB")
# Plot box, horizontally, in pixels of the saved image (the window extent is at figure.dpi)
plot_box = axes.get_window_extent()
plot_left = plot_box.x0 * save_dpi / figure.dpi
plot_right = plot_box.x1 * save_dpi / figure.dpi
# Crop to the non-white pixels (matplotlib's own tight box keeps some white around the drawing)
crop_box = ImageOps.invert(drawing.convert("L")).point(lambda v: 255 if v > 8 else 0).getbbox()
drawing = drawing.crop(crop_box)
drawing_size = 600
scale = drawing_size / max(drawing.size)
drawing = drawing.resize((round(scale * drawing.width), round(scale * drawing.height)),
                         Image.LANCZOS)
# Plot box centre in the scaled, cropped drawing; offsets that put it at the image's centre
plot_centre = scale * (0.5 * (plot_left + plot_right) - crop_box[0])
card = Image.new("RGB", (800, 800), "#ffffff")
card.paste(drawing, (round(400 - plot_centre), (800 - drawing.height) // 2))
card.save(sys.argv[1])
