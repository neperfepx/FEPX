#!/bin/bash
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Make the generated files of the manual's tutorial pages (step 9.11, ROADMAP D95): for every case
# of tutorials/tutorials.txt, run FEPX on the inputs the test suite runs, and write
# to doc/tutorials/generated/<case>/:
#   terminal.txt   the terminal output of the run, shortened (tutorials/tools/tutorial_pages.py excerpt)
#   tree.txt       the tree of the simulation directory
#   *.png          the mesh and field images (neper -V, POV-Ray; the color scale beside each
#                  field, ImageMagick) and the curves (tutorials/tools/tutorial_pages.py plot, Matplotlib,
#                  from the face forces, results/macro, and the volume averages of
#                  neper -S -resmesh stress,strain)
#
# Run by hand, from a clean git tree (the banner of the terminal output names the commit), after
# a Release build, at the end of a step that changes the terminal output, the inputs or the
# results of a tutorial. Needs Neper, POV-Ray, ImageMagick (convert) and Python 3 with NumPy
# and Matplotlib. About 30 min on 4 cores, plus about 1.5 h for the cruciform specimen (tier
# "manual"; left out unless asked for with --case cruciform or --all).
#
# Usage: tutorials/tools/make_tutorial_pages.sh [--build <dir>] [--work <dir>] [--case <case>]... [--all]
#                                     [--figures-only]
#   --build <dir>    the FEPX build (default build-openmpi-Release: Release, Open MPI)
#   --work <dir>     the run directories (default build-tutorial-pages; git-ignored)
#   --case <case>    only this case (repeatable); a restart case needs its prerequisite's run
#   --all            every case, the cruciform specimen included
#   --figures-only   do not run FEPX: make the files from the runs already in the work directory
#
# The pages include the files with literalinclude and figure directives; the docs test
# (tutorials/tools/tutorial_pages.py check) checks that every file a page names exists.
set -euo pipefail

# --- Paths and options -------------------------------------------------------------------------
# repo_dir: the source tree; list_file: the case list; inputs_root: the tutorials' inputs;
# mesh_dir: the meshes of the regression cases; out_root: the generated files.
repo_dir=$(cd "$(dirname "$0")/../.." && pwd)
list_file=${repo_dir}/tutorials/tutorials.txt
inputs_root=${repo_dir}/tutorials
mesh_dir=${repo_dir}/tests/regression/cases/meshes
out_root=${repo_dir}/doc/tutorials/generated
helper="python3 ${repo_dir}/tutorials/tools/tutorial_pages.py"
build_dir=${repo_dir}/build-openmpi-Release
work_root=${repo_dir}/build-tutorial-pages
selected=()
all_cases=0
figures_only=0
while [ $# -gt 0 ]; do
  case $1 in
    --build) build_dir=$(cd "$2" && pwd); shift 2 ;;
    --work) mkdir -p "$2"; work_root=$(cd "$2" && pwd); shift 2 ;;
    --case) selected+=("$2"); shift 2 ;;
    --all) all_cases=1; shift ;;
    --figures-only) figures_only=1; shift ;;
    *) echo "make_tutorial_pages.sh: unknown argument $1" >&2; exit 2 ;;
  esac
done
fepx=${build_dir}/fepx

# --- Checks of the tools -----------------------------------------------------------------------
for tool in neper povray convert python3 mpirun; do
  command -v ${tool} > /dev/null || { echo "make_tutorial_pages.sh: ${tool} not found" >&2; exit 1; }
done
python3 -c "import numpy, matplotlib" || { echo "make_tutorial_pages.sh: NumPy and Matplotlib needed" >&2; exit 1; }
[ ${figures_only} -eq 1 ] || [ -x "${fepx}" ] || { echo "make_tutorial_pages.sh: no ${fepx}" >&2; exit 1; }
if [ ${figures_only} -eq 0 ] && [ -n "$(git -C "${repo_dir}" status --porcelain --untracked-files=no)" ]; then
  echo "make_tutorial_pages.sh: warning: the git tree has changes (the banner will say -dirty)" >&2
fi

# --- The case list ------------------------------------------------------------------------------
# One associative entry per column, keyed by the case name; case_order keeps the list's order.
declare -A case_inputs case_mesh case_ori case_requires case_ranks case_tier
case_order=()
while read -r name inputs mesh ori requires ranks tier; do
  case_inputs[${name}]=${inputs}; case_mesh[${name}]=${mesh}; case_ori[${name}]=${ori}
  case_requires[${name}]=${requires}; case_ranks[${name}]=${ranks}; case_tier[${name}]=${tier}
  case_order+=("${name}")
done < <(grep -v -E '^[[:space:]]*(#|$)' "${list_file}")
if [ ${#selected[@]} -eq 0 ]; then
  for name in "${case_order[@]}"; do
    if [ ${all_cases} -eq 1 ] || [ "${case_tier[${name}]}" = "test" ]; then
      selected+=("${name}")
    fi
  done
fi

# --- Running a case ------------------------------------------------------------------------------
# run_case <case>: a fresh run directory work/<case> with the case's inputs (a restart case starts
# from a copy of its prerequisite's directory, where the page's reader runs it), the run, its
# terminal output in terminal.log (standard error, e.g. the MPI library's messages, in
# terminal.err).
run_case() {
  local name=$1
  local work=${work_root}/${name}
  local inputs=${inputs_root}/${case_inputs[${name}]}
  rm -rf "${work}"
  if [ "${case_requires[${name}]}" != "-" ]; then
    cp -r "${work_root}/${case_requires[${name}]}" "${work}"
    rm -f "${work}/terminal.log"
    # Neper's post-processing of the first run (results/mesh, results/elsets) is not FEPX's:
    # neper -S cannot extend it to the new steps, so it is made again after this run.
    rm -rf "${work}/simulation.sim/results/mesh" "${work}/simulation.sim/results/elsets"
  else
    mkdir -p "${work}"
    cp "${mesh_dir}/${case_mesh[${name}]}" "${work}/simulation.msh"
  fi
  cp "${inputs}/simulation.cfg" "${work}/simulation.cfg"
  [ -f "${inputs}/simulation.phase" ] && cp "${inputs}/simulation.phase" "${work}/"
  case ${case_ori[${name}]} in
    -) ;;
    inputs:*) cp "${inputs}/${case_ori[${name}]#inputs:}" "${work}/simulation.ori" ;;
    *) cp "${mesh_dir}/${case_ori[${name}]}" "${work}/simulation.ori" ;;
  esac
  echo "== ${name}: running on ${case_ranks[${name}]} ranks"
  (cd "${work}" && mpirun -np "${case_ranks[${name}]}" "${fepx}" > terminal.log 2> terminal.err) || {
    echo "make_tutorial_pages.sh: ${name} failed, see ${work}/terminal.log" >&2; exit 1; }
}

# post_process <case>: the volume averages of stress and strain (neper -S -resmesh), when printed.
post_process() {
  local work=${work_root}/$1
  if grep -q '^print stress$' "${work}/simulation.cfg" && grep -q '^print strain$' "${work}/simulation.cfg"; then
    neper -S "${work}/simulation.sim" -resmesh stress,strain > "${work}/neper_s.log" 2>&1
  fi
}

# --- Images --------------------------------------------------------------------------------------
# Common options of the 3D views: the deformed mesh, element edges, the domain edges, the view.
view_options=(-showelt1d all -dataelt1drad 0.004 -dataelt3dedgerad 0.0015
              -dataelt3dedgecol 32:32:32 -cameraangle 13.5 -imagesize 800:800)

# note <case> <command...>: append a command, as typed in the simulation's directory, to the
# case's generated/<case>/commands.txt (the page shows it): arguments with spaces or shell
# characters in double quotes, lines continued with "\" beyond 88 columns.
note() {
  local name=$1
  shift
  local line="$" argument word
  for argument in "$@"; do
    if [[ -z "${argument}" || "${argument}" == *[\ \(\)\[\]\*\?\;\&\|\<\>]* ]]; then
      word="\"${argument}\""
    else
      word=${argument}
    fi
    if [ $(( ${#line} + 1 + ${#word} )) -gt 88 ]; then
      echo "${line} \\" >> "${out_root}/${name}/commands.txt"
      line="   "
    fi
    line+=" ${word}"
  done
  echo "${line}" >> "${out_root}/${name}/commands.txt"
}

# in_work <case> <log> <command...>: run a command in the case's run directory, its output to
# <log> there, and note it.
in_work() {
  local name=$1 log=$2
  shift 2
  note "${name}" "$@"
  (cd "${work_root}/${name}" && "$@" > "${log}" 2>&1) || {
    echo "make_tutorial_pages.sh: ${name}: failed: $* (see ${work_root}/${name}/${log})" >&2; exit 1; }
}

# compose <case> <out_name> [<title>]: the image img_<out_name>.png and its color scale side by
# side, the title above the scale (written here: Neper's scale title can overlap wide tick
# labels), trimmed, with a white margin, as the 8-bit PNG <out_name>.png (flat-shaded elements:
# no visible loss, a tenth of the size), copied to generated/<case>/.
compose() {
  local name=$1 out_name=$2 title=${3:-}
  local image=img_${out_name}
  if [ -f "${work_root}/${name}/${image}-scale3d.png" ]; then
    in_work "${name}" convert_${out_name}.log convert "(" ${image}.png -trim ")" \
      "(" -size 40x1 xc:white ")" "(" -size 240x -background white -font DejaVu-Serif -pointsize 30 \
      caption:"${title}" -bordercolor white -border 0x12 "(" ${image}-scale3d.png -trim ")" \
      -gravity west -append ")" -background white -gravity center +append -bordercolor white \
      -border 20 -resize 600x600 -strip PNG8:${out_name}.png
  else
    in_work "${name}" convert_${out_name}.log convert ${image}.png -trim -bordercolor white \
      -border 20 -resize 600x600 -strip PNG8:${out_name}.png
  fi
  cp "${work_root}/${name}/${out_name}.png" "${out_root}/${name}/"
}

# mesh_image <case> <out_name> <neper options...>: the undeformed mesh, colored as the options
# say (e.g. -dataelsetcol id: by grain).
mesh_image() {
  local name=$1 out_name=$2
  shift 2
  in_work "${name}" neper_v_${out_name}.log neper -V simulation.msh "${view_options[@]}" "$@" \
    -print img_${out_name}
  compose "${name}" "${out_name}"
}

# field_image <case> <out_name> <step> <var> <component> <column> <title> <coo_factor> <scale>
# [options...]: an element field at a step on the deformed mesh (displacements times
# <coo_factor>); <scale> is "<low>:<high>", or "auto" for the 2nd to 98th percentiles of the
# field (tutorials/tools/tutorial_pages.py scale).
field_image() {
  local name=$1 out_name=$2 step=$3 var=$4 component=$5 column=$6 title=$7 factor=$8 scale=$9
  shift 9
  local work=${work_root}/${name}
  if [ "${scale}" = "auto" ]; then
    scale=$(${helper} scale "${work}/simulation.sim/results/elts/${var}/${var}.step${step}" "${column}")
  fi
  in_work "${name}" neper_v_${out_name}.log neper -V simulation.sim -step "${step}" \
    -datanodecoo coo -datanodecoofact "${factor}" "${view_options[@]}" \
    -dataeltcol "${var}${component}" -dataeltscale "${scale}" \
    "$@" -print img_${out_name}
  compose "${name}" "${out_name}" "${title}"
}

# ori_image <case> <out_name> <step>: the mesh at a step colored by element orientation.
ori_image() {
  local name=$1 out_name=$2 step=$3
  in_work "${name}" neper_v_${out_name}.log neper -V simulation.sim -step "${step}" \
    -datanodecoo coo "${view_options[@]}" -dataeltcol ori -print img_${out_name}
  compose "${name}" "${out_name}"
}

# plot <case> <plot case> <simulation directories...>: the curves (tutorials/tools/tutorial_pages.py plot),
# the directories relative to the case's run directory (e.g. simulation.sim,
# ../simple_rve/simulation.sim).
plot() {
  local name=$1 plot_case=$2
  shift 2
  note "${name}" python3 tutorial_pages.py plot "${plot_case}" . "$@"
  (cd "${work_root}/${name}" && ${helper} plot "${plot_case}" "${out_root}/${name}" "$@")
}

# --- The files of each case ----------------------------------------------------------------------
# make_files <case>: the terminal output, the tree, and the case's images and curves.
make_files() {
  local name=$1
  local work=${work_root}/${name}
  local sim=${work}/simulation.sim
  local out=${out_root}/${name}
  rm -rf "${out}"
  mkdir -p "${out}"
  ${helper} excerpt "${work}/terminal.log" "${out}/terminal.txt"
  ${helper} tree "${sim}" "${out}/tree.txt"
  echo "== ${name}: images and curves"
  # The volume averages (computed before, by post_process, for the curves of every case).
  if [ -d "${sim}/results/mesh" ]; then
    note "${name}" neper -S simulation.sim -resmesh stress,strain
  fi
  case ${name} in
    simple_simulation|simple_rve_grip)
      # One color scale for the simple simulation's steps and the grip variant: those of the
      # simple simulation at its last step.
      local simple_elts=${work_root}/simple_simulation/simulation.sim/results/elts
      local strain_scale stress_scale
      strain_scale=$(${helper} scale "${simple_elts}/strain/strain.step2" 3)
      stress_scale=$(${helper} scale "${simple_elts}/stress/stress.step2" 3)
      if [ ${name} = simple_simulation ]; then
        mesh_image ${name} mesh -dataelsetcol id
        for step in 1 2; do
          field_image ${name} strain33_step${step} ${step} strain 33 3 "Strain 33" 10 "${strain_scale}"
          field_image ${name} stress33_step${step} ${step} stress 33 3 "Stress 33 [MPa]" 10 "${stress_scale}"
        done
        plot ${name} simple_simulation simulation.sim
        # The post-processing tutorial's script, on this run: its figures and printout.
        local postprocess_out=${out_root}/postprocessing
        rm -rf "${postprocess_out}"
        mkdir -p "${postprocess_out}"
        cp "${inputs_root}/postprocessing/postprocess.py" "${work}/"
        (cd "${work}" && python3 postprocess.py simulation.sim > postprocess.txt)
        cp "${work}/stress_time.png" "${work}/grains.png" "${work}/postprocess.txt" \
          "${postprocess_out}/"
      else
        field_image ${name} strain33 1 strain 33 3 "Strain 33" 10 "${strain_scale}"
        field_image ${name} stress33 1 stress 33 3 "Stress 33 [MPa]" 10 "${stress_scale}"
      fi ;;
    simple_rve)
      plot ${name} simple_rve ../simple_simulation/simulation.sim simulation.sim \
        ../simple_rve_grip/simulation.sim ;;
    periodic_rve)
      mesh_image ${name} mesh -dataelsetcol id
      field_image ${name} strain11 1 strain 11 1 "Strain 11" 50 auto
      plot ${name} periodic_rve simulation.sim ;;
    multiaxial_csr)
      mesh_image ${name} mesh -dataelsetcol id
      field_image ${name} work_pl 2 work_pl "" 1 "Plastic work [mJ]" 10 auto
      field_image ${name} strain_pl_eq 2 strain_pl_eq "" 1 "Equivalent plastic strain" 10 auto
      plot ${name} multiaxial_csr simulation.sim ;;
    multiaxial_clr)
      mesh_image ${name} phases -dataelsetcol group
      field_image ${name} crss1 4 crss 1 1 "CRSS, system 1 [MPa]" 10 auto
      plot ${name} multiaxial_clr simulation.sim ;;
    restart_cycle1)
      mesh_image ${name} mesh -dataelsetcol id ;;
    restart_cycle2)
      field_image ${name} slip7 4 slip 7 7 "Slip, system 7" 10 auto
      field_image ${name} slip8 4 slip 8 8 "Slip, system 8" 10 auto
      plot ${name} restart simulation.sim ;;
    overwriting)
      ori_image ${name} ori 0
      field_image ${name} defrate_eq 1 defrate_eq "" 1 "Equivalent deformation rate [1/s]" 50 auto
      field_image ${name} work 2 work "" 1 "Work [mJ]" 50 auto
      plot ${name} overwriting simulation.sim ;;
    mixed_targets)
      field_image ${name} stress33_step2 2 stress 33 3 "Stress 33 [MPa]" 10 auto
      plot ${name} mixed_targets simulation.sim ;;
    hcp_basal|hcp_prismatic) ;;
    hcp_pyramidal)
      plot ${name} hcp ../hcp_basal/simulation.sim ../hcp_prismatic/simulation.sim \
        simulation.sim ;;
    cruciform)
      # The specimen seen from above (orthographic), as a thin plate.
      local camera=(-cameraprojection orthographic -cameralookat 0.6:0.6:0.05 -cameraangle 18)
      mesh_image ${name} mesh -dataelsetcol id "${camera[@]}"
      field_image ${name} stress11 1 stress 11 1 "Stress 11 [MPa]" 3 auto "${camera[@]}"
      field_image ${name} stress22 1 stress 22 2 "Stress 22 [MPa]" 3 auto "${camera[@]}"
      field_image ${name} stress12 1 stress 12 6 "Stress 12 [MPa]" 3 auto "${camera[@]}"
      plot ${name} cruciform simulation.sim ;;
    *) echo "make_tutorial_pages.sh: no files defined for ${name}" >&2; exit 1 ;;
  esac
}

# --- Main -----------------------------------------------------------------------------------------
mkdir -p "${work_root}" "${out_root}"
for name in "${selected[@]}"; do
  [ -n "${case_inputs[${name}]+set}" ] || { echo "make_tutorial_pages.sh: unknown case ${name}" >&2; exit 1; }
  if [ ${figures_only} -eq 0 ]; then
    run_case "${name}"
  fi
  post_process "${name}"
done
for name in "${selected[@]}"; do
  make_files "${name}"
done
echo "== done: ${out_root}"
