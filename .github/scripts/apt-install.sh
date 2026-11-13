#!/usr/bin/env bash
# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Install Ubuntu packages on a GitHub Actions runner (used by .github/workflows/ci.yml).
#
# Usage: apt-install.sh <archive-dir> <package>...
#
#   archive-dir  directory of downloaded .deb files; restored from and saved to the actions/cache
#                by the workflow, so that a run with a cache hit downloads (almost) nothing
#   package      packages to install (--no-install-recommends)
#
# Mirrors: the runner's default Azure mirror (azure.archive.ubuntu.com, plain HTTP) has stalled
# repeatedly (installs of 50 MB taking over 10 minutes). apt on the Ubuntu 24.04 runners reads its
# mirror list from /etc/apt/apt-mirrors.txt (a "mirror+file:" source), so the Azure host is
# replaced there as well as in the sources files.

set -euo pipefail

archive_dir=$1
shift

# Use the main Ubuntu archive over HTTPS instead of the Azure HTTP mirror, in every place apt
# reads mirrors from.
for f in /etc/apt/apt-mirrors.txt /etc/apt/sources.list /etc/apt/sources.list.d/*.sources \
         /etc/apt/sources.list.d/*.list; do
  if [ -f "$f" ]; then
    sudo sed -i 's|http://azure.archive.ubuntu.com|https://archive.ubuntu.com|g' "$f"
  fi
done

# Keep downloaded packages in the cached directory (apt needs its partial/ subdirectory).
mkdir -p "$archive_dir/partial"
echo "Cached packages: $(find "$archive_dir" -maxdepth 1 -name '*.deb' | wc -l)"

# Retry and time out instead of hanging on a stalled connection.
apt_opts=(-o Acquire::Retries=3 -o Acquire::http::Timeout=30 -o Acquire::https::Timeout=30
          -o "Dir::Cache::Archives=$archive_dir" -o APT::Keep-Downloaded-Packages=true)

# Update the package lists and install.
sudo apt-get "${apt_opts[@]}" update
sudo apt-get "${apt_opts[@]}" install -y --no-install-recommends "$@"

# Make the archive readable by the cache action (apt wrote it as root); drop apt's lock and
# partial downloads, which must not be cached.
sudo rm -rf "$archive_dir/partial" "$archive_dir/lock"
sudo chown -R "$(id -u):$(id -g)" "$archive_dir"
echo "Packages now in the archive: $(find "$archive_dir" -maxdepth 1 -name '*.deb' | wc -l)"
