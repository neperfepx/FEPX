# FEPX website

Sources of https://fepx.info.

- `website-devel`: working branch, pushed to the private repository neperfepx/FEPX-dev.
  Each push builds a preview (with the manual of FEPX's `devel` branch), available as the
  `website-preview` artifact of the workflow run (see "Viewing the preview").
- `website`: published branch (squash-merged from `website-devel`).  Each push builds the
  website (with the manual of FEPX's `main` branch) and deploys it.

The manual (`doc/`) is not stored here: `make html` copies it from an FEPX branch (committed
files only) before building, together with the files it includes from `tests/`, `tutorials/`
and `COPYING` (in `doc/files/`).

## Building locally

```
pip install -r requirements.txt
make html          # copies doc/ from ../fepx at devel, then builds build/html/
```

Then open `build/html/index.html` (later builds only rebuild the pages that changed).  The FEPX repository and branch can be changed, as in
`make html FEPX_REPO=<path> DOC_REF=main`; `make doc` only copies the manual.

## Viewing the preview

neperfepx/FEPX-dev is private and cannot host the preview.  To download the latest preview
and serve it at http://localhost:8000:

```
./preview.sh [<port>]
```

The first run downloads the preview; later runs reuse it until a new preview is built.
Previews expire after 7 days; to rebuild one:

```
gh workflow run website.yaml --repo neperfepx/FEPX-dev --ref website-devel
```

## Old versions

Old FEPX versions are not stored here: they are the GitHub releases of FEPX (see `downloads.rst`).
