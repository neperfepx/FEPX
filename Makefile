# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = ReadtheDocsSphinxTheme
SOURCEDIR     = .
BUILDDIR      = build

# FEPX repository and branch (or tag, commit) from which the manual is copied.
# Only committed files are copied.  The workflow overrides them (main for the website).
FEPX_REPO     ?= ../fepx
DOC_REF       ?= devel

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile html doc clean

# Copy the manual (doc/) of $(FEPX_REPO) at $(DOC_REF) into doc/, self-contained:
#   1. the manual's files, without its conf.py (the site's applies), its Makefile, its templates and its manual page (doc/man/);
#   2. the files it takes from tests/ and tutorials/, and COPYING (literalinclude, download), in
#      doc/files/ at the same paths, with the references rewritten to them;
#   3. its title page becomes the site's "Documentation" page, without the tutorials (they have
#      their own entry in the site's menu);
#   4. its style sheet is copied to _static/custom.css.
doc:
	@git -C "$(FEPX_REPO)" rev-parse --verify --quiet "$(DOC_REF)^{commit}" > /dev/null 2>&1 \
	  || { echo "No commit $(DOC_REF) in FEPX repository $(FEPX_REPO) (set FEPX_REPO=<path>)"; exit 1; }
	@set -e; \
	rm -rf doc $(BUILDDIR)/doc-src; mkdir -p $(BUILDDIR)/doc-src; \
	git -C "$(FEPX_REPO)" archive "$(DOC_REF)" -- ":(top)doc" | tar -x -C $(BUILDDIR)/doc-src; \
	mv $(BUILDDIR)/doc-src/doc doc; rm -rf $(BUILDDIR)/doc-src doc/conf.py doc/Makefile doc/_templates doc/man; \
	grep -rhoE '(\.\./)+((tests|tutorials)/[^ >`]+|COPYING|VERSIONS)' doc | sed -E 's#^(\.\./)+##' | sort -u | \
	  while read -r f; do \
	    mkdir -p "doc/files/$$(dirname "$$f")"; \
	    git -C "$(FEPX_REPO)" show "$(DOC_REF):$$f" > "doc/files/$$f"; \
	  done; \
	grep -rlE '\.\./(tests/|tutorials/|COPYING|VERSIONS)' doc | xargs -r sed -i -E 's#\.\./(tests/|tutorials/|COPYING|VERSIONS)#files/\1#g'; \
	sed -i -e 's/^\.\. _fepx:$$/.. _documentation:/' \
	  -e 's/^FEPX: Finite Element Polycrystal Plasticity$$/Documentation/' \
	  -e '/^   tutorials\.rst$$/d' doc/index.rst; \
	sed -i '/^Documentation$$/{n;s/^+*$$/+++++++++++++/}' doc/index.rst; \
	if [ -f doc/_static/custom.css ]; then cp doc/_static/custom.css _static/custom.css; fi
	@echo "Copied the manual from $(DOC_REF) ($$(git -C "$(FEPX_REPO)" rev-parse --short "$(DOC_REF)"))"

html: doc
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	cp -rp imgs _static/* $(BUILDDIR)/html
	@echo "The website is in $(BUILDDIR)/html."

clean:
	rm -rf $(BUILDDIR)
