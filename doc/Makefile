# Makefile of the FEPX documentation (Sphinx)
#
# make html    build the documentation in build/html (warnings are errors, as in the CI)
# make clean   remove build/

SPHINXBUILD ?= sphinx-build
SPHINXOPTS  ?= -W --keep-going
BUILDDIR    ?= build

.PHONY: help html clean

help:
	@echo "make html    build the documentation in $(BUILDDIR)/html"
	@echo "make clean   remove $(BUILDDIR)"

html:
	@$(SPHINXBUILD) -M html . "$(BUILDDIR)" $(SPHINXOPTS)
	@echo "The documentation is in $(BUILDDIR)/html."

clean:
	rm -rf "$(BUILDDIR)"
