# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
PYTHON        = python3
SPHINXOPTS    =
SPHINXBUILD   = ${PYTHON} -msphinx
SPHINXPROJ    = Orbits
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@export PYTHONPATH=$(shell pwd)/.. ; $(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
