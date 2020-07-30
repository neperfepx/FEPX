# This file is part of the FEPX software package.
# Copyright (C) 1996-2020 ACME Lab
# See the COPYING file in the top-level directory.

pdf:
	texi2pdf fepx.texi
	cp fepx.pdf ..

info:
	makeinfo fepx.texi
	cp fepx.info ..

html:
	makeinfo --html fepx.texi

clean:
	rm -f *.cp* *.fn* *.ky* *.pg* *.tp* *.vr* *.mv*\
          *.log *.toc *.aux *.dvi *~ *.pdf *.info
install:
	su -c "mkdir -p /usr/local /usr/local/info;\
	cp -f fepx.info /usr/local/info/"
