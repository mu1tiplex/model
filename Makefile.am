SUBDIRS = model tutorial generator

EXTRA_DIST = Makefile.special Doxyfile.in docs_footer.html docs_header.html LICENSE INTRODUCTION.lg docs_mainpageend docs_mainpagestart \
				INTRODUCTION.pdf INTRODUCTION.html README ChangeLog tutorial/ChangeLog model/ChangeLog generator/ChangeLog

CLEANFILES = *.*~ *.dvi *.tex *.log *.aux Doxyfile
MAINTAINERCLEANFILES = doc/latex/* doc/html/* doc/refman.pdf

# These are the special commands for the docs and stuff
include Makefile.special

# regenerate ChangeLog(s) from CVS tree
ChangeLogs:
	cvs2cl.pl -d -b -t --window 86400

# run this to test the lot
test: tutorial/singlemode
	make -C tutorial test

#make sure the docs are included
dist-hook: docs 
	cp -p -r $(srcdir)/doc $(distdir)

AUTOMAKE_OPTIONS = foreign 