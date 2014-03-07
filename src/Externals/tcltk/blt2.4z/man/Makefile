# Generated automatically from Makefile.in by configure.
# ------------------------------------------------------------------------
# 	Makefile for manual page files 
# ------------------------------------------------------------------------

prefix  	= /usr
mandir  	= ${prefix}/man
sectiondir 	= $(mandir)/mann
srcdir		= .
version		= 2.4

instdirs	= $(mandir) $(mandir)/mann $(mandir)/man3

MAN_N		= BLT.n barchart.n beep.n bgexec.n bitmap.n \
		  bltdebug.n busy.n container.n cutbuffer.n \
		  dragdrop.n eps.n graph.n hierbox.n  \
		  hiertable.n htext.n spline.n stripchart.n \
	 	  table.n tabset.n tile.n tree.n treeview.n vector.n \
		  watch.n winop.n

MAN_3 		= Blt_Tree.3 Blt_TreeGetNode.3 \
		  Blt_TreeCreate.3 Blt_TreeGetToken.3 \
		  Blt_TreeCreateNode.3 Blt_TreeName.3 \
		  Blt_TreeDeleteNode.3 Blt_TreeNodeId.3 \
		  Blt_TreeExists.3 Blt_TreeReleaseToken.3

MANPAGES	= $(MAN_N) $(MAN_3)

INSTALL  	= /usr/bin/install -c
INSTALL_DATA  	= ${INSTALL} -m 644
INSTALL_ROOT    = 
SHELL    	= /bin/sh
RM		= rm -rf

all: man.macros $(MANPAGES)

install: mkdirs install-mann install-man3

install-mann: $(MAN_N)
	for i in *.n ; do \
		$(INSTALL_DATA) $$i $(INSTALL_ROOT)$(mandir)/mann; \
	done

install-man3: $(MAN_3)
	for i in *.3 ; do \
		$(INSTALL_DATA) $$i $(INSTALL_ROOT)$(mandir)/man3; \
	done

mkdirs:
	@for i in $(instdirs) ; do \
          if test -d $(INSTALL_ROOT)$$i ; then \
            : ; \
          else \
            echo " mkdir $(INSTALL_ROOT)$$i" ; \
            mkdir $(INSTALL_ROOT)$$i ; \
          fi ; \
	done

.SUFFIXES: .n .mann .3 .man3

.man3.3: $(srcdir)/man.macros
	$(RM) $@
	sed -e "/man\.macros/r $(srcdir)/man.macros" -e '/man\.macros/d'  -e 's/BLT_VERSION/$(version)/' $< > $@

.mann.n: $(srcdir)/man.macros
	$(RM) $@
	sed -e "/man\.macros/r $(srcdir)/man.macros" -e '/man\.macros/d' -e 's/BLT_VERSION/$(version)/' $< > $@

clean:
	$(RM) *.3 *.n

distclean: clean
	$(RM) $(srcdir)/*.bak $(srcdir)/*\~ $(srcdir)/"#"* Makefile

