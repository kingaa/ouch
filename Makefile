REXE = R --vanilla
RCMD = $(REXE) CMD
RCMD_ALT = R --no-save --no-restore CMD
RSCRIPT = Rscript --vanilla
REPODIR = ../www

PDFLATEX = pdflatex
BIBTEX = bibtex
MAKEIDX = makeindex

RM = rm -f
CP = cp
TOUCH = touch
INSTALL = install

PKG = $(shell perl -ne 'print $$1 if /Package:\s+((\w+[-\.]?)+)/;' DESCRIPTION)
VERSION = $(shell perl -ne 'print $$1 if /Version:\s+((\d+[-\.]?)+)/;' DESCRIPTION)
PKGVERS = $(PKG)_$(VERSION)
SOURCE=$(shell ls R/*R src/*.c src/*.h data/*)
CSOURCE=$(shell ls src/*.c)
TESTS=$(shell ls tests/*R)
REVDEPS=surface mvSLOUCH pmc

default:
	@echo $(PKGVERS)

.PHONY: clean win wind tests check

dist manual vignettes: export R_QPDF=qpdf
roxy dist manual vignettes: export R_HOME=$(shell $(REXE) RHOME)
check xcheck xxcheck: export FULL_TESTS=yes
dist revdeps session tests check xcheck xxcheck: export R_KEEP_PKG_SOURCE=yes
revdeps xcheck tests: export R_PROFILE_USER=$(CURDIR)/.Rprofile
revdeps session xxcheck htmldocs vignettes data tests manual: export R_LIBS=$(CURDIR)/library
session: export R_DEFAULT_PACKAGES=datasets,utils,grDevices,graphics,stats,methods,tidyverse,subplex,ouch

htmldocs: inst/doc/*.html

htmlhelp: install
	rsync -avz library/ouch/html/ www/manual
	(cd www/manual;	(cat links.ed && echo w ) | ed - 00Index.html)

vignettes: manual install
	$(MAKE)	-C www/vignettes

news: www/NEWS.html

NEWS: inst/NEWS

inst/NEWS: inst/NEWS.Rd
	$(RCMD) Rdconv -t txt $^ -o $@

www/NEWS.html: inst/NEWS.Rd
	$(RCMD) Rdconv -t html inst/NEWS.Rd -o www/NEWS.html

session: install
	exec $(REXE)

revdeps: install
	mkdir -p library check
	$(REXE) -e "pkgs <- strsplit('$(REVDEPS)',' ')[[1]]; download.packages(pkgs,destdir='library',repos='https://mirrors.nics.utk.edu/cran/')"
	$(RCMD) check --library=library -o check library/*.tar.gz

roxy: $(SOURCE)
	$(REXE) -e "pkgbuild::compile_dll(); devtools::document(roclets=c('rd','collate','namespace'))"

dist: NEWS $(PKGVERS).tar.gz

$(PKGVERS).tar.gz: $(SOURCE) $(TESTS)
	$(RCMD) build --force --no-manual --resave-data --compact-vignettes=both --md5 .

binary: dist
	mkdir -p plib
	$(RCMD) INSTALL --build --library=plib --preclean --clean $(PKGVERS).tar.gz
	rm -rf plib

publish: dist manual news
	$(RSCRIPT) -e 'drat::insertPackage("$(PKGVERS).tar.gz",repodir="$(REPODIR)",action="prune")'
	-$(RSCRIPT) -e 'drat::insertPackage("$(PKGVERS).tgz",repodir="$(REPODIR)",action="prune")'
	-$(RSCRIPT) -e 'drat::insertPackage("$(PKGVERS).zip",repodir="$(REPODIR)",action="prune")'
	$(CP) $(PKG).pdf ../www/manuals

win: dist
	curl -T $(PKGVERS).tar.gz ftp://win-builder.r-project.org/R-release/

wind: dist
	curl -T $(PKGVERS).tar.gz ftp://win-builder.r-project.org/R-devel/

check: dist
	mkdir -p check
	$(RCMD) check --no-stop-on-test-error --library=check -o check $(PKGVERS).tar.gz

qcheck: dist
	mkdir -p check
	$(RCMD) check --library=check -o check --no-vignettes --no-tests $(PKGVERS).tar.gz

qqcheck: dist
	mkdir -p check
	$(RCMD) check --library=check -o check --no-codoc --no-examples --no-vignettes --no-manual --no-tests $(PKGVERS).tar.gz

xcheck: dist
	mkdir -p check library
	$(RCMD_ALT) check --no-stop-on-test-error --as-cran --library=library -o check $(PKGVERS).tar.gz

xxcheck: install xcheck
	mkdir -p check
	$(REXE) -d "valgrind --tool=memcheck --track-origins=yes --leak-check=full" < check/$(PKG).Rcheck/$(PKG)-Ex.R 2>&1 | tee $(PKG)-Ex.Rout

ycheck: dist
	mkdir -p check
	$(RCMD_ALT) check --run-dontrun --run-donttest --as-cran --library=library -o check $(PKGVERS).tar.gz

manual: install $(PKG).pdf

$(PKG).pdf: $(SOURCE)
	$(RCMD) Rd2pdf --no-preview --pdf --force -o $(PKG).pdf .
	$(RSCRIPT) -e "tools::compactPDF(\"$(PKG).pdf\")";

tests: install $(TESTS)
	export R_LIBS
	$(MAKE) -C tests

install: library/$(PKG)

library/$(PKG): dist
	mkdir -p library
	$(RCMD) INSTALL --html --library=library $(PKGVERS).tar.gz

remove:
	if [ -d library ]; then \
		$(RCMD) REMOVE --library=library $(PKG); \
		rmdir library; \
	fi

fresh: clean remove

inst/doc/*.html: install 

%.tex: %.Rnw
	$(RSCRIPT) -e "library(knitr); knit(\"$*.Rnw\")"

%.R: %.Rnw
	$(RSCRIPT) -e "library(knitr); purl(\"$*.Rnw\")"

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.bbl: %.tex
	-$(PDFLATEX) $*
	$(BIBTEX) $*

%.idx: %.tex
	-$(PDFLATEX) $*

%.ind: %.idx
	$(MAKEIDX) $*

%.html: %.Rmd
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\")"

%.html: %.md
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.md\")"

%.R: %.Rmd
	Rscript --vanilla -e "knitr::purl(\"$*.Rmd\",output=\"$*.R\",documentation=2)"

clean:
	$(RM) -r check
	$(RM) src/*.o src/*.so src/symbols.rds www/vignettes/Rplots.*
	$(RM) -r inst/doc/figure inst/doc/cache
	$(RM) -r *-Ex.Rout *-Ex.timings *-Ex.pdf
	$(RM) *.tar.gz $(PKGVERS).zip $(PKGVERS).tgz $(PKG).pdf
