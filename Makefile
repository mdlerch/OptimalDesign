all: clean install

PKG_VERS = $(shell grep -i ^version DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME = $(shell grep -i ^package DESCRIPTION | cut -d : -d \  -f 2)

.PHONY: all install build check clean vignettes

install: $(PKG_NAME)_$(PKG_VERS).tar.gz
	R CMD INSTALL $<

build: $(PKG_NAME)_$(PKG_VERS).tar.gz

check: $(PKG_NAME)_$(PKG_VERS).tar.gz
	R CMD check --as-cran $<

clean:
	-rm $(PKG_NAME)_*.tar.gz
	-rm -rf $(PKG_NAME).Rcheck

$(PKG_NAME)_$(PKG_VERS).tar.gz: DESCRIPTION
	Rscript -e "Rcpp::compileAttributes()"
	R CMD build ../OptimalDesign

vignettes: inst/doc/optimalDesign.pdf

inst/doc/optimalDesign.pdf: vignettes/optimalDesign.Rnw
	Rscript -e "knitr::knit2pdf(\"vignettes/optimalDesign.Rnw\")"
	mv vignettes/optimalDesign.pdf inst/doc/
