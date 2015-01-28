all: install

PKG_VERS = $(shell grep -i ^version DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME = $(shell grep -i ^package DESCRIPTION | cut -d : -d \  -f 2)

.PHONY: all install build check clean

install: $(PKG_NAME)_$(PKG_VERS).tar.gz
	R CMD INSTALL $<

build: $(PKG_NAME)_$(PKG_VERS).tar.gz

check: $(PKG_NAME)_$(PKG_VERS).tar.gz
	R CMD check $<

clean:
	-rm $(PKG_NAME)_$(PKG_VERS).tar.gz
	-rm -rf $(PKG_NAME).Rcheck

$(PKG_NAME)_$(PKG_VERS).tar.gz:
	Rscript -e "Rcpp::compileAttributes()"
	R CMD build ../OptimalDesign

