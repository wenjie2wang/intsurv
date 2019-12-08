objects := DESCRIPTION $(wildcard R/*.R) \
	$(wildcard src/*.cpp) $(wildcard inst/include/*.h) \
	$(wildcard inst/include/intsurv/*.h)
version := $(shell grep "Version" DESCRIPTION | awk '{print $$NF}')
pkg := $(shell grep "Package" DESCRIPTION | awk '{print $$NF}')
tar := $(pkg)_$(version).tar.gz
checkLog := $(pkg).Rcheck/00check.log
# tests := $(wildcard tests/testthat/*.R)
# rmd := vignettes/$(pkg)-intro.Rmd
# vignettes := vignettes/$(pkg)-intro.html

.PHONY: check
check: $(checkLog)

.PHONY: build
build: $(tar)

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)

# .PHONY: preview
# preview: $(vignettes)


$(tar): $(objects)
	@rm -rf src/RcppExports.cpp R/RcppExports.R
	@Rscript -e "library(methods);" \
	-e "Rcpp::compileAttributes()" \
	-e "devtools::document();";
	@$(MAKE) updateTimestamp
	R CMD build .

$(checkLog): $(tar)
	R CMD check $(tar)

.PHONY: check-as-cran
check-as-cran: $(tar)
	R CMD check --as-cran $(tar)

# $(vignettes): $(rmd)
#	Rscript -e "rmarkdown::render('$(rmd)')"

.PHONY: pkgdown
pkgdown:
	@Rscript -e "library(methods); pkgdown::build_site();"


## update copyright year in HEADER, R script and date in DESCRIPTION
.PHONY: updateTimestamp
updateTimestamp:
	@bash misc/update_timestamp.sh

## make tags
.PHONY: TAGS
TAGS:
	Rscript -e "utils::rtags(path = 'R', ofile = 'TAGS')"
	gtags

## do some cleaning
.PHONY: clean
clean:
	@$(RM) -rf *~ */*~ *.Rhistroy src/{*.o,*.so} *.tar.gz *.Rcheck/ .\#*
