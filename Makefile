objects := $(wildcard R/*.R) $(wildcard src/*.[hc]pp) DESCRIPTION
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
	@Rscript -e "library(methods);" \
		-e "devtools::document();" \
		-e "Rcpp::compileAttributes()";
	@$(MAKE) updateTimestamp
	R CMD build .

$(checkLog): $(tar)
# R CMD check --as-cran $(tar)
	R CMD check $(tar)

# $(vignettes): $(rmd)
#	Rscript -e "rmarkdown::render('$(rmd)')"


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
	@rm -rf *~ */*~ *.Rhistroy src/{*.o,*.so} *.tar.gz *.Rcheck/ .\#*
