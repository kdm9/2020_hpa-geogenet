.PHONY: all
all: 00_metadata.html 01_hpa-exploratory.html 02_construct.html 03_ld.html 04_lostruct.html 05_lfmm.html

.PHONY: clean
clean:
	@echo rm -rf *_cache/ *_files/ *.html

%.Rmd: %.R
	jupytext --to Rmd '$^'
	
%.html: %.Rmd
	Rscript -e 'rmarkdown::render("$^", output_format="html_document")'
