.PHONY: all
all: 01_hpa-exploratory.html 02_construct.html

.PHONY: clean
clean:
	@echo rm -rf *_cache/ *_files/ *.html

%.html: %.Rmd
	Rscript -e 'rmarkdown::render("$^", output_format="html_document")'
