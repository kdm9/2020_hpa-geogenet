.PHONY: all
all: 2020-12-14_ggexplore.html

.PHONY: clean
clean:
	@echo rm -rf *_cache/ *_files/ *.html

%.html: %.Rmd
	Rscript -e 'rmarkdown::render("$^", output_format="html_document")'
