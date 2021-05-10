.PHONY: all
nb := $(patsubst %.R,%.ipynb,$(wildcard *.R)) \
	  $(patsubst %.py,%.ipynb,$(wildcard *.py))  \
	  $(patsubst %.sh,%.ipynb,$(wildcard *.sh))
rmd := $(patsubst %.ipynb,%.Rmd,$(nb))
ipynb.html := $(patsubst %.ipynb,%.html,$(nb))

all: $(rmd) $(nb)

.PHONY: clean
clean:
	@echo rm -rf *_cache/ *_files/ *.html

%.Rmd: %.R
	jupytext --to Rmd '$^'

%.Rmd: %.py
	jupytext --to Rmd '$^'

%.Rmd: %.sh
	jupytext --to Rmd '$^'

%.ipynb: %.R
	jupytext --to ipynb --output '$@' '$^'

%.ipynb: %.py
	jupytext --to ipynb --output '$@' '$^'

%.ipynb: %.sh
	jupytext --to ipynb --output '$@' '$^'
	
%.Rmd.html: %.Rmd
	Rscript -e 'rmarkdown::render("$^", output_format="html_document", output_file="$@")'

%.ipynb.html: %.ipynb
	Rscript -e 'rmarkdown::render("$^", output_format="html_document", output_file="$@")'
