.PHONY: all
rmd := $(patsubst %.R,%.Rmd,$(wildcard *.R)) \
	   $(patsubst %.py,%.Rmd,$(wildcard *.py))  \
	   $(patsubst %.sh,%.Rmd,$(wildcard *.sh))

.PRECIOUS: %.Rmd

outmd := $(patsubst %.Rmd,%.out.md,$(rmd))
outhtml := $(patsubst %.Rmd,%.html,$(rmd))


#all: $(rmd)
#md: $(outmd)
#html: $(outhtml)
#
#.PHONY: clean
#clean:
#	@echo rm -rf *_cache/ *_files/ *.html
#
#%.Rmd: %.R
#	jupytext --to Rmd '$^'
#
#%.Rmd: %.py
#	jupytext --to Rmd '$^'
#
#%.Rmd: %.sh
#	jupytext --to Rmd '$^'
#
%.md %.html: %.Rmd
	Rscript -e 'rmarkdown::render("$^", rmarkdown::html_document(keep_md=T), output_file="$@")'

ipynb := $(patsubst %.Rmd,%.ipynb,$(rmd))
ipynb_html := $(patsubst %.ipynb,%.html,$(nb))
nb_html: $(ipynb_html)

.PRECIOUS: %.ipynb

%.ipynb: %.R
	jupytext --to ipynb --execute --output '$@' '$^'

%.ipynb: %.py
	jupytext --to ipynb --execute --output '$@' '$^'

%.ipynb: %.sh
	jupytext --to ipynb --execute --output '$@' '$^'
	
%.html: %.ipynb
	jupyter nbconvert --output $@ --to html --no-prompt $^ 
