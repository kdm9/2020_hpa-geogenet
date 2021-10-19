outhtmls := $(patsubst %.R,%.html,$(wildcard *.R)) \
	   $(patsubst %.py,%.html,$(wildcard *.py))  \
	   $(patsubst %.sh,%.html,$(wildcard *.sh))

.PRECIOUS: %.Rmd

html: $(outhtmls)


.PHONY: clean
clean:
	@echo rm -rf *_cache/ *_files/ *.html

ipynb := $(patsubst %.Rmd,%.ipynb,$(rmd))
ipynb_html := $(patsubst %.ipynb,%.html,$(nb))
nb_html: $(ipynb_html)

.PRECIOUS: %.ipynb

%.html: %.py
	jupytext --to ipynb --execute --output '$@' '$^'

%.html: %.sh
	jupytext --to ipynb --execute --output '$@' '$^'
	
%.html: %.R
	Rscript -e 'rmarkdown::render("$^", output_file="$@")'
