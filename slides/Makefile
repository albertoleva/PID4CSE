LaTeX=PID4CSE-slides

all: compile clean

compile:
	pdflatex $(LaTeX).tex
	pdflatex $(LaTeX).tex

clean: 
	find . -name \*.bbl | xargs rm -f 
	find . -name \*.log | xargs rm -f 
	find . -name \*.lot | xargs rm -f 
	find . -name \*.toc | xargs rm -f 
	find . -name \*.out | xargs rm -f 
	find . -name \*.blg | xargs rm -f 
	find . -name \*.aux | xargs rm -f 
	find . -name \*.ps  | xargs rm -f 
	find . -name \*.dvi | xargs rm -f 
	find . -name \*.snm | xargs rm -f 
	find . -name \*.vrb | xargs rm -f 
	find . -name \*.nav | xargs rm -f
	rm -f *.tex~
	rm -f sec/*.tex~
	rm -rf sec/auto
	rm -rf auto

view:
	open *.pdf
