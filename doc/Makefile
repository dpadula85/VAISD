SHELL := /bin/bash

.PHONY: all view
all:
	pdflatex VAISD.tex &> /dev/null
	bibtex VAISD.aux &> /dev/null
	pdflatex VAISD.tex &> /dev/null
	pdflatex VAISD.tex &> /dev/null
	make clean &> /dev/null

view:
	xdg-open VAISD.pdf

clean: 
	rm -rf *.out *.log *.aux *.bbl *.backup *.blg rsc* *.synctex.gz *~
