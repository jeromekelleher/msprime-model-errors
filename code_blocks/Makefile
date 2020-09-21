code_blocks.pdf: code_blocks.tex
	pdflatex -shell-escape code_blocks.tex

code_blocks.ps: code_blocks.dvi
	dvips code_blocks

code_blocks.dvi: code_blocks.tex code_blocks.bib
	latex code_blocks.tex
clean:
	rm -f *.log *.dvi *.aux
	rm -f *.blg *.bbl
	rm -f *.eps *.[1-9]
	rm -f src/*.mpx *.mpx

mrproper:
	rm -f *.ps *.pdf
