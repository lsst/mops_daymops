ILLUSTRATIONS=illustrations/mopsDiagram.png illustrations/mopsWithinLsst.png

%.png: %.txt
	ditaa $< $@

baseline.pdf: baseline.tex baseline.bib $(ILLUSTRATIONS)
	latexmk -f -pdf baseline

clean: 
	rm -f baseline.pdf *.dvi *.aux *.bbl *.blg *.aux.bak
