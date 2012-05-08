ILLUSTRATIONS=illustrations/mopsDiagram.png illustrations/mopsWithinLsst.png newIllustrations/*.png

%.png: %.txt
	ditaa $< $@

baseline.pdf: baseline.tex baseline.bib overview.tex daymops.tex metrics.tex riskreduction.tex appendices.tex $(ILLUSTRATIONS)
	latexmk -f -pdf baseline

clean: 
	rm -f baseline.pdf *.dvi *.aux *.bbl *.blg *.aux.bak
