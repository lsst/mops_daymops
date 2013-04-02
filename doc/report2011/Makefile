%.png: %.txt
	ditaa $< $@

baseline.pdf: baseline.tex baseline.bib overview.tex daymops_linking.tex metrics_linking.tex precovery_attribution.tex furtherdev.tex appendices.tex larsenbib.tex
	latexmk -f -pdf baseline

clean: 
	rm -f baseline.pdf *.dvi *.aux *.bbl *.blg *.aux.bak
