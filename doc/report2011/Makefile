%.png: %.txt
	ditaa $< $@

baseline.pdf: baseline.tex baseline.bib overview.tex daymops_linking.tex metrics_linking.tex precovery_attribution.tex furtherdev.tex appendices.tex larsenbib.tex
	latexmk -f -pdf baseline

clean: 
	rm -f baseline.pdf *.dvi *.aux *.bbl *.blg *.aux.bak

cover:
	@ echo " "
	@ echo "  Make sure you've updated LDM-156-CoverPage.pdf with "
	@ echo "  the changelog information and re-exported it as "
	@ echo "  LDM-156-CoverPage.pdf."
	@ echo " "

	@ pdfunite LDM-156-CoverPage.pdf baseline.pdf LDM-156.pdf && echo "Done. Docgen in LDM-156.pdf"
