# ---------- Makefile ----------

# Name the input Rmd and desired output file(s)
RMD  = Response-GlobalMetaAnalysisGenetics.Rmd
PDF  = Response-GlobalMetaAnalysisGenetics.pdf          # optional—remove if you don’t build PDF


# Render PDF similarly (needs TeX)
$(PDF): $(RMD)
	Rscript -e "rmarkdown::render('$(RMD)', output_format = 'pdf_document')"

# Convenience target: build everything with `make`
all: $(HTML) $(PDF)

# Remove generated outputs with `make clean`
clean:
	rm -f $(HTML) $(PDF)

.PHONY: all clean
# ------------------------------
