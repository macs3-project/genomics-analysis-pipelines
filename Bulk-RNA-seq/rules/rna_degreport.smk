rule rna_reporthtml:
    input:
        QUANTGENE1 + QUANTGENE2 + ALIGNMENTQC1 + ALIGNMENTQC2
    output:
        report_html = REPORTHTML,
    params:
        tmpdir = OUT_DIR+"/Tmp/",
        sample = SAMPLEFILE,
	cond1 = config["sample1"],
	cond2 = config["sample2"],
	geneinfo = config["annotation"]["geneinfo"],
        MSigDB = config["annotation"]["MSigDB"],
	fdr = config["options"]["fdr"],
        log2fc = config["options"]["log2fc"],
	species = config["options"]["species"],
    shell:
        """
        cp utils/Template_RNAseq_{params.species}.Rmd {params.tmpdir}/
	cd {params.tmpdir}
        ln -s ../Analysis
        ln -s ../QC
	Rscript -e "library(rmarkdown); rmarkdown::render('Template_RNAseq_{params.species}.Rmd', output_format='html_document', output_file='report.html')" --args {params.cond1} {params.cond2} {params.sample} {params.geneinfo} {params.MSigDB} {params.fdr} {params.log2fc}
	mv report.html {output.report_html}
	mv *.tsv ../Analysis/
	mv *.png ../Analysis/
	unlink Analysis
	unlink QC
	rm -f Template_RNAseq_{params.species}.Rmd
	rm -rf {params.tmpdir}/*
        """
