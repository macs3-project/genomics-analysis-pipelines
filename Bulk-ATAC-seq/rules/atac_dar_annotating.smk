# run DAR R script
rule atac_dar:
    input:
        btable = BIN_COUNT_TABLE,
	seqqc = SEQ_QC_SUMMARY,
	peakqc = PEAK_QC_SUMMARY,
	fragstat = FRAG_PNG,
	gssheatmap = GSSHEATMAP,
	gssprofile = GSSPROFILE,
        peaks1 = PEAKS1,
        peaks2 = PEAKS2,
    output:
        report_html = REPORTHTML,
	opendar = DAR_OPEN,
	closedar = DAR_CLOSE,
    params:
        expname = config["outprefix"],
        rmd = DARRMD,
        tmpdir = OUT_DIR+"/Tmp/",
        sample = SAMPLEFILE,
	cond1 = config["sample1"],
        cond2 = config["sample2"],
        fdr = config["options"]["dar_fdr"],
        log2fc = config["options"]["dar_log2fc"],
	gmtfolder = config["annotation"]["MSigDB"],
    shell:
        """
	if [[ -d "{params.tmpdir}" ]]; then rm -rf "{params.tmpdir}"; fi
	mkdir {params.tmpdir}
        cp utils/{params.rmd} {params.tmpdir}
        cd {params.tmpdir}
	Rscript -e "library(rmarkdown); rmarkdown::render('{params.rmd}', output_format='html_document', output_file='report.html', params = list( name ='{params.expname}', seqqc = '../../{input.seqqc}', peakqc = '../../{input.peakqc}', fragstat = '../../{input.fragstat}', gssheatmap = '../../{input.gssheatmap}', gssprofile = '../../{input.gssprofile}', btable = '../../{input.btable}', c1name = '{params.cond1}', c2name = '{params.cond2}', metafile = '../../{params.sample}', c1peaks = '{input.peaks1}', c2peaks = '{input.peaks2}', fdr = {params.fdr}, log2fc = {params.log2fc}, gmtfolder = '{params.gmtfolder}') )"
        mv report.html ../../{output.report_html}
	mv *.bed ../Analysis/
	mv *.tsv ../Analysis/
	cd ../../
        rm -rf {params.tmpdir}
	"""

# homer
rule atac_dar_homer:
    input:
        peak = "{OUT_DIR}/Analysis/{name}.dar.{type}.bed",
    output:
        peaktop500 = temp("{OUT_DIR}/Analysis/{name}.dar.{type}.top500.bed"),
        homeroutput = temp("{OUT_DIR}/Analysis/{name}.dar.{type}.HOMER"),
        homeroutputgz = "{OUT_DIR}/Analysis/{name}.dar.{type}.HOMER.tar.gz",
    params:
        genome = config['options']['species'],
	threads = config['options']['cores'],
    threads:
        config["options"]["cores"],
    shell:
        """
	# column 4 is the -log10pvalue score in the DAR file
        sort -k4nr {input.peak} | cut -f 1,2,3 | head -500 > {output.peaktop500};
        findMotifsGenome.pl {output.peaktop500} {params.genome} {output.homeroutput} -size given -mask -p {params.threads};
        tar -zcf {output.homeroutputgz} {output.homeroutput};
        """
