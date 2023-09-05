# peak annotation with ChIPseeker
rule chip_peakanno:
    input:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
    output:
        pdf = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak_chip_annotation.pdf"
    params:
        gmtfolder = config["annotation"]["GMT"],
    shell:
        "utils/chip_annotate.R {input.peak} {params.gmtfolder};"

# giggle search against CistromeDB
rule chip_giggle_cistrome:
    input:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
    output:
        pdf = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak_giggle_cistrome.pdf"
    params:
        giggleexec = config["annotation"]["giggle"],
        gigglepath = config["annotation"]["giggledb"],
    shell:
        "utils/chip_giggle.R {input.peak} {params.giggleexec} {params.gigglepath};"

# homer
rule chip_homer:
    input:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
    output:
        homeroutput = directory("{OUT_DIR}/Analysis/{name}_HOMER"),
        homeroutputgz = "{OUT_DIR}/Analysis/{name}_HOMER.tar.gz",
    shell:
        """
        sort -k5nr {input.peak} | head -500 | findMotifsGenome.pl - hg38 {output.homeroutput} -size given -mask || echo done! ;
        tar -zcf {output.homeroutputgz} {output.homeroutput};
        """

# run R script on combined peaks
rule report:
    input:
        peakfile = COMBINED_PEAKS, 
	seqqc = SEQ_QC_SUMMARY,
	peakqc = PEAK_QC_SUMMARY,
	gssheatmap = GSSHEATMAP,
	gssprofile = GSSPROFILE,
    output:
        report_html = REPORTHTML,
    params:
        expname = config["outprefix"],
        rmd = REPORTRMD,
        tmpdir = OUT_DIR+"/Tmp/",
	gmtfolder = config["annotation"]["MSigDB"],
        cistromedb = config["annotation"]["giggledb"],
    shell:
        """
	if [[ -d "{params.tmpdir}" ]]; then rm -rf "{params.tmpdir}"; fi
	mkdir {params.tmpdir}
        cp utils/{params.rmd} {params.tmpdir}
        cd {params.tmpdir}
	Rscript -e "library(rmarkdown); rmarkdown::render('{params.rmd}', output_format='html_document', output_file='report.html', params = list( name ='{params.expname}', seqqc = '../../{input.seqqc}', peakqc = '../../{input.peakqc}', gssheatmap = '../../{input.gssheatmap}', gssprofile = '../../{input.gssprofile}', gmtfolder = '{params.gmtfolder}', peakfile = '../../{input.peakfile}', cistromedb = '{params.cistromedb}') )"
        mv report.html ../../{output.report_html}
	mv *.tsv ../Analysis/
	cd ../../
        rm -rf {params.tmpdir}
	"""

