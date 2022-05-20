rule star_map:
    input:
        fastq1 = "{OUT_DIR}/Raw/{samplename}",
        fastq2 = "{OUT_DIR}/Raw/{samplename}",
    output:
        bam1 = "{OUT_DIR}/Alignment/{samplename}.genome.bam",
        bam2 = "{OUT_DIR}/Alignment/{samplename}.transcriptome.bam",
	stat = "{OUT_DIR}/QC/{samplename}.STAR.txt"
    params:
        tmpdir = "{OUT_DIR}/Tmp/{samplename}",
        star_options = " ".join( (starpar_common,starpar_run,starpar_bam,starpar_strand ) )
    shell:
        """
        mkdir -p {params.tmpdir}
	cd {params.tmpdir}
        {star_cmd} --readFilesIn {input.fastq1} {input.fastq2} {params.star_options}
        mv Aligned.sortedByCoord.out.bam {output.bam1} 
        mv Aligned.toTranscriptome.out.bam {output.bam2}
	mv Log.final.out {output.stat}
	cd {OUT_DIR}
	rm -rf {params.tmpdir}
	"""
