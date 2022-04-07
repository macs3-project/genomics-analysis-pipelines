rule rsem_counting:
    input:
        bam = "{OUT_DIR}/Alignment/{samplename}.transcriptome.bam"
    output:
        count1 = "{OUT_DIR}/Analysis/{samplename}_rsem_gene.txt",
        count2 = "{OUT_DIR}/Analysis/{samplename}_rsem_isoform.txt"
    params:
        rsem_index = config["genome"]["rsemindex"],
        tmpdir = "{OUT_DIR}/Tmp/{samplename}",
	rsem_options = " ".join( (rsempar_common, rsempar_run, rsempar_type ) )
    shell:
        """
        mkdir -p {params.tmpdir}
	cd {params.tmpdir}
	{rsem_cmd} {params.rsem_options} {input.bam} {params.rsem_index} Quant > Log.rsem
	mv Quant.genes.results {output.count1}
	mv Quant.isoforms.results {output.count2}
	cd {OUT_DIR}
	rm -rf {params.tmpdir}
	"""
