rule atac_ceas:
    input:
        peak = "Result/Analysis/{fastqid}_peaks.narrowPeak",
        bw = "Result/Analysis/{fastqid}_treat_pileup.bw",
    output:
        ceaspdf = "Result/Analysis/{fastqid}_peaks_ceas.pdf"
    shell:
        ""

rule atac_giggle:
    input:
        peak = "Result/Analysis/{fastqid}_peaks.narrowPeak",
	bw = "Result/Analysis/{fastqid}_treat_pileup.bw",
    output:
        ceaspdf = "Result/Analysis/{fastqid}_peaks_ceas.pdf"
    shell:
        ""
    params:
        name = "{fastqid}",
        chrombed = config["annotation"]["chromBed"],
     	chromlen = config["annotation"]["chromInfo"]
    log:
        "Result/Log/{fastqid}_peak_giggle.log"
    benchmark:
        "Result/Benchmark/{fastqid}_peak_giggle.benchmark"
    shell:
        ""
