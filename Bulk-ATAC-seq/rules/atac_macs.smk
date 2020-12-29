# macs takes clean bam file to call peaks
rule atac_callpeak:
    input:
        bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.clean.bam" 
    output:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
        bdg = "{OUT_DIR}/Analysis/{name}_treat_pileup.bdg",
        xls = "{OUT_DIR}/Analysis/{name}_peaks.xls",
    params:
        name = "{name}",
    log:
        "{OUT_DIR}/Log/{name}_macs3_peak.log"
    benchmark:
        "{OUT_DIR}/Benchmark/{name}_callpeak.benchmark"
    shell:
        "macs3 callpeak --outdir {OUT_DIR}/Analysis -n {params.name} {macs3_option} -t {input.bam};"

# bdg2bw converts bedGraph to bigWig files
rule atac_bdg2bw:
    input:
        bdg = "{OUT_DIR}/Analysis/{bdgname}.bdg",
    output:
        bw = "{OUT_DIR}/Analysis/{bdgname}.bw",
    params:
        chromlen = config["annotation"]["chromInfo"],
    shell:
        "bdg2bw {input.bdg} {params.chromlen};"