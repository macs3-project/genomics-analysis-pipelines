# macs takes clean bam file to call peaks
rule chip_callpeak:
    input:
        bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.clean.bam"
        cbam = expand("%s/Alignment/{control}.sortedByPos.rmdp.clean.bam" % (OUT_DIR), control={cnames})
    output:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
        xls = "{OUT_DIR}/Analysis/{name}_peaks.xls",
        bdg_raw = "{OUT_DIR}/Analysis/{name}_raw.bdg",
	bdg_spmr = "{OUT_DIR}/Analysis/{name}_spmr.bdg",
    params:
        name = "{name}",
    log:
        "{OUT_DIR}/Log/{name}_macs3_callpeak.log"
    benchmark:
        "{OUT_DIR}/Benchmark/{name}_callpeak.benchmark"
    shell:
        "macs3 callpeak --outdir {OUT_DIR}/Analysis -n {params.name} {macs3_callpeak_option} -t {input.bam} -c {input.cbam}; "
	"mv {OUT_DIR}/Analysis/{params.name}_treat_pileup.bdg {output.bdg_spmr}; "
	"mv {OUT_DIR}/Analysis/{params.name}_treat_pileup.bdg {output.bdg_spmr}; "
	"macs3 pileup {macs3_pileup_option} -i {input.bam} -o {output.bdg_raw}; "

# bdg2bw converts bedGraph to bigWig files
rule chip_bdg2bw:
    input:
        bdg = "{OUT_DIR}/Analysis/{bdgname}.bdg",
    output:
        bw = "{OUT_DIR}/Analysis/{bdgname}.bw",
    params:
        chromlen = config["annotation"]["chromInfo"],
    shell:
        "utils/bdg2bw {input.bdg} {params.chromlen};"

