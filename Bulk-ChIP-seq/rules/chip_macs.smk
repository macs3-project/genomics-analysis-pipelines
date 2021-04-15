# macs takes clean bam file to call peaks
rule chip_callpeak:
    input:
        bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.clean.bam",
        cbam = expand("%s/Alignment/{control}.sortedByPos.rmdp.clean.bam" % (OUT_DIR), control=cnames),
    output:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
        xls = "{OUT_DIR}/Analysis/{name}_peaks.xls",
        bdg_raw_t = "{OUT_DIR}/Analysis/{name}_t_raw.bdg",
        bdg_raw_c = "{OUT_DIR}/Analysis/{name}_c_raw.bdg",
    params:
        name = "{name}",
    log:
        "{OUT_DIR}/Log/{name}_macs3_callpeak.log"
    benchmark:
        "{OUT_DIR}/Benchmark/{name}_callpeak.benchmark"
    shell:
        "macs3 callpeak --outdir {OUT_DIR}/Analysis -n {params.name} {macs3_callpeak_option} -t {input.bam} -c {input.cbam}; "
	"mv {OUT_DIR}/Analysis/{params.name}_treat_pileup.bdg {output.bdg_raw_t}; "
	"mv {OUT_DIR}/Analysis/{params.name}_control_lambda.bdg {output.bdg_raw_c}; "

# macs signal tracks
rule chip_signaltracks:
    input:
        bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.clean.bam",
        cbam = expand("%s/Alignment/{control}.sortedByPos.rmdp.clean.bam" % (OUT_DIR), control=cnames),    
        bdg_raw_t = "{OUT_DIR}/Analysis/{name}_t_raw.bdg",
        bdg_raw_c = "{OUT_DIR}/Analysis/{name}_c_raw.bdg",
    output:
        bdg_spmr = "{OUT_DIR}/Analysis/{name}_t_spmr.bdg",
        bdg_logfc = "{OUT_DIR}/Analysis/{name}_logfc.bdg",
        bdg_pscore = "{OUT_DIR}/Analysis/{name}_pscore.bdg",
    params:
        name = "{name}",
    shell:
        "macs3 bdgcmp -t {input.bdg_raw_t} -c {input.bdg_raw_c} -m ppois logfc -o {output.bdg_pscore} {output.bdg_logfc}; "
        "macs3 callpeak --outdir {OUT_DIR}/Analysis -n {params.name}.spmr {macs3_callpeak_option} --SPMR -t {input.bam} -c {input.cbam};"
        "mv {OUT_DIR}/Analysis/{params.name}.spmr_treat_pileup.bdg {output.bdg_spmr}; "
        "rm -f {OUT_DIR}/Analysis/{params.name}.spmr*; "

        
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

