# combined peak call of all replicates
rule chip_callpeak_combined:
    input:
        bam = expand("%s/Alignment/{treat}.sortedByPos.rmdp.clean.bam" % (OUT_DIR), treat=tnames),
        cbam = expand("%s/Alignment/{control}.sortedByPos.rmdp.clean.bam" % (OUT_DIR), control=cnames),
    output:
        peak = COMBINED_PEAKS,
        bdg_raw_t = "%s/Analysis/%s_t_raw.bdg" % (OUT_DIR, config["outprefix"]),
        bdg_raw_c = "%s/Analysis/%s_c_raw.bdg" % (OUT_DIR, config["outprefix"]),
	bdg_logfc = "%s/Analysis/%s_logfc.bdg" % (OUT_DIR, config["outprefix"]),
	bdg_pscore = "%s/Analysis/%s_pscore.bdg" % (OUT_DIR, config["outprefix"]),
	bdg_spmr_t = "%s/Analysis/%s_t_spmr.bdg" % (OUT_DIR, config["outprefix"]),
    params:
        name = config["outprefix"],
    shell:
        """
        macs3 callpeak --outdir {OUT_DIR}/Analysis -n {params.name} --cutoff-analysis {macs3_callpeak_option} -t {input.bam} -c {input.cbam};
        mv {OUT_DIR}/Analysis/{params.name}_treat_pileup.bdg {output.bdg_raw_t};
        mv {OUT_DIR}/Analysis/{params.name}_control_lambda.bdg {output.bdg_raw_c};
	macs3 bdgcmp -t {output.bdg_raw_t} -c {output.bdg_raw_c} -p 0.1 -m ppois logFE -o {output.bdg_pscore} {output.bdg_logfc};
        macs3 callpeak --outdir {OUT_DIR}/Analysis -n {params.name}.spmr {macs3_callpeak_option} --SPMR -t {input.bam} -c {input.cbam};
        mv {OUT_DIR}/Analysis/{params.name}.spmr_treat_pileup.bdg {output.bdg_spmr_t};
        rm -f {OUT_DIR}/Analysis/{params.name}.spmr*;
        """
