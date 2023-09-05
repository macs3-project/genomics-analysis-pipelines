# combined peak call of all replicates
rule chip_callpeak_one_replicate:
    input:
        peak = PEAKS[0],
        bw_raw_t = BIGWIG_RAW[0],
        bw_logfc = BIGWIG_LOGFC[0],
        bw_pscore = BIGWIG_PSCORE[0],
        bw_spmr_t = BIGWIG_SPMR[0],        
    output:
        peak = COMBINED_PEAKS,
        bw_raw_t = COMBINED_BIGWIG_RAW,
        bw_logfc = COMBINED_BIGWIG_LOGFC,
        bw_pscore = COMBINED_BIGWIG_PSCORE,
        bw_spmr_t = COMBINED_BIGWIG_SPMR,
    shell:
        """
        cp {input.peak} {output.peak}
        cp {input.bw_raw_t} {output.bw_raw_t}
        cp {input.bw_logfc} {output.bw_logfc}
        cp {input.bw_pscore} {output.bw_pscore}
        cp {input.bw_spmr_t} {output.bw_spmr_t}
        """
