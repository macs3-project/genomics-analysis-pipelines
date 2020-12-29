# merge to get consensus
rule atac_consensus:
    input:
        peaks1 = PEAKS1,
        peaks2 = PEAKS2,
    output:
        consensus = CONSENSUS,
    params:
        cutoff = 3
    shell:
        "utils/merge_then_call_consensus.sh {output.consensus} {params.cutoff} {input.peaks1} {input.peaks2}; "

