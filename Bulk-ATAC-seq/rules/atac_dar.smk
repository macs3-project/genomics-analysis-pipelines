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

# split consensus peak file into bins
rule atac_binconsensus:
    input:
        consensus = CONSENSUS,
    output:
        binconsensus = BINCONSENSUS,
    params:
        binsize = 100,
    shell:
        "awk -v W={params.binsize} -v OFS=\"\\t\" '{print $1, int($2/W)*W, (int($3/W)+1)*W}' {input.consensus} | bedtools sort -i - | bedtools merge -i - | bedtools makewindows -w {params.binsize} -b - -i winnum > {output.binconsensus}; "

# generate counts for each bigwig
rule atac_bincount:
    input:
        binconsensus = BINCONSENSUS,
        bigwig = "{OUT_DIR}/Analysis/{bdgname}.bw",
    output:
        bincount = "{OUT_DIR}/Analysis/{bdgname}.bincount.txt",
    shell:
        "bigWigAverageOverBed {input.bigwig} {input.binconsensus} {output.bincount}; "
