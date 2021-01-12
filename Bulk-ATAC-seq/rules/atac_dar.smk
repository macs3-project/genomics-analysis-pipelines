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
        "perl -a -n -i -e '$F[3]=~s/^.*?(Peak\d+)$/$1/;print join(\"\\t\",$F[0],$F[1],$F[2],$F[3]),\"\\n\";' {output.consensus}; "

# split consensus peak file into bins
rule atac_binconsensus:
    input:
        consensus = CONSENSUS,
    output:
        binconsensus = BIN_CONSENSUS,
    params:
        binsize = 100,
    shell:
        "awk -v W={params.binsize} -v OFS=\"\\t\" '{{print $1, int($2/W)*W, (int($3/W)+1)*W}}' {input.consensus} | bedtools sort -i - | bedtools merge -i - | cut -f 1,2,3 | awk -v OFS=\"\\t\" '{{n+=1;print $1,$2,$3,\"Peak\"n}}' | bedtools makewindows -w {params.binsize} -b - -i srcwinnum > {output.binconsensus}; "

# generate counts for each bigwig
rule atac_bincount:
    input:
        binconsensus = BIN_CONSENSUS,
        bigwig = "{OUT_DIR}/Analysis/{sample}_raw.bw",
    output:
        bincount = "{OUT_DIR}/Analysis/{sample}.bincount.txt",
    shell:
        "bigWigAverageOverBed {input.bigwig} {input.binconsensus} {output.bincount}.tmp; "
        "echo -e \"bin_id\t{wildcards.sample}\" > {output.bincount};"
        "cut -f 1,6 {output.bincount}.tmp >> {output.bincount}; "
        "rm -f {output.bincount}.tmp; "

# make the big table
rule atac_bincounttable:
    input:
        binconsensus = BIN_CONSENSUS,
        bincount1 = BIN_COUNT1,
        bincount2 = BIN_COUNT2,
    output:
        bincounttable = BIN_COUNT_TABLE,
    shell:
        "echo -e \"bin_id\\tpos\" > {output.bincounttable}; "
        "sort -k4,4 {input.binconsensus} | perl -ane 'print \"$F[3]\\t$F[0]:$F[1]-$F[2]\\n\"' - >> {output.bincounttable}; "
        "for f in {input.bincount1};do sort -k1,1 $f | cut -f 2 - | paste -d\"\\t\" {output.bincounttable} - > {output.bincounttable}.tmp; mv {output.bincounttable}.tmp {output.bincounttable}; done; "
        "for f in {input.bincount2};do sort -k1,1 $f | cut -f 2 - | paste -d\"\\t\" {output.bincounttable} - > {output.bincounttable}.tmp; mv {output.bincounttable}.tmp {output.bincounttable}; done; "
        "rm -f {output.bincounttable}.tmp; "

# run DAR R script
rule atac_dar:
    input:
        bincounttable = BIN_COUNT_TABLE,
    output:
        dar_html = DAR_HTML,
    shell:
        "Rscript -e \"knitr::stitch_rmd(\'pipelineDescriptive.Rmd\')\"; "